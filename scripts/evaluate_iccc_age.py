#!/usr/bin/env python3
"""Compare CTC ICCC major x diagnosis age with a pinned G1 count snapshot.

All numerical work uses Python's standard library. Plotly creates the optional
standalone HTML report. Inputs are read only; --out must not already exist.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime
import hashlib
import html
import io
import json
from pathlib import Path
import platform
import shutil
import sys
import zipfile

SCOPES = ('all', 'index', 'nonindex')
AXES = {'single_year': 'ageRecodeSingle90p', 'grouped_age': 'ageRecode20Groups'}


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def checked(path, expected):
    actual = sha(path)
    if actual != expected:
        raise ValueError(f'Hash mismatch: {path}')
    return {'path': str(path.resolve()), 'sha256': actual, 'size_bytes': path.stat().st_size}


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        yield from csv.DictReader(handle)


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def age_band(age):
    if age == '99':
        return '99'
    if len(age) != 2 or not age.isdigit() or not 0 <= int(age) <= 90:
        raise ValueError(f'Invalid single-year age code: {age!r}')
    value = int(age)
    return '00' if value == 0 else ('01' if value < 5 else f'{value // 5 + 1:02d}')


def summarize_counts(ref, generated):
    """Return total variation plus coverage-explicit conditional comparisons."""
    nr, ng = sum(ref.values()), sum(generated.values())
    if not nr or not ng:
        raise ValueError('Empty comparison scope')
    ra, ga, rm, gm = Counter(), Counter(), Counter(), Counter()
    for (major, age), n in ref.items():
        ra[age] += n
        rm[major] += n
    for (major, age), n in generated.items():
        ga[age] += n
        gm[major] += n
    majors, ages = sorted(rm.keys() | gm.keys()), sorted(ra.keys() | ga.keys())
    cells, by_age, by_major = [], [], []
    for age in ages:
        comparable = ra[age] > 0 and ga[age] > 0
        tv = sum(abs(generated[major, age] / ga[age] - ref[major, age] / ra[age])
                 for major in majors) / 2 if comparable else None
        by_age.append({'age_code': age, 'reference_count': ra[age], 'generated_count': ga[age],
                       'reference_pct': 100 * ra[age] / nr, 'generated_pct': 100 * ga[age] / ng,
                       'iccc_given_age_tv': tv, 'comparable': comparable})
        for major in majors:
            r, g = ref[major, age], generated[major, age]
            rp = r / ra[age] if ra[age] else None
            gp = g / ga[age] if ga[age] else None
            cells.append({'iccc_major': major, 'age_code': age, 'reference_count': r,
                          'generated_count': g, 'reference_joint_pct': 100 * r / nr,
                          'reference_age_total': ra[age], 'generated_age_total': ga[age],
                          'generated_joint_pct': 100 * g / ng,
                          'reference_iccc_given_age_pct': None if rp is None else 100 * rp,
                          'generated_iccc_given_age_pct': None if gp is None else 100 * gp,
                          'conditional_difference_pp': 100 * (gp - rp) if comparable else None,
                          'reference_age_given_iccc_pct': 100 * r / rm[major] if rm[major] else None,
                          'generated_age_given_iccc_pct': 100 * g / gm[major] if gm[major] else None,
                          'zero_observed_reference_cell': r == 0 and g > 0})
    for major in majors:
        comparable = rm[major] > 0 and gm[major] > 0
        tv = sum(abs(generated[major, age] / gm[major] - ref[major, age] / rm[major])
                 for age in ages) / 2 if comparable else None
        by_major.append({'iccc_major': major, 'reference_count': rm[major], 'generated_count': gm[major],
                         'reference_pct': 100 * rm[major] / nr, 'generated_pct': 100 * gm[major] / ng,
                         'age_given_iccc_tv': tv, 'comparable': comparable})
    covered = sum(row['generated_count'] for row in by_age if row['comparable'])
    conditional = sum(row['generated_count'] * row['iccc_given_age_tv']
                      for row in by_age if row['comparable']) / covered if covered else None
    metrics = {'reference_tumor_rows': nr, 'generated_tumor_rows': ng,
               'joint_tv': sum(abs(generated[k] / ng - ref[k] / nr) for k in ref.keys() | generated.keys()) / 2,
               'age_marginal_tv': sum(abs(ga[a] / ng - ra[a] / nr) for a in ages) / 2,
               'iccc_marginal_tv': sum(abs(gm[m] / ng - rm[m] / nr) for m in majors) / 2,
               'iccc_given_age_weighted_tv': conditional,
               'conditional_coverage_generated_rows': covered,
               'conditional_coverage_fraction': covered / ng,
               'generated_rows_in_reference_absent_ages': ng - covered,
               'reference_rows_in_generated_absent_ages': sum(ra[a] for a in ages if not ga[a]),
               'generated_rows_in_zero_observed_reference_cells': sum(generated[k] for k in generated if not ref[k])}
    return metrics, cells, by_age, by_major


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workspace', type=Path, required=True)
    ap.add_argument('--snapshot', type=Path, required=True)
    ap.add_argument('--snapshot-sha256', required=True)
    ap.add_argument('--dataset', action='append', required=True, metavar='LABEL=MANIFEST')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--no-html', action='store_true')
    args = ap.parse_args()
    root, snap, out = args.workspace.resolve(), args.snapshot.resolve(), args.out.resolve()
    if out.exists():
        raise FileExistsError(f'Refusing existing output: {out}')
    manifest_path = snap / 'manifest.json'
    provenance = {'reference_manifest': checked(manifest_path, args.snapshot_sha256), 'reference_files': []}
    sm = json.loads(manifest_path.read_text())
    registry = {x['path']: x for x in sm['files']}

    def reference_file(rel):
        path = snap / rel
        provenance['reference_files'].append(checked(path, registry[rel]['sha256']))
        return path

    rules = json.loads(reference_file('inputs/reference/inputs/rules-v2.3.0.json').read_text())
    hierarchy = next(r['mapping'] for r in rules['rules'] if r['id'] == 'ctc_iccc_2017_hierarchy_is_consistent')
    codebook = json.loads(reference_file('inputs/reference/analysis/codebook.json').read_text())['fields']
    labels = {field: {p['value']: p['description'] for p in codebook[field]['permissible_values']
                      if p['kind'] == 'literal'} for field in ['icccMajor', *AXES.values()]}
    references = {}
    for scope in SCOPES:
        for axis, field in AXES.items():
            name = ('grouped_age_' if axis == 'grouped_age' else '') + f'through2022_{scope}.csv'
            counts = Counter()
            for row in read_csv(reference_file(name)):
                n = int(row['count'])
                if n <= 0:
                    raise ValueError('Nonpositive reference count')
                counts[hierarchy[row['iccc2017ExtendedRecode']]['icccMajor'], row[field]] += n
            if sum(counts.values()) != sm['tables'][name]['count_total']:
                raise ValueError(f'Reference total mismatch: {name}')
            references[scope, axis] = counts
        rolled = Counter()
        for (major, age), n in references[scope, 'single_year'].items():
            rolled[major, age_band(age)] += n
        if rolled != references[scope, 'grouped_age']:
            raise ValueError('Reference age rollup mismatch')
    for axis in AXES:
        if references['index', axis] + references['nonindex', axis] != references['all', axis]:
            raise ValueError('Reference index scope sum mismatch')

    datasets, dataset_provenance, checks = {}, {}, {}
    for specification in args.dataset:
        label, manifest_name = specification.split('=', 1)
        if not label or label == 'Reference' or label in datasets:
            raise ValueError(f'Duplicate/reserved/empty dataset label: {label}')
        mp = Path(manifest_name)
        mp = mp.resolve() if mp.is_absolute() else root / mp
        manifest_hash = sha(mp)
        manifest = json.loads(mp.read_text())
        ctc = [f for f in manifest['files'] if f['source'] == 'CTC']
        if len(ctc) != 1:
            raise ValueError('Expected exactly one CTC manifest entry')
        entry = ctc[0]
        zp = Path(entry['local_path'])
        zp = zp if zp.is_absolute() else root / zp
        payload_provenance = checked(zp, entry['sha256'])
        counts = {(scope, axis): Counter() for scope in SCOPES for axis in AXES}
        ids, total, excluded, hierarchy_errors, band_errors = set(), 0, 0, 0, 0
        with zipfile.ZipFile(zp) as archive:
            members = [name for name in archive.namelist() if name.lower().endswith('.csv')]
            if len(members) != 1:
                raise ValueError('Expected exactly one CTC CSV member')
            with archive.open(members[0]) as handle:
                reader = csv.DictReader(io.TextIOWrapper(handle, encoding='utf-8-sig', newline=''))
                required = {'nccrId', 'yearOfDiagnosis', 'indexCancer', 'icccMajor',
                            'iccc2017ExtendedRecode', 'iccc2017Recode', *AXES.values()}
                if not required <= set(reader.fieldnames):
                    raise ValueError('Missing required CTC columns')
                for row in reader:
                    total += 1
                    ids.add(row['nccrId'])
                    h = hierarchy.get(row['iccc2017ExtendedRecode'])
                    hierarchy_errors += h is None or any(row[k] != h[k] for k in ('icccMajor', 'iccc2017Recode'))
                    band_errors += age_band(row['ageRecodeSingle90p']) != row['ageRecode20Groups']
                    year = row['yearOfDiagnosis']
                    if not (len(year) == 4 and year.isdigit() and 1995 <= int(year) <= 2022):
                        excluded += 1
                        continue
                    scope = 'index' if row['indexCancer'] in ('1', '2') else 'nonindex'
                    for axis, field in AXES.items():
                        for subset in ('all', scope):
                            counts[subset, axis][row['icccMajor'], row[field]] += 1
        if total != entry['rows'] or len(ids) != entry['patients']:
            raise ValueError('Generated count does not match manifest')
        checked(zp, entry['sha256'])
        checked(mp, manifest_hash)
        datasets[label] = counts
        dataset_provenance[label] = {'manifest': checked(mp, manifest_hash), 'ctc': payload_provenance,
                                     'release_id': manifest['release_id'], 'kind': manifest['kind'], 'csv_member': members[0]}
        checks[label] = {'total_rows': total, 'distinct_patients': len(ids), 'excluded_year_rows': excluded,
                         'hierarchy_mismatch_rows': hierarchy_errors, 'age_band_mismatch_rows': band_errors,
                         'manifest_counts_match': True, 'ctc_hash_before_and_after_match': True}

    all_cells, all_ages, all_majors, summaries = [], [], [], []
    for label, counts in datasets.items():
        for scope in SCOPES:
            for axis, field in AXES.items():
                prefix = {'dataset': label, 'scope': scope, 'age_axis': axis}
                metrics, cells, ages, majors = summarize_counts(references[scope, axis], counts[scope, axis])
                summaries.append(prefix | metrics)
                for rows, output in [(cells, all_cells), (ages, all_ages), (majors, all_majors)]:
                    for row in rows:
                        if 'age_code' in row:
                            row['age_label'] = labels[field].get(row['age_code'], row['age_code'] or '(blank)')
                        if 'iccc_major' in row:
                            row['iccc_label'] = labels['icccMajor'].get(row['iccc_major'], row['iccc_major'] or '(blank)')
                        output.append(prefix | row)
    limitations = [
        'Reference is the frozen STAGE_TEST_ONLY v30 file-derived G1 snapshot; clinical representativeness is unverified.',
        'Counts are tumor rows, not unique patients. Year scope is 1995-2022 inclusive for both sources. Unknown/unclassified codes are retained.',
        'Index includes indexCancer 1 or 2 (including benign/borderline); nonindex is the complement. Supplemental comparisons are descriptive, not the r4 sampler target.',
        'R4 preserves r3 ages. Its two-stage conditional sampler uses all-tumor single-year reference counts, with shared patient sex/race and declared fallback for later tumors.',
        'Reference is also a calibration input to r4; this is an in-reference fidelity diagnostic, not held-out validation.',
        'Grouped-age conditional differences can include different single-year mixes within a band. Single-year conditional results are provided separately.',
        'Total variation (TV) is half the sum of absolute probability differences; 0 is identical and 1 is disjoint. It is not a p-value or a row-level error rate.',
        'Conditional TV weights each shared age by generated row count and renormalizes over covered rows; coverage and unsupported ages are reported separately.',
        'A zero observed reference cell is not a biological impossibility. Sparse cells are descriptive; no significance or clinical acceptance threshold is applied.',
        'Evaluation status: Partially evaluated. Hashes, counts, rollups and metrics are checked; clinical validity and release eligibility are not assessed.'
    ]
    out.mkdir(parents=True)
    (out / 'inputs').mkdir()
    shutil.copy2(Path(__file__), out / 'evaluate_iccc_age.py')
    shutil.copy2(manifest_path, out / 'inputs/reference-manifest.json')
    for i, p in enumerate(dataset_provenance.values(), 1):
        shutil.copy2(p['manifest']['path'], out / f'inputs/dataset-{i}-manifest.json')
    write_csv(out / 'cell_comparison.csv', all_cells)
    write_csv(out / 'by_age.csv', all_ages)
    write_csv(out / 'by_iccc_major.csv', all_majors)
    write_csv(out / 'metrics.csv', summaries)
    result = {'schema_version': 1, 'created_at': datetime.now().astimezone().isoformat(),
              'status': 'descriptive_comparison_complete', 'frontier_status': 'Partially evaluated',
              'scope': 'CTC tumor rows; diagnosis years 1995-2022; exact string codes; no imputation',
              'reference_snapshot_id': sm['snapshot_id'], 'provenance': provenance,
              'datasets': dataset_provenance, 'input_checks': checks, 'metrics': summaries,
              'limitations': limitations, 'runtime': {'python': platform.python_version(), 'executable': sys.executable},
              'command': sys.argv, 'script_sha256': sha(Path(__file__))}
    lines = ['# ICCC major group × age at diagnosis', '', 'Comparison of generated data with the frozen G1 reference.', '',
             'Both sides use diagnosis years 1995–2022. Primary unit: tumor row; index/nonindex views are supplemental.', '',
             '| Dataset | Tumor rows | Joint TV (age bands) | Age marginal TV | ICCC mix TV at each single-year age* |',
             '|---|---:|---:|---:|---:|']
    for label in datasets:
        g = next(m for m in summaries if m['dataset'] == label and m['scope'] == 'all' and m['age_axis'] == 'grouped_age')
        s = next(m for m in summaries if m['dataset'] == label and m['scope'] == 'all' and m['age_axis'] == 'single_year')
        lines.append(f"| {label} | {g['generated_tumor_rows']:,} | {g['joint_tv']:.4f} | {g['age_marginal_tv']:.4f} | {s['iccc_given_age_weighted_tv']:.4f} |")
    lines += ['', '*Weighted by generated age frequencies, over shared support. Lower is closer; inspect coverage in metrics.csv.', '',
              f"Reference: {sum(references['all', 'single_year'].values()):,} tumor rows.", '',
              '## Single-year reference coverage', '',
              '| Dataset | Covered generated rows | Rows at ages absent from reference |',
              '|---|---:|---:|']
    for label in datasets:
        s = next(m for m in summaries if m['dataset'] == label and m['scope'] == 'all' and m['age_axis'] == 'single_year')
        lines.append(f"| {label} | {s['conditional_coverage_generated_rows']:,} ({s['conditional_coverage_fraction']:.4%}) | {s['generated_rows_in_reference_absent_ages']:,} |")
    lines += ['', 'These rows remain in joint/age comparisons but cannot receive a reference conditional comparison.', '',
              '## Largest grouped-age conditional differences', '']
    for label in datasets:
        rows = sorted((r for r in all_cells if r['dataset'] == label and r['scope'] == 'all' and r['age_axis'] == 'grouped_age'
                       and r['conditional_difference_pp'] is not None), key=lambda r: abs(r['conditional_difference_pp']), reverse=True)[:6]
        lines += [f'### {label}', '', '| ICCC major | Age | Reference within age | Generated within age | Difference | Reference / generated age rows |', '|---|---|---:|---:|---:|---:|']
        for row in rows:
            lines.append(f"| {row['iccc_major']} {row['iccc_label']} | {row['age_label']} | {row['reference_iccc_given_age_pct']:.2f}% | {row['generated_iccc_given_age_pct']:.2f}% | {row['conditional_difference_pp']:+.2f} pp | {row['reference_age_total']:,} / {row['generated_age_total']:,} |")
        lines.append('')
    lines += ['## Interpretation and limits', ''] + ['- ' + s for s in limitations]
    lines += ['', '## Artifacts', '', '- `report.html`: interactive offline charts.', '- `cell_comparison.csv`: counts, joint percentages, both conditional percentages and differences.',
              '- `by_age.csv`, `by_iccc_major.csv`, `metrics.csv`: all/index/nonindex views at both age resolutions.',
              '- `evaluation.json`: checks, exact input hashes, definitions, code/runtime provenance and limitations.',
              '- `inputs/`: byte-preserved input manifests; payloads remain at their original paths.', '']
    (out / 'summary.md').write_text('\n'.join(lines))
    if not args.no_html:
        render_html(out, datasets, references, labels, summaries, limitations)
        import plotly
        result['runtime']['plotly'] = plotly.__version__
    dump(out / 'evaluation.json', result)
    dump(out / 'completion.json', {'status': 'complete', 'completed_at': datetime.now().astimezone().isoformat(),
                                   'outputs': [{'path': str(p.relative_to(out)), 'sha256': sha(p)}
                                               for p in sorted(out.rglob('*')) if p.is_file()]})
    print(json.dumps({'out': str(out), 'input_checks': checks,
                      'all_tumor_metrics': [r for r in summaries if r['scope'] == 'all']}, indent=2))


def render_html(out, datasets, references, labels, summaries, limitations):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    axis = 'grouped_age'
    pools = {'Reference': references['all', axis], **{k: v['all', axis] for k, v in datasets.items()}}
    majors = sorted({m for c in pools.values() for m, a in c})
    ages = sorted({a for c in pools.values() for m, a in c})
    xlabels = [labels[AXES[axis]].get(a, a).replace(' years', '') for a in ages]
    short = {'01':'I. Leukemias', '02':'II. Lymphomas', '03':'III. CNS', '04':'IV. Neuroblastoma',
             '05':'V. Retinoblastoma', '06':'VI. Renal', '07':'VII. Hepatic', '08':'VIII. Bone',
             '09':'IX. Soft tissue', '10':'X. Germ cell', '11':'XI. Epithelial / melanoma',
             '12':'XII. Other malignant', '99':'99. Unclassified / in situ'}
    ylabels = [short.get(m, m) for m in majors]

    def matrix(counts):
        totals = {a: sum(counts[m, a] for m in majors) for a in ages}
        return [[100 * counts[m, a] / totals[a] if totals[a] else None for a in ages] for m in majors]

    fragments = []
    def add(fig):
        # Age labels such as "01-04" resemble dates to Plotly's auto-detection.
        fig.update_xaxes(type='category', categoryorder='array', categoryarray=xlabels)
        fig.update_layout(template='plotly_white', font={'family':'Arial', 'size':12}, margin={'l':30,'r':30,'t':65,'b':65})
        fragments.append(fig.to_html(full_html=False, include_plotlyjs=True if not fragments else False,
                                      config={'displaylogo':False, 'toImageButtonOptions':{'format':'svg'}}))

    fig = make_subplots(rows=1, cols=len(pools), subplot_titles=list(pools), shared_yaxes=True, horizontal_spacing=.04)
    for i, (label, counts) in enumerate(pools.items(), 1):
        fig.add_trace(go.Heatmap(z=matrix(counts), x=xlabels, y=ylabels, coloraxis='coloraxis',
                                customdata=[[[counts[m,a],sum(counts[k,a] for k in majors)] for a in ages] for m in majors],
                                hovertemplate='%{y}<br>Age %{x}<br>%{z:.2f}% within age<br>%{customdata[0]} / %{customdata[1]} rows<extra></extra>'), row=1, col=i)
    fig.update_yaxes(autorange='reversed')
    fig.update_xaxes(tickangle=-65)
    fig.update_layout(height=590, title='ICCC mix within each diagnosis-age group',
                      coloraxis={'colorscale':'Blues', 'cmin':0, 'cmax':100, 'colorbar':{'title':'% within age'}})
    add(fig)
    refm = matrix(pools['Reference'])
    deltas = {}
    for label in datasets:
        m = matrix(pools[label])
        deltas[label] = [[v - r if v is not None and r is not None else None for v, r in zip(vrow, rrow)]
                         for vrow, rrow in zip(m, refm)]
    bound = max(abs(x) for m in deltas.values() for row in m for x in row if x is not None) or 1
    fig = make_subplots(rows=1, cols=len(datasets), subplot_titles=[f'{k} minus reference' for k in datasets], shared_yaxes=True)
    for i, (label, m) in enumerate(deltas.items(), 1):
        fig.add_trace(go.Heatmap(z=m, x=xlabels, y=ylabels, coloraxis='coloraxis',
                                hovertemplate='%{y}<br>Age %{x}<br>%{z:+.2f} percentage points<extra></extra>'), row=1, col=i)
    fig.update_yaxes(autorange='reversed')
    fig.update_xaxes(tickangle=-65)
    fig.update_layout(height=560, title='Difference in ICCC mix within age groups',
                      coloraxis={'colorscale':'RdBu_r', 'cmin':-bound, 'cmax':bound, 'colorbar':{'title':'pp'}})
    add(fig)
    colors = {'Reference':'#273449', **dict(zip(datasets, ['#d97706','#0891b2','#8b5cf6','#16a34a']))}
    fig = go.Figure()
    for label, counts in pools.items():
        total = sum(counts.values())
        fig.add_trace(go.Bar(name=label, x=xlabels, y=[100 * sum(counts[m, a] for m in majors) / total for a in ages],
                            marker_color=colors.get(label), hovertemplate='Age %{x}<br>%{y:.2f}% of tumor rows<extra>%{fullData.name}</extra>'))
    fig.update_layout(title='Diagnosis-age distribution across all tumor rows', height=420, barmode='group', yaxis_title='% of tumor rows')
    add(fig)
    fig = make_subplots(rows=(len(majors)+2)//3, cols=3, subplot_titles=ylabels, vertical_spacing=.065)
    for j, major in enumerate(majors):
        for label, counts in pools.items():
            total = sum(counts[major, a] for a in ages)
            fig.add_trace(go.Scatter(name=label, x=xlabels, y=[100 * counts[major, a] / total if total else None for a in ages],
                                    mode='lines+markers', line={'color':colors.get(label)}, legendgroup=label, showlegend=j==0,
                                    hovertemplate='Age %{x}<br>%{y:.2f}% within ICCC group<extra>%{fullData.name}</extra>'), row=j//3+1, col=j%3+1)
    fig.update_xaxes(tickangle=-65, nticks=10)
    fig.update_yaxes(title_text='% within ICCC')
    fig.update_layout(height=1500, title='Age at diagnosis within each ICCC major group')
    add(fig)
    metrics_html = '<table><tr><th>Dataset</th><th>Rows</th><th>Joint TV (bands)</th><th>Age TV (bands)</th><th>ICCC | single-year age TV</th></tr>'
    for label in datasets:
        rows = [m for m in summaries if m['dataset']==label and m['scope']=='all']
        g, s = next(m for m in rows if m['age_axis']=='grouped_age'), next(m for m in rows if m['age_axis']=='single_year')
        metrics_html += f"<tr><td>{html.escape(label)}</td><td>{g['generated_tumor_rows']:,}</td><td>{g['joint_tv']:.4f}</td><td>{g['age_marginal_tv']:.4f}</td><td>{s['iccc_given_age_weighted_tv']:.4f}</td></tr>"
    metrics_html += '</table>'
    body = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
    <title>ICCC major × diagnosis age evaluation</title><style>body{font:16px Arial,sans-serif;color:#172033;margin:32px auto;max-width:1500px;padding:0 24px;background:#f9fafb}h1{font-size:34px}p,li{line-height:1.55}table{border-collapse:collapse;width:100%;background:white}th,td{padding:14px;text-align:left;border-bottom:1px solid #dce2e9}th{background:#e8eef4}.chart{background:white;border:1px solid #e2e8f0;border-radius:8px;margin:24px 0;padding:8px}li{margin-bottom:7px}</style>
    <h1>ICCC major group × age at diagnosis</h1><p>Generated r3 release and r4 candidate compared with the frozen G1 reference. Diagnosis years 1995–2022; tumor rows; unknown and unclassified codes retained.</p>
    <p>Lower total variation (TV) means closer agreement: 0 is identical, 1 is disjoint. Conditional TV compares the ICCC mix at each single-year age, weighted by generated age frequencies. This is descriptive agreement with a calibration reference.</p>'''
    body += metrics_html + ''.join('<div class="chart">'+f+'</div>' for f in fragments)
    body += '<h2>Single-year reference coverage</h2><ul>'
    for label in datasets:
        s = next(m for m in summaries if m['dataset']==label and m['scope']=='all' and m['age_axis']=='single_year')
        body += f"<li>{html.escape(label)}: {s['conditional_coverage_generated_rows']:,} of {s['generated_tumor_rows']:,} rows covered ({s['conditional_coverage_fraction']:.4%}); {s['generated_rows_in_reference_absent_ages']:,} rows at ages absent from reference. They remain in joint/age comparisons.</li>"
    body += '</ul>'
    body += '<h2>Interpretation and limits</h2><ul>' + ''.join('<li>'+html.escape(s)+'</li>' for s in limitations) + '</ul>'
    body += '<p>Full code labels, counts, scope-specific metrics and input hashes are available in the accompanying CSV and JSON files. Hover for values; click legend items to toggle a dataset. Charts can be exported as SVG.</p></html>'
    (out / 'report.html').write_text(body)


if __name__ == '__main__':
    main()
