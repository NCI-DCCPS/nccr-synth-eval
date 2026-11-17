#!/usr/bin/env python3
"""Compare ICCC major x diagnosis age across three row-level CSV files.

Python 3.11+. Plotly is required for charts; --no-html uses the standard library.
No workspace, snapshot, release manifest, or patient identifier is required.
Each CSV row counts as one tumor. Source CSVs are opened read-only.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime
import hashlib
import html
import json
from pathlib import Path
import platform
import shutil
import sys

__version__ = '0.1.0'

AXES = {'single_year': 'ageRecodeSingle90p', 'grouped_age': 'ageRecode20Groups'}
# Display labels from the frozen NCCR codebook, not a classification validator.
# Original G1 codebook is hash-bound in its reference snapshot manifest.



LABELS = {'icccMajor': {'01': 'I. Leukemias, Myeloproliferative And Myelodysplastic Diseases',
               '02': 'II. Lymphomas and reticuloendothelial neoplasms',
               '03': 'III. CNS and Miscellaneous Intracranial and Intraspinal Neoplasms',
               '04': 'IV. Neuroblastoma And Other Peripheral Nervous Cell Tumors',
               '05': 'V. Retinoblastoma',
               '06': 'VI. Renal Tumors',
               '07': 'VII. Hepatic Tumors',
               '08': 'VIII. Malignant Bone Tumors',
               '09': 'IX. Soft Tissue And Other Extraosseous Sarcomas',
               '10': 'X. Germ Cell Tumors, Trophoblastic Tumors And Neoplasms Of Gonads',
               '11': 'XI. Other Malignant Epithelial Neoplasms And Malignant Melanomas',
               '12': 'XII. Other And Unspecified Malignant Neoplasms',
               '99': 'Not classified by SEER or in situ',
               '': 'Blank / missing'},
 'ageRecodeSingle90p': {'00': '00 years',
                        '01': '01 years',
                        '02': '02 years',
                        '03': '03 years',
                        '04': '04 years',
                        '05': '05 years',
                        '06': '06 years',
                        '07': '07 years',
                        '08': '08 years',
                        '09': '09 years',
                        '10': '10 years',
                        '11': '11 years',
                        '12': '12 years',
                        '13': '13 years',
                        '14': '14 years',
                        '15': '15 years',
                        '16': '16 years',
                        '17': '17 years',
                        '18': '18 years',
                        '19': '19 years',
                        '20': '20 years',
                        '21': '21 years',
                        '22': '22 years',
                        '23': '23 years',
                        '24': '24 years',
                        '25': '25 years',
                        '26': '26 years',
                        '27': '27 years',
                        '28': '28 years',
                        '29': '29 years',
                        '30': '30 years',
                        '31': '31 years',
                        '32': '32 years',
                        '33': '33 years',
                        '34': '34 years',
                        '35': '35 years',
                        '36': '36 years',
                        '37': '37 years',
                        '38': '38 years',
                        '39': '39 years',
                        '40': '40 years',
                        '41': '41 years',
                        '42': '42 years',
                        '43': '43 years',
                        '44': '44 years',
                        '45': '45 years',
                        '46': '46 years',
                        '47': '47 years',
                        '48': '48 years',
                        '49': '49 years',
                        '50': '50 years',
                        '51': '51 years',
                        '52': '52 years',
                        '53': '53 years',
                        '54': '54 years',
                        '55': '55 years',
                        '56': '56 years',
                        '57': '57 years',
                        '58': '58 years',
                        '59': '59 years',
                        '60': '60 years',
                        '61': '61 years',
                        '62': '62 years',
                        '63': '63 years',
                        '64': '64 years',
                        '65': '65 years',
                        '66': '66 years',
                        '67': '67 years',
                        '68': '68 years',
                        '69': '69 years',
                        '70': '70 years',
                        '71': '71 years',
                        '72': '72 years',
                        '73': '73 years',
                        '74': '74 years',
                        '75': '75 years',
                        '76': '76 years',
                        '77': '77 years',
                        '78': '78 years',
                        '79': '79 years',
                        '80': '80 years',
                        '81': '81 years',
                        '82': '82 years',
                        '83': '83 years',
                        '84': '84 years',
                        '85': '85 years',
                        '86': '86 years',
                        '87': '87 years',
                        '88': '88 years',
                        '89': '89 years',
                        '90': '90+ years',
                        '99': 'Unknown',
                        '': 'Blank / missing'},
 'ageRecode20Groups': {'00': '00 years',
                       '01': '01-04 years',
                       '02': '05-09 years',
                       '03': '10-14 years',
                       '04': '15-19 years',
                       '05': '20-24 years',
                       '06': '25-29 years',
                       '07': '30-34 years',
                       '08': '35-39 years',
                       '09': '40-44 years',
                       '10': '45-49 years',
                       '11': '50-54 years',
                       '12': '55-59 years',
                       '13': '60-64 years',
                       '14': '65-69 years',
                       '15': '70-74 years',
                       '16': '75-79 years',
                       '17': '80-84 years',
                       '18': '85-89 years',
                       '19': '90+ years',
                       '99': 'Unknown',
                       '': 'Blank / missing'}}


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def age_band(age):
    if age == '':
        return ''
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


def format_metric(value):
    return 'N/A' if value is None else f'{value:.4f}'


def render_html(out, datasets, references, labels, summaries, limitations, scope_text):
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
    ylabels = [short.get(m, labels['icccMajor'].get(m, m)) for m in majors]

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
    bound = max((abs(x) for m in deltas.values() for row in m for x in row if x is not None), default=1) or 1
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
        metrics_html += f"<tr><td>{html.escape(label)}</td><td>{g['generated_tumor_rows']:,}</td><td>{g['joint_tv']:.4f}</td><td>{g['age_marginal_tv']:.4f}</td><td>{format_metric(s['iccc_given_age_weighted_tv'])}</td></tr>"
    metrics_html += '</table>'
    body = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
    <title>ICCC major × diagnosis age evaluation</title><style>body{font:16px Arial,sans-serif;color:#172033;margin:32px auto;max-width:1500px;padding:0 24px;background:#f9fafb}h1{font-size:34px}p,li{line-height:1.55}table{border-collapse:collapse;width:100%;background:white}th,td{padding:14px;text-align:left;border-bottom:1px solid #dce2e9}th{background:#e8eef4}.chart{background:white;border:1px solid #e2e8f0;border-radius:8px;margin:24px 0;padding:8px}li{margin-bottom:7px}</style>
    <h1>ICCC major group × age at diagnosis</h1><p>CSV comparison with the supplied reference. SCOPE_PLACEHOLDER; tumor rows; unknown and blank codes retained.</p>
    <p>Lower total variation (TV) means closer agreement: 0 is identical, 1 is disjoint. Conditional TV compares the ICCC mix at each single-year age, weighted by generated age frequencies. This is descriptive agreement with the supplied reference.</p>'''
    body = body.replace('SCOPE_PLACEHOLDER', html.escape(scope_text))
    body += metrics_html + ''.join('<div class="chart">'+f+'</div>' for f in fragments)
    body += '<h2>Single-year reference coverage</h2><ul>'
    for label in datasets:
        s = next(m for m in summaries if m['dataset']==label and m['scope']=='all' and m['age_axis']=='single_year')
        body += f"<li>{html.escape(label)}: {s['conditional_coverage_generated_rows']:,} of {s['generated_tumor_rows']:,} rows covered ({s['conditional_coverage_fraction']:.4%}); {s['generated_rows_in_reference_absent_ages']:,} rows at ages absent from reference. They remain in joint/age comparisons.</li>"
    body += '</ul>'
    body += '<h2>Interpretation and limits</h2><ul>' + ''.join('<li>'+html.escape(s)+'</li>' for s in limitations) + '</ul>'
    body += '<p>Full code labels, counts, scope-specific metrics and input hashes are available in the accompanying CSV and JSON files. Hover for values; click legend items to toggle a dataset. Charts can be exported as SVG.</p></html>'
    (out / 'report.html').write_text(body, encoding='utf-8')


def canonical_code(value, allowed, field, line):
    """Accept padded/unpadded integer codes; retain blanks and count repairs."""
    normalized = value.strip()
    if normalized:
        if not normalized.isascii() or not normalized.isdigit() or len(normalized) > 2:
            raise ValueError(f'CSV line {line}: {field} must be a one/two-digit code, 99, or blank')
        normalized = normalized.zfill(2)
        if normalized not in allowed:
            raise ValueError(f'CSV line {line}: unsupported {field} code {value!r}')
    return normalized


def csv_header(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        fields = next(csv.reader(handle), None)
    if not fields or len(fields) != len(set(fields)):
        raise ValueError(f'{path}: empty or duplicate CSV headers')
    return fields


def count_csv(path, args, include_index):
    before = sha(path)
    header = csv_header(path)
    required = {args.iccc_column, args.age_column}
    filtered = args.min_year is not None or args.max_year is not None
    if filtered:
        required.add(args.year_column)
    if not required <= set(header):
        raise ValueError(f'{path}: missing columns: {sorted(required - set(header))}')
    scopes = ('all', 'index', 'nonindex') if include_index else ('all',)
    counts = {(scope, axis): Counter() for scope in scopes for axis in AXES}
    total = excluded = invalid_year = band_mismatches = 0
    normalizations = Counter()
    major_codes = {f'{n:02d}' for n in range(1, 13)} | {'99'}
    age_codes = {f'{n:02d}' for n in range(91)} | {'99'}
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            total += 1
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'{path}: CSV line {reader.line_num} has a different number of fields than the header')
            if filtered:
                year = row[args.year_column].strip()
                if not year.isascii() or not year.isdigit() or len(year) != 4:
                    excluded += 1
                    invalid_year += 1
                    continue
                if ((args.min_year is not None and int(year) < args.min_year)
                        or (args.max_year is not None and int(year) > args.max_year)):
                    excluded += 1
                    continue
            major = canonical_code(row[args.iccc_column], major_codes, args.iccc_column, reader.line_num)
            age = canonical_code(row[args.age_column], age_codes, args.age_column, reader.line_num)
            normalizations[args.iccc_column] += major != row[args.iccc_column]
            normalizations[args.age_column] += age != row[args.age_column]
            band = age_band(age)
            if args.age_group_column in header:
                existing_band = canonical_code(row[args.age_group_column],
                    {f'{n:02d}' for n in range(20)} | {'99'}, args.age_group_column, reader.line_num)
                band_mismatches += existing_band != band
            selected_scopes = ['all']
            if include_index:
                code = row[args.index_column].strip()
                selected_scopes.append('index' if code in ('1', '2') else 'nonindex')
            for scope in selected_scopes:
                counts[scope, 'single_year'][major, age] += 1
                counts[scope, 'grouped_age'][major, band] += 1
    if before != sha(path):
        raise ValueError(f'{path}: input changed during counting')
    if not sum(counts['all', 'single_year'].values()):
        raise ValueError(f'{path}: no rows remain after filtering')
    if band_mismatches:
        raise ValueError(f'{path}: {band_mismatches} included rows have {args.age_group_column} inconsistent with {args.age_column}')
    return counts, {'path':str(path), 'sha256':before, 'size_bytes':path.stat().st_size,
                    'input_rows':total, 'included_rows':total-excluded, 'excluded_year_rows':excluded,
                    'invalid_or_missing_year_rows':invalid_year, 'code_normalization_counts':dict(normalizations),
                    'age_group_check':'passed' if args.age_group_column in header else 'derived_from_single_year_age',
                    'hash_before_and_after_match':True}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version=f'iccc-age {__version__}')
    parser.add_argument('reference_csv', type=Path, help='Reference: one tumor per CSV row')
    parser.add_argument('generated_a_csv', type=Path)
    parser.add_argument('generated_b_csv', type=Path)
    parser.add_argument('--labels', nargs=2, default=['Generated A', 'Generated B'], metavar=('LABEL_A','LABEL_B'))
    parser.add_argument('--out', type=Path, required=True, help='New output directory (must not exist)')
    parser.add_argument('--iccc-column', default='icccMajor')
    parser.add_argument('--age-column', default='ageRecodeSingle90p', help='Diagnosis age in complete years: 0–90, 90 means 90+, 99 unknown, blank retained')
    parser.add_argument('--age-group-column', default='ageRecode20Groups', help='Optional existing grouped-age column to cross-check')
    parser.add_argument('--year-column', default='yearOfDiagnosis')
    parser.add_argument('--min-year', type=int, help='Optional inclusive lower diagnosis-year limit')
    parser.add_argument('--max-year', type=int, help='Optional inclusive upper diagnosis-year limit')
    parser.add_argument('--index-column', default='indexCancer', help='Optional in all three CSVs: 1/2 index, other values nonindex')
    parser.add_argument('--no-html', action='store_true', help='Write numerical outputs without Plotly')
    args = parser.parse_args(argv)
    paths = [p.resolve() for p in (args.reference_csv,args.generated_a_csv,args.generated_b_csv)]
    out = args.out.resolve()
    if out.exists():
        parser.error(f'Refusing existing output directory: {out}')
    if len(set(args.labels)) != 2 or any(not x.strip() or x == 'Reference' for x in args.labels):
        parser.error('Use two distinct nonempty labels other than Reference')
    if args.min_year is not None and args.max_year is not None and args.min_year > args.max_year:
        parser.error('--min-year cannot exceed --max-year')
    if args.iccc_column == args.age_column:
        parser.error('ICCC and age columns must differ')
    if not args.no_html:
        try:
            import plotly
        except ImportError:
            parser.error('Charts require Plotly: python3 -m pip install plotly (or use --no-html)')
    headers = [csv_header(path) for path in paths]
    include_index = all(args.index_column in h for h in headers)
    pools, provenance = {}, {}
    for label, path in zip(['Reference', *args.labels], paths):
        print(f'Counting {label}: {path}', file=sys.stderr, flush=True)
        pools[label], provenance[label] = count_csv(path, args, include_index)
    references = pools.pop('Reference')
    datasets = pools
    summaries, cells, ages, majors, skipped = [], [], [], [], []
    for label, counts in datasets.items():
        for scope, axis in counts:
            if not sum(references[scope,axis].values()) or not sum(counts[scope,axis].values()):
                skipped.append({'dataset':label,'scope':scope,'age_axis':axis,'reason':'Empty reference or generated subset'})
                continue
            prefix = {'dataset':label,'scope':scope,'age_axis':axis}
            metric, c, a, m = summarize_counts(references[scope,axis], counts[scope,axis])
            summaries.append(prefix | metric)
            for rows, target in [(c,cells),(a,ages),(m,majors)]:
                for row in rows:
                    if 'age_code' in row:
                        row['age_label'] = LABELS[AXES[axis]][row['age_code']]
                    if 'iccc_major' in row:
                        row['iccc_label'] = LABELS['icccMajor'][row['iccc_major']]
                    target.append(prefix | row)
    filtered = args.min_year is not None or args.max_year is not None
    scope_text = (f'Diagnosis years {args.min_year if args.min_year is not None else "earliest"}–{args.max_year if args.max_year is not None else "latest"}'
                  if filtered else 'All input rows; no diagnosis-year filter')
    limitations = [
        scope_text + '. Each included row counts as one tumor, without patient deduplication or weights.',
        'Only two-digit NCCR ICCC major and diagnosis-age encodings are interpreted. Whitespace and one-digit numeric codes are normalized and counted in input checks.',
        'ICCC 99 means unclassified/in situ; age 99 means unknown and 90 means 90+. Blank ICCC/age values are separate categories, not zero or unknown.',
        'Age bands are derived from single-year diagnosis age. If the grouped-age column is present, inconsistent included rows cause an error.',
        'Year filtering excludes out-of-range and invalid/missing diagnosis years; counts are reported for each input.' if filtered else 'Diagnosis-year values are not used because no year limits were requested.',
        'Index/nonindex views use indexCancer 1 or 2 versus other values.' if include_index else 'Index/nonindex views are unavailable because the optional index column is absent from at least one input.',
        'Total variation (TV) is half the sum of absolute probability differences. Zero is identical, one is disjoint; it is not a p-value or a row-level error rate.',
        'Conditional TV weights each shared age by generated row count; rows at reference-absent ages are excluded from that metric and reported separately. N/A means no shared age support.',
        'Grouped-age conditional differences may reflect different single-year mixes within a band; exact-year metrics are also supplied.',
        'Unknown and blank codes are retained. Zero observed reference cells do not establish biological impossibility.',
        'Clinical validity, ICCC hierarchy, reference representativeness, and release eligibility are not assessed. Calibration overlap is not inferred from CSV files.',
        'Evaluation status: Partially evaluated; descriptive comparison with input checks, no clinical acceptance threshold.'
    ]
    out.mkdir(parents=True)
    shutil.copy2(Path(__file__), out / 'compare_iccc_age_csv.py')
    for filename, rows in [('cell_comparison.csv',cells),('by_age.csv',ages),('by_iccc_major.csv',majors),('metrics.csv',summaries)]:
        write_csv(out/filename, rows)
    runtime = {'python':platform.python_version(), 'evaluator':__version__}
    if not args.no_html:
        render_html(out, datasets, references, LABELS, summaries, limitations, scope_text)
        runtime['plotly'] = plotly.__version__
    result = {'schema_version':1,'status':'descriptive_comparison_complete','created_at':datetime.now().astimezone().isoformat(),
              'frontier_status':'Partially evaluated','scope':scope_text,'inputs':provenance,
              'options':{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
              'metrics':summaries,'skipped_subsets':skipped,'limitations':limitations,'runtime':runtime,
              'command':sys.argv,'script_sha256':sha(Path(__file__))}
    dump(out/'evaluation.json',result)
    lines = ['# ICCC major × diagnosis-age CSV comparison','',scope_text+'. One included row = one tumor.','',
             '| Dataset | Included tumors | Joint TV (bands) | Age TV (bands) | ICCC mix TV at single-year age | Conditional coverage |',
             '|---|---:|---:|---:|---:|---:|']
    for label in datasets:
        g = next(m for m in summaries if m['dataset']==label and m['scope']=='all' and m['age_axis']=='grouped_age')
        s = next(m for m in summaries if m['dataset']==label and m['scope']=='all' and m['age_axis']=='single_year')
        lines.append(f"| {label.replace('|', '/')} | {g['generated_tumor_rows']:,} | {g['joint_tv']:.4f} | {g['age_marginal_tv']:.4f} | {format_metric(s['iccc_given_age_weighted_tv'])} | {s['conditional_coverage_fraction']:.4%} |")
    lines += ['', f"Reference rows: {provenance['Reference']['included_rows']:,}.",'', 'Lower TV is closer. Missing conditional support is N/A, not a zero discrepancy.','',
              '## Outputs','', '- `report.html`: four interactive chart sections (unless --no-html).',
              '- `cell_comparison.csv`, `by_age.csv`, `by_iccc_major.csv`, `metrics.csv`: counts and percentages.',
              '- `evaluation.json`: input checks, hashes, options, metrics and limits.','', '## Interpretation and limits','']
    lines += ['- '+text for text in limitations]
    (out/'summary.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    dump(out/'completion.json', {'status':'complete','completed_at':datetime.now().astimezone().isoformat(),
                                'outputs':[{'path':p.name,'sha256':sha(p)} for p in sorted(out.iterdir()) if p.is_file()]})
    print(f'Complete: {out / ("summary.md" if args.no_html else "report.html")}')


def cli():
    """Entry point shared by the installed command and standalone script."""
    try:
        main()
    except (ValueError, OSError, csv.Error) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(cli())
