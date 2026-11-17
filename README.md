# Data Jamboree Evaluation

Workspace for evaluating and comparing outputs from the sibling `gen-*`
experiments.

## Status

ICCC major group × diagnosis-age comparison is implemented. The first evaluation compares
the r3 release and r4 candidate with the frozen G1 reference for diagnosis years
1995–2022. Results are descriptive; clinical validity and release eligibility are
not assessed.

## Layout

- `README.md`: purpose and workspace conventions.
- `scripts/evaluate_iccc_age.py`: pinned snapshot/manifest evaluation.
- `scripts/compare_iccc_age_csv.py`: standalone three-CSV evaluation and the same four chart sections.
- `pyproject.toml`, `INSTALL.md`: installable package and portability guide.
- `examples/`: fictional CSVs for an installation smoke test.
- `dist/`: versioned package downloads and offline bundle.
- `work/`: local evaluation runs, reports and temporary artifacts; ignored by Git.

## Conventions

- Use the shared definitions under `../../reference/` and identify the exact
  rules, contracts and reference snapshots used by each evaluation.
- Reference generator outputs or release manifests by path and hash. Preserve
  existing datasets, manifests and historical audit evidence.
- Keep each evaluation in a new run directory under `work/`, recording input
  hashes, commands, tool versions, results, limitations and checks not performed.
- Structural and statistical checks do not establish clinical validity.
- Release promotion remains a separate, explicit workflow.

See [workspace structure](../../docs/workspace-structure.md) for the broader
project layout. Evaluation is an independent local Git repository. The parent
workspace continues to ignore this folder. Source, tests, examples and documentation
are version-controlled here; generated reports and distributions remain ignored.

## Install in another workspace or on another computer

The three-CSV evaluator is now packaged as `iccc-age-evaluation` version 0.1.0,
with the installed command `iccc-age`.

- [Portable ZIP including runtime dependencies](dist/iccc-age-evaluation-0.1.0-portable.zip)
- [Installation and migration instructions](INSTALL.md)
- [Application wheel](dist/iccc_age_evaluation-0.1.0-py3-none-any.whl)
- [Source distribution](dist/iccc_age_evaluation-0.1.0.tar.gz)
- [Download checksums](dist/SHA256SUMS)

Copy the ZIP to the destination and extract it. From its extracted directory:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links wheels iccc-age-evaluation==0.1.0
.venv/bin/iccc-age --version
```

Python 3.11+ is required; recreate the environment at the destination rather than
copying `.venv`. The bundle installs without internet after Python is available.
Supply the three input CSVs separately. The generated HTML reports also work offline.
See `INSTALL.md` for Windows commands, sample data, input columns, and developer setup.

A fresh macOS/Python 3.12 environment passed offline installation, dependency checks,
eight tests, a run from an unrelated directory with spaces in its path, UTF-8 output,
and chart/metric comparison with the previous implementation. Windows and Linux
have not been tested. The pinned snapshot evaluator remains a separate workspace tool.

## Compare three CSV files directly

Use `scripts/compare_iccc_age_csv.py` with **one reference CSV and two generated
CSVs**, each containing one tumor per row. The script is standalone: it does not
require this workspace, release manifests, frozen snapshots, or patient IDs.
You can copy the single Python file to another folder or computer.

Python 3.11+ is required. Install Plotly once, preferably in a virtual environment:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-iccc-age.txt
```

From `experiments/evaluation/`, run:

```sh
.venv/bin/python scripts/compare_iccc_age_csv.py \
  /path/to/reference.csv \
  /path/to/generated_r3.csv \
  /path/to/generated_r4.csv \
  --labels "r3 release" "r4 candidate" \
  --min-year 1995 --max-year 2022 \
  --out work/my-iccc-age-comparison
```

Open `work/my-iccc-age-comparison/report.html` in a browser. It includes:

1. ICCC mix within each age group: reference and both generated datasets.
2. Generated-minus-reference differences in that mix.
3. The overall diagnosis-age distribution.
4. The age distribution within each ICCC major group.

The report embeds Plotly and works offline. Hover for counts/percentages, toggle
series in legends, and export individual charts as SVG. CSV tables, a Markdown
summary, input SHA-256 hashes, filtering counts, and a completion record accompany
the report. Existing output directories are refused; choose a new path per run.

### CSV columns and defaults

| Column | Requirement | Meaning |
|---|---|---|
| `icccMajor` | Required | Codes 01–12 or 99; 99 is unclassified/in situ. |
| `ageRecodeSingle90p` | Required | Diagnosis age in complete years, 00–90; 90 means 90+, 99 unknown. |
| `yearOfDiagnosis` | Required only when year limits are supplied | Four-digit diagnosis year. |
| `ageRecode20Groups` | Optional | Cross-checked against bands derived from single-year age; inconsistencies fail. |
| `indexCancer` | Optional | If present in all inputs, adds index (1/2) and nonindex comparisons to the CSV tables. |

One-digit numeric codes are zero-padded and surrounding whitespace is removed;
normalization counts are reported. Blank ICCC/age values remain separate
categories. Other nonnumeric values, invalid codes, duplicate headers, and malformed
rows fail clearly. These are **row-level CSV inputs**, not pre-aggregated count tables.

Year filtering applies equally to all three files and reports excluded or invalid
years. **Omitting `--min-year` and `--max-year` compares all input rows.** The
1995–2022 flags above reproduce the prior evaluation's calendar scope.

If column names differ, append options such as:

```sh
--iccc-column ICCCMajor --age-column diagnosis_age --year-column diagnosis_year
```

The mapped age column must still use the diagnosis-age encoding above. Age bands
alone cannot produce the single-year comparison. Use `--help` for all options;
`--no-html` produces numerical results without needing Plotly.

Validation includes seven passing boundary/CLI tests and a full replay against
the original reference plus generated CSVs. All 12 metric rows, all four CSV
tables, and all four chart trace/layout configurations match the retained
snapshot-based results. See [reproduction verification](work/iccc-age-csv-verification-2026-11-16T16-09-00-0500/reproduction-verification.json).

Run the tests with `python3 -B -m unittest discover -s tests -p "test_*.py" -v`.
Clinical validity, ICCC hierarchy and release eligibility are outside this
standalone comparison's scope.

## Latest evaluation

- [Interactive report](work/iccc-major-age-2026-11-16T15-50-16-0500/report.html)
- [Summary](work/iccc-major-age-2026-11-16T15-50-16-0500/summary.md)
- [Numerical verification](work/iccc-major-age-2026-11-16T15-50-16-0500/verification.json)
- [Input provenance and metrics](work/iccc-major-age-2026-11-16T15-50-16-0500/evaluation.json)

Each generated dataset contains 161,469 tumor rows. The comparison reference
contains 1,542,308 tumor rows. The r4 candidate has substantially smaller
age-conditional ICCC discrepancies than r3, while retaining the same age
distribution. Exact-year conditional reference coverage is 99.9709%; 47 generated
rows have ages absent from the reference.

## Reproduce

Run from this directory with Python 3.11+ and Plotly installed. Use `--no-html`
for a standard-library-only numerical report. Choose a new output directory for
every run; existing destinations are refused.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 scripts/evaluate_iccc_age.py \
  --workspace ../.. \
  --snapshot ../../reference/nccr-derived/2026-11-13-stage-test-v30-g1-v1 \
  --snapshot-sha256 64e432745465b40c620b3816dbefb256b73f3c30bcf8f05840ab31a77fb68dda \
  --dataset 'r3 release=releases/jamboree_sdv_rules_v2_3_0_scale_10pct_r3.json' \
  --dataset 'r4 candidate=experiments/gen-sdv/work/runs/jamboree_sdv_rules_v4_0_0_scale_10pct_r4-2026-11-13T23-13-43.154090-0500/candidate.json' \
  --out work/CHOOSE-A-NEW-RUN-ID
```

The retained `verify_evaluation.py` uses pandas and NumPy to independently
recount all three tumor scopes at both age resolutions. It also cross-checks
grouped reference counts against the separately frozen year/index count table.
