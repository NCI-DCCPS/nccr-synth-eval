# ICCC Age Evaluation 0.1.0

Compare one reference CSV and two generated CSVs. Each row represents one tumor.
The `iccc-age` command produces four interactive chart sections, comparison CSVs,
a summary, input hashes and a completion record.

## Move to another workspace or computer

Copy and unzip `iccc-age-evaluation-0.1.0-portable.zip` at the destination.
The bundle includes application source, tests, tiny fictional examples, the
application wheel, and dependency wheels. Supply your three input CSVs separately.
The generator workspace, reference snapshots and old reports are not required.

Install Python 3.11 or newer, with `venv` and `pip` available, before starting.
Create the virtual environment at the destination; do not move an existing `.venv`.
The packaged dependencies are pure-Python wheels. The installation and example
workflow are tested on macOS with Python 3.12; Windows and Linux commands are
provided but have not been executed on those systems.

## Install from the portable bundle (no internet needed)

From the unzipped `iccc-age-evaluation-0.1.0` directory on macOS or Linux:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links wheels iccc-age-evaluation==0.1.0
.venv/bin/iccc-age --version
```

On Windows PowerShell:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --no-index --find-links wheels iccc-age-evaluation==0.1.0
.\.venv\Scripts\iccc-age.exe --version
```

The expected version output is `iccc-age 0.1.0`. Python itself is not included in
the bundle. Python's `venv` support must be installed on the destination system.
The dependency versions and file hashes are recorded in `bundle-manifest.json`.

## Run the included example

```sh
.venv/bin/iccc-age \
  examples/reference.csv examples/generated_a.csv examples/generated_b.csv \
  --labels "Generated A" "Generated B" \
  --min-year 1995 --max-year 2022 \
  --out results/example-01
```

Open `results/example-01/report.html` in a browser. The example has eight tumor
rows in each dataset. Generated B is identical to the reference (zero discrepancy);
Generated A has a grouped joint TV of 0.375. These fictional rows test installation
and chart generation; they are not a clinical reference.

Use `.\.venv\Scripts\iccc-age.exe` in Windows commands. The arguments are the same;
write the command on one line in PowerShell or use PowerShell continuation syntax.

## Run your three CSV files

```sh
.venv/bin/iccc-age \
  /path/to/reference.csv /path/to/generated_a.csv /path/to/generated_b.csv \
  --labels "r3 release" "r4 candidate" \
  --min-year 1995 --max-year 2022 \
  --out results/comparison-01
```

Paths with spaces must be quoted. Input paths and `--out` are resolved from your
current working directory. The installed command can run from any directory;
use the absolute path to `.venv/bin/iccc-age` if you are outside this folder.
Choose a new output directory each time: existing destinations are refused.

Required columns are `icccMajor` and `ageRecodeSingle90p`. `yearOfDiagnosis` is
required only when year limits are specified. Optional `ageRecode20Groups` values
are checked against single-year age. Optional `indexCancer` in all three files
enables supplemental index/nonindex comparisons.

- ICCC major codes: 01–12 or 99 (unclassified/in situ).
- Diagnosis-age codes: 00–90 (90 means 90+) or 99 (unknown).
- Blank ICCC and age values are retained as distinct categories.
- One-digit codes are padded and surrounding whitespace is removed, with counts
  recorded in `evaluation.json`.
- Omit both year limits to compare all rows. The 1995–2022 limits reproduce the
  calendar scope used in the original Data Jamboree evaluation.

For different headers, add `--iccc-column ICCCMajor --age-column diagnosis_age`
and, if needed, `--year-column diagnosis_year`. Mapped columns must use the same
code meanings. Run `iccc-age --help` for all options. `python -m compare_iccc_age_csv`
is also supported when using the environment's Python interpreter.

## Outputs

- `report.html`: ICCC mix within age groups, differences from the reference,
  overall age distribution, and age distribution within each ICCC major group.
- `cell_comparison.csv`, `by_age.csv`, `by_iccc_major.csv`, `metrics.csv`: counts,
  percentages, coverage and total-variation metrics.
- `summary.md`, `evaluation.json`, `completion.json`: summary and provenance.
- `compare_iccc_age_csv.py`: exact source code used for the run.

Charts embed their JavaScript and work offline. No patient-level rows or IDs are
written to the reports; local input file paths are recorded for reproducibility.
This is descriptive reference agreement, not clinical validation or release
certification. Source representativeness and calibration overlap are not inferred.

## Alternative installation with internet access

If you copy only the application wheel, pip can obtain the dependencies:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install ./iccc_age_evaluation-0.1.0-py3-none-any.whl
```

For development, copy the source tree (including `pyproject.toml`, `INSTALL.md`,
`scripts/`, `tests/` and `examples/`) and run `.venv/bin/python -m pip install .`.
Source builds may download the build backend as well as runtime dependencies.
Use the wheel installation above when an offline setup is required.

For the original single-script workflow, copy only
`scripts/compare_iccc_age_csv.py` and `requirements-iccc-age.txt`, install the
requirements in a fresh environment, and run the script with that environment's
Python. The input/output arguments are unchanged. Numerical-only execution with
`--no-html` works without Plotly.

## Verify or rebuild

```sh
.venv/bin/python -m pip check
.venv/bin/python -B -m unittest discover -s tests -p 'test_*.py' -v
```

To build new versioned distributions from source in a development environment:

```sh
.venv/bin/python -m pip install build
.venv/bin/python -m build
```

The 0.1.0 bundle has been tested by installing its wheel and dependencies into a
fresh environment with `--no-index`, running from an unrelated directory with
spaces in its path, checking error exits and UTF-8 output, and comparing chart
data/metrics with the pre-packaging implementation. This verifies the local
installation workflow; other operating systems remain unverified.

The venv and local-wheel installation pattern follows the
[Python Packaging User Guide](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/).
