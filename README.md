# Transcriptomics Differential-Expression Pipeline

A compact, auditable Python workflow for comparing two groups of bulk RNA-seq
count samples. It demonstrates scientific reasoning, statistical hygiene,
reproducibility, and production-minded input validation.

> **Data status:** the built-in demo dataset is synthetic and must not be
> interpreted as a biological result. Every demo output includes metadata that
> identifies it as synthetic.

## What it does

- validates integer-like, non-negative count matrices and sample metadata;
- estimates library size factors with median-ratio normalisation;
- computes log2 fold changes from normalised counts;
- performs two-sided Welch tests on `log2(normalised count + 1)`;
- controls false discoveries with the Benjamini-Hochberg procedure;
- exports a ranked result table and PCA sample coordinates;
- optionally renders a volcano plot and a PCA quality-control plot.

This is an educational portfolio implementation, not a substitute for mature
RNA-seq tools such as DESeq2, edgeR, or limma-voom. A production analysis should
model count dispersion, account for study covariates, and pre-register its
filtering and contrast choices.

## Quick start

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -e ".[dev]"
txde demo --output-dir outputs/demo
pytest
```

## Interactive web application

```bash
streamlit run app.py
```

The app opens with a visibly labelled synthetic demonstration. It also accepts
uploaded count and sample-metadata CSV files, exposes the comparison and
filtering parameters, provides interactive volcano and PCA plots, and exports
results plus machine-readable run metadata.

To run the same application in a container:

```bash
docker build -t txde-app .
docker run --rm -p 8501:8501 txde-app
```

Expected outputs:

- `demo_counts.csv` and `demo_samples.csv` — explicitly labelled synthetic inputs;
- `differential_expression.csv` — ranked statistics and FDR values;
- `pca_coordinates.csv` — two-dimensional sample representation;
- `volcano.png` and `pca.png` — diagnostic figures;
- `run_metadata.json` — parameters, timestamp, software version, and data status.

## Analyse your own data

The count CSV must have gene IDs in the first column and one sample per remaining
column. The metadata CSV requires `sample_id` and `group` columns.

```bash
txde run \
  --counts path/to/counts.csv \
  --samples path/to/samples.csv \
  --group-a control \
  --group-b treated \
  --output-dir outputs/my_analysis
```

## Reproducibility choices

- The demo random seed is fixed by default and recorded.
- Results are sorted deterministically.
- Invalid counts, duplicated identifiers, missing metadata, or undersized groups
  fail early with useful messages.
- Analysis parameters and input paths are written to machine-readable metadata.

## Test coverage

The tests check FDR monotonicity and bounds, size-factor behaviour, recovery of
known synthetic signal, and invalid-input handling.

