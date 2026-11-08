"""Interactive Streamlit application for the transcriptomics workflow."""

from __future__ import annotations

import io
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from txde import __version__, analyse
from txde.demo import generate_demo

st.set_page_config(
    page_title="Transcriptomics DE Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    :root { --accent: #0f766e; --ink: #0f172a; }
    .stApp { color: var(--ink); }
    .block-container { max-width: 1280px; padding-top: 2rem; }
    .demo-banner { border-left: 5px solid #b45309; background: #fffbeb;
      padding: .8rem 1rem; border-radius: .25rem; margin: .8rem 0 1.2rem; }
    .method-note { border-left: 4px solid var(--accent); padding-left: 1rem; }
    [data-testid="stMetricValue"] { font-variant-numeric: tabular-nums; }
    :focus-visible { outline: 3px solid #0f766e !important; outline-offset: 3px; }
    @media (prefers-reduced-motion: reduce) { * { scroll-behavior: auto !important;
      animation-duration: .01ms !important; transition-duration: .01ms !important; } }
    </style>
    """,
    unsafe_allow_html=True,
)


@st.cache_data(show_spinner=False)
def demo_data(seed: int, genes: int, samples_per_group: int):
    return generate_demo(seed, genes, samples_per_group)


@st.cache_data(show_spinner=False)
def run_cached(
    counts_csv: bytes,
    samples_csv: bytes,
    group_a: str,
    group_b: str,
    min_mean_count: float,
):
    counts = pd.read_csv(io.BytesIO(counts_csv), index_col=0)
    samples = pd.read_csv(io.BytesIO(samples_csv))
    return analyse(counts, samples, group_a, group_b, min_mean_count)


def to_csv_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False).encode("utf-8")


def volcano_figure(results: pd.DataFrame, alpha: float) -> go.Figure:
    frame = results.copy()
    frame["status"] = np.where(frame["fdr"] <= alpha, "Passes FDR threshold", "Not significant")
    frame["minus_log10_p"] = -np.log10(frame["p_value"].clip(lower=1e-300))
    figure = px.scatter(
        frame,
        x="log2_fold_change",
        y="minus_log10_p",
        color="status",
        symbol="status",
        hover_name="gene_id",
        hover_data={"p_value": ":.3g", "fdr": ":.3g", "minus_log10_p": False},
        color_discrete_map={"Passes FDR threshold": "#b91c1c", "Not significant": "#64748b"},
        labels={"minus_log10_p": "−log10(p-value)", "log2_fold_change": "log2 fold change"},
    )
    figure.add_vline(x=0, line_dash="dot", line_color="#334155")
    figure.update_layout(legend_title_text="Statistical status", height=520)
    return figure


def pca_figure(pca: pd.DataFrame) -> go.Figure:
    return px.scatter(
        pca,
        x="PC1",
        y="PC2",
        color="group",
        symbol="group",
        text="sample_id",
        hover_name="sample_id",
        height=520,
        labels={"group": "Experimental group"},
    ).update_traces(textposition="top center")


st.title("Transcriptomics Differential-Expression Explorer")
st.caption("A transparent two-group RNA-seq workflow for reproducible scientific analysis")

with st.sidebar:
    st.header("Analysis setup")
    source = st.radio("Data source", ["Synthetic demonstration", "Upload my data"])
    if source == "Synthetic demonstration":
        seed = st.number_input("Random seed", min_value=0, value=2026, step=1)
        genes = st.slider("Number of genes", 100, 3000, 600, 100)
        samples_per_group = st.slider("Samples per group", 2, 20, 6)
        counts, samples, known_truth = demo_data(int(seed), genes, samples_per_group)
        data_status = "SYNTHETIC_DEMO"
    else:
        st.markdown("Count CSV: genes in rows, samples in columns. Metadata requires `sample_id` and `group`.")
        count_upload = st.file_uploader("Count matrix CSV", type="csv")
        sample_upload = st.file_uploader("Sample metadata CSV", type="csv")
        if not count_upload or not sample_upload:
            st.info("Upload both files to configure the comparison.")
            st.stop()
        try:
            counts = pd.read_csv(count_upload, index_col=0)
            samples = pd.read_csv(sample_upload)
        except (OSError, UnicodeError, ValueError, pd.errors.ParserError) as exc:
            st.error(f"The uploaded files could not be read: {exc}")
            st.stop()
        known_truth = set()
        data_status = "USER_SUPPLIED_UNVERIFIED"

    available_groups = sorted(samples["group"].dropna().astype(str).unique()) if "group" in samples else []
    if len(available_groups) < 2:
        st.error("At least two groups are required in the sample metadata.")
        st.stop()
    group_a = st.selectbox("Reference group", available_groups, index=0)
    comparison_options = [group for group in available_groups if group != group_a]
    group_b = st.selectbox("Comparison group", comparison_options, index=0)
    min_mean_count = st.number_input("Minimum mean normalized count", min_value=0.0, value=5.0, step=1.0)
    alpha = st.slider("FDR threshold", 0.001, 0.20, 0.05, 0.001, format="%.3f")

if data_status == "SYNTHETIC_DEMO":
    st.markdown(
        '<div class="demo-banner"><strong>Synthetic demonstration data.</strong> '
        "These counts are simulated and are not biological evidence.</div>",
        unsafe_allow_html=True,
    )
else:
    st.warning("Uploaded data are user supplied and have not been independently verified.")

counts_bytes = counts.to_csv().encode("utf-8")
samples_bytes = samples.to_csv(index=False).encode("utf-8")
try:
    with st.spinner("Normalizing counts and calculating statistics…"):
        results, pca, factors = run_cached(
            counts_bytes, samples_bytes, group_a, group_b, float(min_mean_count)
        )
except (KeyError, TypeError, ValueError, RuntimeError) as exc:
    st.error(f"Analysis stopped: {exc}")
    st.stop()

significant = results["fdr"] <= alpha
metric_columns = st.columns(4)
metric_columns[0].metric("Genes tested", f"{len(results):,}")
metric_columns[1].metric("Passing FDR threshold", f"{int(significant.sum()):,}")
metric_columns[2].metric("Samples", f"{len(pca):,}")
metric_columns[3].metric("Median |log2 FC|", f"{results['log2_fold_change'].abs().median():.2f}")

overview, pca_tab, results_tab, downloads_tab, methods_tab = st.tabs(
    ["Differential expression", "Sample PCA", "Results table", "Downloads", "Methods"]
)

with overview:
    st.plotly_chart(volcano_figure(results, alpha), width="stretch")
    if data_status == "SYNTHETIC_DEMO":
        detected = set(results.loc[significant, "gene_id"])
        true_positive = len(detected & known_truth)
        st.caption(
            f"Demo-only benchmark: {true_positive} of {len(known_truth)} simulated effect genes "
            f"were detected at FDR ≤ {alpha:.3f}."
        )

with pca_tab:
    st.plotly_chart(pca_figure(pca), width="stretch")
    st.dataframe(
        factors.rename_axis("sample_id").rename("size_factor").reset_index(),
        hide_index=True,
        width="stretch",
    )

with results_tab:
    significant_only = st.toggle("Show only genes passing the FDR threshold")
    shown = results.loc[significant] if significant_only else results
    st.dataframe(
        shown,
        hide_index=True,
        width="stretch",
        column_config={
            "p_value": st.column_config.NumberColumn(format="%.3e"),
            "fdr": st.column_config.NumberColumn(format="%.3e"),
            "log2_fold_change": st.column_config.NumberColumn(format="%.3f"),
        },
    )

with downloads_tab:
    metadata = {
        "software": "txde-portfolio-web",
        "version": __version__,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "data_status": data_status,
        "contrast": {"reference": group_a, "comparison": group_b},
        "min_mean_count": min_mean_count,
        "fdr_display_threshold": alpha,
    }
    st.download_button("Download differential-expression results", to_csv_bytes(results), "differential_expression.csv", "text/csv")
    st.download_button("Download PCA coordinates", to_csv_bytes(pca), "pca_coordinates.csv", "text/csv")
    st.download_button("Download run metadata", json.dumps(metadata, indent=2), "run_metadata.json", "application/json")
    if data_status == "SYNTHETIC_DEMO":
        st.download_button("Download synthetic count matrix", counts_bytes, "synthetic_demo_counts.csv", "text/csv")
        st.download_button("Download synthetic sample sheet", samples_bytes, "synthetic_demo_samples.csv", "text/csv")

with methods_tab:
    st.markdown(
        """
        <div class="method-note">
        <strong>Workflow.</strong> Median-ratio library normalization is followed by
        two-sided Welch tests on log2(normalized count + 1). P-values are adjusted
        with the Benjamini–Hochberg procedure. PCA is computed by singular-value
        decomposition of the transformed normalized matrix.
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.warning(
        "This portfolio implementation is educational. Production RNA-seq studies should use a "
        "validated count model such as DESeq2 or edgeR and account for dispersion, covariates, "
        "batch effects, study design, and pre-specified filtering."
    )

