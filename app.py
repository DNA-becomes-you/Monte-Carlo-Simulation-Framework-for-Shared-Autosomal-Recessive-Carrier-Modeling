#!/usr/bin/env python3
"""
Streamlit UI for Monte Carlo simulation framework.
"""

from __future__ import annotations

import json
import random
from io import BytesIO

import streamlit as st

from simulation.population import generate_population
from simulation.statistics import evaluate_population

st.set_page_config(
    page_title="Monte Carlo Carrier Simulation",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("Monte Carlo Simulation Framework for Shared Autosomal Recessive Carrier Modeling")

st.markdown(
    """
    This tool simulates pairwise overlap of pathogenic recessive carrier variants
    in synthetic populations. It is intended for research and methodology validation.
    """
)

with st.sidebar:
    st.header("Parameters")
    users = st.number_input(
        "Number of individuals",
        min_value=50,
        max_value=5000,
        value=500,
        step=50,
    )
    genes = st.number_input(
        "Gene panel size",
        min_value=100,
        max_value=2000,
        value=500,
        step=100,
    )
    avg_carriers = st.slider(
        "Average carriers per individual (λ)",
        min_value=0.5,
        max_value=6.0,
        value=3.0,
        step=0.5,
    )
    population_structure = st.slider(
        "Population structure (founder fraction)",
        min_value=0.0,
        max_value=0.1,
        value=0.02,
        step=0.01,
        format="%.2f",
    )
    seed = st.number_input(
        "Random seed",
        min_value=0,
        value=42,
        step=1,
    )
    n_bootstrap = st.number_input(
        "Bootstrap samples",
        min_value=100,
        max_value=5000,
        value=1000,
        step=100,
    )

    run = st.button("Run Simulation", type="primary")

if run:
    rng = random.Random(seed)
    with st.spinner("Generating population..."):
        population = generate_population(
            users=users,
            panel_size=genes,
            avg_carriers=avg_carriers,
            population_structure=population_structure,
            rng=rng,
        )
    with st.spinner("Computing statistics and bootstrap CIs..."):
        stats = evaluate_population(
            population,
            max_pairs=min(1000, users * (users - 1) // 2),
            n_bootstrap=n_bootstrap,
            rng=rng,
        )

    st.subheader("Results")

    col1, col2, col3 = st.columns(3)
    with col1:
        e = stats["pctEligibleZero"]
        e_lo, e_hi = stats["pctEligibleZeroCI"]
        st.metric(
            "Compatibility % (no shared loci)",
            f"{e:.1f}%",
            f"95% CI: {e_lo:.1f}–{e_hi:.1f}",
        )
    with col2:
        o = stats["pctOverlapOne"]
        o_lo, o_hi = stats["pctOverlapOneCI"]
        st.metric(
            "% Overlap exactly 1",
            f"{o:.1f}%",
            f"95% CI: {o_lo:.1f}–{o_hi:.1f}",
        )
    with col3:
        a = stats["avgIntersection"]
        a_lo, a_hi = stats["avgIntersectionCI"]
        st.metric(
            "Avg intersection size",
            f"{a:.4f}",
            f"95% CI: {a_lo:.4f}–{a_hi:.4f}",
        )

    st.subheader("Intersection distribution")
    dist = stats.get("intersectionDistribution", [])
    if dist:
        import plotly.express as px

        isects = [d["intersection"] for d in dist]
        counts = [d["count"] for d in dist]
        fig = px.bar(
            x=isects,
            y=counts,
            labels={"x": "Intersection size", "y": "Pair count"},
            title="Histogram of pairwise intersection sizes",
        )
        fig.update_layout(showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

    st.subheader("Download results")
    output = {
        "metadata": {
            "users": users,
            "genes": genes,
            "avgCarriers": avg_carriers,
            "populationStructure": population_structure,
            "seed": seed,
            "nBootstrap": n_bootstrap,
        },
        "statistics": {
            "pctEligibleZero": stats["pctEligibleZero"],
            "pctEligibleZeroCI": list(stats["pctEligibleZeroCI"]),
            "pctOverlapOne": stats["pctOverlapOne"],
            "pctOverlapOneCI": list(stats["pctOverlapOneCI"]),
            "avgIntersection": stats["avgIntersection"],
            "avgIntersectionCI": list(stats["avgIntersectionCI"]),
            "pairCount": stats["pairCount"],
        },
    }
    buf = BytesIO()
    buf.write(json.dumps(output, indent=2).encode())
    buf.seek(0)
    st.download_button(
        "Download JSON",
        data=buf,
        file_name="simulation_results.json",
        mime="application/json",
    )

st.sidebar.markdown("---")
st.sidebar.caption(
    "This tool models shared autosomal recessive pathogenic carrier overlap. "
    "It is not a trait optimization or polygenic scoring tool. For research only."
)
