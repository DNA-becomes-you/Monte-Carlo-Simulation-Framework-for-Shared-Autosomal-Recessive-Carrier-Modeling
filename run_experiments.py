#!/usr/bin/env python3
"""
CLI entrypoint for Monte Carlo simulation experiments.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from simulation.experiments import run_grid


def _parse_list(value: str, conv: type) -> list:
    """Parse comma-separated string to list. Single value also works."""
    s = str(value).strip()
    if "," in s:
        return [conv(x.strip()) for x in s.split(",")]
    return [conv(s)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Monte Carlo carrier overlap experiments"
    )
    parser.add_argument(
        "--users",
        type=int,
        default=200,
        help="Number of individuals per experiment (default: 200)",
    )
    parser.add_argument(
        "--genes",
        type=str,
        default="500",
        help="Comma-separated gene panel sizes (default: 500)",
    )
    parser.add_argument(
        "--avg-carriers",
        type=str,
        default="3",
        help="Comma-separated avg-carriers values (default: 3)",
    )
    parser.add_argument(
        "--population-structures",
        type=str,
        default="0.01,0.02,0.03,0.05,0.1",
        help="Comma-separated population-structure values",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("experiment_results.json"),
        help="Output JSON path",
    )
    return parser.parse_args()


def print_table(experiments: list[dict]) -> None:
    """Print formatted results table."""
    header = (
        "Genes | AvgCarriers | Structure | Eligible% (CI)           | "
        "%Overlap1 (CI)        | AvgInt (CI)"
    )
    print(header)
    print("-" * len(header))
    for exp in experiments:
        g = exp["genes"]
        ac = exp["avgCarriers"]
        ps = exp["populationStructure"]
        e = exp["pctEligibleZero"]
        e_lo, e_hi = exp["pctEligibleZeroCI"]
        o = exp["pctOverlapOne"]
        o_lo, o_hi = exp["pctOverlapOneCI"]
        a = exp["avgIntersection"]
        a_lo, a_hi = exp["avgIntersectionCI"]
        row = (
            f"{g:<5} | {ac:<10.1f} | {ps:<9.2g} | "
            f"{e:.1f}% ({e_lo:.1f}–{e_hi:.1f})      | "
            f"{o:.1f}% ({o_lo:.1f}–{o_hi:.1f})   | "
            f"{a:.4f} ({a_lo:.4f}–{a_hi:.4f})"
        )
        print(row)


def main() -> None:
    args = parse_args()
    genes_list = _parse_list(args.genes, int)
    avg_carriers_list = _parse_list(args.avg_carriers, float)
    structures = _parse_list(args.population_structures, float)

    experiments = run_grid(
        genes_list=genes_list,
        avg_carriers_list=avg_carriers_list,
        structures_list=structures,
        users=args.users,
        seed=args.seed,
    )

    output_experiments = []
    for exp in experiments:
        output_experiments.append({
            "genes": exp["genes"],
            "avgCarriers": exp["avgCarriers"],
            "populationStructure": exp["populationStructure"],
            "pctEligibleZero": exp["pctEligibleZero"],
            "pctEligibleZeroCI": list(exp["pctEligibleZeroCI"]),
            "pctOverlapOne": exp["pctOverlapOne"],
            "pctOverlapOneCI": list(exp["pctOverlapOneCI"]),
            "avgIntersection": exp["avgIntersection"],
            "avgIntersectionCI": list(exp["avgIntersectionCI"]),
        })

    output = {
        "metadata": {
            "users": args.users,
            "genesSweep": genes_list,
            "avgCarriersSweep": avg_carriers_list,
            "populationStructuresSweep": structures,
            "seed": args.seed,
        },
        "experiments": output_experiments,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    print(f"Results saved to {args.output}\n")
    print_table(experiments)


if __name__ == "__main__":
    main()
