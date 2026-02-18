#!/usr/bin/env python3
"""
DNABY DATE — Recessive Carrier Experiment Runner (Zone 2)

Runs parameter sweeps (genes × avg_carriers × population-structure).
Includes bootstrap 95% confidence intervals.
Models carrier overlap in known recessive genes.
No database usage. Standard library only. Deterministic with --seed.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

# Reuse population logic from generator (in-memory)
from generate_population import (
    generate_locus_weights,
    generate_user_vector,
)


def _parse_list(value: str, conv: type) -> list:
    """Parse comma-separated string to list. Single value also works."""
    s = str(value).strip()
    if "," in s:
        return [conv(x.strip()) for x in s.split(",")]
    return [conv(s)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run carrier compatibility experiments"
    )
    parser.add_argument(
        "--users",
        type=int,
        default=200,
        help="Number of users per experiment (default: 200)",
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
        "--variance",
        type=float,
        default=1.5,
        help="Variance around avg-carriers (default: 1.5)",
    )
    parser.add_argument(
        "--population-structures",
        type=str,
        default="0.01,0.02,0.03,0.05,0.1",
        help="Comma-separated population-structure values (default: 0.01,0.02,0.03,0.05,0.1)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for determinism (default: 42)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("fixtures/synthetic/experiment_results.json"),
        help="Output JSON path",
    )
    return parser.parse_args()


def compute_intersection_size(a: list[int], b: list[int]) -> int:
    """
    Two-pointer intersection size. O(n + m).
    Requires a and b to be sorted.
    """
    i, j = 0, 0
    count = 0
    while i < len(a) and j < len(b):
        if a[i] == b[j]:
            count += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    return count


def bootstrap_ci(
    values: list[float],
    rng: random.Random,
    n_bootstrap: int = 1000,
) -> tuple[float, float]:
    """
    Bootstrap 95% CI for the mean.
    Returns (lower_95, upper_95).
    """
    if not values:
        return (0.0, 0.0)
    n = len(values)
    means: list[float] = []
    for _ in range(n_bootstrap):
        sample = rng.choices(values, k=n)
        means.append(sum(sample) / n)
    means.sort()
    lo = int(0.025 * n_bootstrap)
    hi = int(0.975 * n_bootstrap)
    return (means[lo], means[hi])


def generate_population(
    users: int,
    panel_size: int,
    avg_carriers: float,
    variance: float,
    population_structure: float,
    rng: random.Random,
) -> list[dict]:
    """Generate synthetic carrier population in memory."""
    locus_weights = generate_locus_weights(panel_size, population_structure, rng)
    population: list[dict] = []
    for i in range(users):
        user = generate_user_vector(
            user_index=i,
            panel_size=panel_size,
            avg_carriers=avg_carriers,
            variance=variance,
            locus_weights=locus_weights,
            rng=rng,
        )
        population.append(user)
    return population


def compute_pair_metrics(
    population: list[dict],
    max_pairs: int = 50000,
    rng: random.Random | None = None,
) -> tuple[list[dict], dict]:
    """
    Compute per-pair metrics and aggregate stats.
    Returns (pair_results, summary).
    pair_results: list of {intersection_size, eligible, overlap_one}
    """
    n = len(population)
    total_pairs = n * (n - 1) // 2

    if total_pairs == 0:
        return ([], {
            "pctEligibleZero": 0.0,
            "pctOverlapOne": 0.0,
            "avgIntersection": 0.0,
        })

    all_pairs: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            all_pairs.append((i, j))

    if len(all_pairs) <= max_pairs:
        pairs = all_pairs
    else:
        r = rng if rng is not None else random.Random()
        pairs = r.sample(all_pairs, max_pairs)

    pair_results: list[dict] = []
    for i, j in pairs:
        a = population[i]["pathogenicLoci"]
        b = population[j]["pathogenicLoci"]
        sz = compute_intersection_size(a, b)
        pair_results.append({
            "intersection_size": sz,
            "eligible": sz == 0,
            "overlap_one": sz == 1,
        })

    pair_count = len(pair_results)
    eligible_count = sum(1 for p in pair_results if p["eligible"])
    overlap_one_count = sum(1 for p in pair_results if p["overlap_one"])
    total_intersection = sum(p["intersection_size"] for p in pair_results)

    summary = {
        "pctEligibleZero": 100 * eligible_count / pair_count,
        "pctOverlapOne": 100 * overlap_one_count / pair_count,
        "avgIntersection": total_intersection / pair_count,
    }
    return (pair_results, summary)


def run_experiment(
    genes: int,
    avg_carriers: float,
    population_structure: float,
    users: int,
    variance: float,
    seed: int,
    scenario_index: int,
) -> dict:
    """Run one experiment for a given parameter combination."""
    scenario_seed = seed + scenario_index
    rng = random.Random(scenario_seed)
    population = generate_population(
        users, genes, avg_carriers, variance, population_structure, rng
    )
    pair_results, summary = compute_pair_metrics(
        population, max_pairs=500, rng=rng
    )

    # Bootstrap CIs (use same RNG for determinism)
    eligible_vals = [1.0 if p["eligible"] else 0.0 for p in pair_results]
    overlap_one_vals = [1.0 if p["overlap_one"] else 0.0 for p in pair_results]
    intersection_vals = [float(p["intersection_size"]) for p in pair_results]

    ci_eligible = bootstrap_ci(eligible_vals, rng)
    ci_overlap_one = bootstrap_ci(overlap_one_vals, rng)
    ci_avg_int = bootstrap_ci(intersection_vals, rng)

    # Convert eligible/overlap CIs to percentages
    ci_eligible_pct = (100 * ci_eligible[0], 100 * ci_eligible[1])
    ci_overlap_one_pct = (100 * ci_overlap_one[0], 100 * ci_overlap_one[1])

    return {
        "genes": genes,
        "avgCarriers": avg_carriers,
        "populationStructure": population_structure,
        "pctEligibleZero": round(summary["pctEligibleZero"], 2),
        "pctEligibleZeroCI": (round(ci_eligible_pct[0], 2), round(ci_eligible_pct[1], 2)),
        "pctOverlapOne": round(summary["pctOverlapOne"], 2),
        "pctOverlapOneCI": (round(ci_overlap_one_pct[0], 2), round(ci_overlap_one_pct[1], 2)),
        "avgIntersection": round(summary["avgIntersection"], 4),
        "avgIntersectionCI": (round(ci_avg_int[0], 4), round(ci_avg_int[1], 4)),
    }


def print_table(experiments: list[dict]) -> None:
    """Print comparison table with bootstrap CIs."""
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

    experiments: list[dict] = []
    scenario_index = 0

    for genes in genes_list:
        for avg_carriers in avg_carriers_list:
            for ps in structures:
                exp = run_experiment(
                    genes=genes,
                    avg_carriers=avg_carriers,
                    population_structure=ps,
                    users=args.users,
                    variance=args.variance,
                    seed=args.seed,
                    scenario_index=scenario_index,
                )
                experiments.append(exp)
                scenario_index += 1

    # JSON-serializable output (tuples -> lists for CI)
    output_experiments = []
    for exp in experiments:
        out = {
            "genes": exp["genes"],
            "avgCarriers": exp["avgCarriers"],
            "populationStructure": exp["populationStructure"],
            "pctEligibleZero": exp["pctEligibleZero"],
            "pctEligibleZeroCI": list(exp["pctEligibleZeroCI"]),
            "pctOverlapOne": exp["pctOverlapOne"],
            "pctOverlapOneCI": list(exp["pctOverlapOneCI"]),
            "avgIntersection": exp["avgIntersection"],
            "avgIntersectionCI": list(exp["avgIntersectionCI"]),
        }
        output_experiments.append(out)

    output = {
        "metadata": {
            "users": args.users,
            "genesSweep": genes_list,
            "avgCarriersSweep": avg_carriers_list,
            "populationStructuresSweep": structures,
            "variance": args.variance,
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
