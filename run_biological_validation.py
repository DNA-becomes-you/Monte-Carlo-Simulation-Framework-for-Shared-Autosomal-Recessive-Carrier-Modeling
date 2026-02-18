#!/usr/bin/env python3
"""
DNABY DATE — Biological Validation Test Script (Zone 2)

Tests multiple biological scenarios: gene panel size × carrier load.
Validates compatibility distribution across parameter combinations.
No database usage. Standard library only. Deterministic with --seed.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

# Reuse population logic from generator (in-memory)
from generate_population import (
    generate_shared_pool,
    generate_user_vector,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run biological validation across parameter combinations"
    )
    parser.add_argument(
        "--users",
        type=int,
        default=500,
        help="Number of users per scenario (default: 500)",
    )
    parser.add_argument(
        "--genes",
        type=str,
        default="300,500,1000",
        help="Comma-separated gene panel sizes (default: 300,500,1000)",
    )
    parser.add_argument(
        "--avg-carriers",
        type=str,
        default="2,3,4,5",
        help="Comma-separated avg-carriers values (default: 2,3,4,5)",
    )
    parser.add_argument(
        "--variance",
        type=float,
        default=1.5,
        help="Variance around avg-carriers (default: 1.5)",
    )
    parser.add_argument(
        "--population-structure",
        type=float,
        default=0.02,
        help="Population structure fraction (default: 0.02)",
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
        default=Path("fixtures/synthetic/biological_validation.json"),
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


def generate_population(
    users: int,
    panel_size: int,
    avg_carriers: float,
    variance: float,
    population_structure: float,
    rng: random.Random,
) -> list[dict]:
    """Generate synthetic carrier population in memory."""
    shared_pool = generate_shared_pool(panel_size, population_structure, rng)
    population: list[dict] = []
    for i in range(users):
        user = generate_user_vector(
            user_index=i,
            panel_size=panel_size,
            avg_carriers=avg_carriers,
            variance=variance,
            shared_pool=shared_pool,
            population_structure=population_structure,
            rng=rng,
        )
        population.append(user)
    return population


def compute_metrics(
    population: list[dict],
    max_pairs: int = 1000,
    rng: random.Random | None = None,
) -> dict:
    """
    Compute eligibility and overlap metrics.
    pctEligibleZero: intersection == 0
    pctEligibleOne: intersection <= 1
    pctOverlapOne: intersection == 1
    """
    n = len(population)
    total_pairs = n * (n - 1) // 2

    if total_pairs == 0:
        return {
            "pctEligibleZero": 0.0,
            "pctEligibleOne": 0.0,
            "pctOverlapOne": 0.0,
            "avgIntersection": 0.0,
            "maxIntersection": 0,
        }

    all_pairs: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            all_pairs.append((i, j))

    if len(all_pairs) <= max_pairs:
        pairs = all_pairs
    else:
        r = rng if rng is not None else random.Random()
        pairs = r.sample(all_pairs, max_pairs)

    eligible_zero = 0
    eligible_one = 0
    overlap_one = 0
    total_intersection = 0
    max_intersection = 0

    for i, j in pairs:
        a = population[i]["pathogenicLoci"]
        b = population[j]["pathogenicLoci"]
        sz = compute_intersection_size(a, b)
        if sz == 0:
            eligible_zero += 1
            eligible_one += 1
        elif sz == 1:
            eligible_one += 1
            overlap_one += 1
        total_intersection += sz
        max_intersection = max(max_intersection, sz)

    pair_count = len(pairs)
    return {
        "pctEligibleZero": 100 * eligible_zero / pair_count,
        "pctEligibleOne": 100 * eligible_one / pair_count,
        "pctOverlapOne": 100 * overlap_one / pair_count,
        "avgIntersection": total_intersection / pair_count,
        "maxIntersection": max_intersection,
    }


def run_scenario(
    genes: int,
    avg_carriers: float,
    users: int,
    variance: float,
    population_structure: float,
    base_seed: int,
) -> dict:
    """Run one scenario (genes × avg_carriers) with deterministic seed."""
    # Derive unique seed per combination for determinism
    scenario_seed = base_seed + genes * 1000 + int(avg_carriers * 100)
    rng = random.Random(scenario_seed)
    population = generate_population(
        users, genes, avg_carriers, variance, population_structure, rng
    )
    metrics = compute_metrics(population, max_pairs=1000, rng=rng)
    return {
        "genes": genes,
        "avgCarriers": avg_carriers,
        "pctEligibleZero": round(metrics["pctEligibleZero"], 2),
        "pctEligibleOne": round(metrics["pctEligibleOne"], 2),
        "pctOverlapOne": round(metrics["pctOverlapOne"], 2),
        "avgIntersection": round(metrics["avgIntersection"], 4),
        "maxIntersection": metrics["maxIntersection"],
    }


def print_table(results: list[dict]) -> None:
    """Print comparison table."""
    header = (
        f"{'Genes':<8} | {'AvgCarriers':<10} | {'Eligible(0)':<12} | "
        f"{'Eligible(<=1)':<14} | {'%Overlap1':<10} | {'AvgInt':<8}"
    )
    print(header)
    print("-" * len(header))
    for r in results:
        row = (
            f"{r['genes']:<8} | "
            f"{r['avgCarriers']:<10.1f} | "
            f"{r['pctEligibleZero']:<12.2f} | "
            f"{r['pctEligibleOne']:<14.2f} | "
            f"{r['pctOverlapOne']:<10.2f} | "
            f"{r['avgIntersection']:<8.4f}"
        )
        print(row)


def main() -> None:
    args = parse_args()

    genes_list = [int(x.strip()) for x in args.genes.split(",")]
    avg_carriers_list = [float(x.strip()) for x in args.avg_carriers.split(",")]

    results: list[dict] = []

    for genes in genes_list:
        for avg_carriers in avg_carriers_list:
            res = run_scenario(
                genes=genes,
                avg_carriers=avg_carriers,
                users=args.users,
                variance=args.variance,
                population_structure=args.population_structure,
                base_seed=args.seed,
            )
            results.append(res)

    output = {
        "metadata": {
            "users": args.users,
            "variance": args.variance,
            "populationStructure": args.population_structure,
            "seed": args.seed,
        },
        "results": results,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    print(f"Results saved to {args.output}\n")
    print_table(results)


if __name__ == "__main__":
    main()
