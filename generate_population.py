#!/usr/bin/env python3
"""
DNABY DATE — Synthetic Recessive Carrier Population Generator (Zone 2)

Generates synthetic pathogenic carrier vectors for carrier compatibility simulation.
Models rare autosomal recessive variants (1–5 carriers per individual).
No database access. No PHI. Output is JSON only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate synthetic recessive carrier vectors"
    )
    parser.add_argument(
        "--users",
        type=int,
        required=True,
        help="Number of users to generate",
    )
    parser.add_argument(
        "--genes",
        type=int,
        default=500,
        help="Gene panel size (known autosomal recessive loci) (default: 500)",
    )
    parser.add_argument(
        "--avg-carriers",
        type=float,
        default=3.0,
        help="Average pathogenic variants per individual (default: 3)",
    )
    parser.add_argument(
        "--variance",
        type=float,
        default=1.5,
        help="Variance (kept for CLI compatibility; carrier count uses Poisson, not variance)",
    )
    parser.add_argument(
        "--population-structure",
        type=float,
        default=0.02,
        help="Fraction of loci common within subpopulations (default: 0.02)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for determinism",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output JSON file path",
    )
    return parser.parse_args()


def sample_poisson(lmbda: float, rng: random.Random) -> int:
    """
    Sample from Poisson(λ). For λ <= 0 returns 0.
    Recessive pathogenic variants are rare independent events;
    Poisson models count of rare independent events (mean = variance = λ).
    """
    if lmbda <= 0:
        return 0
    L = math.exp(-lmbda)
    k = 0
    p = 1.0
    while p > L:
        k += 1
        p *= rng.random()
    return max(0, k - 1)


def hash_loci(loci: list[int]) -> str:
    """SHA256 of sorted array (JSON-encoded)."""
    payload = json.dumps(sorted(loci)).encode()
    return hashlib.sha256(payload).hexdigest()


def generate_locus_weights(
    panel_size: int,
    population_structure: float,
    rng: random.Random,
) -> list[float]:
    """
    Generate per-locus selection weights.
    A small subset of loci have elevated weight (founder effect).
    """
    weights = [1.0] * panel_size

    # Number of founder loci (scales with population_structure)
    n_founder = max(1, int(panel_size * population_structure))
    founder_loci = rng.sample(range(panel_size), n_founder)

    # Increase weight for founder loci
    for idx in founder_loci:
        weights[idx] = 10.0  # 10x frequency boost (founder effect)

    return weights


def generate_user_vector(
    user_index: int,
    panel_size: int,
    avg_carriers: float,
    variance: float,
    locus_weights: list[float],
    rng: random.Random,
) -> dict:
    """
    Generate a single user's pathogenic carrier loci.
    Uses weighted sampling to model founder effects.
    variance is kept for API compatibility but does not affect carrier count.
    """
    # 1) carrier_count ~ Poisson(avg_carriers)
    # Recessive pathogenic variants are rare independent events; Poisson
    # models count of rare independent events (mean = variance = λ).
    carrier_count = sample_poisson(avg_carriers, rng)
    carrier_count = min(carrier_count, panel_size)

    if carrier_count == 0:
        return {
            "userIndex": user_index,
            "pathogenicLoci": [],
            "hash": hash_loci([]),
        }

    # 2) Weighted sampling without replacement
    available_indices = list(range(panel_size))
    available_weights = locus_weights.copy()
    selected = []

    for _ in range(carrier_count):
        if not available_indices:
            break

        total_weight = sum(available_weights)
        r = rng.random() * total_weight
        cumulative = 0.0

        for i, (idx, w) in enumerate(zip(available_indices, available_weights)):
            cumulative += w
            if cumulative >= r:
                selected.append(idx)
                available_indices.pop(i)
                available_weights.pop(i)
                break

    pathogenic_loci = sorted(selected)

    # Validation
    assert len(pathogenic_loci) == len(set(pathogenic_loci)), "Duplicate loci"
    assert all(0 <= x < panel_size for x in pathogenic_loci), "Out-of-range"

    return {
        "userIndex": user_index,
        "pathogenicLoci": pathogenic_loci,
        "hash": hash_loci(pathogenic_loci),
    }


def intersection_size(a: list[int], b: list[int]) -> int:
    """Two-pointer O(n+m) intersection size. Arrays must be sorted."""
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


def is_eligible(a: list[int], b: list[int]) -> bool:
    """Eligible = no shared pathogenic loci (empty intersection)."""
    return intersection_size(a, b) == 0


def main() -> None:
    args = parse_args()

    rng = random.Random(args.seed) if args.seed is not None else random.Random()
    panel_size = args.genes
    locus_weights = generate_locus_weights(
        panel_size,
        args.population_structure,
        rng,
    )

    population: list[dict] = []
    for i in range(args.users):
        user = generate_user_vector(
            user_index=i,
            panel_size=panel_size,
            avg_carriers=args.avg_carriers,
            variance=args.variance,
            locus_weights=locus_weights,
            rng=rng,
        )
        population.append(user)

    output = {
        "metadata": {
            "users": args.users,
            "genePanelSize": panel_size,
            "avgCarriers": args.avg_carriers,
            "variance": args.variance,
            "populationStructure": args.population_structure,
            "seed": args.seed,
        },
        "population": population,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    print(f"Wrote {args.users} users to {args.output}")

    # Distribution summary: sample first 100 pairs
    n_pairs = min(100, args.users * (args.users - 1) // 2)
    eligible_count = 0
    overlap_one_count = 0
    total_intersection = 0
    pair_count = 0

    for i in range(args.users):
        for j in range(i + 1, args.users):
            if pair_count >= n_pairs:
                break
            a = population[i]["pathogenicLoci"]
            b = population[j]["pathogenicLoci"]
            sz = intersection_size(a, b)
            if sz == 0:
                eligible_count += 1
            elif sz == 1:
                overlap_one_count += 1
            total_intersection += sz
            pair_count += 1
        if pair_count >= n_pairs:
            break

    if pair_count > 0:
        pct_eligible = 100 * eligible_count / pair_count
        pct_overlap_one = 100 * overlap_one_count / pair_count
        avg_intersection = total_intersection / pair_count
        print(f"\n--- Distribution summary (sample of {pair_count} pairs) ---")
        print(f"  % eligible (no shared loci): {pct_eligible:.1f}%")
        print(f"  % with exactly 1 overlap:    {pct_overlap_one:.1f}%")
        print(f"  Average intersection size:   {avg_intersection:.2f}")


if __name__ == "__main__":
    main()
