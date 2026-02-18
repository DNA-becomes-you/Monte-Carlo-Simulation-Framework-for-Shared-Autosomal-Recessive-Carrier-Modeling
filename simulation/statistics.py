"""
Bootstrap confidence intervals and population evaluation.
"""

from __future__ import annotations

import random

from simulation.intersection import intersection_size


def bootstrap_confidence_interval(
    values: list[float],
    rng: random.Random,
    n_bootstrap: int = 1000,
) -> tuple[float, float]:
    """
    Bootstrap 95% CI for the mean (percentile method).
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


def evaluate_population(
    population: list[dict],
    max_pairs: int = 500,
    n_bootstrap: int = 1000,
    rng: random.Random | None = None,
) -> dict:
    """
    Sample pairs, compute metrics, and bootstrap 95% CIs.
    Returns structured dict with point estimates and CIs.
    """
    n = len(population)
    total_pairs = n * (n - 1) // 2

    if total_pairs == 0:
        return {
            "pctEligibleZero": 0.0,
            "pctEligibleZeroCI": (0.0, 0.0),
            "pctOverlapOne": 0.0,
            "pctOverlapOneCI": (0.0, 0.0),
            "avgIntersection": 0.0,
            "avgIntersectionCI": (0.0, 0.0),
            "pairCount": 0,
            "intersectionDistribution": [],
        }

    all_pairs: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            all_pairs.append((i, j))

    r = rng if rng is not None else random.Random()
    if len(all_pairs) <= max_pairs:
        pairs = all_pairs
    else:
        pairs = r.sample(all_pairs, max_pairs)

    pair_results: list[dict] = []
    for i, j in pairs:
        a = population[i]["pathogenicLoci"]
        b = population[j]["pathogenicLoci"]
        sz = intersection_size(a, b)
        pair_results.append({
            "intersection_size": sz,
            "eligible": sz == 0,
            "overlap_one": sz == 1,
        })

    pair_count = len(pair_results)
    eligible_count = sum(1 for p in pair_results if p["eligible"])
    overlap_one_count = sum(1 for p in pair_results if p["overlap_one"])
    total_intersection = sum(p["intersection_size"] for p in pair_results)

    pct_eligible = 100 * eligible_count / pair_count
    pct_overlap_one = 100 * overlap_one_count / pair_count
    avg_intersection = total_intersection / pair_count

    eligible_vals = [1.0 if p["eligible"] else 0.0 for p in pair_results]
    overlap_one_vals = [1.0 if p["overlap_one"] else 0.0 for p in pair_results]
    intersection_vals = [float(p["intersection_size"]) for p in pair_results]

    ci_eligible = bootstrap_confidence_interval(eligible_vals, r, n_bootstrap)
    ci_overlap = bootstrap_confidence_interval(overlap_one_vals, r, n_bootstrap)
    ci_avg = bootstrap_confidence_interval(intersection_vals, r, n_bootstrap)

    # Build intersection distribution (counts per value)
    dist: dict[int, int] = {}
    for p in pair_results:
        sz = p["intersection_size"]
        dist[sz] = dist.get(sz, 0) + 1
    intersection_distribution = [
        {"intersection": k, "count": v} for k, v in sorted(dist.items())
    ]

    return {
        "pctEligibleZero": pct_eligible,
        "pctEligibleZeroCI": (100 * ci_eligible[0], 100 * ci_eligible[1]),
        "pctOverlapOne": pct_overlap_one,
        "pctOverlapOneCI": (100 * ci_overlap[0], 100 * ci_overlap[1]),
        "avgIntersection": avg_intersection,
        "avgIntersectionCI": ci_avg,
        "pairCount": pair_count,
        "intersectionDistribution": intersection_distribution,
    }
