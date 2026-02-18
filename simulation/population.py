"""
Population generation for shared autosomal recessive carrier modeling.

Poisson carrier counts. Weighted founder-effect locus sampling.
Pure functions, RNG-injected, fully deterministic.
"""

from __future__ import annotations

import hashlib
import json
import math
import random


def sample_poisson(lmbda: float, rng: random.Random) -> int:
    """
    Sample from Poisson(λ). For λ <= 0 returns 0.

    Recessive pathogenic variants are rare independent events; Poisson
    models count of rare independent events (mean = variance = λ).
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
    n_founder = max(1, int(panel_size * population_structure))
    n_founder = min(n_founder, panel_size)
    founder_loci = rng.sample(range(panel_size), n_founder)
    for idx in founder_loci:
        weights[idx] = 10.0
    return weights


def generate_user_vector(
    user_index: int,
    panel_size: int,
    avg_carriers: float,
    locus_weights: list[float],
    rng: random.Random,
) -> dict:
    """
    Generate a single individual's pathogenic carrier loci.
    Carrier count ~ Poisson(avg_carriers). Weighted sampling without replacement.
    """
    carrier_count = sample_poisson(avg_carriers, rng)
    carrier_count = min(carrier_count, panel_size)

    if carrier_count == 0:
        return {
            "userIndex": user_index,
            "pathogenicLoci": [],
            "hash": hash_loci([]),
        }

    available_indices = list(range(panel_size))
    available_weights = locus_weights.copy()
    selected: list[int] = []

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
    assert len(pathogenic_loci) == len(set(pathogenic_loci))
    assert all(0 <= x < panel_size for x in pathogenic_loci)

    return {
        "userIndex": user_index,
        "pathogenicLoci": pathogenic_loci,
        "hash": hash_loci(pathogenic_loci),
    }


def generate_population(
    users: int,
    panel_size: int,
    avg_carriers: float,
    population_structure: float,
    rng: random.Random,
) -> list[dict]:
    """Generate synthetic carrier population."""
    locus_weights = generate_locus_weights(panel_size, population_structure, rng)
    population: list[dict] = []
    for i in range(users):
        user = generate_user_vector(
            user_index=i,
            panel_size=panel_size,
            avg_carriers=avg_carriers,
            locus_weights=locus_weights,
            rng=rng,
        )
        population.append(user)
    return population
