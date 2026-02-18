"""
Parameter sweep and experiment orchestration.
"""

from __future__ import annotations

import random

from simulation.population import generate_population
from simulation.statistics import evaluate_population


def run_experiment(
    genes: int,
    avg_carriers: float,
    population_structure: float,
    users: int,
    seed: int,
    scenario_index: int,
    max_pairs: int = 500,
    n_bootstrap: int = 1000,
) -> dict:
    """
    Run one experiment. Returns structured result dict.
    Deterministic via seed + scenario_index.
    """
    scenario_seed = seed + scenario_index
    rng = random.Random(scenario_seed)
    population = generate_population(
        users=users,
        panel_size=genes,
        avg_carriers=avg_carriers,
        population_structure=population_structure,
        rng=rng,
    )
    stats = evaluate_population(
        population,
        max_pairs=max_pairs,
        n_bootstrap=n_bootstrap,
        rng=rng,
    )
    return {
        "genes": genes,
        "avgCarriers": avg_carriers,
        "populationStructure": population_structure,
        "pctEligibleZero": round(stats["pctEligibleZero"], 2),
        "pctEligibleZeroCI": (
            round(stats["pctEligibleZeroCI"][0], 2),
            round(stats["pctEligibleZeroCI"][1], 2),
        ),
        "pctOverlapOne": round(stats["pctOverlapOne"], 2),
        "pctOverlapOneCI": (
            round(stats["pctOverlapOneCI"][0], 2),
            round(stats["pctOverlapOneCI"][1], 2),
        ),
        "avgIntersection": round(stats["avgIntersection"], 4),
        "avgIntersectionCI": (
            round(stats["avgIntersectionCI"][0], 4),
            round(stats["avgIntersectionCI"][1], 4),
        ),
    }


def run_grid(
    genes_list: list[int],
    avg_carriers_list: list[float],
    structures_list: list[float],
    users: int,
    seed: int,
    max_pairs: int = 500,
    n_bootstrap: int = 1000,
) -> list[dict]:
    """
    Run full parameter grid. Returns list of experiment results.
    """
    results: list[dict] = []
    scenario_index = 0
    for genes in genes_list:
        for avg_carriers in avg_carriers_list:
            for ps in structures_list:
                exp = run_experiment(
                    genes=genes,
                    avg_carriers=avg_carriers,
                    population_structure=ps,
                    users=users,
                    seed=seed,
                    scenario_index=scenario_index,
                    max_pairs=max_pairs,
                    n_bootstrap=n_bootstrap,
                )
                results.append(exp)
                scenario_index += 1
    return results
