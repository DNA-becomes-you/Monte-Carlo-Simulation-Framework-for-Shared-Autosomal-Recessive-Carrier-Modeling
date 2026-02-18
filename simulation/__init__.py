"""
Monte Carlo Simulation Framework for Shared Autosomal Recessive Carrier Modeling.

Pure scientific simulation. No product references.
"""

from simulation.population import (
    generate_population,
    generate_user_vector,
    generate_locus_weights,
    sample_poisson,
)
from simulation.intersection import intersection_size, is_eligible
from simulation.statistics import bootstrap_confidence_interval, evaluate_population
from simulation.experiments import run_experiment, run_grid

__all__ = [
    "sample_poisson",
    "generate_locus_weights",
    "generate_user_vector",
    "generate_population",
    "intersection_size",
    "is_eligible",
    "bootstrap_confidence_interval",
    "evaluate_population",
    "run_experiment",
    "run_grid",
]
