# Monte Carlo Simulation Framework for Shared Autosomal Recessive Carrier Modeling

A standalone, reproducible simulation framework for modeling pairwise overlap of pathogenic autosomal recessive carrier variants in synthetic populations.

## Overview

This framework simulates the probability that two randomly selected individuals share pathogenic variants at the same autosomal recessive locus—a scenario that confers elevated reproductive risk (25% affected offspring when both partners are carriers). The simulation is designed for:

- **Methodology validation** of carrier overlap models
- **Reproducible experiments** with deterministic seeding
- **Parameter sensitivity analysis** across gene panel size, carrier load, and population structure

### Biological Rationale

- **Autosomal recessive inheritance**: Carriers (heterozygotes) are typically asymptomatic. When both partners carry a pathogenic variant at the *same* gene, each pregnancy has a 25% risk of an affected child.
- **Poisson carrier counts**: Pathogenic recessive variants are rare and arise independently across loci. The number of carrier loci per individual is modeled as Poisson(λ), where λ is the average carrier count.
- **Rare event overlap**: Because variants are rare, shared loci between unrelated pairs are uncommon. The framework estimates compatibility (zero shared loci) and the distribution of overlap sizes.
- **Bootstrap CIs**: Uncertainty in point estimates is quantified via bootstrap percentile intervals.

## Model Assumptions

- **Carrier count**: Poisson(λ) per individual.
- **Independent loci**: No linkage disequilibrium; loci are sampled independently (with weighted founder effects).
- **Founder effect**: A fraction of loci have elevated selection weights to model subpopulation structure.
- **Compatibility**: A pair is compatible if the intersection of their carrier locus sets is empty.

## Installation

```bash
pip install -r requirements.txt
```

## Reproducing Manuscript Results

### CLI

```bash
python run_experiments.py \
  --users 3000 \
  --genes 500 \
  --avg-carriers 2,3,4,5 \
  --population-structures 0.0 \
  --seed 42 \
  --output results.json
```

### Example output

```
Genes | AvgCarriers | Structure | Eligible% (CI)           | %Overlap1 (CI)        | AvgInt (CI)
500   | 2.0         | 0         | 96.2% (94.1–97.8)      | 2.8% (1.9–3.9)   | 0.0241 (0.0182–0.0312)
...
```

## Running the Streamlit UI

```bash
streamlit run app.py
```

The UI allows interactive parameter adjustment, runs simulations on demand, displays results with 95% CIs, and provides a histogram of intersection sizes. JSON results can be downloaded.

## Project Structure

```
.
├── simulation/
│   ├── population.py    # Poisson sampling, locus weights, population generation
│   ├── intersection.py  # O(n+m) intersection and eligibility
│   ├── statistics.py    # Bootstrap CIs, population evaluation
│   ├── experiments.py   # run_experiment, run_grid
│   └── __init__.py
├── app.py               # Streamlit UI
├── run_experiments.py   # CLI entrypoint
├── requirements.txt
├── README.md
├── LICENSE
└── CITATION.cff
```

## Ethical Framing

This framework models **shared pathogenic recessive carrier overlap** for research and methodology validation. It:

- Does **not** perform trait optimization
- Does **not** perform polygenic scoring
- Does **not** select for or against non-pathogenic variants
- Is intended for **academic and research purposes** only

Interpretation of results in clinical or reproductive decision-making requires appropriate professional guidance.

## License

MIT.
