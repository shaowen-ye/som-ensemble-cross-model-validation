<img width="432" height="28" alt="image" src="https://github.com/user-attachments/assets/d5e636c3-0756-4b90-9f19-7ad3c3e5575f" /># SOM ensemble and cross-model validation

This repository provides the analysis-ready data and R workflow associated with the manuscript:

> **When are self-organizing map classes defensible? An ensemble and cross-model validation workflow for environmental monitoring data**
> — manuscript prepared by Shaowen Ye et al. 

The workflow evaluates whether self-organizing map (SOM) structures persist across map sizes, random initializations, grouped resampling, alternative clustering algorithms, monitoring campaigns, variable modules, and missing-data conditions. It is designed to distinguish reproducible continuous gradients from potentially fragile hard partitions.

## Contents

- `data/analysis_events.csv`: analysis-ready event table used by the public workflow.
- `data/data_dictionary.csv`: variable descriptions and units.
- `R/`: data preparation, SOM ensemble, cross-model validation, transfer, sensitivity, figure, and table scripts.
- `run_all.R`: single entry point for the complete analysis.
- `outputs/figures/`: associated analytical figures.
- `outputs/tables/`: key statistical tables generated directly from the analyses.
- `outputs/results/sessionInfo.txt`: R session information from the validated release run.

## Analysis design

The primary analysis uses 391 events with complete observations for nine indicators. The public event table retains 514 events and their original missing-data structure so that the variable-module and native-missingness sensitivity analyses can be reproduced from the same input.

The public workflow includes:

1. log10 transformation of chlorophyll-a and nutrient variables, followed by z standardization;
2. 250 SOM fits from five hexagonal grids and 50 prespecified seeds;
3. candidate partitions for `k = 2–8`;
4. controlled comparison with K-means, Ward.D2, and Gaussian mixture models;
5. 200 grouped resamples using the same training blocks for all four model families;
6. consensus support, entropy, ARI, AMI, clusterwise Jaccard, and auxiliary internal indices;
7. leave-one-campaign-out transfer and prespecified sensitivity analyses;
8. scripted generation of associated analytical figures and key statistical tables.

Transformations and scaling in resampling and transfer analyses are estimated from training observations only.

## Reproduce the analysis

The validated release used R 4.5.3. From the repository root:

```r
install.packages("renv")
renv::restore()
```

Then run:

```bash
Rscript run_all.R
```

The complete ensemble and resampling workflow is computationally intensive. Generated model objects are stored locally under `outputs/results/models/` and are not tracked by Git.

To verify the released data and tracked outputs without refitting the models:

```bash
Rscript tests/validate_release.R
```

## Reproducibility notes

- All random seeds and model budgets are defined in `R/00_config.R`.
- All model families use the same grouped resampling indices.
- The workflow records model failures rather than silently discarding them.
- Numerical labels are used for hard partitions; the repository does not treat consensus labels as ecological truth.
- Input identifiers and resampling blocks are pseudonymized.

## License

R code is released under the MIT License. The analysis-ready data are released under the Creative Commons Attribution 4.0 International License; see `DATA_LICENSE.md`.

## Maintainer

Shaowen Ye
GitHub: [Shaowen-Ye](https://github.com/Shaowen-Ye)
