# Reproducibility scope

This repository is a public, self-contained release of the R workflow used for the associated manuscript. The statistical implementation was copied from the manuscript analysis directory and reorganized only where necessary to provide a portable input interface and relative paths.

## Preserved without analytical changes

- nine analysis variables and their transformations;
- five SOM grid sizes, 50 SOM seeds, and 1,000 training iterations;
- candidate partitions from `k = 2` to `k = 8`;
- K-means, Ward.D2, and Gaussian mixture model comparators;
- 200 grouped resamples and shared resampling indices;
- consensus, ARI, AMI, clusterwise Jaccard, support, entropy, and internal indices;
- leave-one-campaign-out transfer analysis;
- variable-module, transformation, missingness, and campaign-exclusion sensitivities;
- plotting code for manuscript Figs. 3–6 and supplementary figures.

## Public-release adaptations

- the 514-event analysis table uses English column names and pseudonymized event and block identifiers;
- internal source locators, coordinates, project paths, and governance fields were removed;
- all file access uses repository-relative paths;
- code and assets for Figs. 1–2 and the graphical abstract were excluded;
- a display-only script that reformatted secondary supplementary tables was excluded.

## Validation completed for this release

- all 514/390/440/459/438/454 analytical cohort counts were reproduced, after the
  chlorophyll-a quality-control rule described in `README.md`;
- the 250 SOM objects and their `k = 2–8` assignments matched the manuscript analysis objects;
- grouped resampling indices and sample-level consensus results matched;
- continuous-gradient summaries matched;
- all 33 result tables shared with the manuscript analysis directory were byte-identical;
- the tracked PNG files are the exact manuscript versions of Figs. 3–6 and Figs. S1–S5.

The complete workflow records, rather than silently removes, model failures. No failures occurred in the validated 200-resample run or in the leave-one-campaign-out transfer analysis.

## First-round revision (version 1.2.0)

- Scripts `10_r1_som_topology.R` to `15_r1_make_tables.R` were added and the figure script was
  extended (Fig. 4D legend, Fig. 5D medians and means, Fig. 6 with matched references, Figs S6-S7).
- Refitting the 200 grouped-resample SOMs reproduces the stored partitions at every k
  (`TableS_r1_resample_som_refit_check.csv`), and refitting the leave-one-survey-out models
  reproduces the stored transfer results exactly.
- On the 390 complete events the native missing-distance SOM reproduces all 50 corresponding
  primary maps, with identical partitions at k = 2-8 (`TableS_r1_missing_distance_identity_check.csv`).
- Transfer counts use only held-out sets whose own reference contains more than one cluster;
  `tests/validate_release.R` recomputes them from `TableS_r1_transfer_vs_matched_omission.csv`.
- Random block omissions list blocks in order of first appearance, and every jittered plotting layer
  has a fixed seed, so the revision tables, `results_facts_r1.json` and all figures are byte-identical
  to those in the manuscript analysis directory (checked before release).
