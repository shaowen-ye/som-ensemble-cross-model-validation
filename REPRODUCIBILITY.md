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

- all 514/391/440/459/439/454 cohort counts were reproduced;
- the 250 SOM objects and their `k = 2–8` assignments matched the manuscript analysis objects;
- grouped resampling indices and sample-level consensus results matched;
- continuous-gradient summaries matched;
- all 33 result tables shared with the manuscript analysis directory were byte-identical;
- the tracked PNG files are the exact manuscript versions of Figs. 3–6 and Figs. S1–S5.

The complete workflow records, rather than silently removes, model failures. No failures occurred in the validated 200-resample run or in the leave-one-campaign-out transfer analysis.
