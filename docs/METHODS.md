# Methods and estimands

## Primary spatial summary

The reference tumor population is the author's exact `tumor cells` annotation; the immune population is `CD8+ T cells`. Four debris/ambiguous labels are excluded as listed in the source audit. Distances are Euclidean distances from a CD8 centroid to its nearest eligible tumor centroid, within the same region.

At a preselected radius of 80 native pixels, **access fraction** is the proportion of eligible CD8 cells within that distance. **Available-site fraction** is the corresponding proportion among all eligible non-tumor cells. Their ratio describes CD8 proximity relative to the available non-tumor positions. A centroid distance is not a histological boundary distance or evidence of physical contact. Regions require at least 20 tumor cells, 10 CD8 cells and 20 eligible non-tumor cells. Seventeen insufficient regions are retained with an explicit undefined estimand.

Two lower-tail random-labelling nulls test unusually low CD8 access. Both keep tumor locations and actual cell coordinates fixed. The global null allocates the observed number of CD8 labels without replacement among eligible non-tumor positions. The block null also fixes CD8 abundance in 400-pixel blocks anchored at the region minimum coordinate. It asks a more local conditional question and can remove broad-scale spatial differences. Neither assumes cells are uniform in a rectangular field.

Each uses 199 randomizations. P-values include the observed-case correction `(1 + lower-tail count) / 200`; the minimum is 0.005. Benjamini–Hochberg adjustment is applied separately to the two families of 123 eligible regional tests. These exploratory, dependent regional summaries are not a clinical diagnostic test. Sparse effects may be undetectable at this resolution.

Regional access ratios are averaged within patients; patients are then equally weighted. The descriptive mean interval resamples the 35 patient means 2,000 times with seed 17. It does not correct selection of patients or prove generalization to all colorectal cancers.

## Sensitivity and inferred compartments

Each region has seven retained perturbations: radii 40/80/120 px; centroid jitter with SD 3/8 px; random 80% cell retention; and inclusion of mixed tumor/immune annotations as tumor. Each stochastic condition uses one predetermined regional seed. These are scoped sensitivity checks, not a fitted segmentation-error distribution.

The displayed tumor-rich territory is the tumor fraction among the nearest 20 eligible cells. A grid point is supported only within 80 px of an observed eligible centroid. Fractions at least 0.5 define tumor-rich territory for the conditional movement model. This new compartment inference is distinct from the author's neighborhood labels and from the primary centroid-distance estimand.

## Protein reconstruction

The target is `log1p` of author-processed CD8 fluorescence intensity. Train patients have ID modulo 5 in 1/2/3 (21 patients); development patients have remainder 4 (seven); reporting patients have remainder 0 (seven). All regions belonging to a patient remain together.

The marker ridge model uses log-transformed CD3, CD4, CD45, CD2, CD5, CD7, FOXP3, CD20, CD68, Cytokeratin, aSMA and CD45RO. It standardizes predictors using training means/SDs, fits centered ridge regression, and clips negative predicted log intensities to zero. CD8, author cell types, author neighborhoods, patient identity and coordinates are excluded from predictor features. Labels only define the common analysis exclusion mask.

Development compares ridge penalties 1/10/100 and spatial-neighbor counts 5/15/30 using equal-patient mean RMSE over both mask patterns, seed 17. Selected values are penalty 100 and 30 neighbors. The simple baseline imputes the training median log intensity. The spatial method uses inverse-squared-distance weighting among observed cells in the same region. It may use the unmasked channel values in a reporting region; the cross-patient marker model does not. These are explicitly different acquisition assumptions.

Masks hide floor(30% of eligible cells), either uniformly at random or in a contiguous x-ranked band. Each reporting case uses seed 101, 107 or 113. Metrics use hidden cells only for channel RMSE/MAE. There are 28 reporting regions × 2 patterns × 3 seeds × 3 methods = 504 cases. Reporting aggregates average cases within each patient before averaging patients.

For a downstream preservation check, marker-high cells exceed the training 90th percentile of log CD8 among eligible training cells. Their proximity fraction is compared before and after imputation among non-tumor cells; unmasked values remain unchanged. This illustrative signal threshold is not a validated CD8 T-cell phenotype gate. It is distinct from the reference author-annotated CD8 analysis.

A ridge residual interval adds the training residual 5th/95th percentiles to predictions. Its empirical reporting coverage is measured, not guaranteed. It is not an independently calibrated conditional uncertainty model.

## Conditional movement

Walkers begin at observed supported coordinates in the leftmost 20% of positions with tumor fraction below 0.2. A proposed 8 px step mixes a random direction with a 0.5 bias toward the nearest tumor centroid. Proposals reflect at the observed bounding box; unsupported voids are rejected. Every proposed step into a tumor-rich position, including movement within that area, is accepted with the assumed permeability probability.

Each run contains 64 walkers and 120 steps. “Entry” means ever entering the inferred tumor-rich territory. “Centroid contact” means ever coming within 20 px of a tumor centroid; it is not observed membrane contact. Neither step units nor permeability are calibrated to migration measurements. All paths may be exported; the display shows at most 32 paths for legibility.

The frozen movement design uses the first two lexicographic regions of patients 5/20/35, permeability 0.1/0.5/1.0 and seeds 101/107/113/127/131: 90 runs. Same-seed proposals are reused across permeability values, while rejected proposals can change later positions. These are conditional scenarios, not independent biological experiments.
