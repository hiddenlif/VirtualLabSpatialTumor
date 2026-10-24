# Investigation: recovering measurements does not guarantee recovering a conclusion

The first completed local investigation measures CD8-to-tumor proximity in published colorectal tissue and tests how missing protein measurements affect a downstream spatial summary. It supports inspectable regional description and evaluated reconstruction. It does not establish functional immune exclusion or predict immune-cell migration.

## Reference cohort

The complete source contains 258,385 cells from 140 regions and 35 patients. After retaining but excluding 17,831 debris/ambiguous annotations, 240,554 cells are eligible. Seventeen regions lack the minimum tumor/CD8 counts; all remain in the evidence bank. All 35 patients retain at least one eligible region.

At the planned 80 px radius, the equal-patient mean CD8 access ratio is **1.029**, with a descriptive patient bootstrap interval of **0.969–1.089**. It includes 1. Patient means range from about 0.50 to 1.59, so the cohort mean should not replace inspection of individual regions.

No region passed the adjusted 0.05 depletion threshold under either the global or block conditional null. With 199 replicates, minimum attainable p is 0.005; adjustment across 123 regions limits discovery of sparse effects. The design was retained rather than tuned after the result. Zero discoveries is not proof that immune exclusion is absent. These centroids and conditional nulls do not establish an independent histological exclusion diagnosis.

## Reconstruction on seven separate reporting patients

All methods use the same masks and cases. Rows below average within patients before averaging the seven patients.

| Missingness | Method | RMSE, log1p CD8 | Marker-high proximity absolute error, percentage points |
|---|---|---:|---:|
| Random 30% | Training median | 2.318 | 2.293 |
| Random 30% | Spatial interpolation | 2.393 | 2.221 |
| Random 30% | Marker ridge regression | 1.957 | 2.341 |
| Contiguous band | Training median | 2.335 | 4.659 |
| Contiguous band | Spatial interpolation | 2.212 | 4.675 |
| Contiguous band | Marker ridge regression | 1.976 | 4.274 |

Ridge reduced intensity error in both missingness conditions. Its downstream proximity error was slightly worse under random masking and lower under band masking. This mixed result does not justify treating lower channel RMSE as evidence of a preserved spatial conclusion. The spatial method also illustrates that local interpolation can help with one missingness pattern and perform poorly with another.

Training-residual 90% intervals covered **90.4%** of reporting hidden values under random masks and **90.0%** under band masks after equal-patient averaging. These observed coverage values are not guarantees, and do not establish conditional calibration or transfer to new acquisition settings.

## A retained failure

In **P35 / reg070_A**, the band mask with seed 101 hides 533 cells. The true marker-high proximity fraction is **71.4%**; all three reconstructions yield **92.3%**, an error of **20.9 percentage points**. Ridge has the best channel RMSE of the three in this case (2.044 log units), yet the same large downstream error. The original measured values remain accessible for validation; no case was removed because it contradicts the overall channel-error advantage.

The post-run failure inventory also includes **P35 / reg070_B**, band seed 113: marker-high proximity shifts from **46.7% to 28.6%**. These are illustrative intensity-defined populations, not independent CD8 T-cell identity calls.

## Movement is a separate conditional experiment

The 90 seeded movement runs use six preselected measured regions. Across 30 runs at each assumed permeability, mean tumor-rich entry fractions are **0.16%, 0.31% and 0.31%** for 0.1, 0.5 and 1.0. Mean proximity within 20 px of a tumor centroid is **95.05%** at all three values.

Thus this particular walk model and region selection yield many centroid-proximity events and very few entries into the inferred tumor-rich territory. Changing the permeability assumption does not materially change the proximity endpoint here. This is a limitation of the constructed scenario and its endpoint definitions, not evidence about stromal permeability, T-cell infiltration, physical contact, or therapy effectiveness.

## Decision and next evidence

Keep measured observations, author annotations, new compartment inference, reconstructed channel values and simulated paths distinct. Use reconstruction only with an explicit downstream preservation check. Retain the median and spatial baselines alongside the fitted marker model, including cases where the simpler method is competitive.

A stronger next scientific update would link a justified image subset and segmentation masks to these exact cell IDs, establish coordinate scaling and boundary interpretation, and evaluate uncertainty against independent image/phenotype annotations. Movement calibration would separately require observed trajectories and tissue-barrier measurements. These extensions are outside this first local release; no unsupported patient-specific 3D reconstruction is supplied.

Evidence: `reports/study-v1/` retains the plan-bound model, all reference regions, 980 sensitivity rows, all 504 reconstruction cases and all 90 movement outputs. `reports/secondary-analysis-v1.json` labels the post-run descriptive failure/movement summaries. The new computation reproduced eight substantive study artifacts exactly; the summary matched apart from its execution timestamp.
