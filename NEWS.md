# MRKCFC 0.2.0

This is the release that produced the results of the accompanying manuscript.

* `MKCFC()` fits the multivariate K-centers model. `robust = TRUE` with the
  default `mfpca_scatter = "sign"` is the proposed method, with Kendall FPCA in
  the first stage and the spatial-sign scatter in the second; `robust = TRUE`
  with `mfpca_scatter = "moment"` keeps the second-moment matrix in the second
  stage; `robust = FALSE` is the non-robust comparator. `init.labels` supplies a
  starting partition directly, and `basis_common = TRUE` shares one set of
  eigenfunctions across classes for ablation.
* `MFPCA_k()`, `MKCFC_init()` and `predicted_yhat()` expose the per-class
  two-stage estimator, the starting partition and the class-specific
  reconstruction that the loop iterates over.
* `KFPCA()` is the Kendall FPCA runtime, derived from KFPCA 2.0.
* `select_K()`, `ch_index()`, `silhouette_index()` and `whiten_scores()` choose
  the number of classes on the multivariate scores of one global fit, with a
  robust Mahalanobis whitening; `silhouette_score()` is the curve-based
  silhouette.
* `rre()`, `robust_normalize()`, `back_transform()`, `gen_scores()`,
  `vm_intermediate()`, `spatial_median()`, `sscm()`, `robust_stage2()`,
  `principal_angles()` and `stage2_diagnostic()` are the metrics, preprocessing
  and simulation helpers of the manuscript.

The version number continues from the package this code was split from. The
univariate method is distributed separately as `RKCFC`.
