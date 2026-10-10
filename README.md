# MRKCFC

Multivariate robust K-centers functional clustering for several longitudinal
variables that are observed sparsely and irregularly and may be non-Gaussian.

The method clusters subjects on the joint course of several markers. Each class
is described by a multivariate functional principal component analysis in two
stages, after Happ and Greven (2018): a marginal principal component analysis of
each variable, then an eigendecomposition of the stacked scores, which gives the
multivariate eigenfunctions. Subjects are reassigned to the class whose model
reconstructs their observed trajectories best, and the two steps alternate, as
in the K-centers clustering of Chiou and Li (2007).

Both stages can be estimated from second moments, and both then follow the few
subjects whose values run far above the rest, which in medical longitudinal
data is common. MRKCFC makes both stages robust. The first stage uses Kendall's
tau functional principal component analysis (Zhong et al., 2022), whose
eigenspace agrees with that of the covariance operator under an elliptical
process without a moment assumption on the scores. The second stage replaces
the sample second-moment matrix of the stacked scores by their spatial-sign
covariance matrix, whose eigenvectors agree with those of the scatter matrix
under elliptical symmetry while each subject enters through its direction only;
the spread along each direction is then recovered by a robust scale. A robust
first stage with a second-moment second stage loses at the second stage what it
gained at the first, and the package keeps that construction available as a
comparator.

## Installation

```r
# install.packages("remotes")
remotes::install_github("BingbingFan/MRKCFC")

# the release that produced the results of the accompanying manuscript
remotes::install_github("BingbingFan/MRKCFC@v0.2.0")
```

## Quick start

```r
library(MRKCFC)
data(pbc2)

# Two markers, subjects with at least two visits on both.
vars <- c("serBilir", "albumin")
d <- pbc2[stats::complete.cases(pbc2[, c("id", "year", vars)]), ]
keep <- names(which(table(d$id) >= 2))
d <- d[d$id %in% keep, ]

# Time rescaled to [0, 1], so that the basis and the bandwidths below are on
# the same footing as the observation grid.
rng <- range(d$year)
d$t <- (d$year - rng[1]) / diff(rng)

ids <- unique(d$id)
Ly <- lapply(vars, function(v) lapply(ids, function(i) d[[v]][d$id == i]))
Lt <- rep(list(lapply(ids, function(i) d$t[d$id == i])), length(vars))

optns <- list(
  interval = c(0, 1), nK = 2, bw = 0.1, bwK = "GCV", bwmean = "GCV",
  nRegGrid = 51,
  basis = fda::create.bspline.basis(c(0, 1), nbasis = 13, norder = 4))

fit <- MKCFC(Ly = Ly, Lt = Lt, K = 2, robust = TRUE, optns = optns,
             mfpca_scatter = "sign", maxIter = 20, parallel = FALSE,
             trace = 0)
table(fit$res_clust$g)
fit$conv
```

`Ly` and `Lt` are lists with one element per variable, each a list with one
vector per subject. The data-frame form (`data`, `id`, `ys`, `timevar`) is
accepted instead.

## One function, three models

| Model | `robust` | `mfpca_scatter` | First stage | Second stage |
|---|---|---|---|---|
| MRKCFC-S, the proposed method | `TRUE` | `"sign"` | Kendall's tau FPCA | spatial-sign scatter |
| MRKCFC-M | `TRUE` | `"moment"` | Kendall's tau FPCA | second-moment matrix |
| MKCFC | `FALSE` | `"moment"` | ordinary FPCA | second-moment matrix |

The two branches take different option lists. With `robust = TRUE`: `interval`,
`nK`, `bw`, `bwK`, `bwmean`, `nRegGrid` and `basis`. With `robust = FALSE`:
`dataType`, `methodMuCovEst`, `nRegGrid`, `methodSelectK` and `useBinnedData`,
as taken by fdapace.

## Starting values and convergence

The loop converges when the class labels repeat an earlier iteration exactly
(`conv.prop = 0`). Under heavy tails a few boundary subjects can alternate
indefinitely while the partition is otherwise settled, so read `conv` together
with `iter`: `iter == maxIter` means the cap bound the fit. The loop also stops
when any class would fall below ten curves, since a covariance surface cannot
be estimated from fewer, and the K-means start (`init.method = "kmeans"`) can
propose such a class on heavy-tailed data. The recommended practice is several
random starts (`init.method = "random"`, or a partition passed through
`init.labels`), keeping the fit with the smallest total squared reconstruction
error over subjects; `predicted_yhat()` returns the reconstruction of one
subject under a fitted class, and `rre()` the relative reconstruction error of
a whole partition.

## Choosing the number of classes

`select_K()` evaluates the silhouette and the Calinski-Harabasz index on the
multivariate scores of one global fit (`MFPCA_k()` with a single class over all
subjects), with a Euclidean or a robust Mahalanobis metric. The scores are the
one coordinate system in which sparsely observed subjects can be compared, and
the criteria compare partitions of that fixed space. Neither criterion is
defined at `K = 1`, so they answer how many classes there are given at least
two. `silhouette_score()` is the curve-based alternative, which measures
distances to class mean functions and therefore favours separation in level.

## What is robust, and what is not

| Component | Estimator | Robust |
|---|---|---|
| Marginal covariance surfaces and eigenfunctions | Kendall's tau operator | yes |
| Second-stage rotation (multivariate eigenfunctions) | spatial-sign scatter | yes |
| Spread along each multivariate direction | interquartile range | yes |
| Mean functions | local linear least squares | no |
| Principal component scores | least squares | no |
| Subject-level pre-smoothing | Nadaraya-Watson | no |
| Class-assignment rule | squared Euclidean residuals | no |

MRKCFC is robust to non-Gaussian principal component scores, including scores
without finite variance. It is not protected against gross outliers in the
observed values, and the least-squares score step is ill-conditioned for a
subject whose number of visits is close to the number of retained components.

## Also included

The metrics, preprocessing and simulation machinery of the accompanying
manuscript ship here, so that the version producing the published numbers is
the version that installs. The scripts of the simulation study and of the
real-data analysis are in the companion repository
[MRKCFC_Simulation](https://github.com/BingbingFan/MRKCFC_Simulation).

- `rre()`, the relative reconstruction error of a partition, 0 for a perfect
  reconstruction and 1 for one no better than the mean of the observations.
- `robust_normalize()` and `back_transform()`, median and interquartile range
  rescaling of a variable and its inverse.
- `gen_scores()` and `vm_intermediate()`, six score distributions ordered by how
  many moments exist, each calibrated to a common interquartile range.
- `select_K()`, `ch_index()`, `silhouette_index()` and `whiten_scores()`,
  score-space cluster-validity criteria.
- `spatial_median()`, `sscm()`, `robust_stage2()`, `principal_angles()` and
  `stage2_diagnostic()`, the second-stage estimator and its diagnostics.

## Related

The univariate method, for one longitudinal variable, is the package
[RKCFC](https://github.com/BingbingFan/RKCFC).

## Attribution

The Kendall FPCA runtime is derived from KFPCA version 2.0 by Rou Zhong and
Jingxiao Zhang (<https://github.com/BingbingFan/KFPCA>), incorporated under
GPL (>= 3) and identified in the source. MRKCFC-specific adaptations are by
Bingbing Fan.

## References

Chiou, J.-M. and Li, P.-L. (2007). Functional clustering and identifying
substructures of longitudinal data. *Journal of the Royal Statistical Society,
Series B* **69**(4), 679-699.

Happ, C. and Greven, S. (2018). Multivariate functional principal component
analysis for data observed on different (dimensional) domains. *Journal of the
American Statistical Association* **113**(522), 649-659.

Zhong, R., Liu, S., Li, H. and Zhang, J. (2022). Robust functional principal
component analysis for non-Gaussian longitudinal data. *Journal of Multivariate
Analysis* **189**, 104864.

Visuri, S., Koivunen, V. and Oja, H. (2000). Sign and rank covariance matrices.
*Journal of Statistical Planning and Inference* **91**(2), 557-575.

Yao, F., Muller, H.-G. and Wang, J.-L. (2005). Functional data analysis for
sparse longitudinal data. *Journal of the American Statistical Association*
**100**(470), 577-590.

## License

GPL (>= 3).
