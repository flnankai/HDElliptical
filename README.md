# HDElliptical

`HDElliptical` is the companion R package for *High-Dimensional Data Analysis
for Elliptically Symmetric Distributions*. Version 0.1.0 implements Chapters
1--7, with Rcpp/RcppArmadillo kernels for the computationally intensive
parts and direct links between the book, the source papers, the public API,
and regression tests.

## Current status

The Chapter 1 foundation layer is implemented and tested:

- spatial signs, the modified-Weiszfeld spatial median, and empirical spatial
  ranks;
- spatial sign covariance, multivariate spatial Kendall, and spatial-rank
  covariance matrices;
- Tyler's trace-normalized shape estimator and the angular central Gaussian
  log-likelihood;
- the joint Hettmansperger-Randles location-shape estimator; and
- simulation from spherical and elliptically symmetric distributions.

The first Chapter 2 block is also implemented: one- and two-sample Hotelling
tests, fixed-dimensional spatial sign/rank tests, and the high-dimensional
Srivastava--Du, Park--Ayyala, Bai--Saranadasa, Chen--Qin, and
Srivastava--Katayama--Kano mean tests, plus the Cai--Liu--Xia
precision-adjusted maximum test, the Feng--Zou--Wang--Zhu scale-invariant
Behrens--Fisher test, the common-covariance Composite T2 test, and the
Wang--Peng--Li high-dimensional raw spatial-sign test and Feng--Sun
scalar-invariant one-sample sign test, together with the Feng--Liu--Ma
one-sample INST, the Feng--Zou--Wang and Li--Wang--Zou two-sample sign tests,
the Feng--Zhang--Liu high-dimensional two-sample spatial-rank test, and the
Huang--Liu--Zhou--Feng two-sample inverse-norm sign test (tINST), together
with the Xu--Lin--Wei--Pan analytical adaptive sum-of-powers test (aSPU), the
scaled-spatial-median spatial-sign MAX/MAXSUM procedures, the Yan--Zhao--Feng
general-\(m\) weighted MAX/MAXSUM procedures, and the Zhang--Feng adaptive
rank max/sum/Cauchy tests.
The late Chapter 2 block is now implemented as well: the Zhang--Zhou--Guo
one-sample normal-reference test with paired and same-unit linear-hypothesis
wrappers, the Zhang--Zhu--Zhang scale-invariant two-sample normal-reference
test, the Wang--Xu approximate-randomization test, the Feng--Wang PDQ
two-sample spatial-sign test, the Zhao--Feng strong-correlation one-sample
spatial-sign test, and the Feng--Zhou--Wang fixed-ridge and Cauchy-aggregated
elliptical regularized Hotelling tests.
The formula-complete generic weighted layer is also public through
`generic_weighted_hr_location()` and `oracle_weighted_sign_sum_test()`. The
oracle statistic requires an explicit `null_sd`, or supplied `nu2` and
`trace_R2`, before reporting a calibrated p-value; no feasible plug-in is
invented from incomplete prose.
The aSPU interface distinguishes the primary paper's raw mean-difference
powers from the book's scale-invariant studentized-coordinate variant. CLX exposes
distinct oracle, supplied feasible, and adaptive-threshold precision paths.
The Zhang--Feng sum and combined rank paths require either a supplied long-run
variance or an explicitly requested Ouyang--Parzen lag; the package does not
invent an automatic bandwidth absent from the source paper.
The compiled kernels expose the raw formula components and reject singular or
invalid variance estimates instead of silently applying a ridge or
pseudoinverse.

The only Chapter 2 item not implemented is the structured-correlation review
at `ch2_location.tex:2438-2451`. It is explicitly blocked pending source: a
legally accessible primary full text has not been obtained, and the book gives
no executable statistic or calibration formulas. No public API is invented
from that prose.

Chapter 3 is now implemented as 38 callable APIs. It covers classical and
high-dimensional covariance/sphericity tests, covariance thresholding and
POET, the convex-l1 EC2 covariance program, Gaussian off-diagonal graphical
lasso and Gaussian CLIME, elliptical sign/rank/adaptive sphericity, SSCM
proportionality, certificate-checked SCLIME/SGLASSO precision estimation,
tensor elliptical graphs, the Ollila shrinkage/BASIC/BASICS/linear-pooling line,
high-dimensional HR, and elliptical factor reconstruction/precision. Solver
functions return feasibility, positive-definiteness, KKT, and convergence
diagnostics and do not silently add a ridge, eigenvalue floor, pseudoinverse,
or nearest-PD repair.
The distinct SCIO and scaled-lasso precision reviews, and the adaptive/MC+
extensions of EC2, remain review-only because the short book passages do not
close their additional solver and tuning contracts.

Chapter 4 is now implemented as 33 callable APIs. It covers unconditional and
conditional alpha tests, mutual-fund FDR, classical/dense-max/spatial-sign/ERHT
change-point tests and ERHT-WBS segmentation, Gaussian and robust white-noise
tests, the radial-directional diagnostic, weighted-sign alpha completion,
book Gaussian-alpha and conditional-Wald benchmarks, and
Gaussian/panel/vector/exact rank-U independence tests. The separate
mutual-independence max--sum studentization review remains source-blocked; the
available sources do not uniquely determine its calibration, and the exact
rank-U vector method is not presented as a substitute. Randomization, wild
bootstrap, permutation, and Gaussian-process draws appear only where they are
intrinsic to a method's calibration; the package contains no paper-specific
simulation grids, size/power replications, tables, scripts, or data.

Chapter 5 is now implemented as 23 callable APIs on a common certified
binary-classifier fit and prediction contract. It covers oracle and classical
Gaussian/elliptical LDA and QDA, independence/FAIR/thresholded and direct sparse
linear rules, sparse Gaussian and spatial-sign QDA, certified precision
plug-ins, and robust GQDA/HR-GQDA. Fits record score orientation and scale,
priors, tuning, feasibility/KKT and determinant/SPD certificates; invalid fits
cannot be used for prediction and receive no hidden numerical repair.

Chapter 6 is now implemented as 16 callable APIs across four modules:
classical PCA/CCA/Bartlett testing with robust factor subspaces, RTS factors
and factor-number selection; spatial-sign, Kendall and generalized-sign PCA;
TPM/Fantope/PMD sparse PCA and CCA; and the primary metric-lasso `sscca()`
interface, kept distinct from the book's `sign_whitened_sparse_cca()`.
Centers, divisors, ranks, eigengaps, ties, supports, l1 boundaries, ridge
operators, KKT certificates and failure states are explicit; there is no
hidden floor, pseudoinverse or numerical repair.

Chapter 7 is now implemented as 14 exported APIs backed by 19 registered
compiled kernels. It covers Lloyd K-means, full-covariance Gaussian-mixture
EM, sparse K-means and K-medians, CHIME, IF-PCA, K-spatial medians, SM-SSCM,
Sparse-SM and its tau/K selectors, plus `semc_fit()`, its prediction method,
and `semc_select_k_gap()`. Initialization, ties, empty clusters, cycles, seeds,
solver certificates and returned final states are explicit.
CHIME canonicalizes component labels so `omega <= 1/2` after initialization
and every update, and the spatial kernels preserve finite close-large,
subnormal and opposite-extreme geometry.

SEMC is no longer source-blocked. It is an independent rewrite against the
primary source and the pinned
`flnankai/GEMcluster@10fce04fe690fe274dd5d237cfcd3d5c6a4139f6` MIT release.
Primary-paper and official-software-only controls are distinguished explicitly;
formula and pinned-software fixtures are tested without copying author-code
expressions or reproducing a paper simulation.

The fresh Chapter 1--7 gate contains 176 exports, 3 S3 registrations, and 140
Rcpp attributes matched one-for-one by 140 native registrations. Roxygen
generated 223 Rd files covering 169 public topics. Across 56 test files there
are 688 named `test_that` blocks and 5,781 static `expect_*` calls; the isolated
installed-package run reported 6,267 passes and zero failures, errors, warnings,
or skips. Two independent generation passes produced the same 226-file hash
manifest, and all 226 generated files matched the synchronized release tree.
The final 0.1.0 source tarball built from this release state completed both
release checks. Standard `R CMD check --no-manual` reported 0 ERROR, 0 WARNING,
and 1 NOTE; remote-enabled `R CMD check --as-cran --no-manual` reported
0 ERROR, 0 WARNING, and 2 NOTEs. Both ran examples, tests, and vignettes; the
as-CRAN run also executed `donttest` examples. Both checks share the conservative
Windows DLL scan for linked `_exit`, `abort`, and `exit` symbols; an exact
direct-call scan of the R and `src/` sources found zero calls. The additional
as-CRAN NOTE identifies this as a new CRAN submission.

The public repository and issue tracker were bootstrapped on 2026-08-16 and
now resolve at `https://github.com/flnankai/HDElliptical`. The final as-CRAN
run kept incoming remote checks enabled and completed without errors or
warnings. The annotated `v0.1.0` tag and GitHub release were created only after
the Windows, macOS, and Linux Actions jobs passed and the post-bootstrap
tarball, check logs, and SHA-256 manifest were regenerated.

The normative method ledger has 87 rows: 83 implemented, 2 review-only, and 2
source-blocked. The two blocked rows are the Chapter 2 structured-correlation
review and Chapter 4 mutual-independence studentization; SEMC is implemented.

The dated [method coverage](development/METHOD_COVERAGE.md),
[machine-readable traceability](development/METHOD_TRACEABILITY.csv),
[book errata](development/BOOK_ERRATA.md), and
[validation snapshot](development/VALIDATION_20260814.md) keep the book
formula, primary-paper correction, public API, compiled kernel, and test
evidence traceable for release and manuscript review.

The installed-package [deterministic benchmark harness](inst/benchmarks/README.md)
runs one representative callable workflow for each of Chapters 1--7. Its
[machine-local report](output/benchmarks/BENCHMARK-REPORT.md) and
[CSV output](output/benchmarks/HDElliptical-benchmarks.csv) record timings and
deterministic fingerprints. They are software benchmarks, not paper-specific
simulation, size/power, or empirical-reproduction studies.

## Installation

```r
# Current GitHub version:
remotes::install_github("flnankai/HDElliptical")
```

For a local checkout:

```r
install.packages(".", repos = NULL, type = "source")
```

## Quick start

```r
library(HDElliptical)

set.seed(20260814)
x <- relliptical(
  n = 200,
  location = c(0, 0, 0),
  shape = matrix(c(2, .5, 0, .5, 1, .2, 0, .2, .7), 3),
  radial = function(n) abs(rt(n, df = 3))
)

spatial_median(x)
sscm(x)
spatial_kendall(x)
tyler_shape(x)
hr_estimator(x)

# Chapter 2 examples
srivastava_du_one_sample_test(x)
park_ayyala_one_sample_test(x)

y <- relliptical(
  n = 200,
  location = c(0.1, 0, 0),
  shape = matrix(c(2, .5, 0, .5, 1, .2, 0, .2, .7), 3),
  radial = function(n) abs(rt(n, df = 3))
)
paired_y <- x + matrix(rnorm(length(x), sd = 0.2), nrow(x))
cai_liu_xia_two_sample_test(x, y)
feng_zou_wang_zhu_two_sample_test(x[1:10, ], y[1:11, ])
composite_t2_two_sample_test(x[1:8, 1:3], y[1:8, 1:3])
wang_peng_li_one_sample_test(x)
feng_sun_one_sample_test(x[1:10, ], tol = 1e-6)
inst_one_sample_test(x[1:10, ], tol = 1e-6)
feng_zou_wang_two_sample_sign_test(x[1:10, ], y[1:11, ], tol = 1e-6)
li_wang_zou_two_sample_sign_test(x[1:10, ], y[1:11, ], tol = 1e-6)
feng_zhang_liu_spatial_rank_test(x[1:10, ], y[1:11, ], tol = 1e-6)
tinst_two_sample_test(x, y)
xu_lin_wei_pan_aspu_test(x, y)
scaled_spatial_median(x)
spatial_sign_max_test(x)
spatial_sign_maxsum_test(x[1:10, ], tol = 1e-6)
weighted_scaled_spatial_median(x, m = 0)
yan_zhao_feng_weighted_max_test(x, m = -1)
yan_zhao_feng_weighted_maxsum_test(x[1:10, ], m = -1, tol = 1e-6)
zhang_feng_rank_one_sample_test(x, component = "max")
zhang_feng_rank_two_sample_test(x, y, component = "max")
zhang_zhou_guo_one_sample_test(x)
zhang_zhou_guo_paired_test(x, paired_y)
zhang_zhou_guo_linear_hypothesis_test(x, rbind(c(1, -1, 0)))
zhang_zhu_zhang_two_sample_test(x, y)
wang_xu_approx_randomization_test(
  x[1:8, ], y[1:8, ], calibration = "exact"
)
feng_wang_pdq_two_sample_test(
  x[1:10, ], y[1:11, ], B = 19, seed = 2605
)
zhao_feng_strongcorr_sign_test(x[1:12, ], B = 19, seed = 2601)
elliptical_regularized_hotelling_test(x, rho = 0.5)
elliptical_regularized_hotelling_cauchy_test(
  x, rho = c(0.2, 0.5, 1)
)

# Chapter 3 examples
john_sphericity_test(x)
zou_peng_feng_wang_sphericity_test(x)
ollila_raninen_shrinkage_covariance(x, sphericity = "ell3")
elliptical_factor_number(x, max_factors = 1)
spatial_sign_poet(x, factors = 1, threshold = 0.1)

# Chapter 5 examples
training <- rbind(x, y)
group <- factor(rep(c("class1", "class2"), each = nrow(x)))
lda_fit <- classical_lda_classifier(training, group)
predict(lda_fit, training[1:4, ], type = "class")
```

Every function treats rows as observations and columns as variables. Exact
zero directions follow the convention `U(0) = 0`; estimators expose
convergence and zero-residual diagnostics instead of silently discarding
observations.

## Reproducibility policy

Each method is accepted into the package only after it has:

1. a source location and bibliographic key;
2. a documented callable API;
3. numerical or invariance tests;
4. a compiled implementation when it materially reduces runtime; and
5. a check against the original paper whenever the book omits an implementable
   quantity or contains a suspected typographical error.

The package is under active development and is not yet ready for scientific
production use.

The package targets reusable method implementations. It does not ship the
paper-specific Monte Carlo grids, simulation tables, or power-study scripts
from the many source articles; compact deterministic and fixed-seed checks are
kept only to verify formulas, invariances, numerical stability, and API
contracts.

Randomization and wild-bootstrap draws inside
`wang_xu_approx_randomization_test()`, `feng_wang_pdq_two_sample_test()`, and
`zhao_feng_strongcorr_sign_test()` are intrinsic parts of those tests' null
calibrations. They are callable-method computations, not replications of the
papers' simulation studies. No paper size/power simulation scenarios or
result tables are included.

Chapter 7 selector permutations and IF-PCA empirical-null generation are also
intrinsic method calibrations. They use supplied plans or explicit local seeds
and restore the caller's RNG state; they are not paper simulation replications.
