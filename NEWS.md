# HDElliptical 0.1.1

- Replaced README links to build-excluded audit and benchmark files with
  public GitHub URLs so they remain resolvable from the CRAN source package.
- Added an author-year citation and public source URL to the package
  description for CRAN submission metadata.

# HDElliptical 0.1.0

## Chapter 1 foundations

- Added Rcpp implementations of spatial signs, the modified-Weiszfeld spatial
  median, empirical spatial ranks, the multivariate spatial Kendall matrix,
  Tyler's fixed-point estimator, and the Hettmansperger-Randles iteration.
- Added R interfaces for SSCM and spatial-rank covariance matrices, ACG
  likelihood evaluation, shape normalization, and elliptical simulation.
- Added tests of defining equations and translation, scale, orthogonal, and
  affine equivariance properties.
- Hardened directional and shape calculations against overflow and underflow,
  including finite inputs near the limits of double precision.
- Changed iterative convergence certificates to use the defining estimating
  equation residuals, preventing false convergence near an observation.
- Defined the one-dimensional HR interface as the ordinary sample median with
  unit shape and an explicit not-applicable joint-shape diagnostic.

## Chapter 2 location tests

- Added exact one- and two-sample Hotelling tests and fixed-dimensional
  spatial sign, signed-rank, and pooled spatial-rank tests.
- Added Rcpp implementations of the Srivastava--Du, Park--Ayyala,
  Bai--Saranadasa, Chen--Qin, and Srivastava--Katayama--Kano
  high-dimensional mean tests.
- Added the Cai--Liu--Xia precision-adjusted maximum test with deliberately
  separate oracle, supplied-feasible, and adaptive-threshold precision paths.
  The feasible path uses transformed within-group empirical variances rather
  than the diagonal of an estimated precision matrix.
- Added the Wang--Peng--Li one-sample high-dimensional spatial-sign test with
  its literal leave-two-out cross-validation variance estimator. Exact zero
  signs follow `U(0) = 0`; invalid variance estimates are not repaired.
- Added the Feng--Sun scalar-invariant one-sample spatial-sign test with
  literal pair-specific leave-two-out location/diagonal fits and the primary
  paper's feasible trace calibration omitted from the short book presentation.
- Added the Huang--Liu--Zhou--Feng two-sample inverse-norm sign test (tINST)
  with the paper's feasible cross statistic, observation-specific leave-one-
  out nuisance fits, published variance decomposition, and explicit iteration
  stability versus score-residual diagnostics.
- Added the Feng--Zou--Wang--Zhu scale-invariant high-dimensional
  Behrens--Fisher test with the primary paper's asymptotic centering and exact
  leave-four-out / two-plus-two leaveout trace estimators.
- Added the common-covariance Composite T2 test with paper-greedy correlation
  blocks, a full two-plus-two leaveout statistic, and the original first-group
  leave-four-out trace calibration.
- Added the Feng--Liu--Ma one-sample INST procedure with null-centered
  endpoints, delete-two nuisance fits, and the supplement's direct
  `2 * n^(-4)` feasible variance.
- Added the Feng--Zou--Wang two-sample multivariate-sign test with crossed
  leave-one-out locations, full/leave-one-out diagonal HR fits, and the
  primary paper's Proposition 2 feasible trace and variance calibration.
- Added the Li--Wang--Zou simpler two-sample spatial-sign test with two
  full-sample diagonal HR fits, the exact feasible bias correction, and the
  published distinction between ordered-pair trace denominators and the
  outer squared-sample-size variance factors.
- Added the Feng--Zhang--Liu high-dimensional two-sample spatial-rank test
  with its feasible leave-four/leave-two trace calibration. The API exposes
  both a coherent geometric scale identification and the article's literal
  trace-normalized finite-sample definition.
- Added the Xu--Lin--Wei--Pan analytical aSPU test. Its default finite-power
  scores reproduce the primary paper's raw sample-mean differences, while an
  explicit `book_studentized` path implements the book's scale-invariant
  standardized-coordinate variant; both use deterministic odd/even Gaussian
  and maximum extreme-value calibration without paper-specific simulations.
- Added the scaled spatial median, the spatial-sign MAX test, and its MAXSUM
  Cauchy combination with the existing Feng--Sun SUM component. The API
  exposes fixed-point stability and estimating-equation residuals separately.
- Added the Yan--Zhao--Feng weighted scaled spatial median and general-`m`
  weighted MAX/MAXSUM tests. The `m = -1` SUM path is locked component by
  component to the existing INST implementation, and the published half-power,
  endpoint-index, and radial-moment typos are corrected explicitly.
- Added Zhang--Feng one- and two-sample adaptive rank tests with max, squared-
  rank sum, and equal-weight Cauchy components. Sum calibration requires a
  supplied positive long-run variance or an explicitly selected
  Ouyang--Parzen lag; there is no hidden bandwidth or variance floor.
- Added `generic_weighted_hr_location()` and the deliberately oracle
  `oracle_weighted_sign_sum_test()` for the formula-complete generic weighted
  layer. Calibration is reported only from an explicit `null_sd`, or supplied
  `nu2` and `trace_R2`; no feasible p-value is guessed.
- Added `zhang_zhou_guo_one_sample_test()`,
  `zhang_zhou_guo_paired_test()`, and
  `zhang_zhou_guo_linear_hypothesis_test()` for normal-reference one-sample
  means, rowwise paired differences, and same-unit linear hypotheses. All three
  wrappers reuse the compiled centred-U-statistic and unbiased-trace kernel,
  and deterministic tests lock the appendix `v + 4` denominator, wrapper
  transformations, primal/dual agreement, invariances, and no-repair errors.
- Added `zhang_zhu_zhang_two_sample_test()`, the Zhang--Zhu--Zhang
  scale-invariant two-sample normal-reference test, with crossed covariance
  weights, bias-corrected trace calibration, and
  distinct published-threshold and unadjusted degrees-of-freedom paths.
  Deterministic formula, primal/dual, invariance, extreme-tail, and boundary
  checks cover the compiled implementation.
- Added `wang_xu_approx_randomization_test()` with the Wang--Xu full-sample CQ
  statistic and adjacent within-group half-difference reference sample. The
  compiled kernel supports exhaustive small-sample sign enumeration and a
  reproducible counter-based Monte Carlo mode; deterministic tests lock odd-
  row discard, tie handling, exact tails without +1, Monte Carlo tails with
  +1, seed isolation, and no-repair boundaries.
- Added `feng_wang_pdq_two_sample_test()`, the Feng--Wang PDQ two-sample
  spatial-sign test, with empirical
  pairwise-difference U-quantiles, standardized spatial medians, feasible
  bias/variance matrices, and Rademacher wild-bootstrap calibration. Literal
  small-sample and fixed-counter tests verify the compiled statistic,
  bootstrap, dual finite-`B` decision conventions, invariances, zero
  directions, extreme units, and degeneracy contracts.
- Added `zhao_feng_strongcorr_sign_test()`, the Zhao--Feng strong-correlation
  one-sample spatial-sign test, with Rademacher or Gaussian wild-bootstrap
  multipliers. The compiled pair-sum
  kernel and deterministic tests preserve the distinct null and fitted
  centerings, multiplier-dependent diagonal, common-`tau` cancellation,
  finite-`B` decisions, RNG contract, invariances, and zero/no-repair cases.
- Added `elliptical_regularized_hotelling_test()` and
  `elliptical_regularized_hotelling_cauchy_test()` for fixed-ridge and
  Cauchy-grid Feng--Zhou--Wang elliptical regularized Hotelling testing. One
  economy-SVD compiled kernel evaluates the feasible
  companion-matrix centre and variance for every ridge; deterministic tests
  cover literal margins, the `p > n` complement route, small positive ridges,
  stable signed-log Cauchy aggregation, invariances, extreme scales, and
  explicit failures. The paper's unspecified simulation-program Bartlett
  modification is not guessed.
- Marked `ch2_location.tex:2438-2451` structured-correlation combinations as
  the sole blocked Chapter 2 item. A legally accessible primary full text has
  not been obtained and the book contains no executable formulas, so no API
  is fabricated pending source access.
- Retained randomization and wild bootstrap only where they are intrinsic
  test calibrations. These callable calibrations are not paper simulation
  replications; no paper-specific Monte Carlo grids, size/power scenarios,
  tables, reproduction scripts, or data were added.
- Stabilized CLX group centering and adaptive covariance moments under large
  common translations and representable extreme global scales.
- Kept the original common-covariance Bai--Saranadasa statistic separate from
  the unequal-covariance Chen--Qin U-statistic, correcting their conflation in
  the current book draft.
- Added stable primal/dual Gram calculations, leave-out kernels, scaled
  diagnostic components, and independent formula/invariance/degeneracy tests.

## Chapter 3 matrix estimation and testing

- Added six classical covariance, sphericity, and spatial sign/rank shape
  tests, with the primary Mauchly, John, and Nagao finite-sample factors.
- Added Bickel--Levina, Rothman--Levina--Zhu, and Cai--Liu covariance
  thresholding plus supplied-factor POET, together with five
  high-dimensional Gaussian/light-tail covariance tests.
- Added four elliptical sign/rank/adaptive sphericity tests and two
  SSCM/spatial-rank proportionality tests with their feasible bias and
  variance calibration.
- Added certificate-checked spatial-sign SCLIME and SGLASSO precision
  estimation and exact support thresholding. Failed feasibility, KKT, or
  positive-definiteness checks are returned as failures rather than repaired.
- Added `ec2_covariance()` for the formula-complete convex-l1 EC2 branch,
  `gaussian_graphical_lasso()` with an off-diagonal l1 penalty, and
  `clime_precision()` with primary smaller-absolute-value symmetrization. Each
  reports strict feasibility/KKT or primal-dual certificates without hidden
  ridge, floor, pseudoinverse, or SPD repair. Adaptive/MC+ EC2, SCIO, and
  scaled-lasso remain review-only because their distinct contracts are not
  closed by the short review passages.
- Added tensor spatial-sign precision estimation and graph thresholding with
  column-major Kolda unfolding, the published pilot switch, and modewise
  objective/KKT diagnostics.
- Added the Ollila--Raninen shrinkage covariance estimator, regularized SSCM,
  BASIC, BASICS, and linear covariance pooling with self-contained quadrature,
  bounded inversion, and pooling optimisation.
- Added the primary joint high-dimensional HR fixed point, which bands the
  standardised-sign SSCM within every iteration and reports map, score,
  convergence, conditioning, and pilot diagnostics.
- Added elliptical factor-number selection, spatial-sign POET, one-step Tyler
  POET, and sparse factor precision reconstruction with explicit factor count
  and tuning constants.
- At the earlier Chapter 1--3 gate, generation produced 107 Rd files.
  That source package passed all examples, 35 test files, both
  vignettes, and `R CMD check --no-manual` with zero errors and warnings; the
  only NOTE is the existing Windows linked-symbol scan.
- Added no paper-specific simulation grids, size/power scenarios, result
  tables, reproduction scripts, or simulation data. Tests are deterministic
  formula, invariance, solver-certificate, and numerical-boundary checks.

## Chapter 4 other high-dimensional tests

- Added 33 callable APIs spanning unconditional and conditional alpha tests,
  mutual-fund FDR, change-point testing and ERHT-WBS segmentation, white-noise
  tests, a radial-directional diagnostic, and Gaussian/panel/vector
  independence testing.
- Completed the formula-closed Chapter 4 layer with
  `weighted_spatial_sign_alpha_oracle_test()`,
  `zhao_chen_zi_inst_alpha_test()`, `book_gaussian_alpha_cauchy_test()`,
  `conditional_factor_wald_test()`, and
  `wang_liu_feng_vector_u_independence_test()`. The separate mutual-independence
  max--sum studentization remains source-blocked rather than being replaced by
  the adjacent exact rank-U vector method.
- Corrected the book-level DMS, spatial-sign change-point, alpha-test,
  conditional-residual, white-noise variance, panel-Gumbel, and serial-panel
  formula conflations against the primary sources; every software choice is
  exposed in returned diagnostics.
- Added strict nuisance-fit, covariance, determinant, variance, solver, and
  seed contracts. Invalid inputs or uncertified fits fail explicitly rather
  than receiving an implicit ridge, floor, pseudoinverse, or nearest-PD repair.
- Retained bootstrap, randomization, permutation, and Gaussian-process draws
  only where they define a method's own calibration. No paper simulation
  design, size/power grid, empirical table, reproduction script, or data set
  was added.
- At the earlier Chapter 1--4 gate, generation produced 113 exports and 147
  Rd files. That source tarball passed 41 test files, 426 named test blocks, all examples, both vignettes,
  and `R CMD check --no-manual` with zero errors and warnings; the only NOTE is
  the existing Windows linked-symbol scan.

## Chapter 5 classification

- Added a common `hd_classifier_fit` contract and 23 callable APIs covering
  oracle/classical LDA and QDA, independence and FAIR rules, thresholded and
  direct sparse LDA, sparse Gaussian/spatial-sign QDA, certified plug-in rules,
  and robust GQDA/HR-GQDA.
- Added strict vector- and matrix-Dantzig, lasso, and interaction-loss solver
  certificates. Tuning paths, priors, score orientation/scale, determinant
  signs, SPD status, feasibility, duality, and KKT residuals are retained in
  every fit; failed fits are not silently repaired or made predictable.
- Corrected the book's unequal-prior elliptical shortcut, DSDA coding/objective,
  Jiang versus Cai--Zhang attribution, explicit ridge-operator semantics,
  SSQDA trace computation, and HR-QDA/GQDA distinction against primary sources.
- Added deterministic formula, prediction, invariance, tuning-isolation,
  solver-certificate, and degeneracy tests. At the earlier Chapter 1--5 gate,
  the source tarball passed 45 test files, 495 named test blocks, 4,772
  expectations, all 170 Rd files, examples and both vignettes with zero errors
  and warnings; the only
  NOTE is the existing Windows linked-symbol scan.
- No paper-specific classification simulation grid, error-rate table,
  reproduction script, benchmark data, or guessed optimizer was added.

## Chapter 6 dimension reduction

- Added 16 callable APIs across four modules: classical PCA/CCA/Bartlett
  testing and robust factor/subspace/number methods; spatial-sign, Kendall and
  generalized-sign PCA; TPM/Fantope/PMD sparse PCA and CCA; and the primary
  metric-lasso `sscca()` kept distinct from the book's whitened-PMD variant.
- Made centers, divisors, ranks, eigengaps, ties, supports, l1 boundaries and
  ridge operators explicit. KKT, PSD and convergence certificates are
  returned; invalid fits are never promoted or silently repaired.
- Separated generalized-sign median/raw-MAD and h-order cutoff families and
  recorded the radial population scope. Preserved the primary SSCCA metric
  diagonals and BIC program rather than silently whitening them.
- Added 77 named deterministic test blocks with 539 static source `expect_*`
  call occurrences. They verify formulas, invariances, solvers and boundaries;
  no paper Monte Carlo design, tuning result, table, script or data was added.

## Chapter 7 clustering

- Added Lloyd K-means, full-covariance Gaussian-mixture EM, sparse K-means and
  its selector, sparse K-medians, CHIME, IF-PCA, K-spatial medians, SM-SSCM,
  Sparse-SM and its tau/K selectors, plus SEMC fitting/prediction and Gap-LSE
  selection: 14 exported APIs and 19 registered compiled kernels.
- Canonicalized CHIME labels after initialization and every update so
  `omega <= 1/2`, with coherent mean, responsibility, beta, score, class and
  history orientation. The full supplied lambda path and exact lasso KKT
  certificates are preserved without a hidden ridge.
- Verified ordered-pair BCSS as `2 * (TSS - WCSS)`, the s=1 soft-threshold
  boundary, full-covariance GMM rank/SPD checks, and explicit tie, start,
  farthest-empty, cycle and returned-final-state contracts.
- Hardened robust spatial geometry with raw-first close-large differences,
  overflow-safe scaled accumulation, a Cholesky metric root, and finite
  subnormal/opposite-DBL_MAX paths.
- Kept selector permutations and IF-PCA null generation local, fixed-seed and
  caller-RNG preserving. These are intrinsic calibrations, not paper
  simulations.
- Added `semc_fit()`, `predict.semc_fit()`, and `semc_select_k_gap()` as an
  independent rewrite against the primary source and pinned
  `flnankai/GEMcluster@10fce04fe690fe274dd5d237cfcd3d5c6a4139f6` MIT release.
  Primary-paper and software-only controls remain explicit and separate; fixed
  formula/software fixtures contain no copied author-code expression or paper
  simulation.
- The fresh unified gate has 176 exports, 3 S3 registrations, 140 Rcpp
  attributes matched by 140 native registrations, 223 Rd files, and 169 public
  topics. Its 56 test files contain 688 `test_that` blocks and 5,781 static
  `expect_*` calls; the installed run reported 6,267 passes and zero failures,
  errors, warnings, or skips. Two generation passes produced an identical
  226-file hash manifest and all 226 files matched the release tree.
- The final 87-row ledger records 83 implemented, 2 review-only, and 2
  source-blocked method rows. The blocked rows are the Chapter 2 structured-
  correlation review and Chapter 4 mutual-independence studentization; SEMC is
  implemented.
- The final 0.1.0 source tarball built from this release state passed standard
  `R CMD check --no-manual` with 0 ERROR, 0 WARNING, and 1 NOTE, and
  remote-enabled `R CMD check --as-cran --no-manual` with 0 ERROR, 0 WARNING,
  and 2 NOTEs. Both ran examples, tests, and vignettes, and the as-CRAN check
  also ran `donttest`. Both checks share the conservative DLL linked-symbol
  scan for `_exit`, `abort`, and `exit`; exact R/`src` direct-call scanning
  found zero. The additional as-CRAN NOTE identifies a new CRAN submission.
- The public `flnankai/HDElliptical` repository and issue tracker were
  bootstrapped before the final tag. The post-bootstrap as-CRAN check kept
  incoming remote checks enabled, and the Windows, macOS, and Linux Actions
  jobs passed before the annotated `v0.1.0` tag, release, and regenerated
  release assets were published.
- Added a deterministic installed-package benchmark harness and checked-in CSV/
  Markdown report with one representative workflow per chapter. These are
  machine-local software timings and fingerprints, not paper simulations,
  size/power experiments, or empirical reproductions.
