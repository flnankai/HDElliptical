# Validation snapshot: 2026-08-14

This record freezes the verification evidence for the package state built from
the 2026-08-14 book snapshot. It is not a substitute for continuous
integration; it documents the independent formula checks and numerical stress
tests that must remain green as later chapters are added.

## Build environment

- R 4.5.2 on Windows 11 x64 (UCRT).
- GCC 14.3.0 and C++17 through Rtools.
- `Rcpp` and `RcppArmadillo` compiled from a clean ASCII staging directory so
  that the Unicode parent path could not affect tool discovery.
- The current fresh source compiled, installed, and loaded from an isolated
  ASCII staging directory. The package source retains both HTML vignettes for
  the final release check.
- `Rcpp::compileAttributes()` found 140 attributes, and the generated registry
  contains the same 140 native entries. The namespace contains 176 exports and
  3 S3 registrations.
- Roxygen generated 223 Rd files covering 169 public topics; all Rd files
  parsed strictly.
- The current suite contains 56 test files, 688 named `test_that` blocks, and
  5,781 static `expect_*` calls. The installed-package reporter recorded 6,267
  passes with zero failures, errors, warnings, or skips.
- Two independent generation passes produced the same 226-file hash manifest,
  and all 226 generated files matched the synchronized release tree.
- The JSS manuscript executed its embedded R examples, resolved its BibTeX
  references with the official JSS style, and produced the visually inspected
  22-page release PDF recorded by the JSS build materials.
- The final 0.1.0 source tarball built from this release state completed
  standard `R CMD check --no-manual` with 0 ERROR, 0 WARNING, and 1 NOTE, and
  remote-enabled `R CMD check --as-cran --no-manual` with 0 ERROR, 0 WARNING,
  and 2 NOTEs. Both checks ran examples, tests, and vignettes, and the as-CRAN
  run also ran `donttest` examples. Both share the conservative Windows DLL
  scan for linked `_exit`, `abort`, and `exit`; exact scanning found zero direct
  calls in the R and `src/` sources. The additional as-CRAN NOTE identifies a
  new CRAN submission.
- The post-bootstrap as-CRAN run kept incoming remote checks enabled. The public
  GitHub repository and issue URLs resolved, and all incoming checks completed
  before the annotated `v0.1.0` tag and release were created. A console-only
  Bioconductor index timeout did not persist: dependency checking completed
  successfully.

The machine's startup locale requested `C.UTF-8`, which Windows R could not
set. Checks were therefore run with `LC_ALL=C`; the earlier locale messages
were environmental startup output rather than package failures.

## Chapter 1 independent QA

- Spatial signs returned the correct unit direction for both
  `c(1.7e308, 1.7e308)` and two minimum positive subnormal values.
- Modified Weiszfeld regression cases covered near-observation false
  convergence, subnormal separation, weights near the double limit, and global
  scales from `1e-200` to `1e200`.
- Pairwise ranks and Kendall matrices covered subnormal differences, opposite
  finite extremes, and cancellation between nearby values of order `1e308`.
- Tyler shape passed general nonsymmetric affine transformations, scales from
  `1e-200` to `1e200`, and an external `ICSNP::tyler.shape` comparison.
- The Hettmansperger--Randles estimator passed general and condition-number
  stress affine transformations, scale tests, default-initializer collision,
  one-dimensional conventions, invalid controls, and an external
  `SpatialNP` comparison. When an over-strict tolerance falls below its
  floating-point accuracy floor, it correctly returns `converged = FALSE`.
- ACG log likelihood agreed with independent log-sum-exp and ordinary direct
  formula calculations, including a shape condition number near `1e300`.
- Elliptical generators covered zero positive-semidefinite shape, invalid
  fractional dimensions, and deterministic radial behavior.

Known non-blocking limit: a shape such as `diag(c(1e308, 1e-308))` spans a
condition number near `1e616`, beyond the usable dynamic range of double
precision after normalization. The package rejects or loses that tiny
direction; the documented and tested range extends through condition numbers
near `1e300`.

## Chapter 2 formula and regression evidence

These are software formula and regression checks, not replications of paper
Monte Carlo designs, tables, data sets, size studies, or power studies. At the
earlier Chapter 1--2 gate, the tree had 26 test suites/files, 232 named test
blocks, and 2,536 passing
expectations: Chapter 1 contributes 1 file, 12 blocks, and 104 expectations;
Chapter 2 contributes 25 files, 220 blocks, and 2,432 expectations. The unified
run had zero failures, errors, warnings, or skips:

| Module | Focused expectations | Independent reference |
|---|---:|---|
| Hotelling one/two sample | 41 | direct covariance solve, univariate t identities |
| Fixed-\(p\) sign/signed-rank/pooled-rank | 57 | literal sign and ordered-pair R loops |
| Srivastava--Du / Bai--Saranadasa | 109 | direct matrix formula plus 100 randomized fixtures |
| Park--Ayyala / Chen--Qin | 53 | literal leave-out loops plus 50 randomized fixtures per method |
| Srivastava--Katayama--Kano | 69 | direct matrix formula plus 100 randomized fixtures |
| Cai--Liu--Xia | 139 | oracle, feasible, and adaptive paper equations, including stable large-translation and extreme-scale references |
| Common-covariance Composite \(T^2\) | 74 | literal two-plus-two leaveout statistic, greedy blocks, and first-group leave-four trace |
| Xu--Lin--Wei--Pan aSPU | 129 | primary raw and book-studentized scores, Isserlis moments, Miwa/Gumbel calibration, and PSD contracts |
| Scaled spatial-sign MAX/MAXSUM | 113 | literal joint location/diagonal equations, complete MAX factors, Gumbel inversion, Feng--Sun reuse, and stable Cauchy endpoints |
| Yan--Zhao--Feng weighted MAX/MAXSUM | 145 | general-\(m\) equations and direct feasible variance; \(m=-1\) locks every endpoint, leaveout fit, kernel, trace, variance, and statistic to INST |
| Zhang--Feng adaptive rank tests | 108 | exact signed-rank/WMW moments, max/sum/Cauchy components, supplied and explicit Ouyang--Parzen long-run variance paths, tie/order/model contracts |
| Wang--Peng--Li | 76 | literal equation (7) leave-two-out loops and equation (8) agreement on its valid domain |
| Feng--Sun | 80 | literal pair-specific leave-two-out numerator and feasible ordered-pair trace/variance |
| Feng--Liu--Ma INST | 101 | literal null-centered endpoints, leave-two nuisance fits, and direct supplement variance |
| Feng--Zou--Wang sign test | 103 | literal feasible crossed-location statistic and Proposition 2 trace calibration |
| Feng--Zou--Wang--Zhu Behrens--Fisher | 99 | literal leave-four / two-plus-two trace formulae and corrected bias factors |
| Li--Wang--Zou SST | 127 | full-sample crossed locations, bias, ordered traces, and squared-sample-size variance factors |
| Feng--Zhang--Liu spatial rank | 91 | both scale identifications, feasible trace calibration, and corrected variance factor |
| Huang--Liu--Zhou--Feng tINST | 80 | literal feasible cross statistic, leave-one-out nuisance fits, and published variance decomposition |
| Zhang--Zhou--Guo normal-reference one-sample/paired/linear-hypothesis tests | 121 | centred U-statistic rather than motivating oracle, exact Gaussian/Wishart trace corrections, three-cumulant chi-square matching, and primal/dual identities |
| Zhang--Zhu--Zhang normal-reference scale-invariant test | 102 | crossed covariance weights, Eqs. (24)--(25), both paper-adjusted and unadjusted degrees of freedom, and scaled-chi-square right tail |
| Wang--Xu approximate randomization | 101 | full-sample observed CQ, adjacent half-differences, exhaustive sign references, plus-one Monte Carlo convention, and discarded odd rows |
| Feng--Wang PDQ spatial-sign test | 115 | exact coordinatewise U-quantiles, crossed spatial medians, plug-in matrices, centred Rademacher bootstrap, and separate finite-\(B\) decisions |
| Zhao--Feng strong-correlation spatial-sign test | 109 | null-centred observed signs, fitted bootstrap signs, raw pair sums, Rademacher/Gaussian multipliers, zero-sign identity, and both finite-\(B\) decisions |
| Feng--Zhou--Wang ERHT/ERHT--CC | 90 | direct feasible companion-matrix center and variance, full-matrix/SVD references, stable ridge quadratic and Cauchy aggregation, with no guessed Bartlett mapping |

The final Chapter 2 source also adds `generic_weighted_hr_location()` and the
oracle `oracle_weighted_sign_sum_test()`. Their tests reconstruct the arbitrary
weight callback, weighted location equation, diagonal HR update, and quadratic
score. Calibration is accepted only when `null_sd`, or `nu2` with `trace_R2`,
is supplied; a feasible p-value is not manufactured from the review prose.

There is no paper-simulation replication in these counts. Random fixtures are
deterministic formula QA, and stochastic checks use recorded fixed seeds. The
only bootstrap or randomization exercised is calibration intrinsic to the
published method itself (including Wang--Xu, PDQ, and Zhao--Feng); those tests
verify exact small-space enumeration or fixed-seed algorithmic formulas rather
than reproducing paper scenarios, tables, empirical sizes, or powers.

The randomized maximum observed absolute or scaled discrepancies were about
`2.84e-13` for Srivastava--Du, `5.11e-15` for Bai--Saranadasa,
`8.4e-15` for Park--Ayyala, `7.8e-16` for Chen--Qin, and `1.07e-14`
for Srivastava--Katayama--Kano. CLX tests separately reproduce the pooled
covariance, entrywise variance, thresholds, eigenvalue adjustment, oracle and
feasible denominators, extreme-value CDF/quantile inversion, invariances,
large-translation stability, representable scales through `1e150` and
`1e-150`, and failure contracts. WPL and tINST add literal leave-out formula
checks and explicit no-repair/nonconvergence contracts. The two new max--sum
families add log-tail Cauchy tests through `log(p) = -1000`, exact endpoint
rules, full radial-moment factors, and explicit update-stability versus
estimating-equation diagnostics. The adaptive rank tests use deterministic
small-sample score references and never infer an unreported automatic HAC
bandwidth.

Extreme-scale CQ tests retain finite internal `*.scaled` components at
`1e150` and `1e-150`; raw fourth-order components may honestly overflow or
underflow, while the reported statistic remains exactly reconstructible from
the scaled numerator and variance.

## Chapter 3 formula and regression evidence

At the earlier Chapter 1--3 gate, the tree had 35 test files and 332 named `test_that`
blocks. The focused evidence totals 3,555 passing expectations: Chapter 1
contributes 104, Chapter 2 contributes 2,432, and Chapter 3 contributes 1,019.
The final Chapter 3 source tarball completed a fresh C++17 build and
installation, all examples, all tests, all 107 Rd topics, both vignettes, and
vignette rebuilding with zero errors and zero warnings. The only NOTE was the
known Windows linked-symbol scan described above.

| Module | Focused expectations | Independent reference |
|---|---:|---|
| Classical covariance/sphericity/sign-rank tests | 56 | Wishart/F/chi-square formulae, exact score moments, corrected Mauchly/John/Nagao factors, invariances and singularity boundaries |
| Gaussian/light-tail covariance estimators | 56 | direct threshold maps, product variability, all four RLZ rules, adaptive thresholds and supplied-factor POET identities |
| Gaussian EC2/off-diagonal glasso/CLIME completion | included in final unified suite | convex-l1 EC2 equality/subgradient/spectral-dual/complementarity certificates; graphical-lasso SPD/full-KKT checks; CLIME per-column primal/dual/gap certificates and primary symmetrization |
| High-dimensional Gaussian/light-tail tests | 68 | direct trace and leaveout U-statistics, corrected centring/standardisation, primal/dual calculations and label exchange |
| Elliptical sign/rank/adaptive sphericity | 104 | literal bias choices, ordered quadruples, SSCM maximum/Gumbel and stable truncated-Cauchy combination |
| SSCM proportionality and certified precision | 174 | literal ordered formulae; SCLIME feasibility/duality/symmetrisation; SGLASSO PD and KKT certificates |
| Tensor elliptical graphical model | 184 | Kolda unfolding, pilot switch, modewise scatter/objective, penalty scaling, SPD/descent/KKT and threshold support |
| Ollila shrinkage/BASIC/BASICS/linear pooling | 165 | direct robust moments and shrinkage, self-contained quadrature/inversion, pooling QP and boundary failures without official-code repairs |
| High-dimensional HR | 157 | primary joint banded fixed point, map/score residuals, pilot handling, trace-p identification, scale/permutation and no-repair cases |
| Elliptical factor reconstruction and precision | 55 | threshold/factor selectors, POET residual reconstruction, one-step Tyler, CLIME/GLASSO and Woodbury identities |

These are deterministic formula, invariance, solver-certificate, and numerical
boundary checks. The final Chapter 3 surface contains 38 exports, including the
formula-complete convex-l1 EC2, Gaussian off-diagonal graphical-lasso, and
Gaussian CLIME APIs. Adaptive/MC+ EC2, SCIO, and scaled-lasso remain
review-only; no duplicate or uncertified solver is exposed. The checks do not
reproduce any source paper's Monte Carlo design, simulation grid, empirical
size/power table, or data analysis.

## Chapter 4 formula and regression evidence

At the earlier Chapter 1--4 gate, the tree had 41 test files and 426 named `test_that`
blocks. Focused evidence totals 4,327 passing expectations: Chapter 1
contributes 104, Chapter 2 contributes 2,432, Chapter 3 contributes 1,019, and
Chapter 4 contributes 772. The final Chapter 4 source tarball completed a fresh
C++17 build and installation, all examples, all tests, all 147 Rd topics, both
vignettes, and vignette rebuilding with zero errors and zero warnings. The
only NOTE was the known Windows linked-symbol scan described above.

| Module | Focused expectations | Independent reference |
|---|---:|---|
| Unconditional alpha tests | 123 | exact GRS/Pesaran--Yamagata/MAX formulae, primary combinations, robust sign/radial corrections, projection and failure contracts |
| Conditional alpha and FDR | 122 | literal restricted/unrestricted residual constructions, feasible sum/max/CSS/CSM components, sieve ordering, FDR thresholds, RNG isolation |
| Change point and ERHT-WBS | 198 | direct CUSUM/DMS/sign/ERHT scan references, corrected pivots, endpoint nuisance fits, WBS geometry, intrinsic calibration and no-repair paths |
| White noise | 108 | portmanteau baseline, feasible FLM inner-product sum/max, spatial-sign variance, primary rank maxima and unsupported-family boundaries |
| Radial-directional diagnostic | 64 | direct radial/directional score and null calibration, invariances, extreme-scale and degenerate-radius cases |
| Independence | 157 | Wilks/Bartlett, Pesaran CD, panel sum/max, serial covariance thresholding, rank-permutation variance, Fisher/min-p and design boundaries |
| Formula-closed Chapter 4 completion | included in final unified suite | weighted-sign alpha oracle and feasible inverse-norm endpoint; book Gaussian-alpha and supplied conditional-Wald benchmarks; exact Hoeffding D/BKR R/tau-star rank-U kernels with workload and permutation contracts |

Intrinsic method calibration is tested with fixed seeds, exact small reference
spaces, or direct deterministic reconstructions. The final Chapter 4 surface
contains 33 exports. Its mutual-independence max--sum studentization row remains
source-blocked because the available sources do not close the calibration; the
implemented exact rank-U vector method is not presented as an adjacent-paper
substitution. No Chapter 4 test reproduces a paper's simulation scenario,
empirical size/power table, or data analysis.

## Chapter 5 formula, solver and regression evidence

At the earlier Chapter 1--5 gate, the tree had 45 test files and 495 named `test_that`
blocks. The final installed-package run reported 4,772 passing expectations:
Chapter 1 contributes 104, Chapter 2 contributes 2,432, Chapter 3 contributes
1,019, Chapter 4 contributes 772, and Chapter 5 contributes 445. The source
tarball completed a fresh C++17 build and installation, all examples, all 170
Rd topics, both vignettes, and vignette rebuilding with zero errors and zero
warnings. The only NOTE was the known Windows linked-symbol scan.

| Module | Focused expectations | Independent reference |
|---|---:|---|
| Shared contract, oracle and classical classifiers | 103 | direct Gaussian/elliptical scores, priors, label and feature contracts, canonical/twice-log-likelihood orientation, SPD and non-Gaussian counterexamples |
| Sparse and spatial-sign linear classifiers | 178 | FAIR/Shao thresholds, exact DSDA coding, vector-Dantzig primal/dual/KKT, explicit ridge operators, certified plug-in fits and leakage-free tuning |
| Sparse Gaussian and spatial-sign QDA | 115 | Li--Shao thresholding, Jiang losses/intercept, matrix-free SDAR/SSQDA Dantzig constraints, triple-U trace identity, determinant sign and post-symmetry feasibility |
| GQDA and HR-GQDA | 49 | literal primary inequality for positive/negative/zero log-determinant contrasts, exhaustive breakpoint selection, Gaussian endpoint identities and certified robust adapters |

All classifier tests are deterministic formula, prediction, invariance,
solver-certificate, tuning-isolation and failure-contract checks. They do not
reproduce a paper simulation design, error-rate table, tuning grid or benchmark
data analysis.

## Chapter 6 formula, solver and regression evidence

| Module | Independent reference and contract |
|---|---|
| Classical PCA/CCA and factor methods | Explicit centering and divisors; direct SVD/eigen and Wilks--Bartlett references; view exchange; repeated-root projector comparison; square-root-p loadings; OLS factor scores; MKER/MKTCR, Kendall-ratio and eigengap boundaries |
| Robust spectral PCA | Certified Chapter 1 spatial signs and Kendall matrices; generalized radial transforms at all endpoints; median/raw-MAD versus h-order cutoff families; LTS support, tie and zero rules; projector rather than basis-vector comparison; no hidden eigenvalue floor |
| Sparse PCA/CCA | Direct TPM truncation and ties; Fantope PSD, trace, primal/dual and residual certificates; PMD l1 boundaries, deflation, paired-view invariance, and explicit nonconvergence/failure states |
| Primary and book SSCCA | The primary metric-lasso/BIC program retains p-scaled metric blocks, non-unit diagonals and KKT certificates; the distinct book whitened-PMD interface requires an explicit ridge and never substitutes a pseudoinverse or silent whitening repair |

The four Chapter 6 implementation modules expose 16 public APIs. Their four
test files contain 77 named `test_that` blocks and 539 static source
occurrences of `expect_*`. These are source counts, not a claim
about a separately measured runtime assertion total. They cover deterministic
formula, invariance, rank, support, certificate and failure contracts. No
paper tuning grid, Monte Carlo factor scenario, empirical table or benchmark
data analysis is reproduced.

## Chapter 7 formula, solver and regression evidence

| Module | Static `expect_*` calls | Independent reference and contract |
|---|---:|---|
| CHIME and IF-PCA | 136 | Exact E/M and lasso KKT reconstruction; stable logistic evaluation; canonical component exchange after initialization and every update so `omega <= 1/2`; full supplied lambda path; paper/software KS scaling, HCT null contracts, deterministic ties, local seeds and caller-RNG restoration |
| Classical and sparse clustering | 150 | Direct WCSS, full-covariance likelihood and E/M references; ordered-pair BCSS `2 * (TSS - WCSS)`; s=1 soft-threshold KKT boundary; deterministic initialization/ties; explicit empty repair, rank/SPD, zero-objective, Gap and RNG contracts |
| Spatial clustering and selectors | 152 | Unsquared spatial-median center objective versus squared or metric assignment; `U(0)`; SSCM ridge/inverse certificates; active, empty, cycle and final-state reconstruction; Gap/BWDM selectors; close-large, subnormal, opposite-DBL_MAX and non-diagonal Cholesky metric fixtures |
| SEMC fitting, prediction, and Gap-LSE selection | 93 | primary formula plus pinned official-software fixtures; strict Tyler/POET/glasso SPD and solver certificates; paper hard-label versus software soft-delta dispersion; prediction invalidation and provenance contracts |

Chapter 7 now has 14 exported APIs plus the prediction S3 method and 19
registered compiled kernels. Its four test files contain 76 named blocks and
531 static `expect_*` calls; runtime outcomes are included in the unified 6,267-
pass installed-suite result below. Strict `sourceCpp` gates passed and shared
Chapter 7 source/test hashes matched the staging snapshot.

SEMC is no longer source-blocked. `semc_fit()`, `predict.semc_fit()`, and
`semc_select_k_gap()` are an independent rewrite against `FengZhuang2026SEMC`
and the pinned
`flnankai/GEMcluster@10fce04fe690fe274dd5d237cfcd3d5c6a4139f6` MIT release.
Primary-paper and software-only contracts remain distinct; fixed formula and
official-software fixtures contain no copied author-code expression or paper
simulation. Chapter 7 selector permutations and IF-PCA null draws are likewise
intrinsic calibrations with fixed inputs or explicit local seeds, not paper
simulation replications.

## Final Chapter 1--7 integration outcome (2026-08-15)

The fresh ASCII staging run contains 176 exports, 3 S3 registrations, and 140
Rcpp attributes matched one-for-one by 140 native registrations. Roxygen
produced 223 Rd files covering 169 public topics. The source built, installed,
and loaded successfully; all Rd files parsed.

Across 56 test files, the tree contains 688 named `test_that` blocks and 5,781
static `expect_*` calls. The installed-package reporter recorded 6,267 passes
with zero failures, errors, warnings, or skips. Two independent generation runs
produced the same 226-file hash manifest, and the synchronized release tree
matched 226/226 generated files.

The normative 87-row traceability ledger records 83 implemented, 2
review-only, and 2 source-blocked rows. SCIO/scaled-lasso and unaudited
classifier families remain review-only. The two source-blocked rows are the
Chapter 2 structured-correlation review and Chapter 4 mutual-independence
studentization. SEMC is implemented. On the exact tarball identified above,
the standard check completed with 0 ERROR, 0 WARNING, and 1 NOTE; the as-CRAN
check completed with 0 ERROR, 0 WARNING, and 2 NOTEs. The shared linked-symbol
NOTE, the additional new-submission NOTE, and the zero-direct-call result are
recorded above. The final as-CRAN check kept incoming remote checks enabled
after the repository and issue tracker were published.

## Deterministic software benchmarks

The installed-package [benchmark harness](../inst/benchmarks/README.md) runs one
representative callable workflow for every book chapter. The checked-in
[Markdown report](../output/benchmarks/BENCHMARK-REPORT.md) and
[CSV output](../output/benchmarks/HDElliptical-benchmarks.csv) record
machine-local timings and deterministic result fingerprints. Repeated
fingerprints matched. This is practical software evidence, not a paper
simulation, size/power study, empirical application, or performance guarantee.

## Release gate

Every subsequent integration must rerun, from a source tarball:

1. generated Rcpp registration and roxygen documentation;
2. all testthat suites;
3. examples and both vignettes;
4. namespace load/unload checks;
5. a clean source build with no compiled objects left in `src/`; and
6. formula-specific randomized references for any changed statistical kernel.
