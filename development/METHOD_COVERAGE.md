# Method coverage

This file is the implementation ledger for the book dated 2026-08-14. A method
is marked **implemented** only when its public API, documentation, compiled
kernel where appropriate, and regression tests are all present.
Package status in this ledger is current through 2026-08-15.

The normative row-level inventory is `METHOD_TRACEABILITY.csv`; it records an
explicit status, manuscript range or primary-only boundary, primary key/URL,
public API, exact native symbol or `PURE_R`, test file, and Rd topic for every
method family. Run `Rscript development/validate-method-traceability.R` from
the package root (or the analogous workspace command) to validate it. Aggregate
counts later in this file are the preceding full integration snapshot and are
intentionally not recomputed by this traceability-only edit.

| Chapter | Scope | Status | Public API or next unit |
|---|---|---|---|
| 1 | Directional foundations | Implemented; independent final verification complete | `spatial_sign()`, `spatial_median()`, `spatial_rank()`, `sscm()`, `spatial_kendall()`, `spatial_rank_covariance()`, `tyler_shape()`, `acg_loglik()`, `hr_estimator()`, `relliptical()`, `rspherical()`, `normalize_shape()` |
| 2 | Location estimation and testing | Implemented except for one source-blocked passage; all formula-complete methods have callable APIs and deterministic checks | `hotelling_one_sample_test()`, `hotelling_two_sample_test()`, `spatial_sign_test()`, `spatial_signed_rank_test()`, `spatial_rank_test()`, `srivastava_du_one_sample_test()`, `park_ayyala_one_sample_test()`, `bai_saranadasa_two_sample_test()`, `chen_qin_two_sample_test()`, `srivastava_katayama_kano_two_sample_test()`, `cai_liu_xia_two_sample_test()`, `feng_zou_wang_zhu_two_sample_test()`, `composite_t2_two_sample_test()`, `xu_lin_wei_pan_aspu_test()`, `wang_peng_li_one_sample_test()`, `feng_sun_one_sample_test()`, `inst_one_sample_test()`, `scaled_spatial_median()`, `spatial_sign_max_test()`, `spatial_sign_maxsum_test()`, `weighted_scaled_spatial_median()`, `yan_zhao_feng_weighted_max_test()`, `yan_zhao_feng_weighted_maxsum_test()`, `feng_zou_wang_two_sample_sign_test()`, `li_wang_zou_two_sample_sign_test()`, `feng_zhang_liu_spatial_rank_test()`, `tinst_two_sample_test()`, `zhang_feng_rank_one_sample_test()`, `zhang_feng_rank_two_sample_test()`, `zhang_zhou_guo_one_sample_test()`, `zhang_zhou_guo_paired_test()`, `zhang_zhou_guo_linear_hypothesis_test()`, `zhang_zhu_zhang_two_sample_test()`, `wang_xu_approx_randomization_test()`, `feng_wang_pdq_two_sample_test()`, `zhao_feng_strongcorr_sign_test()`, `elliptical_regularized_hotelling_test()`, `elliptical_regularized_hotelling_cauchy_test()`, `generic_weighted_hr_location()`, `oracle_weighted_sign_sum_test()`; blocked pending source: `ch2_location.tex:2438-2451` structured-correlation review (no API invented) |
| 3 | Matrix estimation and testing | Implemented; 35 exported APIs, deterministic formula checks, and clean source-package check complete | Classical and high-dimensional covariance tests, certified Gaussian EC2/graphical-lasso/CLIME estimation, elliptical sphericity/proportionality, robust precision, tensor graphs, SSCM shrinkage, HDHR, and factor reconstruction; see Chapter 3 traceability below |
| 4 | Other high-dimensional tests | Implemented; 28 exported APIs, deterministic formula checks, intrinsic-calibration checks, and a clean source-package check complete | Alpha/FDR and conditional-alpha procedures, formula-complete weighted-alpha/Wald completions, change point, white noise, radial-directional diagnostics, and simple plus exact degenerate vector independence; see Chapter 4 traceability below |
| 5 | Classification | Implemented; 23 exported APIs, a common certified binary-classifier contract, deterministic solver/formula checks, and clean source-package check complete | Oracle and plug-in LDA/QDA, FAIR and sparse linear rules, sparse QDA, spatial-sign classifiers, and robust GQDA/HR-GQDA; see Chapter 5 traceability below |
| 6 | Dimension reduction | Implemented; deterministic formula and solver-certificate verification complete | 16 APIs for classical PCA/CCA, robust factors, sign/Kendall/generalized-sign PCA, sparse PCA/CCA, and distinct primary/book SSCCA interfaces; see Chapter 6 traceability below |
| 7 | Clustering | Implemented, including SEMC; deterministic formula, selector, RNG, solver-certificate, and numerical-boundary verification complete | Lloyd/GMM, sparse methods, CHIME/IF-PCA, K-spatial/SM-SSCM/Sparse-SM selectors, and `semc_fit()`/`predict.semc_fit()`/`semc_select_k_gap()`; see Chapter 7 traceability below |

## Chapter 1 traceability

| Book lines | Method | Citation key | R API | Compiled kernel | Verification |
|---|---|---|---|---|---|
| `ch1_foundations.tex:499-542` | Spatial sign | `Oja2010` | `spatial_sign()` | `cpp_spatial_sign` | Zero convention; positive-scale invariance |
| `ch1_foundations.tex:547-590` | Spatial median | `Oja2010` | `spatial_median()` | `cpp_spatial_median` | Estimating equation; translation equivariance; coincident-point handling |
| `ch1_foundations.tex:592-624` | Empirical spatial rank | `Oja2010` | `spatial_rank()` | `cpp_spatial_rank` | Direct all-pairs R calculation |
| `ch1_foundations.tex:629-703` | SSCM | `VisuriOjaKoivunen2000` | `sscm()` | Shared sign kernel | Unit trace; symmetry; orthogonal equivariance |
| `ch1_foundations.tex:706-778` | Multivariate spatial Kendall | `HanLiu2018ECA; FanLiuWang2018` | `spatial_kendall()` | `cpp_spatial_kendall` | Exact U-statistic; translation invariance; unit trace |
| `ch1_foundations.tex:781-793` | Rank covariance | `VisuriOjaKoivunen2000` | `spatial_rank_covariance()` | Shared rank kernel | Direct rank outer-product calculation |
| `ch1_foundations.tex:804-1055` | Tyler shape | `Tyler1987` | `tyler_shape()` | `cpp_tyler_shape` | Fixed-point residual; trace normalization; scale and affine equivariance |
| `ch1_foundations.tex:948-1017` | ACG likelihood | `Tyler1987` | `acg_loglik()` | R/BLAS | Direction and shape-scale invariance |
| `ch1_foundations.tex:1058-1242` | Hettmansperger-Randles estimator | `HettmanspergerRandles2002` | `hr_estimator()` | `cpp_hr_estimator` | Joint equation residuals; trace normalization |
| `ch1_foundations.tex:74-214` | Elliptical stochastic representation | foundational | `relliptical()`, `rspherical()` | R/BLAS | Dimension, singular-shape, and fixed-radius tests |

## Chapter 2 traceability

| Book lines | Method | Citation key | R API | Compiled kernel | Verification |
|---|---|---|---|---|---|
| `ch2_location.tex:146-168` | One-sample Hotelling \(T^2\) | `Hotelling1931` | `hotelling_one_sample_test()` | `cpp_hotelling_one_sample` | Exact F calibration; univariate t-test identity; affine invariance; singular-covariance rejection |
| `ch2_location.tex:197-231` | Two-sample Hotelling \(T^2\) | `Hotelling1931` | `hotelling_two_sample_test()` | `cpp_hotelling_two_sample` | Pooled-covariance formula; univariate pooled-t identity; affine and sample-swap invariance; dimension/rank rejection |
| `ch2_location.tex:315-330` | One-sample classical spatial-sign test | `MottonenOja1995; Oja2010` | `spatial_sign_test()` | `cpp_spatial_sign_test` | Independent R sign calculation; translation/scale/orthogonal invariance; zero convention and singular-moment rejection |
| `ch2_location.tex:332-361` | One-sample classical spatial signed-rank test | `MottonenOja1995; Oja2010` | `spatial_signed_rank_test()` | `cpp_spatial_signed_rank_test` | Independent ordered-pair R calculation; Hoeffding factor \(1/4\); translation/scale/orthogonal invariance; zero-pair and degeneracy checks |
| `ch2_location.tex:370-393` | Two-sample pooled spatial-rank test | `MottonenOja1995; Oja2010` | `spatial_rank_test()` | `cpp_spatial_rank_test` | Independent pooled-rank R calculation with sample-covariance denominator \(N-1\); old \(N\)-denominator correction; group-swap and similarity invariance; extreme-value, self/tie, variable/rank checks |
| `ch2_location.tex:409-422` | Srivastava--Du one-sample diagonal test | `SrivastavaDu2008` | `srivastava_du_one_sample_test()` | `cpp_srivastava_du_one_sample` | Direct covariance/correlation formula; primal/dual Gram agreement; coordinatewise rescaling and extreme global-scale checks; invalid variance rejection |
| `ch2_location.tex:423-442` | Park--Ayyala one-sample leave-two-out test | `ParkAyyala2013` | `park_ayyala_one_sample_test()` | `cpp_park_ayyala_one_sample` | Literal ordered-pair leave-two-out reference; finite-sample correction in numerator and variance; permutation/diagonal-scale/translation tests; degenerate leave-out variance rejection |
| `ch2_location.tex:443-457` | Original common-covariance Bai--Saranadasa test | `BaiSaranadasa1996` | `bai_saranadasa_two_sample_test()` | `cpp_bai_saranadasa_two_sample` | Pooled-covariance formula kept distinct from CQ; primal/dual Gram agreement; group swap, orthogonal and global-scale checks |
| `ch2_location.tex:458-519` | Chen--Qin unequal-covariance U-statistic | `ChenQin2010` | `chen_qin_two_sample_test()` | `cpp_chen_qin_two_sample` | Direct U-statistic and literal leave-out variance references; stable scaled components; group swap/rotation/translation/extreme-scale checks; non-positive variance rejection |
| `ch2_location.tex:520-627` | Srivastava--Katayama--Kano two-sample diagonal test | `SrivastavaKatayamaKano2013` | `srivastava_katayama_kano_two_sample_test()` | `cpp_skk_two_sample` | Original-paper/corrigendum scaling; direct matrix reference; primal/dual Gram agreement; group swap and coordinatewise-scale invariance; degenerate denominator rejection |
| `ch2_location.tex:787-849` | Cai--Liu--Xia precision-adjusted maximum test | `CaiLiuXia2014; CaiLiu2011AdaptiveThreshold` | `cai_liu_xia_two_sample_test()` | `cpp_clx_two_sample`; `cpp_clx_adaptive_precision` | Separate oracle/supplied-feasible/adaptive-threshold paths; feasible transformed-variance denominator with divisors \(n_k\); pooled adaptive threshold with divisor \(N\); extreme-value inversion, equivariance, indefinite-feasible, eigen-floor and failure-contract tests. CLIME is not presented as an additional package backend |
| `ch2_location.tex:628-777` | Scale-invariant high-dimensional Behrens--Fisher test | `FengZouWangZhu2015BF` | `feng_zou_wang_zhu_two_sample_test()` | `cpp_fzwz_bf_two_sample` | Primary-paper asymptotic centering; corrected single power of \(\gamma\) in the group-2 bias term; exact leave-four-out and 2+2 leaveout trace estimators; local dual-Gram computation; formula, equivariance, minimum-sample and no-repair tests |
| `ch2_location.tex:778-784` | Common-covariance Composite \(T^2\) test | `FengZouWangZhu2017CT2` | `composite_t2_two_sample_test()` | `cpp_composite_t2_two_sample` | Every 2+2 leaveout term rebuilds the pooled covariance, greedy correlation blocks, and block solves; group-1 leave-four-out trace calibration with denominator \(2P_{n_1}^4\); formula, label-asymmetry, equivariance, extreme-scale and strict-SPD tests |
| `ch2_location.tex:855-925` | Analytical adaptive sum-of-powers test (aSPU) | `XuLinWeiPan2016` | `xu_lin_wei_pan_aspu_test()` | `cpp_aspu_standardize`; `cpp_aspu_power_moments` | Explicit primary-paper raw and book-studentized score scales; Isserlis moments for arbitrary Gaussian covariance; deterministic Miwa odd/even families and Gumbel maximum; covariance-path, PSD, scale-definition, extreme-scale and input-contract tests |
| `ch2_location.tex:1215-1415` | Generic weighted HR location and oracle weighted-sign sum statistic | `FengLiuMa2021INST; YanFengZhang2025InverseNormMaxsum` | `generic_weighted_hr_location()`, `oracle_weighted_sign_sum_test()` | `cpp_ch2_generic_weighted_initial`; `cpp_ch2_generic_weighted_geometry`; `cpp_ch2_generic_weighted_step`; `cpp_ch2_generic_weighted_quadratic` | Literal weighted location and unweighted diagonal HR updates; arbitrary scalar callback and power endpoints; the test remains explicitly oracle and is calibrated only by supplied `null_sd` or supplied `nu2` plus `trace_R2`; no feasible plug-in is guessed |
| `ch2_location.tex:1422-1550` | Wang--Peng--Li raw spatial-sign test | `WangPengLi2015` | `wang_peng_li_one_sample_test()` | `cpp_wang_peng_li_one_sample` | Literal primary-paper equation (7) leave-two-out variance; `U(0)=0` retained; equation (8) shortcut deliberately avoided when zero signs occur; extreme-scale, overflow-fallback, invariance, formula and degeneracy tests |
| `ch2_location.tex:1551-1638` | Feng--Sun scalar-invariant one-sample spatial-sign test | `FengSun2016` | `feng_sun_one_sample_test()` | `cpp_feng_sun_one_sample` | Literal pair-specific leave-two-out joint location/diagonal fits; numerator uses diagonal only while feasible trace uses location and diagonal; ordered-pair variance; equation-residual, invariance, extreme-scale and no-repair tests |
| `ch2_location.tex:1639-1765` | One-sample inverse-norm sign test (INST) | `FengLiuMa2021INST` | `inst_one_sample_test()` | `cpp_inst_one_sample` | Endpoint signs remain centered at the null value; diagonal nuisance fits are leave-two-out unweighted joint location/scale fits; direct feasible variance uses the supplement's \(2n^{-4}\) ordered-pair formula; oracle factorization is diagnostic only; formula, invariance, nonconvergence, zero-residual, negative-variance and overflow-fallback tests |
| `ch2_location.tex:927-1324,1766-1892,1931-1993` | Scaled spatial median and spatial-sign MAX/MAXSUM | `LiuFengZhaoWang2025MaxsumLocation` | `scaled_spatial_median()`, `spatial_sign_max_test()`, `spatial_sign_maxsum_test()` | `cpp_scaled_spatial_median` | Joint location/diagonal fixed-point equations; explicit update-stability versus score-residual diagnostics; exact \(\widehat\zeta_1^2\), \(1-n^{-1/2}\), Gumbel, and stable signed-log Cauchy factors; Feng--Sun SUM reuse; formula, equivariance, extreme-tail, zero-radius and nonconvergence tests |
| `ch2_location.tex:927-1324,1893-1996` | General-\(m\) weighted scaled median and weighted MAX/MAXSUM | `YanFengZhang2025InverseNormMaxsum` | `weighted_scaled_spatial_median()`, `yan_zhao_feng_weighted_max_test()`, `yan_zhao_feng_weighted_maxsum_test()` | `cpp_weighted_scaled_spatial_median`; `cpp_yzf_weighted_max`; `cpp_yzf_weighted_maxsum` | Estimating-equation-derived \(D^{1/2}\) update; sample radial moments; direct leave-two feasible SUM variance; \(m=-1\) componentwise identity with INST; stable Gumbel/Cauchy, equivariance, zero-radius and strict/no-repair tests |
| `ch2_location.tex:2003-2126` | Two-sample multivariate-sign test | `FengZouWang2016JASA` | `feng_zou_wang_two_sample_sign_test()` | `cpp_feng_zou_wang_two_sample_sign` | Full and leave-one-out diagonal HR fits; crossed-location numerator; Proposition 2 feasible \(c_k\), three trace terms and upper-tail variance calibration; corrected bridge matrices and first powers of the \(c\)-ratios; formula, invariance, extreme-scale and no-repair tests |
| `ch2_location.tex:2127-2241` | Simpler bias-corrected two-sample spatial-sign test | `LiWangZou2016SimpleTwoSample` | `li_wang_zou_two_sample_sign_test()` | `cpp_li_wang_zou_two_sample_sign` | Two full-sample diagonal HR fits; crossed-location statistic; exact feasible bias and three trace plug-ins; variance uses outer \(n_k^2\) factors; formula, label swap, invariance, extreme-scale, convergence and no-repair tests |
| `ch2_location.tex:2242-2351` | High-dimensional two-sample spatial-rank test | `FengZhangLiu2020SpatialRank` | `feng_zhang_liu_spatial_rank_test()` | `cpp_feng_zhang_liu_spatial_rank` | Feasible leave-four/leave-two trace calibration; corrected second \(p^2\) factor; explicit geometric versus literal paper-trace scale identification; ordered-factor, equivariance, extreme-scale, convergence and no-repair tests |
| `ch2_location.tex:2352-2411` | Two-sample inverse-norm sign test (tINST) | `HuangLiuZhouFeng2023TwoSampleINST` | `tinst_two_sample_test()` | `cpp_tinst_two_sample` | Primary-paper feasible cross statistic; all nuisance fits are observation-wise leave-one-out; corrected group-size, diagonal-index and group-2 trace typos; full raw variance decomposition and iteration/score diagnostics; no variance repair |
| `ch2_location.tex:2412-2437` | Adaptive rank max, squared-rank sum, and Cauchy tests | `ZhangFeng2024AdaptiveMean; OuyangLiuTongXu2022Rank` | `zhang_feng_rank_one_sample_test()`, `zhang_feng_rank_two_sample_test()` | `cpp_zhang_feng_one_sample_scores`; `cpp_zhang_feng_two_sample_scores`; `cpp_zhang_feng_parzen_tau` | Literal signed-rank/Wilcoxon--Mann--Whitney score moments; Gumbel max; supplied or explicitly lagged Ouyang--Parzen long-run variance; stable Cauchy; tie, symmetry/pure-shift, coordinate-order and no-floor contracts |
| `ch2_location.tex:2438-2451` | Structured-correlation combinations | `LiuZhaoFengWang2025StructuredCorr` | None: blocked pending source; no API invented | None | A legally accessible primary full text has not been obtained and the book gives only a review paragraph, without an executable statistic or calibration formula; implementation is blocked rather than inferred |
| `ch2_location.tex:2593-2621` | One-sample normal-reference mean test, paired test, and same-unit linear hypothesis | `ZhangZhouGuo2022` | `zhang_zhou_guo_one_sample_test()`, `zhang_zhou_guo_paired_test()`, `zhang_zhou_guo_linear_hypothesis_test()` | `cpp_zhang_zhou_guo_one_sample` | Deterministic literal checks of the centred U-statistic, three unbiased traces including the appendix \(v+4\) denominator, matched chi-square calibration, primal/dual agreement, exact wrapper transformations, invariances, extreme scales, and no-repair failures |
| `ch2_location.tex:2622-2641` | Normal-reference scale-invariant two-sample test | `ZhangZhuZhang2023` | `zhang_zhu_zhang_two_sample_test()` | `cpp_zhang_zhu_zhang_two_sample` | Deterministic literal checks of crossed covariance weights, bias-corrected squared traces, Welch--Satterthwaite and paper-threshold degrees of freedom, primal/dual agreement, invariances, tail underflow, and degenerate-calibration rejection |
| `ch2_location.tex:2642-2666` | Approximate randomization for the high-dimensional Behrens--Fisher problem | `WangXu2022` | `wang_xu_approx_randomization_test()` | `cpp_wang_xu_approx_randomization` | Deterministic literal checks of the full-sample CQ statistic, adjacent half-differences and odd-row discard, exhaustive small-sample sign enumeration, counter-based Monte Carlo signs, \(\geq\) ties, exact/no-+1 versus Monte Carlo/+1 tails, seed isolation, invariances, and boundaries |
| `ch2_location.tex:2670-3083` | Pairwise-difference-quantile two-sample spatial-sign test | `FengWang2026PDQ` | `feng_wang_pdq_two_sample_test()` | `cpp_feng_wang_pdq_two_sample` | Deterministic literal checks of empirical U-quantiles, standardised spatial medians, \(K_1/K_2/K_3\), diagonal-deleted bias and variance, counter-based Rademacher bootstrap, paper critical and auxiliary +1 decisions, invariances, extreme units, zero directions, and no-repair failures |
| `ch2_location.tex:3084-3206` | Strong-correlation one-sample spatial-sign refinement | `ZhaoFeng2026NoteOneSample` | `zhao_feng_strongcorr_sign_test()` | `cpp_zhao_feng_strongcorr_sign_bootstrap` | Deterministic literal Rademacher/Gaussian pair-form checks; null-versus-fitted centring, diagonal multiplier term, cancellation of the common \(\tau\) factor, paper critical and auxiliary +1 decisions, seed isolation, invariances, zero signs, convergence, extreme scales, and no-repair boundaries |
| `ch2_location.tex:3207-3741` | Elliptical regularized Hotelling fixed-ridge and Cauchy-grid tests | `FengZhouWang2026ERHT` | `elliptical_regularized_hotelling_test()`, `elliptical_regularized_hotelling_cauchy_test()` | `cpp_erht_grid` | Deterministic literal checks of the feasible companion-matrix centre/variance, fixed-ridge and grid-margin identity, economy-SVD row/null-space routes, \(p>n\) complement calibration, small positive ridges, stable signed-log Cauchy aggregation, invariances, extreme scales, and explicit unsupported/degenerate failures |

### Chapter 2 completion status

The rows below retain the late-chapter implementation queue as a release
checklist. Every formula-complete method in the current Chapter 2 draft is
implemented. The structured-correlation review at lines 2438--2451 is the
only blocked item: a legally accessible primary full text has not been
obtained and the book supplies no executable formulas, so no API is invented.
The formula-complete generic weighted equations at lines 1215--1415 are now
exposed through `generic_weighted_hr_location()` and the deliberately oracle
`oracle_weighted_sign_sum_test()`; absent supplied population calibration the
latter returns the raw score without manufacturing a p-value.

| Book lines | Method family | Citation key | Status / dependency |
|---|---|---|---|
| `ch2_location.tex:778-784` | Common-covariance Composite \(T^2\) test | `FengZouWangZhu2017CT2` | Implemented, exported, and source-package checked |
| `ch2_location.tex:855-925` | Adaptive sum-of-powers (aSPU) test | `XuLinWeiPan2016` | Implemented, exported, and source-package checked; dual primary-paper raw and book-studentized APIs, no paper-specific simulations |
| `ch2_location.tex:927-1324` | Scaled/weighted spatial-median equations and Bahadur support layer | chapter support | Implemented and exported as `scaled_spatial_median()` and `weighted_scaled_spatial_median()` with distinct equation/update diagnostics |
| `ch2_location.tex:1639-1765` | Weighted signs and one-sample INST | `FengLiuMa2021INST` | Implemented, exported, and source-package checked |
| `ch2_location.tex:1766-1996` | Spatial-sign / inverse-norm max and max--sum tests | `LiuFengZhaoWang2025MaxsumLocation; YanFengZhang2025InverseNormMaxsum` | Implemented, exported, and source-package checked; no paper-specific simulations |
| `ch2_location.tex:2003-2126` | Two-sample multivariate-sign test | `FengZouWang2016JASA` | Implemented, exported, and source-package checked |
| `ch2_location.tex:2127-2241` | Simplified bias-corrected two-sample sign test | `LiWangZou2016SimpleTwoSample` | Implemented, exported, and source-package checked |
| `ch2_location.tex:2242-2351` | High-dimensional two-sample spatial-rank test | `FengZhangLiu2020SpatialRank` | Implemented, exported, and source-package checked; distinct from fixed-\(p\) pooled-rank API |
| `ch2_location.tex:2412-2437` | Adaptive rank combinations | `ZhangFeng2024AdaptiveMean` | Implemented, exported, and source-package checked; no invented general rank-\(L_q\) family or automatic long-run-variance bandwidth |
| `ch2_location.tex:2438-2451` | Structured-correlation combinations | `LiuZhaoFengWang2025StructuredCorr` | **Blocked / pending source.** A legally accessible primary full text has not been obtained, and the book contains no statistic, estimator, calibration, or rejection formula; no public API is fabricated from the review prose |
| `ch2_location.tex:2593-2621` | One-sample, paired, and same-unit linear-hypothesis normal-reference tests | `ZhangZhouGuo2022` | Implemented as three wrappers over one audited compiled kernel; documented and deterministically verified |
| `ch2_location.tex:2622-2641` | Scale-invariant two-sample normal-reference test | `ZhangZhuZhang2023` | Implemented, documented, and deterministically verified with distinct published and unadjusted degrees-of-freedom paths |
| `ch2_location.tex:2642-2666` | Approximate-randomization calibration | `WangXu2022` | Implemented, documented, and deterministically verified with exact enumeration and reproducible Monte Carlo modes |
| `ch2_location.tex:2670-3083` | Pairwise-difference-quantile spatial-sign test with wild bootstrap | `FengWang2026PDQ` | Implemented, documented, and deterministically verified; intrinsic wild bootstrap retained as the test calibration |
| `ch2_location.tex:3084-3206` | Strong-correlation one-sample sign refinement | `ZhaoFeng2026NoteOneSample` | Implemented, documented, and deterministically verified with Rademacher and Gaussian wild-bootstrap multipliers |
| `ch2_location.tex:3207-3741` | Elliptical regularized Hotelling test and ridge aggregation | `FengZhouWang2026ERHT` | Implemented, documented, and deterministically verified for fixed ridges and analytic Cauchy aggregation; the unspecified simulation-program Bartlett modification is not guessed |

The Wang--Xu randomization and the Feng--Wang and Zhao--Feng wild bootstraps
are intrinsic reference-distribution calibrations of their callable tests.
They are not replications of a paper's Monte Carlo size/power study. The
package contains no paper simulation grids, power tables, or reproduction
scripts for these methods.

## Chapter 3 traceability

| Book lines | Method family | Citation key | R API | Compiled kernel | Verification |
|---|---|---|---|---|---|
| `ch3_matrix.tex:228-405` | Classical covariance likelihood, sphericity, and sign/rank shape tests | `Mauchly1940; John1971; Nagao1973; HallinPaindaveine2006` | `gaussian_covariance_lrt()`, `mauchly_sphericity_test()`, `john_sphericity_test()`, `nagao_identity_test()`, `hallin_paindaveine_shape_test()`, `spatial_sign_sphericity_test()` | `cpp_ch3_gaussian_sign_geometry`; `cpp_ch3_gaussian_rank_shape` | Literal Wishart and rank-score references; corrected finite-sample factors; affine/orthogonal invariances; singular and minimum-df failures |
| `ch3_matrix.tex:408-445,535-549` | Gaussian/light-tail covariance thresholding and POET | `BickelLevina2008Cov; RothmanLevinaZhu2009; CaiLiu2011AdaptiveThreshold; FanLiaoMincheva2013` | `bickel_levina_covariance_threshold()`, `rothman_levina_zhu_covariance_threshold()`, `cai_liu_adaptive_covariance_threshold()`, `poet_covariance()` | `cpp_ch3_gaussian_threshold_matrix`; `cpp_ch3_gaussian_product_variability` | Primary divisor conventions; generalized thresholds; supplied factor count; symmetry, scaling, PSD and no-repair checks |
| `ch3_matrix.tex:446-483` | EC2 sparse covariance, convex off-diagonal-lasso branch | `LiuWangZhao2014EC2` | `ec2_covariance()` | `cpp_ch3gp_ec2_l1` | Primary correlation-scale unit-diagonal program, explicit eigenvalue constraint, and equality/fixed-point/subgradient/spectral-dual/complementarity certificates; adaptive and MC+ branches remain review-only |
| `ch3_matrix.tex:484-492` | Gaussian off-diagonal graphical lasso | `YuanLin2007; FriedmanHastieTibshirani2008` | `gaussian_graphical_lasso()` | `cpp_ch3gp_offdiag_glasso` | SPD-preserving majorization/backtracking with full diagonal and off-diagonal KKT checks; no ridge or pseudoinverse |
| `ch3_matrix.tex:493-522` | Gaussian CLIME | `CaiLiuLuo2011CLIME` | `clime_precision()` | `cpp_ch3gp_clime` | Every raw column has primal/dual/stationarity/gap certificates; smaller-absolute symmetrisation is reported without inventing an SPD guarantee |
| `ch3_matrix.tex:523-534` | SCIO and scaled-lasso precision reviews | `LiuLuo2015SCIO; SunZhang2012ScaledLasso` | None: review-only | None | The short review does not close their distinct solver/tuning contracts; no duplicate API is inferred |
| `ch3_matrix.tex:550-690` | High-dimensional Gaussian/light-tail covariance tests | `WangYao2013; ChenZhangZhong2010; FisherSunGallagher2010; LiChen2012` | `wang_yao_corrected_lrt()`, `wang_yao_corrected_john_test()`, `chen_zhang_zhong_covariance_test()`, `fisher_sun_gallagher_sphericity_test()`, `li_chen_covariance_test()` | `cpp_ch3_gaussian_trace_u_statistics`; `cpp_ch3_gaussian_li_chen` | Direct trace/U-statistic references, finite-sample centring, primal/dual identities, label exchange, extreme scales, and invalid-variance failures |
| `ch3_matrix.tex:691-1059` | Elliptical sign/rank and adaptive sphericity | `ZouPengFengWang2014Sphericity; FengLiu2017RankSphericity; ZhaoYangZhangFengWang2026AdaptiveSphericity` | `zou_peng_feng_wang_sphericity_test()`, `feng_liu_rank_sphericity_test()`, `zhao_yang_zhang_feng_wang_sign_max_test()`, `zhao_yang_zhang_feng_wang_adaptive_sphericity_test()` | `cpp_elliptical_sphericity_sign_core`; `cpp_elliptical_sphericity_rank_core` | Literal bias choices, ordered quadruple factors, SSCM maximum/Gumbel calibration, stable Cauchy combination, invariance and zero-radius/no-repair boundaries |
| `ch3_matrix.tex:1060-1316` | SSCM equality and spatial-rank proportionality tests | `ChengLiuPengZhangZheng2019SSCM; FengZhangLiu2022Proportionality` | `cheng_sscm_equality_test()`, `feng_spatial_rank_proportionality_test()` | `cpp_ch3pp_sscm_components`; `cpp_ch3pp_spatial_rank_components` | Independent ordered-pair/four-index formulae, feasible bias/variance, group exchange and similarity invariance, strict convergence and invalid-calibration contracts |
| `ch3_matrix.tex:1317-1463` | Spatial-sign precision estimation and support thresholding | `LuFeng2025Precision` | `spatial_sign_precision()`, `threshold_spatial_sign_precision()` | `cpp_ch3pp_sclime`; `cpp_ch3pp_sglasso` | SCLIME primal/dual/feasibility certificates and smaller-absolute-value symmetrisation; full-matrix SGLASSO KKT/PD certificates; exact support threshold and no hidden ridge/floor |
| `ch3_matrix.tex:1464-1769` | Tensor elliptical spatial-sign graphical model | `LiuLuZhouFengWang2025TensorEGM` | `tensor_spatial_sign_precision()`, `threshold_tensor_spatial_sign_precision()` | `cpp_ch3teg_mode_crossproducts`; `cpp_ch3teg_whitened_mode_scatter`; `cpp_ch3teg_offdiag_glasso` | Column-major/Kolda unfolding, pilot switch, modewise objective and penalty scaling, SPD/descent/KKT certificates, threshold signs, dimension/permutation/equivariance checks |
| `ch3_matrix.tex:1770-1841` | SSCM shrinkage, BASIC/BASICS, and linear pooling | `OllilaRaninen2019Shrinkage; RaninenOllila2022BASIC; RaninenTylerOllila2022LinearPooling` | `ollila_raninen_shrinkage_covariance()`, `regularized_spatial_sign_covariance()`, `basic_shape()`, `basics_shape()`, `linear_pool_covariance()` | `cpp_ollila_sample_moments`; `cpp_ollila_sign_components`; `cpp_ollila_basic_inverse`; `cpp_ollila_pool_qp` | Direct moment/shrinkage/pooling formulae, self-contained quadrature and bounded inversion, scale/eigenvalue properties, theoretical-boundary and no-extrapolation/no-floor contracts |
| `ch3_matrix.tex:1842-1968` | High-dimensional Hettmansperger--Randles estimator | `YanFengZhang2025HR` | `high_dimensional_hr()` | `cpp_ch3_hdhr_fit` | Primary joint Algorithm 2 fixed point, per-iteration banding before congruence, trace-p normalisation, final direct inverse, pilot certificates, map/score/convergence and no-repair checks |
| `ch3_matrix.tex:1969-2123` | Elliptical factor number, POET reconstruction, precision, and one-step Tyler refinement | `XuMaWangFeng2026EllipticalFactor` | `elliptical_factor_number()`, `spatial_sign_poet()`, `poet_tme()`, `elliptical_factor_precision()` | `cpp_ch3ef_tyler_one_step`; shared threshold and certified SCLIME/SGLASSO kernels | Explicit threshold constant and factor cap, ER/GR selectors, one-step Tyler identity, raw-residual POET/precision reconstruction, Woodbury formula, invariance, extreme-scale and failure-contract tests |

Chapter 3 contributes 35 exported functions across nine deterministic test
suites. The formula-complete convex EC2 l1 branch, Gaussian off-diagonal graphical
lasso, and Gaussian CLIME now have distinct certified public APIs. Adaptive or
MC+ EC2, SCIO, and scaled-lasso remain review-only because the brief book text
does not close their additional solver and tuning contracts.
No paper simulation design, tuning grid, size/power table, or reproduction
data is included.

## Chapter 4 traceability

| Method family | R API | Compiled kernel | Verification and boundary |
|---|---|---|---|
| Unconditional alpha tests | `grs_alpha_test()`, `pesaran_yamagata_alpha_test()`, `feng_lan_liu_ma_alpha_max_test()`, `gaussian_alpha_combination_test()`, `liu_feng_ma_spatial_sign_alpha_test()`, `zhao_feng_wang_wang_robust_alpha_test()` | `cpp_ch4_alpha_ols`; `cpp_ch4_lfm_spatial_sign_core`; `cpp_scaled_spatial_median` | Exact regression projection and finite-sample factors; primary Bonferroni/Fisher or truncated-Cauchy combinations; literal sign/trace/radial corrections |
| Weighted alpha completion | `weighted_spatial_sign_alpha_oracle_test()`, `zhao_chen_zi_inst_alpha_test()` | `cpp_ch4_completion_standardized_radii`; `cpp_ch4_completion_weighted_alpha_q`; shared `cpp_ch4_alpha_ols` and `cpp_ch4_lfm_spatial_sign_core` | The general K formula is explicitly oracle; the inverse-norm endpoint is the primary feasible method with audited leave-two-out trace and no radial floor |
| Book Gaussian alpha and conditional Wald benchmarks | `book_gaussian_alpha_cauchy_test()`, `conditional_factor_wald_test()` | `cpp_ch4_alpha_ols`; `PURE_R` Wald solve | Book Gaussian benchmark requires supplied `trace_R2`; the fixed-dimensional Wald benchmark requires supplied delta/covariance and strict Cholesky |
| Conditional alpha and mutual-fund FDR | `conditional_alpha_sieve_design()`, `conditional_alpha_sieve_fit()`, `wang_zhao_feng_wang_mutual_fund_fdr()`, `ma_lan_su_tsai_conditional_alpha_sum_test()`, `ma_feng_wang_bao_conditional_alpha_test()`, `zhao_conditional_spatial_sign_sum_test()`, `zhao_wang_conditional_spatial_sign_test()` | `cpp_ch4_afc_project`; `cpp_ch4_afc_spatial_kendall`; `cpp_ch4_afc_light_components`; `cpp_ch4_afc_css_components`; `cpp_scaled_spatial_median` | Restricted/unrestricted projection contracts, explicit sieve or supplied nuisance fits, feasible sum/max/CSS/CSM calibration and FDR thresholding |
| Change-point tests and segmentation | `classical_cusum_test()`, `wang_feng_dms_test()`, `spatial_sign_change_point_test()`, `erht_change_point_test()`, `erht_wbs()` | `ch4_cp_cusum_cpp`; `ch4_cp_dms_moments_cpp`; `ch4_cp_scaled_hr_cpp`; `ch4_cp_spatial_median_cpp`; `ch4_cp_ordered_pair_square_sum_cpp`; `ch4_cp_erht_moments_cpp` | Corrected DMS dense/max pivots, endpoint sign nuisance estimates, ERHT ridge scan and WBS recursion |
| White-noise tests | `white_noise_portmanteau_test()`, `feng_liu_ma_white_noise_test()`, `zhao_chen_wang_spatial_sign_white_noise_test()`, `chen_song_feng_rank_white_noise_test()` | `cpp_ch4wn_flm_core`; `cpp_ch4wn_spatial_sign_core`; `cpp_ch4wn_rank_max_core` | Feasible inner-product SUM/MAX, corrected sign variance, and primary-backed rank maxima; unsupported rank-sum combinations are not invented |
| Radial-directional diagnostic | `zhang_feng_radial_directional_test()` | `cpp_zhang_feng_radial_directional` | Literal radial-direction score, feasible centring/scaling, extreme-value tail and explicit degeneracy contracts |
| Gaussian, panel, serial-panel, and simple vector independence | `gaussian_wilks_independence_test()`, `pesaran_cd_test()`, `feng_jiang_liu_xiong_panel_independence_test()`, `wang_liu_feng_ma_serial_panel_test()`, `wang_liu_feng_vector_independence_test()` | `cpp_ch4ind_wilks_core`; `cpp_ch4ind_pairwise_core`; `cpp_ch4ind_serial_panel_core`; `cpp_ch4ind_rank_vector_core` | Cholesky Wilks/Bartlett, panel contracts, and Spearman/Kendall permutation variance |
| Exact degenerate rank-U vector independence | `wang_liu_feng_vector_u_independence_test()` | `cpp_ch4_completion_vector_u_core` | Exact Hoeffding D, BKR R and tau-star symmetrized kernels with an explicit workload guard and intrinsic permutation variance |

Chapter 4 contributes 28 exported functions across six deterministic test
suites. Functional/temporal change-point prose and unsupported rank white-noise
sum combinations remain review-only. Closed-source mutual-independence
studentisation is source-blocked and remains unimplemented because the
available source does not uniquely determine a practical calibrated method.
The distinct vector-independence primary does determine exact Hoeffding D,
BKR R and tau-star kernels; those are implemented without pretending that they
solve the still-blocked mutual-independence problem. Bootstrap,
permutation, or Gaussian-process draws are retained only when they are part of
the callable method's own reference distribution. They are not paper-specific
Monte Carlo size/power replications.

## Chapter 5 traceability

| Book lines | Method family | R API | Compiled kernel | Verification and boundary |
|---|---|---|---|---|
| `ch5_classification.tex:64-383` | Shared classifier contract and oracle/classical Gaussian or elliptical LDA/QDA | `hd_classifier_fit()`, `gaussian_lda_oracle()`, `gaussian_qda_oracle()`, `elliptical_oracle_classifier()`, `classical_lda_classifier()`, `classical_qda_classifier()` | `cpp_ch5_classical_linear_scores`; `cpp_ch5_classical_quadratic_scores`; `cpp_ch5_classical_mahalanobis_pairs`; `cpp_ch5_classical_two_class_moments` | Exact base-R scores and priors; canonical versus twice-log-likelihood score scales; factor levels, feature names and tie orientation; a non-Gaussian unequal-prior counterexample; singular/invalid fit rejection |
| `ch5_classification.tex:384-521` | Independence, FAIR, thresholded LDA and direct sparse linear discrimination | `independence_classifier()`, `fair_classifier()`, `shao_threshold_lda()`, `lpd_classifier()`, `dsda_classifier()`, `sparse_plugin_lda()` | `cpp_ch5_classical_fair_criterion`; `cpp_ch5_classical_hard_threshold_covariance`; `cpp_c5lin_dantzig`; `cpp_c5lin_dsda` | Literal variance divisors, FAIR ranking, strict threshold equality, primary DSDA labels/objective, vector-Dantzig primal/dual/KKT certificates, stratified leakage-free tuning, and explicit ridge/operator contracts |
| `ch5_classification.tex:591-955,1060-1127` | Spatial-sign sparse LDA and certified precision plug-in LDA | `sslda()`, `spatial_sign_precision_lda()` | `cpp_c5lin_dantzig`; `cpp_ch3pp_sclime`; `cpp_ch3pp_sglasso` | Classwise spatial medians/SSCMs, pooled sign operator, Dantzig feasibility and duality, valid-fit adapters, midpoint/equal-prior contract, invariance and no-repair failures |
| `ch5_classification.tex:522-590,956-1059` | Sparse Gaussian and spatial-sign QDA | `li_shao_sparse_qda()`, `jiang_da_qda()`, `sdar_qda()`, `ssqda()`, `sparse_plugin_qda()` | `cpp_c5qda_operator`; `cpp_c5qda_li_shao`; `cpp_c5qda_jiang_matrix`; `cpp_c5qda_jiang_vector`; `cpp_c5qda_dantzig_matrix`; `cpp_c5qda_dantzig_vector`; `cpp_c5qda_ssqda_moments`; `cpp_c5qda_signed_logdet` | Li--Shao three thresholds and LOO bisection; Jiang penalised interaction/linear losses and exact training-error intercept; Cai--Zhang SDAR attribution; SSQDA trace U-statistic identity; post-symmetry feasibility, positive determinant sign, equal-prior and no-repair contracts |
| `ch5_classification.tex:189-243,1060-1185` | Generalized QDA and robust HR-GQDA | `gqda_classifier()`, `robust_gqda()`, `hr_gqda()`, `hr_qda()` | `cpp_ch5_gqda_components`; `cpp_ch3_hdhr_fit` | Direct primary inequality over signed/zero log-determinant contrasts; deterministic breakpoint/end-point selection; `c=0/1` identities; class exchange and coordinate equivariance; certified robust scale/precision adapters; alias semantics |

Chapter 5 contributes 23 exports across four deterministic test suites, with
69 named test blocks and 445 passing expectations. The package score convention
is positive for class 1 and every fit records its primary orientation, score
scale, priors, tuning path, solver certificates and failure stage. Invalid or
uncertified fits are not predictable and receive no hidden ridge, pseudoinverse,
determinant absolute value, eigenvalue floor, or nearest-PD repair.

Scout, penalised Fisher discrimination, sparse optimal scoring, CODA, and
specific fixed-dimensional M/MVE/MCD/S/SD estimators remain review-only until
their distinct primary algorithms and solver contracts are audited. The generic
plug-in adapters consume already certified covariance/precision fits rather
than fabricating duplicate optimisers. No paper classification simulation,
cross-validation grid, empirical error table, or benchmark data is included.

## Chapter 6 traceability

| Book lines | Method family | R API | Compiled kernel | Verification and boundary |
|---|---|---|---|---|
| `ch6_pca_factor.tex:57-310,1018-1207` | Classical PCA/CCA, Bartlett testing, robust factor subspaces, RTS factors, and factor-number selection | `classical_pca()`, `classical_cca()`, `cca_bartlett_test()`, `robust_factor_subspace()`, `rts_factor()`, `kendall_factor_number()`, cross-chapter `elliptical_factor_number()` | `cpp_ch6cf_center_moments`; `cpp_ch6cf_cca_moments`; `cpp_ch6cf_rts_components`; certified Chapter 1 sign/Kendall/median kernels | Explicit divisor and centering; MKER/MKTCR remains distinct from the Chapter 3 spatial-sign ER/GR selector |
| `ch6_pca_factor.tex:385-419` | Approximate factor models and POET | cross-chapter `poet_covariance()` | `cpp_ch3_gaussian_threshold_matrix`; `cpp_ch3_gaussian_product_variability` | Chapter 6 explicitly points back to the Chapter 3 POET covariance implementation; it is not an omitted Chapter 6 API |
| `ch6_pca_factor.tex:448-818,942-1017` | Spatial-sign, Kendall, and generalized-spatial-sign PCA | `spatial_sign_pca()`, `kendall_pca()`, `generalized_sign_pca()` | `cpp_ch6rs_radial_transform` plus certified Chapter 1 sign/Kendall kernels | Zero and tie rules; divisor conventions; radial endpoints; median/MAD versus h-order cutoffs; LTS support; projector comparison; no hidden eigenvalue floor |
| `ch6_pca_factor.tex:345-447,819-941,1208-1342` | Sparse PCA and sparse CCA | `truncated_power_pca()`, `sparse_spatial_sign_pca()`, `fantope_pca()`, `pmd_sparse_pca()`, `pmd_sparse_cca()` | `cpp_ch6spc_tpm`; `cpp_ch6spc_fantope`; `cpp_ch6spc_pmd_pca`; `cpp_ch6spc_pmd_cca` | Exact ties; PSD and Fantope residual certificates; l1 boundaries; deflation and paired-view contracts; explicit failure rather than repair |
| `ch6_pca_factor.tex:1208-1342` plus the cited primary | Primary metric-lasso SSCCA versus the distinct book whitened-PMD summary | `sscca()`, `sign_whitened_sparse_cca()` | `cpp_ch6_sscca_primary`; `cpp_ch6_sign_whitened_pmd` | p-scaled metric blocks and non-unit diagonals; BIC and KKT certificates; explicit ridge; no pseudoinverse or whitening repair |

Chapter 6 exposes 16 reusable APIs across four implementation modules. Its
four test files contain 77 named blocks and 539 static `expect_*` call
occurrences covering formulas, invariances, ranks and solver certificates.
No paper tuning grid, Monte Carlo design, factor table, benchmark script or
source-article data is included. The cross-chapter POET and spatial-sign
factor-number rows are validated against the Chapter 3 tests and Rd topics in
`METHOD_TRACEABILITY.csv`, rather than being counted as missing Chapter 6 code.

## Chapter 7 traceability

| Book lines | Method family | R API | Compiled kernel | Verification and boundary |
|---|---|---|---|---|
| `ch7_clustering.tex:9-102` | Lloyd K-means and full-covariance Gaussian-mixture EM | `lloyd_kmeans()`, `gaussian_mixture_em()` | `cpp_ch7cs_lloyd_core`; `cpp_ch7cs_gmm_core` | WCSS and E/M identities; deterministic ties, max-min starts and seeds; explicit empty repair; equivariance; likelihood, rank and SPD certificates; no hidden ridge |
| `ch7_clustering.tex:103-138,467-512` | Sparse K-means, its selector, and sparse K-medians | `sparse_kmeans()`, `sparse_kmeans_select_s()`, `sparse_kmedian()` | `cpp_ch7cs_sparse_kmeans_core`; `cpp_ch7cs_sparse_kmedian_core` | Ordered-pair BCSS equals `2 * (TSS - WCSS)`; soft-threshold and s=1 KKT boundaries; Gap calibration, RNG isolation, and zero-objective failures |
| `ch7_clustering.tex:139-385` | CHIME and IF-PCA | `chime_clustering()`, `if_pca()` | `cpp_ch7_chime_fit`; `cpp_ch7_if_scores`; `cpp_ch7_if_null_scores`; `cpp_ch7_if_deterministic_kmeans` | Exact E/M, KKT and stable-logistic checks; canonical `omega <= 1/2` at every stage; full supplied lambda path; KS/HCT scale, ties, local seeds and caller-RNG restoration |
| `ch7_clustering.tex:386-670` | K-spatial medians, SM-SSCM, Sparse-SM and the permutation tau selector | `k_spatial_median()`, `sm_sscm()`, `sparse_k_spatial_median()`, `sparse_sm_select_tau()` | `cpp_ch7sc_assign_euclidean`; `cpp_ch7sc_sscm_metric`; `cpp_ch7sc_assign_metric`; `cpp_ch7sc_feature_scores`; `cpp_ch7sc_geometry` | Unsquared center versus squared assignment geometry; `U(0)`; active, empty, cycle and final-state contracts; book-backed permutation Gap tau selection |
| Primary only: `ZhaoZhuangFeng2026SparseKSM` (no Chapter 7 line attribution) | Sparse-SM BWDM K selector | `sparse_sm_select_k()` | `cpp_ch7sc_geometry` plus certified Sparse-SM kernels | This selector is mapped only to the primary paper; it is not the distinct SEMC Gap-LSE rule at book lines 816--834 |
| `ch7_clustering.tex:671-929` | Semiparametric elliptical mixture clustering (SEMC) | `semc_fit()`, `predict.semc_fit()`, `semc_select_k_gap()` | `cpp_ch7_semc_delta`; `cpp_ch7_semc_softmax`; `cpp_ch7_semc_weighted_sign_scatter`; `cpp_ch7_semc_weighted_tyler`; `cpp_ch7_semc_offdiag_glasso`; `cpp_ch7_semc_weighted_kde` | Formula-complete primary plus pinned official-software contracts; all `shape = tyler/poet/glasso` paths require strict SPD/certificates; invalid fits are not predictable; paper hard-label Gap-LSE and explicit software soft-delta dispersion remain distinct |

SEMC is no longer source-blocked. Its implementation is an independent rewrite against `FengZhuang2026SEMC` and `flnankai/GEMcluster@10fce04fe690fe274dd5d237cfcd3d5c6a4139f6` (MIT), with fixed formula/official fixtures, public examples, provenance diagnostics, and no copied author-code expression or paper simulation.

Chapter 7 exposes 12 reusable APIs backed by 13 compiled kernels. Its 66 named
test blocks execute 446 assertions in fresh installed-package focused runs.
Selector permutations and IF-PCA empirical-null draws are intrinsic method
calibrations with fixed inputs or explicit local seeds and caller-RNG
isolation; they are not paper simulations. No paper grids, tables, scripts or
data were added.

The signed Chapter 1--7 integration contains 176 exports, 140 registered native
routines, and 223 Rd topics. All 56 test files, 688 named blocks and 6,267
runtime assertions passed with zero failure, error, warning, or skip. All 169
standard example topics and the newly added `donttest` examples passed. The
release tarball check is rerun only after the final release-document and JSS
synchronization, so its outcome is not inferred from this source-level gate.

## Source cautions carried into implementation

- The unused `chapters/appendix_notation.tex` is not part of `main.tex` or the
  current PDF and is excluded from coverage.
- Several later-chapter LaTeX escapes are missing (`ref`/`neq`), and the tINST
  draft contains index inconsistencies. Later implementations must be checked
  against the original papers and numerical tables before being marked
  implemented.
- The displayed fixed-\(p\) signed-rank statistic in the book omits the
  Hoeffding-projection factor \(1/4\). `spatial_signed_rank_test()` uses
  \(Q_{SR}=n\bar R^T B_R^{-1}\bar R/4\), because
  \(\operatorname{Var}\{\sqrt n\bar R\}\to4B_R\).
- Review-only mentions without a self-contained algorithm are catalogued but
  will not automatically become public functions; established solver packages
  will be used where reimplementing a generic optimizer would reduce numerical
  reliability.




