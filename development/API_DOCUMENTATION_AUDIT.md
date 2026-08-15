# Public API Documentation Audit

## Scope and method

This audit covers executable examples and reference metadata for exported public APIs. It compares the checked-in `NAMESPACE` and existing `man` index with a fresh ASCII-path staging run of `Rcpp::compileAttributes()` and `roxygen2::roxygenise()`. The shared `NAMESPACE`, `man`, native registration, tests, C++, and function bodies were not edited by this documentation task.

A direct-call pass means that the generated `\examples{}` block contains a syntactic call to the exported alias. One compatibility exception was approved: `hr_qda` shares `hr_gqda.Rd`, whose example directly calls `hr_gqda()` and exercises the shared implementation. A reference pass means the export-owning Rd topic has a non-empty `\references{}` block. Book-only or package-adapter entries are explicitly labeled and are not presented as independent primary methods.

## Snapshot counts

| Metric | Checked-in/current | Fresh generated | Result |
|---|---:|---:|---|
| Exported aliases | 164 | 176 | 12 source-only exports recovered by fresh generation |
| Public export-owning Rd topics | 157 | 169 | 12 source-only public topics recovered |
| All generated Rd files | not rewritten in shared tree | 223 | all parsed |
| Alias-level direct-call gaps from initial audit | 28 | 1 | only approved shared-topic alias `hr_qda` remains |
| Topic-level reference gaps from initial audit | 35 | 0 | closed |
| Generated topics with examples executed | - | 169 | zero failures |
| Newly added direct-call fixtures executed with `run.donttest = TRUE` | - | 27 | zero failures |

## Direct-call traceability

| Exported alias | R source | Initial direct call | Fresh status | Action/evidence |
|---|---|---:|---|---|
| `acg_loglik` | `chapter1-foundations.R` | no | pass | fixed deterministic matrix call |
| `bai_saranadasa_two_sample_test` | `chapter2-quadratic-tests.R` | no | pass | fixed seeded two-sample fixture |
| `bickel_levina_covariance_threshold` | `chapter3-gaussian-estimators.R` | no | pass | fixed seeded covariance fixture |
| `cai_liu_adaptive_covariance_threshold` | `chapter3-gaussian-estimators.R` | no | pass | fixed seeded covariance fixture |
| `chen_qin_two_sample_test` | `chapter2-leaveout-tests.R` | no | pass | fixed seeded two-sample fixture |
| `chen_song_feng_rank_white_noise_test` | `chapter4-white-noise.R` | no | pass | deterministic short series call |
| `chen_zhang_zhong_covariance_test` | `chapter3-gaussian-highdim.R` | no | pass | fixed seeded covariance fixture |
| `feng_liu_ma_white_noise_test` | `chapter4-white-noise.R` | no | pass | deterministic short series call |
| `fisher_sun_gallagher_sphericity_test` | `chapter3-gaussian-highdim.R` | no | pass | fixed seeded matrix call |
| `hallin_paindaveine_shape_test` | `chapter3-gaussian-classical.R` | no | pass | fixed seeded shape fixture |
| `hr_qda` | `chapter5-gqda.R` | no | shared-pass | compatibility alias shares `hr_gqda.Rd`; `hr_gqda()` is called directly |
| `john_sphericity_test` | `chapter3-gaussian-classical.R` | no | pass | fixed seeded matrix call |
| `li_chen_covariance_test` | `chapter3-gaussian-highdim.R` | no | pass | fixed seeded two-sample covariance fixture |
| `mauchly_sphericity_test` | `chapter3-gaussian-classical.R` | no | pass | fixed seeded matrix call |
| `nagao_identity_test` | `chapter3-gaussian-classical.R` | no | pass | fixed seeded matrix call |
| `normalize_shape` | `chapter1-foundations.R` | no | pass | deterministic diagonal-matrix call |
| `park_ayyala_one_sample_test` | `chapter2-leaveout-tests.R` | no | pass | fixed seeded one-sample fixture |
| `poet_covariance` | `chapter3-gaussian-estimators.R` | no | pass | fixed seeded covariance fixture |
| `rothman_levina_zhu_covariance_threshold` | `chapter3-gaussian-estimators.R` | no | pass | fixed seeded covariance fixture |
| `rspherical` | `chapter1-foundations.R` | no | pass | explicit seed and dimensions |
| `sparse_sm_select_k` | `chapter7-spatial-clustering.R` | no | pass | fixed data, supplied permutation, explicit tuning; `donttest` executed |
| `sparse_sm_select_tau` | `chapter7-spatial-clustering.R` | no | pass | fixed data, supplied permutations, explicit tuning; `donttest` executed |
| `spatial_sign_sphericity_test` | `chapter3-gaussian-classical.R` | no | pass | fixed seeded matrix and explicit center |
| `srivastava_du_one_sample_test` | `chapter2-quadratic-tests.R` | no | pass | fixed seeded one-sample fixture |
| `wang_yao_corrected_john_test` | `chapter3-gaussian-highdim.R` | no | pass | fixed seeded matrix call |
| `wang_yao_corrected_lrt` | `chapter3-gaussian-highdim.R` | no | pass | fixed seeded matrix call |
| `white_noise_portmanteau_test` | `chapter4-white-noise.R` | no | pass | deterministic sinusoidal series |
| `zhao_chen_wang_spatial_sign_white_noise_test` | `chapter4-white-noise.R` | no | pass | deterministic short series call |

## Reference traceability

| Public topic | R source | Fresh status | Evidence used |
|---|---|---|---|
| `acg_loglik` | `chapter1-foundations.R` | pass | Tyler (1987), DOI 10.1214/aos/1176350263 |
| `classical_lda_classifier` | `chapter5-classical.R` | pass | Anderson (2003), standard Gaussian discriminant formulation |
| `classical_qda_classifier` | `chapter5-classical.R` | pass | Anderson (2003), standard Gaussian discriminant formulation |
| `conditional_alpha_sieve_design` | `chapter4-alpha-fdr-conditional.R` | pass | Ma et al. (2020), DOI 10.1080/07350015.2018.1482758 |
| `conditional_alpha_sieve_fit` | `chapter4-alpha-fdr-conditional.R` | pass | Ma et al. (2020), DOI 10.1080/07350015.2018.1482758 |
| `elliptical_oracle_classifier` | `chapter5-classical.R` | pass | Fang and Anderson (1990); Wakaki (1994) |
| `gaussian_covariance_lrt` | `chapter3-gaussian-classical.R` | pass | Anderson (2003) |
| `gaussian_lda_oracle` | `chapter5-classical.R` | pass | Anderson (2003) |
| `gaussian_mixture_em` | `chapter7-classical-sparse.R` | pass | book manuscript Chapter 7; standard E/M construction, not a new primary method |
| `gaussian_qda_oracle` | `chapter5-classical.R` | pass | Anderson (2003) |
| `hd_classifier_fit` | `chapter5-classical.R` | pass | book manuscript Chapter 5; package infrastructure topic |
| `independence_classifier` | `chapter5-classical.R` | pass | Bickel and Levina (2004) |
| `john_sphericity_test` | `chapter3-gaussian-classical.R` | pass | John (1971) |
| `lloyd_kmeans` | `chapter7-classical-sparse.R` | pass | book manuscript Chapter 7 deterministic contract; MacQueen (1967) labeled foundational |
| `nagao_identity_test` | `chapter3-gaussian-classical.R` | pass | Nagao (1973) |
| `normalize_shape` | `chapter1-foundations.R` | pass | Fang and Anderson (1990) |
| `relliptical` | `chapter1-foundations.R` | pass | Fang and Anderson (1990) |
| `robust_gqda` | `chapter5-gqda.R` | pass | Bose et al. (2015), DOI 10.1016/j.patcog.2015.02.016 |
| `rspherical` | `chapter1-foundations.R` | pass | Fang and Anderson (1990) |
| `sign_whitened_sparse_cca` | `chapter6-sscca.R` | pass | book manuscript Chapter 6; distinct Qian-Liu-Feng method linked as arXiv:2504.13018 |
| `sm_sscm` | `chapter7-spatial-clustering.R` | pass | Zhao-Zhuang-Feng, arXiv:2605.00598 |
| `sparse_k_spatial_median` | `chapter7-spatial-clustering.R` | pass | Zhao-Zhuang-Feng, arXiv:2605.00598 |
| `sparse_kmeans` | `chapter7-classical-sparse.R` | pass | Witten and Tibshirani (2010), DOI 10.1198/jasa.2010.tm09415 |
| `sparse_kmeans_select_s` | `chapter7-classical-sparse.R` | pass | Witten and Tibshirani (2010), DOI 10.1198/jasa.2010.tm09415 |
| `sparse_kmedian` | `chapter7-classical-sparse.R` | pass | book manuscript Chapter 7; explicitly book-defined baseline |
| `sparse_plugin_lda` | `chapter5-linear.R` | pass | Anderson (2003) plus book manuscript Chapter 5; package adapter labeled |
| `sparse_plugin_qda` | `chapter5-sparse-qda.R` | pass | Anderson (2003) plus book manuscript Chapter 5; package adapter labeled |
| `sparse_sm_select_k` | `chapter7-spatial-clustering.R` | pass | Zhao-Zhuang-Feng, arXiv:2605.00598 |
| `sparse_sm_select_tau` | `chapter7-spatial-clustering.R` | pass | Zhao-Zhuang-Feng, arXiv:2605.00598 |
| `spatial_rank` | `chapter1-foundations.R` | pass | Oja (2010), DOI 10.1007/978-1-4419-0468-3 |
| `spatial_rank_covariance` | `chapter1-foundations.R` | pass | Oja (2010), DOI 10.1007/978-1-4419-0468-3 |
| `spatial_sign` | `chapter1-foundations.R` | pass | Oja (2010), DOI 10.1007/978-1-4419-0468-3 |
| `spatial_sign_sphericity_test` | `chapter3-gaussian-classical.R` | pass | Hallin and Paindaveine (2006) |
| `wang_yao_corrected_john_test` | `chapter3-gaussian-highdim.R` | pass | Wang and Yao primary reference already used by the paired correction topic |
| `white_noise_portmanteau_test` | `chapter4-white-noise.R` | pass | Box and Pierce (1970); Ljung and Box (1978) |

## Fresh-only exports recovered by generation

| Export | Owning source group | Direct example | Reference |
|---|---|---|---|
| `generic_weighted_hr_location` | Chapter 2 generic weighted | pass | pass |
| `oracle_weighted_sign_sum_test` | Chapter 2 generic weighted | pass | pass |
| `ec2_covariance` | Chapter 3 Gaussian precision | pass | pass |
| `gaussian_graphical_lasso` | Chapter 3 Gaussian precision | pass | pass |
| `clime_precision` | Chapter 3 Gaussian precision | pass | pass |
| `weighted_spatial_sign_alpha_oracle_test` | Chapter 4 completion | pass | pass |
| `zhao_chen_zi_inst_alpha_test` | Chapter 4 completion | pass | pass |
| `book_gaussian_alpha_cauchy_test` | Chapter 4 completion | pass | pass |
| `conditional_factor_wald_test` | Chapter 4 completion | pass | pass |
| `wang_liu_feng_vector_u_independence_test` | Chapter 4 completion | pass | pass |
| `semc_fit` | Chapter 7 SEMC (excluded from this edit scope) | pass | pass |
| `semc_select_k_gap` | Chapter 7 SEMC (excluded from this edit scope) | pass | pass |

## Reproducible gates

- Fresh staging path used for this audit: `C:/Users/flnankai/AppData/Local/Temp/HDElliptical_api_docs_20260815_175351304/HDElliptical`.
- `Rcpp::compileAttributes()` completed.
- `roxygen2::roxygenise(..., roclets = c("rd", "namespace"), load_code = roxygen2::load_source)` completed.
- All 223 generated Rd files passed `tools::parse_Rd()`.
- The isolated staging package compiled, installed, and loaded successfully under R 4.5.2.
- All 169 generated topics with standard examples executed with zero failures.
- All 27 newly added direct-call fixtures executed individually; `run.donttest = TRUE` exercised both Sparse--SM selection examples; zero failures.
- All 15 documentation-edited R files parsed, and their parsed expressions matched the pre-documentation baseline 15/15, confirming comment-only changes.
- Shared-source hashes matched the fresh staging copy 15/15 at the end of the gate.
- Source-only roxygen emitted in-package unresolved-link ordering diagnostics; every public export still resolved to a generated Rd topic, and these diagnostics caused no generation, parse, install, or example failure.
