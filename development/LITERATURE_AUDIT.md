# Literature and Code Audit

审计快照：2026-08-14 书稿；本文档服务于 `HDElliptical` 的方法溯源、API
规划、数值复核和软件许可决策。这里记录的是正文实际引用和优先实现的核心论文，
不是对高维椭圆分布文献的完整综述。

## 审计范围与引用闭包

本次只读审计使用以下书稿文件：

- PDF：`High-Dimensional-Data-Analysis-for-Elliptically-Symmetric-Distributions_complete_book_latest_arxiv_20260814.pdf`；
- 主源文件：`book-source/book_latest_arxiv_20260814_source/main.tex`；
- 参考文献数据库：`book-source/book_latest_arxiv_20260814_source/references.bib`；
- 编译参考文献：`book-source/book_latest_arxiv_20260814_source/main.bbl`；
- 正文：`book-source/book_latest_arxiv_20260814_source/chapters/ch1_foundations.tex`
  至 `ch7_clustering.tex`；
- 概率附录：`book-source/book_latest_arxiv_20260814_source/chapters/appendix_probability.tex`。

正文逐章去重后的 citation key 数量为：

| 章节 | 去重 citation keys | 主要实现范围 |
|---|---:|---|
| 1 | 12 | spatial sign/rank、SSCM、Kendall、Tyler、HR |
| 2 | 39 | 位置估计与一/两样本检验、max-sum、PDQ、ERHT |
| 3 | 46 | 球形性、比例性、精度矩阵、张量图、因子矩阵 |
| 4 | 44 | alpha/FDR、变点、白噪声、独立性、椭圆拟合优度 |
| 5 | 25 | LDA/QDA、SSLDA、SSQDA |
| 6 | 28 | PCA、因子模型、稀疏 CCA |
| 7 | 29 | 稀疏/稳健聚类、CHIME、IF-PCA、SEMC |

七章正文并集为 191 个 key。概率附录引用 4 个 key，其中 `Tropp2012` 与第 6 章
重叠，因此被编译书稿（七章加概率附录）的并集为 194 个 key；它们全部能在
`references.bib` 中解析，且与 `main.bbl` 的 194 个 `\bibitem` 一致。
`references.bib` 共 200 条，以下 6 条未被当前书稿引用：

```text
FengSun2015Note
FopMurphyScrucca2019
WangFengLiuZhou2021Directional
WittenTibshirani2011
ZhouPanShen2009
ZouYinFengWang2014NPMLECP
```

逐 key 的机器可读证据位于 `development/BIBLIOGRAPHY_TRACEABILITY.csv`，包括章节、
BibTeX/BBL 存在性、Pandoc 标准化作者/年份/题名、条目类型、稳定标识状态和
`METHOD_TRACEABILITY.csv` 行映射。证据快照日期为 2026-08-15；从工作区根目录运行
`Rscript HDElliptical/development/validate-bibliography-traceability.R` 可复核
194/200/6 闭包、UTF-8、严格 CSV 往返及 Pandoc BibTeX/CSV 门禁。

## 第 1 章：方向性基础、位置与散布

| Citation key | 方法或论文 | 稳定标识 |
|---|---|---|
| `Oja2010` | spatial sign/rank、spatial median 的系统性参考 | [DOI](https://doi.org/10.1007/978-1-4419-0468-3) |
| `Tyler1987` | Tyler distribution-free M-scatter / ACG shape | [DOI](https://doi.org/10.1214/aos/1176350263) |
| `HettmanspergerRandles2002` | affine-equivariant HR location/shape | [DOI](https://doi.org/10.1093/biomet/89.4.851) |
| `VisuriOjaKoivunen2000` | sign and rank covariance matrices | 书中期刊元数据 |
| `TaskinenKankainenOja2012` | spatial-sign PCA | [DOI](https://doi.org/10.1016/j.spl.2012.01.001) |
| `HanLiu2018ECA` | elliptical component analysis / multivariate Kendall | [DOI](https://doi.org/10.1080/01621459.2016.1246366), [arXiv:1310.3561](https://arxiv.org/abs/1310.3561) |
| `FanLiuWang2018` | elliptical factor covariance estimation | [DOI](https://doi.org/10.1214/17-AOS1588) |
| `YanFengZhang2025HR` | high-dimensional HR estimator | [arXiv:2505.01669](https://arxiv.org/abs/2505.01669) |

实现优先级：本章是所有后续章节的公共内核。先稳定 `U(0)=0`、中心化、
shape normalization 和退化样本约定，再扩展检验或分类器。`BOOK_ERRATA.md`
记录了 HR 仿射证明中的尺度遗漏；实现不能机械照抄该证明等式。

## 第 2 章：高维位置估计与检验

| Citation key | 方法或论文 | 稳定标识 |
|---|---|---|
| `BaiSaranadasa1996` | common-covariance Bai--Saranadasa test | [期刊原文页](https://www3.stat.sinica.edu.tw/statistica/j6n2/j6n21/j6n21.htm) |
| `ChenQin2010` | heteroscedastic two-sample U-statistic | [arXiv:1002.4547](https://arxiv.org/abs/1002.4547) |
| `SrivastavaKatayamaKano2013` | diagonal, scale-invariant two-sample test | [DOI](https://doi.org/10.1016/j.jmva.2012.08.014), [corrigendum](https://doi.org/10.1016/j.jmva.2013.04.016) |
| `CaiLiuXia2014` | precision-adjusted max test | [DOI](https://doi.org/10.1111/rssb.12034), [作者 PDF](https://faculty.wharton.upenn.edu/wp-content/uploads/2014/06/Two_Sample_Test_of_High_Dimensional_Means_Under_Dependence.pdf) |
| `WangPengLi2015` | high-dimensional nonparametric mean-vector test | [DOI](https://doi.org/10.1080/01621459.2014.988215) |
| `FengSun2016` | one-sample spatial-sign test | [DOI](https://doi.org/10.1214/16-EJS1176) |
| `FengZouWang2016JASA` | two-sample multivariate-sign test | [DOI](https://doi.org/10.1080/01621459.2015.1035380) |
| `LiWangZou2016SimpleTwoSample` | simplified two-sample spatial-sign test | [DOI](https://doi.org/10.1016/j.jmva.2016.04.004) |
| `FengLiuMa2021INST` | inverse-norm sign test | [DOI](https://doi.org/10.1080/07350015.2020.1736084) |
| `FengZhangLiu2020SpatialRank` | high-dimensional two-sample spatial-rank test | [DOI](https://doi.org/10.1016/j.csda.2019.106889) |
| `HuangLiuZhouFeng2023TwoSampleINST` | two-sample INST/tINST | [DOI](https://doi.org/10.1002/cjs.11731) |
| `ZhangFeng2024AdaptiveMean` | adaptive rank-based mean tests | [DOI](https://doi.org/10.1016/j.spl.2024.110226) |
| `LiuZhaoFengWang2025StructuredCorr` | structured-correlation adaptive location test | [DOI](https://doi.org/10.1007/s42952-025-00305-7) |
| `LiuFengZhaoWang2025MaxsumLocation` | spatial-sign max-sum test | [DOI](https://doi.org/10.5705/ss.202024.0051) |
| `YanFengZhang2025InverseNormMaxsum` | inverse-norm weighted max-sum | [arXiv:2501.14168](https://arxiv.org/abs/2501.14168) |
| `ZhaoFeng2026NoteOneSample` | one-sample spatial-sign refinement | [arXiv:2601.08736](https://arxiv.org/abs/2601.08736) |
| `FengWang2026PDQ` | projection/directional-quantile two-sample method | [arXiv:2605.03265](https://arxiv.org/abs/2605.03265) |
| `FengZhouWang2026ERHT` | elliptical regularized Hotelling test | [arXiv:2606.25942](https://arxiv.org/abs/2606.25942) |

实现注意：BS、CQ、SKK、CLX 和 fixed-
\(p\) signed-rank 的书稿公式不能直接转录，具体差异及固定实现口径见
`BOOK_ERRATA.md`。书稿 `ch2_location.tex:1215-1415` 的一般 weighted-HR
方程与 oracle weighted-sign statistic 已形成
`generic_weighted_hr_location()` 和 `oracle_weighted_sign_sum_test()`；后者只在
用户 supplied `null_sd` 或 `nu2`+`trace_R2` 时校准，不把 special-case primary 的
leave-out 方法外推为任意 K 的通用 feasible procedure。

## 第 3 章：矩阵估计与检验

| Citation key | 方法或论文 | 稳定标识 |
|---|---|---|
| `ZouPengFengWang2014Sphericity` | multivariate-sign sphericity test | [DOI](https://doi.org/10.1093/biomet/ast040) |
| `FengLiu2017RankSphericity` | high-dimensional rank sphericity tests | [DOI](https://doi.org/10.1016/j.jmva.2017.01.003), [arXiv:1502.04558](https://arxiv.org/abs/1502.04558) |
| `ZhaoYangZhangFengWang2026AdaptiveSphericity` | adaptive sphericity tests | [DOI](https://doi.org/10.1016/j.jmva.2026.105634) |
| `ChengLiuPengZhangZheng2019SSCM` | equality of two SSCMs | [DOI](https://doi.org/10.1111/sjos.12350) |
| `FengZhangLiu2022Proportionality` | covariance proportionality test | [DOI](https://doi.org/10.1080/24754269.2021.1984373) |
| `OllilaRaninen2019Shrinkage` | optimal covariance shrinkage under elliptical sampling | [DOI](https://doi.org/10.1109/TSP.2019.2908144) |
| `RaninenOllila2022BASIC` | bias-adjusted sign covariance | 书中期刊元数据 |
| `RaninenTylerOllila2022LinearPooling` | linear covariance pooling | 书中期刊元数据 |
| `LuFeng2025Precision` | SCLIME/SGLASSO robust sparse precision | [arXiv:2503.03575](https://arxiv.org/abs/2503.03575) |
| `LiuLuZhouFengWang2025TensorEGM` | tensor elliptical graphic model | [arXiv:2508.00333](https://arxiv.org/abs/2508.00333) |
| `LiuWangZhao2014EC2` | EC2 sparse covariance with eigenvalue constraints | [DOI](https://doi.org/10.1080/10618600.2013.782818) |
| `YuanLin2007` / `FriedmanHastieTibshirani2008` | Gaussian graphical lasso | [Yuan--Lin DOI](https://doi.org/10.1093/biomet/asm018), [graphical-lasso DOI](https://doi.org/10.1093/biostatistics/kxm045) |
| `CaiLiuLuo2011CLIME` | CLIME sparse precision | [DOI](https://doi.org/10.1198/jasa.2011.tm10155) |
| `LiuLuo2015SCIO` | SCIO precision review | [DOI](https://doi.org/10.1016/j.jmva.2014.11.005) |
| `SunZhang2012ScaledLasso` | scaled-lasso precision review | [arXiv:1202.2723](https://arxiv.org/abs/1202.2723), [JMLR 14 (2013)](https://www.jmlr.org/papers/v14/sun13a.html) |
| `XuMaWangFeng2026EllipticalFactor` | SSCM-POET / POET-TME | [arXiv:2512.19325](https://arxiv.org/abs/2512.19325) |

实现状态：`ec2_covariance()`、`gaussian_graphical_lasso()` 和
`clime_precision()` 已分别绑定 EC2 correlation program、off-diagonal
graphical lasso 和 CLIME column programs，并返回可审核的 feasibility/KKT
证书。它们不是 SCLIME/SGLASSO 的别名。Adaptive/MC+ EC2、SCIO 和
scaled-lasso 仍是 review-only，因为其不同的 penalty、tuning 与 solver contract
没有被短 review 唯一确定。

## 第 4 章：其他高维检验

### Alpha 与 FDR

| Citation key | 方法或论文 | 稳定标识 |
|---|---|---|
| `LanFengLuo2018AssetPricing` | high-dimensional linear asset-pricing test | 书中期刊元数据 |
| `FengLanLiuMa2022AlphaSparse` | sparse-alternative alpha test | [DOI](https://doi.org/10.1016/j.jeconom.2021.07.011) |
| `LiuFengMa2023HeavyAlpha` | heavy-tailed alpha test | [DOI](https://doi.org/10.5705/ss.202021.0134) |
| `ZhaoFengWangWang2024RobustAlpha` | robust high-dimensional alpha test | [DOI](https://doi.org/10.1111/obes.70080) |
| `MaFengWang2025DependentAlpha` | dependent-observation alpha test | [arXiv:2401.14052](https://arxiv.org/abs/2401.14052) |
| `ZhaoMaFeng2026LqAlpha` | \(L_q\)-norm alpha test | [arXiv:2603.29764](https://arxiv.org/abs/2603.29764) |
| `WangZhaoFengWang2025MutualFundFDR` | robust mutual-fund selection with FDR | 书中期刊元数据 |
| `ZhaoChenZi2022INSTAlpha` | weighted spatial-sign / inverse-norm alpha tests | [DOI](https://doi.org/10.1002/sta4.490) |
| `LiYang2011ConditionalFactor` / `AngKristensen2012ConditionalFactor` | fixed-dimensional conditional-factor Wald benchmark | 书中期刊元数据 |

### 变点、白噪声和独立性

| Citation key | 方法或论文 | 稳定标识 |
|---|---|---|
| `WangFeng2023JRSSBChangePoint` | data-adaptive high-dimensional change point | [DOI](https://doi.org/10.1093/jrsssb/qkad048), [arXiv:2205.00709](https://arxiv.org/abs/2205.00709) |
| `YuFengZhu2025FunctionalCP` | functional change-point detection | [arXiv:2506.15143](https://arxiv.org/abs/2506.15143) |
| `LiuFengPengWang2025SpatialSignCP` | spatial-sign change-point inference | [arXiv:2504.19306](https://arxiv.org/abs/2504.19306) |
| `WangLiuFeng2025TemporalCP` | temporally dependent change-point inference | [arXiv:2511.01487](https://arxiv.org/abs/2511.01487) |
| `SongWenFeng2026ERHTCP` | ERHT change-point detection | [arXiv:2607.22162](https://arxiv.org/abs/2607.22162) |
| `FengLiuMa2024WhiteNoise` | high-dimensional white-noise test | [DOI](https://doi.org/10.5705/ss.202023.0300) |
| `ChenSongFeng2025RankWhiteNoise` | rank-based white-noise tests | [DOI](https://doi.org/10.5705/ss.202022.0382) |
| `LongJiangLiuXiong2022PanelIndep` | panel max-sum independence | [DOI](https://doi.org/10.1214/21-AOS2142) |
| `WangLiuFengMa2024MutualIndep` | rank max-sum mutual independence | [DOI](https://doi.org/10.1016/j.jeconom.2023.105578) |
| `WangLiuFengMa2024FisherPanel` | Fisher combination for panel independence | [DOI](https://doi.org/10.5705/ss.202023.0348) |
| `WangLiuFeng2026VectorIndep` | independence between high-dimensional vectors | [DOI](https://doi.org/10.1111/sjos.70063) |
| `ZhangFeng2026RadialDirectional` | elliptical GOF by radial--directional dependence | [arXiv:2605.03592](https://arxiv.org/abs/2605.03592) |

实现边界：Chapter 4 completion 的五个入口分别绑定 general-K oracle formula、
Zhao--Chen--Zi inverse-norm feasible endpoint、book-only Gaussian alpha benchmark、
supplied-estimate conditional Wald benchmark，以及 Wang--Liu--Feng exact
D/BKR/tau-star vector-U primary。后者不替代 `WangLiuFengMa2024MutualIndep`
中仍无法从合法来源闭合的 mutual-independence studentisation。

## 第 5 章：稳健稀疏判别

| Citation key | 方法或论文 | 稳定标识 |
|---|---|---|
| `CaiLiu2011LDA` | direct sparse LDA benchmark | 书中期刊元数据 |
| `MaiZouYuan2012` | direct sparse discriminant analysis | 书中期刊元数据 |
| `JiangWangLeng2018` | direct sparse QDA benchmark | 书中期刊元数据 |
| `CaiZhang2021QDA` | convex sparse QDA benchmark | 书中期刊元数据 |
| `ZhuangFeng2025SSLDA` | spatial-sign direct sparse LDA | [arXiv:2504.11117](https://arxiv.org/abs/2504.11117) |
| `ShenFeng2025SSQDA` | spatial-sign direct sparse QDA | [arXiv:2504.11187](https://arxiv.org/abs/2504.11187) |
| `LuFeng2025Precision` | precision plug-in used by robust classifiers | [arXiv:2503.03575](https://arxiv.org/abs/2503.03575) |
| `YanFengZhang2025HR` | HR location/shape plug-in | [arXiv:2505.01669](https://arxiv.org/abs/2505.01669) |

## 第 6 章：PCA、因子模型与稀疏 CCA

| Citation key | 方法或论文 | 稳定标识 |
|---|---|---|
| `TaskinenKankainenOja2012` | spatial-sign PCA | [DOI](https://doi.org/10.1016/j.spl.2012.01.001) |
| `HanLiu2018ECA` | Kendall/ECA PCA | [DOI](https://doi.org/10.1080/01621459.2016.1246366), [arXiv:1310.3561](https://arxiv.org/abs/1310.3561) |
| `ZhaoWangFeng2026SPCA` | high-dimensional spatial-sign PCA | [arXiv:2409.13267](https://arxiv.org/abs/2409.13267) |
| Leyder--Raymaekers--Verdonck (2024) | published generalized spherical PCA primary | [DOI](https://doi.org/10.1007/s11222-024-10413-9) |
| `RaymaekersRousseeuw2019` | foundational generalized spatial-sign covariance/PCA | [DOI](https://doi.org/10.1016/j.jmva.2018.11.010) |
| `WangWangFeng2026GSPCA` | book-attributed high-dimensional extension | 原稿不可取得；不作为实现证据 |
| `FanLiuWang2018` | elliptical factor covariance estimation | [DOI](https://doi.org/10.1214/17-AOS1588) |
| `HeKongYuZhang2022FactorNoMoments` | factor analysis without moment constraints | 书中期刊元数据 |
| `XuMaWangFeng2026EllipticalFactor` | robust matrix/factor reconstruction | [arXiv:2512.19325](https://arxiv.org/abs/2512.19325) |
| `QianLiuFeng2025SparseCCA` | sparse CCA under elliptical symmetry | [arXiv:2504.13018](https://arxiv.org/abs/2504.13018) |

跨章映射：`ch6_pca_factor.tex:385-419` 的 POET 由 Chapter 3
`poet_covariance()` 执行；`ch6_pca_factor.tex:1018-1131` 的 spatial-sign
ER/GR factor-number route 映射到 `elliptical_factor_number()`。
`kendall_factor_number()` 是后续 Kendall MKER/MKTCR 的不同方法，不能用它
掩盖前一个映射，也无需另造 Chapter 6 wrapper。

## 第 7 章：高维稳健聚类

| Citation key | 方法或论文 | 稳定标识 |
|---|---|---|
| `WittenTibshirani2010Clustering` | sparse K-means | [DOI](https://doi.org/10.1198/jasa.2010.tm09415) |
| `CaiMaZhang2019` | CHIME | [DOI](https://doi.org/10.1214/18-AOS1711) |
| `JinWang2016IFPCA` | influential-features PCA clustering | [DOI](https://doi.org/10.1214/15-AOS1423) |
| `JinKeWang2017` | clustering phase transitions | 书中期刊元数据 |
| `ZhaoZhuangFeng2026SparseKSM` | sparse K-spatial-median clustering | [arXiv:2605.00598](https://arxiv.org/abs/2605.00598) |
| `FengZhuang2026SEMC` | semiparametric elliptical mixture clustering | [arXiv:2605.08995](https://arxiv.org/abs/2605.08995) |

SEMC 已从旧 source-blocked 状态转为 implemented。论文公式与官方软件默认值
分层处理：`semc_fit()`/`predict.semc_fit()`/`semc_select_k_gap()` 保留 paper
hard-label Gap-LSE 为默认，同时只通过显式参数开放 software soft-delta、初始化、
KDE/spline、damping、factor/glasso grids 等 contract。六个 Rcpp 内核、固定 official
fixture、公式/RNG/边界测试和两个公开入口的 examples 均可追踪；没有加入论文 simulation。

## 一手代码与许可证

以下链接只指向作者仓库、作者主页或 CRAN；未使用第三方代码聚合站作为代码
来源。

| 方法 | 一手来源 | 技术/许可 | 对本包的决定 |
|---|---|---|---|
| ERHTCP | [flnankai/ERHTCP](https://github.com/flnankai/ERHTCP) | R、Rcpp、RcppArmadillo；GPL-3-or-later | 可作为数值 oracle；MIT 主包不得复制 GPL 源码 |
| SEMC | [flnankai/GEMcluster@10fce04](https://github.com/flnankai/GEMcluster/tree/10fce04fe690fe274dd5d237cfcd3d5c6a4139f6) | Rcpp/RcppArmadillo、C++17；SPDX `MIT`; copyright 2026 Dan Zhuang and Long Feng | 已完成独立重写；只把 pinned 软件作为 contract/fixture oracle，不复制任何可版权源码表达 |
| temporally dependent CP | [flnankai/HDChangePoint](https://github.com/flnankai/HDChangePoint) | R；仓库论文题名/作者与书目不完全一致 | 核实 manuscript 版本后再映射 API |
| CHIME | [作者论文页](https://drjingma.com/papers/cai-chime), [drjingma/gmm](https://github.com/drjingma/gmm) | MATLAB；未发现明确许可证 | 仅作算法和数值 oracle，不复制源码 |
| IF-PCA | [作者软件页](https://www.stat.cmu.edu/~jiashun/Research/software/HCClustering/Tutorial.html) | MATLAB；未发现明确许可证 | 仅作验证 oracle，独立 Rcpp 实现 |
| sparse K-means | [CRAN `sparcl`](https://CRAN.R-project.org/package=sparcl) | GPL-2 | 可选比较依赖或测试 oracle，不复制进 MIT 源树 |
| sign/rank foundations | [CRAN `ICSNP`](https://CRAN.R-project.org/package=ICSNP) | GPL >= 2 | 回归测试 oracle |
| spatial nonparametrics | [CRAN `SpatialNP`](https://CRAN.R-project.org/package=SpatialNP) | GPL-2 | 回归测试 oracle |
| M-scatter | [CRAN `fastM`](https://CRAN.R-project.org/package=fastM) | Rcpp/RcppArmadillo；GPL >= 2 | Tyler/M-scatter 性能与数值 oracle |

### SEMC pinned provenance

- 审计 commit：`10fce04fe690fe274dd5d237cfcd3d5c6a4139f6`。
- 官方仓库 license：SPDX `MIT`；精确版权元数据为 year `2026`，holders
  `Dan Zhuang and Long Feng`。
- `R/chapter7-semc.R`、`src/chapter7_semc.cpp` 与
  `tests/testthat/test-chapter7-semc.R` 是完全独立重写；没有复制、翻译式复制或逐行
  改写任何 GEMcluster 的可版权源码表达，也没有复制 `huge`/GPL 源码。
- 官方仓库只用于确定非唯一的软件 contract 与冻结 fixture；paper-only 与
  software-only 选择均在 API 参数和 diagnostics 中显式区分。

作者账号下另有
[alphaInstability](https://github.com/flnankai/alphaInstability) 和
[ALtest](https://github.com/flnankai/ALtest)，但它们不在当前 194 个正文引用 key 的
闭包内，不能仅因存在仓库就列为本书已覆盖方法。

在已核查的一手论文页、arXiv 页面和作者仓库中，尚未定位到 ECA、SSLDA、
SSQDA、SPCA、SCLIME/SGLASSO、tensor EGM、high-dimensional HR、sparse CCA
和 sparse K-spatial-median 的直接代码链接。“未定位”不等于代码不存在；新增来源
前应重新检查作者主页和论文版本。

## 许可证策略

`HDElliptical` 当前在 `DESCRIPTION` 中声明 `MIT + file LICENSE`。因此采用以下
固定策略：

1. 可以直接整合 MIT/BSD/Apache 兼容代码，但必须保留原版权和许可证文本，并在
   对应源码头和 `inst/COPYRIGHTS`（建立后）中注明来源。
2. GPL 源码只能作为外部依赖、可选 Suggests 或数值 oracle；不得复制、翻译式复制
   或做逐行 Rcpp 改写后放入 MIT 主包。
3. 如果决定直接合并 ERHTCP 或 GPL CRAN 源码，应先把整个包改为 GPL-3-compatible，
   这是项目级决策，不能由单个方法实现隐式触发。
4. 未明确授权的 CHIME/IF-PCA 作者代码只能用于理解输入输出和验证数值结果；实现
   必须从论文公式独立完成，并保留独立测试记录。
5. 论文中的数学算法可以独立实现；测试中比较公开软件输出时，记录软件版本、参数、
   随机种子和容差，不将其源码纳入发布 tarball。
6. SEMC 采用完全独立重写，因此不把上游源码并入发布物；仍记录 pinned commit、MIT
   SPDX 与版权人/年份，以便 `inst/COPYRIGHTS` 和第三方 provenance 审计。

## 书目异常与待核验项

1. **SPCA 重复。** `ZhaoWangFeng2025SSPCA` 与
   `ZhaoWangFeng2026SPCA` 的作者、题名和 arXiv:2409.13267 完全相同，分别在
   第 2/5 章和第 6 章使用，PDF 因此将同一论文列为 2024a/2024b。包内引用应只保留
   一个规范 key；稳定 DOI 为
   [10.48550/arXiv.2409.13267](https://doi.org/10.48550/arXiv.2409.13267)。
2. **作者和题名缺失。** `MaFengWang2025DependentAlpha` 当前少列 Bao Jigang，
   且题名少 “High Dimensional”。公开版本为 Ma, Feng, Wang and Bao，
   *Testing Alpha in High Dimensional Linear Factor Pricing Models with Dependent
   Observations*，[arXiv:2401.14052](https://arxiv.org/abs/2401.14052)。
3. **key 中作者错误。** `YanFengZhang2025InverseNormMaxsum` 的实际作者为 Yan,
   Zhao and Feng，没有 Zhang；建议规范为 `YanZhaoFeng2025InverseNormMaxsum`。
4. **缺年份（已解决）。** `Zhao2023ConditionalAlpha` 当前已有 `year = {2023}`
   字段；追踪表直接读取 BibTeX/Pandoc 元数据，不从 key 推断年份。
5. **key 年份与记录年份不一致。** `XuMaWangFeng2026EllipticalFactor` 的 BibTeX
   年份为 2025；两个 SPCA key 的 BibTeX 年份均为 2024。不要从 key 自动生成
   `citation()` 年份。
6. **“Accepted” 状态过期。** 下列条目已有 DOI，应核对并补充卷期页码或 article
   number：`LiuFengZhaoWang2025MaxsumLocation`、
   `ZhaoYangZhangFengWang2026AdaptiveSphericity`、
   `FengLiuMa2024WhiteNoise`、`WangLiuFengMa2024FisherPanel`、
   `ZhaoFengWangWang2024RobustAlpha`、`WangLiuFeng2026VectorIndep`。
7. **online-first 与卷期年。** ECA、INST、two-sample INST、sphericity、SSCM、
   proportionality 等论文的在线年早于印刷卷期年。包的 citation 数据应分别保存
   `year` 和 `published_online`，不应无条件用 Crossref 年份覆盖卷期年。
8. **新核验的预印本与剩余空缺。** `FengWang2026PDQ` 为
   [arXiv:2605.03265](https://arxiv.org/abs/2605.03265)，
   `MaFengWang2026TimeVaryingAlpha` 为 [arXiv:2604.13772](https://arxiv.org/abs/2604.13772)，
   `ZhaoWang2024ConditionalMaxAlpha` 为 [arXiv:2604.12252](https://arxiv.org/abs/2604.12252)；
   最后一项公开论文年份为 2026，与 key 和 BibTeX 的 2024 不一致。
   `WangWangFeng2026GSPCA` 仅是书内 `@unpublished` 归因；工作区不存在该原稿，
   也未找到 DOI、arXiv 或出版社页面，因此不得以 `local_source` 标记。
9. **字段缺失。** 当前 200 条 BibTeX 均没有 `doi`、`url` 或 `eprint` 字段；JSS
   稿件和包级 `CITATION` 建立前应以本审计中的稳定标识为起点逐项补齐，而不是从
   citation key 猜测。
10. **方法主文献独立证据与元数据纠错。**
    `PRIMARY_METHOD_VERIFICATION.csv` 逐项绑定 52 条独立核验的 primary
    evidence、MTR 行与论文算法位置；另以一行记录不可取得的 GSPCA 书内归因。
    验证器现在要求全部 113 条 `method_primary` 书目行具有可独立审计的稳定标识；
    当前汇总为 122 条 `stable_identifier`、71 条 `book_metadata_only` 和一条
    `no_stable_identifier`。出版社/机构仓储还确认：
    `BosePalSahaRayNayak2015` 的作者应为 Smarajit Bose、Amita Pal、
    Rita SahaRay、Jitadeepa Nayak；`FisherSunGallagher2010` 的第三作者应为
    Colin M. Gallagher；`ParkAyyala2013` 的作者应为 Junyong Park 和
    Deepak Nag Ayyala；`OllilaRaninen2019Shrinkage` 的正式 DOI 是
    [10.1109/TSP.2019.2908144](https://doi.org/10.1109/TSP.2019.2908144)，
    不是包文档曾使用的 `10.1109/TSP.2019.2906691`。

## 维护规则

- 一个公开函数至少绑定一个正文公式/算法位置和一个规范 citation key。
- 书稿与原论文冲突时，以原论文、正式 corrigendum 和可复现数值结果为准，并在
  `BOOK_ERRATA.md` 记录差异。
- “已实现”要求 API、文档、测试和必要的 compiled kernel 同时存在；仅有论文公式
  或外部仓库不算完成。
- 新发现的作者代码必须记录仓库 URL、commit、许可证和是否允许复制；缺少许可证
  时默认不复制。



