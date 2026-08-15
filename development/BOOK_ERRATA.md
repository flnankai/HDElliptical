# Book Errata Relevant to Implementation

状态：已确认；基准书稿为 2026-08-14 的
`book-source/book_latest_arxiv_20260814_source`。本文档只记录会改变算法、统计量、
索引、归一化或软件文档含义的勘误。行号对应这一快照；后续书稿重排时以公式标签
为稳定定位符。

本轮没有修改书稿正文。`HDElliptical` 的实现应遵循下面的“实现决策”，并在测试中
固定相应约定。

本包的目标是提供可直接调用的方法函数，而不是复现各论文的模拟研究。因此勘误的
关闭证据限于公式对照、确定性数值夹具、不变性和退化输入测试；不要求、也不收录
Monte Carlo 场景、功效表或论文模拟数据。

## 概览

| ID | 来源位置 | 问题 | 实现影响 |
|---|---|---|---|
| `CH1-HR-01` | `ch1_foundations.tex:1164-1187` | HR 仿射证明漏掉 \(c^{-1/2}\) | 证明修正；估计方程不变 |
| `CH1-HR-02` | `ch1_foundations.tex:1240` | HR location 步被误称为 Newton update | 应称 Weiszfeld-type/fixed-point update |
| `CH1-TRACE-01` | `ch1_foundations.tex:648-696,726-727,877-916,1099-1138` | unit-trace/shape 命题漏掉无零半径或无 ties 条件 | 补 `P(xi=0)=0`；trace 扣除零方向占比 |
| `CH2-SR-01` | `ch2_location.tex:333-352` | signed-rank Wald 统计量漏 \(1/4\) | 统计量必须乘 \(1/4\) |
| `CH2-BS-CQ-01` | `ch2_location.tex:443-472` | 书中 BS numerator 实为 CQ numerator | BS、CQ 必须作为两个 API 实现 |
| `CH2-SKK-01` | `ch2_location.tex:531-571` | 缺 \(p^{-1/2}\)，\(F_k\) 与 \(c_{p,n}\) 口径有误 | 按原文及 corrigendum 实现 |
| `CH2-CLX-01` | `ch2_location.tex:793-847` | feasible CLX 的标准化描述不完整 | 用 transformed empirical variance |
| `CH2-FZWZ-BF-01` | `ch2_location.tex:670-685` | group-2 bias 项多写一重 \(\gamma\)，且把渐近中心化称为 exact mean | 按原论文单次 \(\gamma\) 与渐近中心化实现 |
| `CH2-CT2-01` | `ch2_location.tex:778-784` | 把 Composite \(T^2\) 误述为 BF/多个检验组合 | common covariance 下按变量相关分块的 quadratic test |
| `CH2-ASPU-01` | `ch2_location.tex:855-925` | 把候选幂称为偶数且把校准概括为 bootstrap/simulation | 原文允许任意正整数幂，并给出 odd/even/max 三组的联合渐近校准 |
| `CH2-FS-01` | `ch2_location.tex:1551-1638` | Feng--Sun 实际检验所需 feasible trace/variance 被省略 | 补原文 leave-two-out trace 与 \(\widehat\sigma_n^2\) |
| `CH2-INST-01` | `ch2_location.tex:1416-1418,1639-1765` | 泛化叙述会误导 INST endpoint 的中心化和 nuisance 角色 | endpoint 始终减 \(\mu_0\)；leave-two-out location 只参与估计 \(D_{ij}\) |
| `CH2-FZW-SIGN-01` | `ch2_location.tex:2003-2126` | \(A_1,A_2\) 的 \(c\)-比率与对角桥矩阵写错，且遗漏 feasible calibration | 按原文 Proposition 2 的一次比率、两组 \(D^{-1/2}\) 桥与三项 trace 实现 |
| `CH2-LWZ-01` | `ch2_location.tex:2127-2241` | 只泛称 feasible plug-in，未列 trace 与 variance 的不同样本量因子 | trace 用 \(n_k(n_k-1)\)，variance 外层用 \(n_k^2\) |
| `CH2-FZL-RANK-01` | `ch2_location.tex:2242-2351` | 省略 feasible traces；原文第二方差项漏 \(p^2\)，且尺度识别与有限样本不变性冲突 | 补 leaveout 校准、恢复 \(p^2\)，显式区分 geometric/paper-trace |
| `CH2-SS-MAX-01` | `ch2_location.tex:1128-1135,1766-1993` | 逆半径矩条件漏标准化，最大特征值条件漏 \(\sqrt{\cdot}\)，且未交代迭代无一般收敛理论 | 按 accepted primary 的 scaled median、MAX 与 Cauchy MAXSUM 实现 |
| `CH2-WMAX-01` | `ch2_location.tex:1079-1085,1234-1243,1893-1996` | weighted MAX/SUM 的来源元数据及原论文若干方向、下标和中心化式有冲突 | 从估计方程推导 general-\(m\) 更新，并按直接 feasible variance 实现 |
| `CH2-GENERIC-WEIGHTED-01` | `ch2_location.tex:1215-1415` | 一般权函数的估计方程和 oracle statistic 可执行，但正文未给任意 K 的 feasible plug-in 校准 | 发布一般 weighted-HR 估计与显式 oracle test；p 值只来自 supplied `null_sd` 或 `nu2`+`trace_R2` |
| `CH2-ZF-RANK-01` | `ch2_location.tex:2412-2451` | 把一般 rank-\(L_q\) 与 min-p 误归给 Zhang--Feng，且遗漏长期方差校准边界 | 只实现论文的 squared-rank sum、max 与 Cauchy；\(\tau^2\) 必须显式给定或显式选 Ouyang--Parzen |
| `CH2-STRUCT-01` | `ch2_location.tex:2438-2447` | Liu--Zhao--Feng--Wang 小节无公式，且 robust elliptical 与 banded 表述没有可核实的 primary 依据 | 合法全文缺失时不发布猜测性实现 |
| `CH2-ZZG-NR-01` | `ch2_location.tex:2479-2480,2595-2618` | 把 motivating oracle \(n\|\bar X\|^2\) 当成可行检验起点 | 实际 API 使用减去 \(\operatorname{tr}(S)\) 的 centered U-statistic |
| `CH2-ZZZ-NRSI-01` | `ch2_location.tex:2622-2640` | 混合 raw-\(L_2\)、F-type 与 primary scale-invariant 统计量 | 按对角标准化统计量及 \(\chi_d^2/d\) 参考实现 |
| `CH2-WX-01` | `ch2_location.tex:2642-2666` | 把 Wang--Xu 的组内半差 Rademacher 校准误述为中心化后的标签重分组 | observed CQ 用完整样本；参考分布用相邻半差与随机符号 |
| `CH2-PDQ-01` | `ch2_location.tex:2684-2689,3028-3077` | radial positivity 被放宽为 \(\xi\ge0\)，且理论只给条件分位数、未定义 finite-\(B\) p 值 | 坚持正半径条件；分别报告 paper critical rule 与显式辅助 MC p 值 |
| `CH2-ZF-STRONG-01` | `ch2_location.tex:3093-3199` | primary 的 \(\kappa_4\)、未定义 \(S\) 及 \(U(0)\) 证明恒等式彼此不自洽 | 不替换符号；pair sum 使用允许零 sign 的一般恒等式 |
| `CH2-ERHT-01` | `ch2_location.tex:3207-3712` | 正式 feasible 公式可核实，但数值段提到的 Bartlett center mapping 未定义 | 只实现直接 feasible center/variance，不猜 Bartlett 修正 |
| `CH3-CLASSICAL-01` | `ch3_matrix.tex:287-405` | Mauchly、John、Nagao 的有限样本乘子有误 | 按 primary 的 Wishart df 与 1/2 因子实现 |
| `CH3-GAUSSIAN-01` | `ch3_matrix.tex:408-690` | covariance divisor、RLZ rule、Wang--Yao attribution 与 Li--Chen SE 标记混杂 | 各 API 按各自 primary 口径实现 |
| `CH3-GAUSSIAN-PREC-01` | `ch3_matrix.tex:446-534` | EC2 的 primary correlation program 被缩写为 covariance display，且 graphical lasso/CLIME solver certificates 与 SCIO/scaled-lasso 边界未闭合 | 实现 EC2 l1、off-diagonal glasso、CLIME；adaptive/MC+ EC2、SCIO、scaled-lasso 保持 review-only |
| `CH3-SPHERICITY-01` | `ch3_matrix.tex:691-1059` | sign bias plug-in 与 adaptive Cauchy 外层定义不完整 | 补 feasible bias；明确 score 与 combined p-value |
| `CH3-PROP-PREC-01` | `ch3_matrix.tex:1060-1463` | proportionality 的 trace normalization 未同步缩放；precision solver 契约被省略 | 恢复 primary scaling，并要求 feasibility/KKT certificates |
| `CH3-TENSOR-01` | `ch3_matrix.tex:1464-1769` | primary 一处 Kronecker 次序与 threshold 绝对值排版不一致 | 依 vec/Kolda objective 与 sign-consistency theorem 实现 |
| `CH3-OLLILA-01` | `ch3_matrix.tex:1770-1841` | 作者代码含未写入理论定义的数值修补 | 实现 published estimator，不复制 clipping/extrapolation/ridge |
| `CH3-HDHR-01` | `ch3_matrix.tex:1842-1968` | 书稿把联立 banded fixed point 改写成 fixed-pilot 后处理 | 实现 primary Algorithm 2；不做第二次 band |
| `CH3-FACTOR-01` | `ch3_matrix.tex:1969-2123` | factor 数、threshold 常数与 precision solver 均非可自动识别量 | API 要求显式给定；不从 simulation tuning 猜默认值 |
| `CH4-COMPLETION-01` | `ch4_other_tests.tex:237-384,538-617,1010-1088,2849-2988` | 三个 formula-complete alpha/Wald benchmarks 及 vector degenerate-U primary 曾被遗漏或误列 review-only | 发布 5 个窄口径 API；oracle/book/primary 范围与 mutual-independence blocked 边界分开 |
| `CH6-XCHAPTER-MAP-01` | `ch6_pca_factor.tex:385-419,1018-1207` | POET 与 spatial-sign factor-number 已在 Chapter 3 实现，但 Chapter 6 ledger 未做跨章映射 | 显式映射 `poet_covariance()` 与 `elliptical_factor_number()`；不与 Kendall MKER/MKTCR 混同 |
| `CH6-GSPCA-01` | `ch6_pca_factor.tex:574-661` | The book conflates the 2024 median/raw-MAD and 2019 h-order cutoffs, while its compact population display suppresses the general radial variable and conditions | Keep cutoff conventions separate, report radial scope, and do not overstate an eigenspace claim |
| `CH6-SSCCA-01` | `ch6_pca_factor.tex:1208-1342` | The book's whitened constrained-PMD summary is not the primary Qian--Liu--Feng metric-lasso/BIC program | Publish separate named APIs; no attribution conflation or hidden whitening repair |
| `CH7-CHIME-ID-01` | `ch7_clustering.tex:145-235` | The model imposes `omega <= 1/2`, but the printed iterate does not state the required post-update relabeling | Canonically exchange labels after initialization and every update, including means, responsibilities and beta orientation |
| `CH7-SPATIAL-CONTRACT-01` | `ch7_clustering.tex:390-568` | Unsquared center fitting, squared/metric assignment, `U(0)`, active thresholds, empty/cycle/final-state and selector rules are not one monotone objective | Expose each contract and never apply an unrequested silent repair |
| `CH7-SEMC-CONTRACT-01` | `ch7_clustering.tex:671-929` | 论文公式与官方软件的初始化、generator、tuning 和 selector contracts 必须分层记录 | `semc_fit()`/prediction/Gap-LSE 已实现；software-only controls 显式标注，三种 shape 均要求 SPD/certificate；旧 source-blocked 结论撤销 |
| `CH5-ELLIP-PRIOR-01` | `ch5_classification.tex:189-243` | A common elliptical generator/shape does not yield the displayed log-prior linear shortcut for unequal priors | Use the exact generator score; the shortcut is restricted to equal priors or Gaussian exponential generators |
| `CH5-QDA-SCALE-01` | `ch5_classification.tex:147-187` | The direct QDA display is twice the canonical Gaussian log-likelihood ratio | Preserve the decision boundary but record `score_scale` explicitly |
| `CH5-DSDA-01` | `ch5_classification.tex:485-521` | DSDA labels/objective are not the primary paper's coding and scale | Use class-1 (-n/n_1), class-2 (+n/n_2), the primary SSE scale, and a separate classification midpoint |
| `CH5-SDAR-ATTR-01` | `ch5_classification.tex:561-590` | The displayed Dantzig QDA constraints are attributed to Jiang rather than Cai--Zhang SDAR | Expose distinct Jiang penalised-loss and Cai--Zhang Dantzig APIs |
| `CH5-RIDGE-01` | `ch5_classification.tex:424-521,824-955` | Practical LPD/SSLDA ridge terms alter the constraint operator but are described like numerical initialisation | Ridge is an explicit fitted tuning parameter; no hidden stabilisation |
| `CH5-SSQDA-01` | `ch5_classification.tex:956-1059` | The triple ordered trace statistic is computationally presented and unequal-prior factors conflict in the primary | Use its (O(np)) sample-variance identity and restrict the first release to equal priors |
| `CH5-HRGQDA-01` | `ch5_classification.tex:1128-1185` | The book's Gaussian plug-in HR-QDA is not the cited primary method | Implement certified HR-GQDA with an explicit (c\in[0,1]) threshold |
| `CH5-GQDA-ALG-01` | `ch5_classification.tex:1128-1185` | The primary sorting shortcut reverses logic for a negative log-determinant contrast and has a separation-endpoint typo | Select (c) by evaluating the original inequality at all legal breakpoints/endpoints |
| `TINST-SRC-01` | `uotwo_cjs_r.tex:251-253,348` | 三处组别索引错误 | leave-one-out 循环按 \(n_k,D_k\) |
| `SLIDES-SSLDA-01` | `slides.tex:1039-1042` | 目标写成 \(V_0^{-1}\mu_d\) | 目标按约束取 \(V_0\mu_d\) |
| `TEX-REF-01` | ch2/ch3/ch4 三处 | `\ref` 被转义为 CR/断行加 `ef` | 不传播损坏交叉引用 |
| `TEX-NEQ-01` | `ch5_classification.tex:213,234` | `\neq` 丢失反斜杠 | 语义按“不等于”处理 |

## CH1-HR-01：HR 仿射证明的尺度因子

**来源。**
`chapters/ch1_foundations.tex:1141-1205`，命题
`prop:hr-affine-equivariance`；重点为
`eq:hr-affine-equivariance-shape` 及源码第 1164--1187 行。

书稿定义

\[
\widehat V_{\mathrm{HR}}^{Y}
  = c B\widehat V_{\mathrm{HR}}B^\top,
\qquad
c=\frac{p}{\operatorname{tr}(B\widehat V_{\mathrm{HR}}B^\top)},
\]

以及

\[
O=(B\widehat V_{\mathrm{HR}}B^\top)^{-1/2}
  B\widehat V_{\mathrm{HR}}^{1/2}.
\]

第 1185--1187 行当前把白化残差直接写成
\(O\widehat V_{\mathrm{HR}}^{-1/2}(X_i-\widehat\mu_{\mathrm{HR}})\)。
正确等式是

\[
(\widehat V_{\mathrm{HR}}^{Y})^{-1/2}
(Y_i-\widehat\mu_{\mathrm{HR}}^{Y})
=c^{-1/2}O\widehat V_{\mathrm{HR}}^{-1/2}
(X_i-\widehat\mu_{\mathrm{HR}}).
\]

因为 \(c>0\) 且空间符号满足 \(U(az)=U(z)\)（\(a>0\)），下一步的符号
等式仍然成立，因此命题结论不受影响。

**实现决策。** `hr_estimator()` 继续在每次 shape 更新后做
`trace(V) = p` normalization，并直接对标准化残差取方向；不人为把
\(c^{-1/2}\) 塞入估计方程。测试必须覆盖一般非正交 \(B\) 下的位置仿射等变、
shape 按 trace 归一后的仿射等变，以及正比例缩放不改变方向。

## CH1-HR-02：location update 不是 Newton 步

**来源。** `chapters/ch1_foundations.tex:1218-1240`，尤其是
`eq:hr-iter-location` 后的第 1240 行。正文把该 location update 称为
“one-step Newton update”。固定当前 shape \(V\) 后，显示的更新实际是
Mahalanobis 空间中位数方程的 Weiszfeld-type/fixed-point 步：

\[
\mu^{+}=\mu+V^{1/2}
\frac{\sum_i U\{V^{-1/2}(X_i-\mu)\}}
     {\sum_i 1/\|V^{-1/2}(X_i-\mu)\|}.
\]

真正的 Newton 步需要估计 score 的矩阵 Jacobian，其中包含
\(\sum_i\{I-U_iU_i^\top\}/\|\varepsilon_i\|\)；它一般不是上式的标量
分母。

**实现决策。** `hr_estimator()` 保留书中显示的 fixed-point/Weiszfeld 更新，
但软件文档和算法说明不把它称为 Newton 方法。收敛以两条估计方程残差判断，
而不是用“Newton step 很小”作为替代。

## CH1-TRACE-01：样本 sign/Kendall 矩阵的 unit-trace 条件

**来源。** `chapters/ch1_foundations.tex:509,648-696,726-727,877-916,
1099-1138`。第 509 行的 sign 椭圆命题明确要求
\(\Pr(\xi=0)=0\)，但后续 SSCM、Tyler shape 恢复和 HR population 命题
没有继承这一条件。这些位置把 SSCM 或 multivariate Kendall matrix 的
trace 直接写成 1，但本书同时采用 `U(0)=0`。因此每个零残差或 tied pair
的外积 trace 为 0，而不是 1。样本 SSCM 的精确关系为
\[
\operatorname{tr}(\widehat S_U)=1-\frac{n_0}{n},
\]
其中 \(n_0\) 是零中心化残差数；样本 Kendall matrix 的精确关系为
\[
\operatorname{tr}(\widehat K)
=1-\frac{N_0}{\binom n2},
\]
其中 \(N_0\) 是 tied pairs 数。总体 unit-trace 结论同样要求中心或相等
事件没有原子质量。若 \(\Pr(\xi=0)>0\)，Tyler 在中心点的二次型比值是
`0/0`，HR 证明中的 \(U(RU)=U\) 也不成立，且 shape 方程的 trace 小于
\(p\)。所以相关 population 命题应显式继承
\(\Pr(\xi=0)=0\)；Kendall 结论另需
\(\Pr(X=\widetilde X)=0\)。

**实现决策。** `sscm()` 返回 `n_zero`，`spatial_kendall()` 返回
`n_pairs` 和 `n_zero_pairs`，帮助页明确上述 trace 关系。实现保留
`U(0)=0`，不删除 tied observations，也不把矩阵强行重新归一到 unit trace。

## CH2-SR-01：one-sample spatial signed-rank 的 \(1/4\)

**来源。** `chapters/ch2_location.tex:332-360`：

- `eq:ch2-lowdim-signed-rank` 定义 \(R_i\) 和 \(\bar R\)；
- `eq:ch2-lowdim-signed-rank-stat` 定义 \(Q_{\mathrm{SR}}\)；
- `eq:ch2-fixedp-sign-chi` 声明卡方极限。

当前书稿写作

\[
Q_{\mathrm{SR}}
=n\bar R^\top\widehat B_R^{-1}\bar R,
\qquad
\widehat B_R=n^{-1}\sum_iR_iR_i^\top.
\]

但 signed-rank 是二阶 U-statistic；其一阶 Hoeffding 投影使
\(\sqrt n\bar R\) 的渐近协方差为 \(4B_R\)，而不是 \(B_R\)。在保持书中
\(\widehat B_R\) 定义不变时，正确 Wald 统计量为

\[
Q_{\mathrm{SR}}
=\frac n4\bar R^\top\widehat B_R^{-1}\bar R.
\]

**实现决策。** 对外函数使用显式的 `n / 4` 版本，并返回所用
`covariance_factor = 4`。不要同时把 \(\widehat B_R\) 乘 4 再把统计量乘
\(1/4\)，否则会重复修正。测试以直接 U-statistic/Hoeffding 投影计算和固定
\(p\) 模拟卡方分位数为准。

## CH2-BS-CQ-01：Bai--Saranadasa 与 Chen--Qin 被写成同一 numerator

**来源。**

- `chapters/ch2_location.tex:443-456`，`eq:ch2-BS-stat`；
- `chapters/ch2_location.tex:458-517`，`eq:ch2-CQ-stat`、
  `eq:ch2-CQ-variance`；
- Bai--Saranadasa 原文：[`BaiSaranadasa1996`](https://www3.stat.sinica.edu.tw/statistica/j6n2/j6n21/j6n21.htm)；
- Chen--Qin：[`ChenQin2010`](https://arxiv.org/abs/1002.4547)，原文 Eq. (2.2)
  也明确复述了 BS numerator。

书中的 BS 式为

\[
\|\bar X_1-\bar X_2\|^2
-\frac{\operatorname{tr}(\widehat\Sigma_1)}{n_1}
-\frac{\operatorname{tr}(\widehat\Sigma_2)}{n_2}.
\]

当 \(\widehat\Sigma_k\) 是分母 \(n_k-1\) 的无偏组内样本协方差时，这个式子
与书中 `eq:ch2-CQ-stat` 的 U-statistic 代数完全相同，因此当前正文把两个方法
合并了。

BS 的原始设定是共同协方差 \(\Sigma_1=\Sigma_2=\Sigma\)，其 numerator 应为

\[
M_{\mathrm{BS}}
=\|\bar X_1-\bar X_2\|^2
-\frac{n_1+n_2}{n_1n_2}\operatorname{tr}(S_p),
\]

其中

\[
S_p=\frac{(n_1-1)S_1+(n_2-1)S_2}{n_1+n_2-2}
\]

是共同协方差的无偏 pooled estimator。CQ 则保留书中 cross-product
U-statistic，并允许 \(\Sigma_1\ne\Sigma_2\)。

**实现决策。** 建立两个独立函数/内部 kernel：

- BS 明确要求 common covariance，使用 pooled covariance 和 BS 自己的方差估计；
- CQ 使用三项 U-statistic 和 heteroscedastic trace-variance estimator；
- 明确记录所有 covariance denominator，禁止依赖 R 函数默认值而不加说明；
- 单元测试验证 CQ numerator 与“分别减无偏组内 trace”的代数恒等式，同时验证
  不平衡样本下它不等于 pooled BS numerator。

## CH2-SKK-01：SKK 的缩放、\(F_k\) 和可行校正

**来源。** `chapters/ch2_location.tex:520-615`，尤其是
`eq:ch2-SKK-hatq`、`eq:ch2-SKK-FG`、`eq:ch2-SKK-varhat` 和
`eq:ch2-SKK-stat`。一手来源为
[`SrivastavaKatayamaKano2013`](https://doi.org/10.1016/j.jmva.2012.08.014)
及其正式 [corrigendum](https://doi.org/10.1016/j.jmva.2013.04.016)。

已确认三处实现级差异：

1. 当前 `eq:ch2-SKK-hatq` 写成 \(\widehat Q-p\)。原文标准化量是
   \[
   \widehat q_{\mathrm{SKK}}
   =\frac{\widehat Q-p}{\sqrt p}.
   \]
2. `eq:ch2-SKK-FG` 的首项必须明确为矩阵平方在 trace 内：
   \[
   \operatorname{tr}\{(\widehat D^{-1}\widehat\Sigma_k)^2\},
   \]
   不能解释成
   \(\{\operatorname{tr}(\widehat D^{-1}\widehat\Sigma_k)\}^2\)。bias
   correction 的分母是 \(n_k-1\)，不是当前书稿中的 \(n_k\)。因此采用
   \[
   \widehat F_k=\frac1p\left[
   \operatorname{tr}\{(\widehat D^{-1}\widehat\Sigma_k)^2\}
   -\frac{\{\operatorname{tr}(\widehat D^{-1}\widehat\Sigma_k)\}^2}{n_k-1}
   \right].
   \]
3. `eq:ch2-SKK-stat` 中的有限样本因子应由样本相关矩阵
   \(\widehat R\) 计算，即
   \[
   \widehat c_{p,n}=1+\frac{\operatorname{tr}(\widehat R^2)}{p^{3/2}},
   \]
   不能在可行统计量中使用未知总体 \(R\)。

**实现决策。** 代码内部对象分别命名为 `Q_raw`、`q_scaled`、`F1`、`F2`、
`G`、`R_hat` 和 `c_hat`，防止把 \(Q-p\) 与 \((Q-p)/\sqrt p\) 混用。
默认实现采用正式 corrigendum 口径，并在返回值中保留这些中间量，便于与论文表格
逐项复核。

## CH2-CLX-01：feasible CLX 不能用 \(\widehat\Omega_{jj}\) 直接标准化

**来源。**

- 书稿 `chapters/ch2_location.tex:787-847`；
- `eq:ch2-CLX-zbar`、`eq:ch2-CLX-stat`、`eq:ch2-CLX-null`；
- Cai--Liu--Xia [作者 PDF](https://faculty.wharton.upenn.edu/wp-content/uploads/2014/06/Two_Sample_Test_of_High_Dimensional_Means_Under_Dependence.pdf)，
  oracle Eq. (2) 及 feasible Eq. (6)--(7)。

书中 oracle 式是正确的：已知 \(\Omega\) 时，\(\Omega X\) 的协方差是
\(\Omega\)，所以 denominator 是 \(\omega_{jj}\)。问题在第 845--847 行的
feasible 描述过于简化；不能仅把所有 \(\Omega\) 替换成 \(\widehat\Omega\)，再用
\(\widehat\Omega_{jj}\) 标准化。

原论文的可行统计量先令

\[
\widehat Z=\widehat\Omega(\bar X-\bar Y),
\]

然后从变换后的两组观测分别计算经验协方差
\(\widehat\Omega^{(1)}\)、\(\widehat\Omega^{(2)}\)，其组内 denominator 是
\(n_1\)、\(n_2\)，并取

\[
\widehat\omega^{(0)}_{jj}
=\frac{n_1}{n_1+n_2}\widehat\omega^{(1)}_{jj}
+\frac{n_2}{n_1+n_2}\widehat\omega^{(2)}_{jj}.
\]

最终

\[
M_{\widehat\Omega}
=\frac{n_1n_2}{n_1+n_2}
 \max_j\frac{\widehat Z_j^2}{\widehat\omega^{(0)}_{jj}}.
\]

用于估计 \(\widehat\Omega\) 的 pooled covariance 按论文 Eq. (6)--(7) 使用
denominator \(n_1+n_2\)，不是 R 的无偏协方差默认 denominator；这是算法定义的
一部分。

**实现决策。** oracle API 接受已知 \(\Omega\) 并用 `diag(Omega)`；feasible API
必须显式形成 `X %*% t(Omega_hat)`、`Y %*% t(Omega_hat)` 的组内经验方差并按上式
合并，绝不以 `diag(Omega_hat)` 代替。测试用直接循环实现复核 Eq. (6)--(7)，并
分别覆盖 oracle 和 estimated-precision 路径。

## CH2-FZWZ-BF-01：Behrens--Fisher bias 与可行 trace 口径

**来源。** `chapters/ch2_location.tex:670-685`，以及 Feng、Zou、Wang、Zhu
(2015) 的 published statistic 与 feasible plug-in 公式。

书稿的 \(b_1\) 第二项写成与 \(\gamma^2\widehat\sigma_{2k}^4\) 成正比；原论文
及其 feasible plug-in 均只有一重 \(\gamma\)：

\[
\widehat b_1=\sum_k\left\{
\frac{2\widehat\sigma_{1k}^4}{n_1(n_1-1)D_k^2}+
\frac{2\gamma\widehat\sigma_{2k}^4}{n_2(n_2-1)D_k^2}
\right\}.
\]

此外，正文称相应展开为 exact mean，但原文在显式主项后保留
\(o\{\sqrt{\operatorname{var}(T_n)}\}\)，所以它是渐近零假设中心化，不是有限样本
精确期望。书稿还没有列出实际检验所需的 feasible trace 公式：两个组内 trace
使用 leave-four-out 边际方差和分母 \(2P_{n_s}^4\)，交叉 trace 每组各删两个观测并
使用分母 \(4P_{n_1}^2P_{n_2}^2\)。

**实现决策。** `feng_zou_wang_zhu_two_sample_test()` 使用一重 \(\gamma\)，把返回量
明确称为 asymptotic null centering，并逐式实现原文的 leave-four-out 与 2+2
leaveout trace 估计器；因此每组至少需要六个观测。任何非正边际组合方差或最终
方差均报错，不以 ridge、绝对值或 floor 修补。

## CH2-CT2-01：Composite \(T^2\) 不是 Behrens--Fisher 组合检验

**来源。** `chapters/ch2_location.tex:778-784`；Feng、Zou、Wang、Zhu (2017)，
doi:10.5705/ss.202015.0199。

原论文的正式名称是 Composite \(T^2\) test，并在两样本部分假设两组共享未知
common covariance \(\Sigma\)。它既不是 unequal-covariance Behrens--Fisher
程序，也不是将 BF、CQ、SKK 等现有检验的 p 值按相关性组合。方法按绝对样本相关
把变量分成小块，对每块使用 Hotelling quadratic form；总和等价于用 block-
diagonal pooled-covariance inverse。

可行 numerator 是论文 Eq. (3.1) 的 2+2 leaveout U-statistic：每个
\(i_1\ne i_2,j_1\ne j_2\) 组合都从两组各删两行，重新估计 pooled covariance、
重新分块并解 block covariance。零假设方差的 trace 按原文“for simplicity”只用
第一组的 leave-four-out estimator，因此有限样本标准化对组标签并不完全对称。
论文 Remark 1 的实用分块规则是 greedy absolute-correlation algorithm，而不是
组合爆炸的全局 subset search。

**实现决策。** 对外函数命名为 `composite_t2_two_sample_test()`，明确
`covariance.model = "common"` 和 first-group trace calibration；默认采用论文
推荐的 \(K=2\) 与 deterministic greedy tie-break。每个 leaveout block 必须严格
正定；不以 ridge、广义逆或方差 floor 改写方法。

## CH2-ASPU-01: the aSPU powers and calibration are misstated

**Source.** `chapters/ch2_location.tex:855-925`; Xu, Lin, Wei and Pan
(2016), doi:10.1093/biomet/asw029.

The sentence before `eq:ch2-Xu-SPU` calls \(\gamma\) an even integer, but the
same section recommends \(\Gamma=\{1,2,3,4,5,6,\infty\}\).  The primary paper
defines the finite candidates for any positive integer \(\gamma\).  Odd powers
are directional and use a two-sided Gaussian calibration; even powers are
nonnegative aggregate statistics and use an upper-tail Gaussian calibration.

There is also a substantive scale-definition mismatch.  The book defines
\(T_{\mathrm{SPU}}(\gamma)\) from inverse-variance-standardized coordinates
\(W_j\), whereas the primary paper and the authors' `highmean` implementation
raise the unstandardized sample-mean differences to the power \(\gamma\).
Consequently, the raw-paper statistic uses the covariance matrix of the mean
difference in its Gaussian moment calculation; the book variant uses the
correlation matrix of the studentized coordinates.  These two finite-sample
procedures are not numerically interchangeable.

The primary paper also derives an analytic calibration for the adaptive test.
Finite odd-power statistics have a joint Gaussian limit, finite even-power
statistics have a separate joint Gaussian limit, and the standardized maximum
has the stated extreme-value limit.  The odd, even, and maximum groups are
asymptotically independent.  With all three groups present, the final
calibration is therefore

\[
  P_{\mathrm{aSPU}}
  =1-\{1-\min(P_O,P_E,P_\infty)\}^{3},
\]

where \(P_O\) and \(P_E\) are multivariate-normal tail probabilities and
\(P_\infty\) is the extreme-value tail probability.  Empirical permutation or
parametric-bootstrap calibration is an optional alternative in the authors'
software, not a requirement for the primary analytic procedure.

**Implementation decision.** The practical package API will expose the
primary-paper raw-difference definition and the book's studentized-coordinate
variant as two explicitly labelled score scales, with the primary-paper scale
as the default.  Both use the analytic odd/even/maximum calibration and expose
the component statistics, moments, covariance or correlation matrices,
covariance regularization, and any positive-semidefinite adjustment.  The
package will not include the source paper's Monte Carlo designs, simulation
tables, or power-study scripts.

## CH2-GENERIC-WEIGHTED-01：一般权函数可执行，但校准边界是 oracle

**来源。** `chapters/ch2_location.tex:1215-1415`；
`FengLiuMa2021INST` 与 `YanFengZhang2025InverseNormMaxsum`。

书稿给出了任意径向权函数 \(K\) 的 weighted location estimating equation、
保留 unweighted HR diagonal equation 的迭代更新，以及
\(T_n(K)\) 和
\(\sigma_{n,K}^2=2\nu_{2,K}^2\operatorname{tr}(R^2)/
\{n(n-1)p^2\}\)。这些式子足以形成一般估计器和 oracle score。
但正文只说明 published feasible procedures 需要各自的 leave-out
replacement；它没有给出一个对任意用户函数 \(K\) 都成立的通用 plug-in 校准。

**实现决策。** `generic_weighted_hr_location()` 逐式实现 weighted location
与 unweighted diagonal update，支持固定 endpoint、power family 和逐标量 R
callback；zero radius、非正 denominator、非正 diagonal 或未收敛均直接失败。
`oracle_weighted_sign_sum_test()` 要求 supplied \(\theta,D\)，且只有显式
`null_sd` 或成对给出的 `nu2`、`trace_R2` 才返回标准化统计量和 p 值；
否则只返回 raw score 并标记 `calibrated = FALSE`。这不是 source-blocked，
也不虚构一般 feasible leave-out estimator。

## CH2-FS-01：Feng--Sun 检验的 feasible 方差被省略

**来源。** `chapters/ch2_location.tex:1551-1638`；Feng 和 Sun (2016)，
doi:10.1214/16-EJS1176。

书稿给出 \(T_{SS}\) 和 oracle \(\sigma_n\)，但没有列出将其变成可调用检验所需的
原论文 feasible calibration。对每个 ordered pair \(i\ne j\)，先从删去 \(i,j\)
的样本联合估计 \((\widehat\theta_{ij},\widehat D_{ij})\)，然后使用

\[
\widehat{\operatorname{tr}(R^2)}=
\frac{p^2}{n(n-1)}\sum_{i\ne j}
\left[
U\{\widehat D_{ij}^{-1/2}(X_i-\widehat\theta_{ij})\}^{\!\top}
U\{\widehat D_{ij}^{-1/2}(X_j-\widehat\theta_{ij})\}
\right]^2,
\]

以及

\[
\widehat\sigma_n^2=
\frac{2\widehat{\operatorname{tr}(R^2)}}{n(n-1)p^2}.
\]

正文的 local-alternative 分母还加入了 refined signal-dependent 项；原论文
Theorem 2 在其条件下把该项视为可忽略并直接按 \(\sigma_n\) 标准化。因此这个
refinement 不能替代零假设下实际 p 值所需的 feasible 方差。

**实现决策。** `feng_sun_one_sample_test()` 的 numerator 只使用 pair-specific
\(\widehat D_{ij}\) 和零假设中心 \(\mu_0\)，而 feasible trace 同时使用
\(\widehat\theta_{ij}\) 与 \(\widehat D_{ij}\)。主 p 值按上述原文零假设方差给出；
不把 local-alternative oracle 项混入校准，也不修补非正方差或未收敛的 leaveout
拟合。

## CH2-FZW-SIGN-01: two-sample sign bridge matrices and feasible calibration

**Source.** `chapters/ch2_location.tex:2003-2126`; Feng, Zou and Wang
(2016), doi:10.1080/01621459.2015.1035380.

The matrices displayed around lines 2042--2046 are not the matrices in the
primary paper.  The correct population quantities are

\[
A_1=\frac{c_2}{c_1}\Sigma_1^{1/2}D_1^{-1/2}D_2^{-1/2}\Sigma_1^{1/2},
\qquad
A_2=\frac{c_1}{c_2}\Sigma_2^{1/2}D_2^{-1/2}D_1^{-1/2}\Sigma_2^{1/2}.
\]

Thus each \(c\)-ratio appears to the first power, not squared, and each
matrix uses a bridge between the two groups' inverse square-root diagonals;
the own-group factor is not \(D_k^{-1}\).

The book section also omits the feasible calibration in Proposition 2.  It
requires the two full-sample diagonal fits, every leave-one-out direction,
\(\widehat c_k\), two ordered within-group trace estimates with the bridge
diagonals, the cross-group trace estimate, and their three variance
contributions.  The resulting \(R_n/\widehat\sigma_n\) is calibrated by the
upper standard-normal tail.  The local-alternative oracle variance is not the
variance used for the reported p-value.  In addition, the original null
theorem uses conditions C1--C4; the book should not import the separate TS5
local-alternative condition into the null calibration.

**Implementation decision.** `feng_zou_wang_two_sample_sign_test()` implements
the crossed leave-one-out numerator and Proposition 2 feasible variance with
the first-power ratios and two-group bridge diagonals.  It exposes every fit,
trace, and variance contribution.  It does not add the optional bootstrap or
the source paper's simulation study, and it does not repair unconverged fits,
non-positive scales, or non-positive variance.

## CH2-INST-01: one-sample INST endpoint centering and feasible variance

**Source.** `chapters/ch2_location.tex:1416-1418,1639-1765`; Feng, Liu and
Ma (2021), doi:10.1080/07350015.2020.1736084, and its official supplement
doi:10.6084/m9.figshare.11914095.v2.

The broad discussion near lines 1416--1418 suggests replacing both
\((\theta,D)\) by a leave-out weighted estimator in every feasible weighted-
sign statistic.  For one-sample INST this is too broad.  For pair \((i,j)\),
the delete-two joint unweighted location/diagonal fit supplies
\(\widehat D_{ij}\), but the two endpoint residuals in the statistic are
\(X_i-\mu_0\) and \(X_j-\mu_0\), not residuals about
\(\widehat\theta_{ij}\).  The leaveout location is used only while estimating
the nuisance diagonal scale.

The directly feasible null variance in Supplement S.3 has the exact factor
\(2n^{-4}\) multiplying its ordered-pair sum.  It must not be replaced by
\(2\{n(n-1)\}^{-2}\), nor by the oracle factorization involving
\(\nu_2,c_0\), and \(\operatorname{tr}(R^2)\).  Those factorized quantities are
useful diagnostics but do not replace the supplement's feasible variance.

**Implementation decision.** `inst_one_sample_test()` keeps every endpoint
centered at `mu`, uses the delete-two location only inside the diagonal fit,
and calibrates with the direct \(2n^{-4}\) variance.  The factorized variance
is returned under diagnostics only.  Zero training residuals, unstable fits,
and non-positive direct variance fail explicitly; no ridge, weight clipping,
absolute-value repair, or variance floor is applied.

## CH2-LWZ-01：Li--Wang--Zou feasible plug-in 的样本量因子

**来源。** `chapters/ch2_location.tex:2127-2241`；Li、Wang、Zou (2016)，
doi:10.1016/j.jmva.2016.04.004。

书稿的 oracle statistic、bias 和 variance 结构与原文一致，但只说各未知量可由
ratio-consistent plug-in 替换，没有列出使方法真正可调用的有限样本公式。原文
组内 squared-trace 估计器对 ordered distinct pairs 求和，分母是
\(n_k(n_k-1)\)；把 trace 代回最终 feasible variance 时，外层因子却是
\(2/(n_k^2p^2)\)，不是 \(2/\{n_k(n_k-1)p^2\}\)。两者不能合并或约掉。
交叉 trace 使用全部 \(n_1n_2\) 个 pairs。科学备择是双侧位置不等，但标准化的
二次型统计量按右尾正态校准。

**实现决策。** `li_wang_zou_two_sample_sign_test()` 使用两组 full-sample
diagonal HR fits，不引入书稿未要求的 leaveout；逐式返回 bias、三个 trace 与
三个 variance 项，并用原文的 \(n_k^2\) 外层因子。非收敛、零半径或非正方差均
明确失败，不做 ridge、绝对值或 floor 修补。

## CH2-FZL-RANK-01：两样本空间秩的 feasible 校准与尺度识别

**来源。** `chapters/ch2_location.tex:2242-2351`；Feng、Zhang、Liu (2020)，
doi:10.1016/j.csda.2019.106889；作者源 arXiv:1506.08315。

书稿只展示 oracle variance，省略原文实际检验所需的 feasible trace 估计器：
两个组内 trace 使用 mutually distinct ordered quadruples 和 group-specific
leave-four-out diagonal fits，交叉 trace 使用组内 ordered pairs 与 pooled
leave-two-out diagonals。作者 TeX 在 feasible variance 的第二组项漏写一重
\(p^2\)；这与 oracle 的组对称性和 trace estimator 的量纲冲突，软件中恢复该
因子。

原文还对每一次 diagonal fit 强制 \(\operatorname{tr}(D)=p\)，同时宣称统计量
具有精确逐坐标尺度不变性。由于 pooled diagonal 会混合不同 leaveout fit，逐 fit
trace normalization 产生不同的全局乘子，有限样本下两项陈述并不相容。单位几何
均值识别则在逐坐标重标度后只产生共同乘子，从而保持论文宣称的不变性。

**实现决策。** `feng_zhang_liu_spatial_rank_test()` 默认
`scale_identification = "geometric"`，并保留 `"paper_trace"` 作为逐式审计
模式；diagnostics 明示选择。两种口径均实现同一 leaveout trace 结构并恢复第二
组 \(p^2\)。正式零假设理论只覆盖 common scatter；unequal-scatter 经验结果不被
写成理论保证。

## CH2-SS-MAX-01：scaled spatial median、MAX 与 MAXSUM 的条件

**来源。** `chapters/ch2_location.tex:1128-1135,1766-1993`；Liu、Feng、
Zhao 和 Wang 的 accepted Statistica Sinica 论文，
doi:10.5705/ss.202024.0051；作者接受稿与 arXiv:2402.01381。

书稿关于逆半径矩的早期条件把 \(E(R^{-k})\) 本身写成一致有界。primary 的条件是

\[
  E\{(R/\sqrt p)^{-k}\}=p^{k/2}E(R^{-k})=O(1),
\]

否则在高维半径通常为 \(\sqrt p\) 的模型下，书稿版本并不是论文使用的标准化条件。
书稿第 1852--1855 行还把渐近独立条件写成
\(\lambda_{\max}(R)\ll p(\log p)^{-1}\)；accepted primary 的条件是

\[
  \lambda_{\max}(R)\ll \sqrt p\, (\log p)^{-1}.
\]

这不是排版上的无关差异：漏掉平方根会显著放宽 MAX 与 SUM 渐近独立所需的相关性
限制。

full-sample nuisance estimator 解

\[
 n^{-1}\sum_iU\{D^{-1/2}(X_i-\theta)\}=0,
 \qquad
 p\,\operatorname{diag}\!left(n^{-1}\sum_iU_iU_i^\top\right)=I_p,
\]

并以 Weiszfeld/diagonal fixed-point 步求解。primary 明确说明该递推尚无一般的收敛、
解存在或唯一性证明；书稿应继承这一限制，而不能把一个数值上稳定的末次迭代自动
称为估计方程的解。MAX 的实际统计量包含
\(\widehat\zeta_1^2=\{n^{-1}\sum_i\widehat r_i^{-1}\}^2\)、
有限样本因子 \(1-n^{-1/2}\) 和 Gumbel 中心化
\(-2\log p+\log\log p\)。MAXSUM 的另一分量是 Feng--Sun 的可行 SUM 检验，
正式组合规则是等权 Cauchy，而不是模拟校准。

**实现决策。** `scaled_spatial_median()`、`spatial_sign_max_test()` 和
`spatial_sign_maxsum_test()` 逐式实现上述估计方程、MAX 和 Cauchy 组合；SUM 分量
复用 `feng_sun_one_sample_test()`。所有迭代返回更新量与方程残差，严格模式下未稳定
即失败；零半径、非正尺度和非正可行方差不得用扰动、ridge 或 floor 修补。Cauchy 尾部
使用 signed-log/`atan2` 计算，避免极端 p 值的 `tan()` 溢出；精确的对向端点
\((0,1)\) 或 \((1,0)\) 没有唯一连续延拓，必须显式报错。软件不包含原论文模拟。

## CH2-WMAX-01：Yan--Zhao--Feng weighted MAX/SUM 的公式口径

**来源。** `chapters/ch2_location.tex:1079-1085,1234-1243,1893-1996`；
arXiv:2501.14168。该预印本作者实际为 Guowei Yan、Ping Zhao、Long Feng；当前
书稿的 citation key/作者写成 `YanFengZhang...`，应改为 Yan--Zhao--Feng。

设 \(e_i=D^{-1/2}(X_i-\theta)\)、\(r_i=\|e_i\|\)、\(U_i=U(e_i)\)，
且 \(m\le1\)。位置与尺度估计方程是

\[
 \sum_i r_i^mU_i=0,
 \qquad
 p\,\operatorname{diag}\!\left(n^{-1}\sum_iU_iU_i^\top\right)=I_p.
\]

因此相应 fixed-point 步必须为

\[
 \theta^+=\theta+D^{1/2}
   \frac{\sum_i r_i^mU_i}{\sum_i r_i^{m-1}}
 =\frac{\sum_i r_i^{m-1}X_i}{\sum_i r_i^{m-1}},
 \qquad
 D^+=pD\,\operatorname{diag}\!\left(n^{-1}\sum_iU_i^2\right).
\]

预印本算法的两个显示步分别误放了一个 \(D^{-1/2}\)，会破坏单位等变性并且不能解
其自身的估计方程；书稿第 1079--1085 和 1234--1243 行给出的方向反而是正确修复。
主文使用 \(\zeta_k=E(R^k)\)，附录一处写成逆幂是 typo。MAX 的径向倍率应为
\(\widehat\zeta_{m-1}^2/\widehat\zeta_{2m}\)。书稿第 1964--1969 行讨论 MAX/SUM
独立性时漏了 MAX 的 \(-2\log p+\log\log p\) 中心化；第 1972--1975 行的版本正确。

SUM 中每个 pair 的第二个方向必须是 \(U_j\)，而不是再次写 \(U_i\)；delete-two
训练集必须同时排除 \(i,j\)。实际 p 值使用直接可行方差

\[
 2n^{-4}\sum_{i\ne j}r_i^{2m}r_j^{2m}
 \{(U_i-\widetilde\mu_{ij})^\top U_j\}
 \{(U_j-\widetilde\mu_{ij})^\top U_i\},
\]

其中 \(\widetilde\mu_{ij}\) 是 delete-two 样本上、仍以零假设中心化的**未加权**
spatial signs 的均值。它不是 fitted location，也不随 \(m\) 改成 weighted mean。

**实现决策。** `weighted_scaled_spatial_median()`、
`yan_zhao_feng_weighted_max_test()` 和
`yan_zhao_feng_weighted_maxsum_test()` 从估计方程实现 general-\(m\) 递推，
而不复制预印本的半幂 typo；MAX 使用样本径向矩，SUM 使用上述直接可行方差。
当 \(m=-1\) 时 SUM 必须与 `inst_one_sample_test()` 数值一致。所有负幂所需的零半径、
未收敛 fit 和非正方差均显式失败，不做权重截断、ridge 或 floor；不收录模拟复现。

## CH2-ZF-RANK-01：Zhang--Feng rank 方法不是一般 \(L_q\) family

**来源。** `chapters/ch2_location.tex:2412-2451`；Zhang 和 Feng (2024)，
doi:10.1016/j.spl.2024.110226，arXiv:2401.00255；SUM 分量来源 Ouyang 等
(2022)，doi:10.1016/j.csda.2022.107495。

primary 只组合三个明确对象：平方 signed-rank/Mann--Whitney 的 SUM 分量、新提出的
MAX 分量，以及二者的等权 Cauchy 组合。它没有提出书稿第 2421--2435 行所述的一般
\(T_q=\sum_j|W_j|^q\) family 或 q-grid；论文结论反而把 multiple-\(L_q\) rank
U-statistics 列为 future work。`min(p_S,p_M)` 只在功效比较中出现，不是论文定义的
另一种正式组合校准；正式 adaptive test 是 Cauchy 组合。

论文定义坐标序列的长期方差

\[
 \tau^2=1+2\sum_{k\ge1}\gamma(k),
\]

但没有提供 estimator、kernel、截断阶数或自动 bandwidth。不能仅凭该论文悄悄选一个
默认值。被引用的 Ouyang 等方法给出 Parzen 加权 autocovariance estimator，但其
`lag = 5` 只是数值研究选择，不是普适默认。Zhang--Feng 文中的实现相关 typo 还包括：
autocovariance denominator 写成 \(p+k\) 而应为 \(p-k\)；相关上界应针对
\(i<j\) 的 \(|\sigma_{ij}|\)，不能包含 \(\sigma_{ii}=1\)；集合条件缺
\(|C_p|/p\to0\) 的箭头；p 值应为 \(1-\Phi(T_S)\) 和 \(1-G(T_M)\)，不是
`1-Phi^{-1}` 或 `1-G^{-1}`。

理论范围也比“重尾稳健均值检验”更窄：一样本要求连续、逐坐标对称的 location family；
两样本要求共同连续分布的纯平移模型。长期方差理论把坐标按给定顺序视为 stationary
strong-mixing sequence，所以任意列重排会改变 HAC 校准；ties/zeros 也不属于原理论。

**实现决策。** `zhang_feng_rank_one_sample_test()` 和
`zhang_feng_rank_two_sample_test()` 只公开 `max`、`sum`、`combined` 三个分量。
`max` 可直接校准；`sum`/`combined` 要求用户给定正的 `tau_sq`，或显式选择
`tau_method = "ouyang_parzen"` 并显式给出 `lag`。实现返回 autocovariances、Parzen
weights 和 \(\widehat\tau^2\)，对非正值直接失败，不 floor；连续性/ties、坐标顺序和
模型边界写入帮助页。测试只做逐式枚举/公式、不变性和退化契约，不复现论文模拟。

## CH2-STRUCT-01：structured-correlation 小节没有可实现公式

**来源。** `chapters/ch2_location.tex:2438-2447`；Liu、Zhao、Feng 和 Wang
(2025)，doi:10.1007/s42952-025-00305-7。合法来源审计覆盖 Springer 正式页面、
KCI 一页预览、RISS/ScienceON 跳转、Crossref、OpenAlex、Semantic Scholar、作者主页
和机构站点。

书稿只有“dense statistic + sparse statistic + combination”的概括，没有定义线性
precision 模型、结构基矩阵、precision estimator、可行 max/sum 统计量或组合 p 值。
第 2440--2441 行把 primary 明示的 linear precision structure 扩写为“linear or
banded”，第 2446--2447 行又称整个构造位于 robust elliptical setting；公开摘要和
附录只能核实原始两样本均值、线性 precision、max statistic、已有 sum statistic 和
Cauchy combination，不能核实 banded、rank、spatial-sign 或一般椭圆稳健构造。

Springer 主文为 subscription content，且没有可用 SharedIt；KCI 合法预览只有标题、
摘要和引言开头，RISS 只回链 KCI，开放获取索引均标记 closed，也没有登记 supplement
或 author manuscript。公开附录虽出现

\[
 \widetilde M_n=\frac{n_1n_2}{n_1+n_2}
   \max_i\widetilde W_i^2,
 \qquad
 \widetilde{\boldsymbol W}
   =\widehat\Omega^{1/2}(\bar X_1-\bar X_2),
\]

及相应 Gumbel 极限，却没有定义实际使用的
\(\widehat{\widehat\Omega}\)、\(y_{n-2}\)、\(T_n\)、\(\mu_0\)、\(\sigma_0\)
或 finite-sample/Cauchy 规则。不能用其引用的 Yang--Zheng--Li sum test 反推这些遗漏。

**实现决策。** 当前不增加公开 API，也不从相邻论文拼接实现。只有在合法取得正式
全文或作者接受稿并逐式核实上述定义后才关闭此项；同时应删除或限定书稿中的
“banded”与“robust elliptical”断言。

## CH2-ZZG-NR-01：ZZG oracle quadratic 不是可行检验统计量

**来源。** `chapters/ch2_location.tex:2479-2480,2595-2618`；Zhang、Zhou 和 Guo
(2022)，doi:10.1007/s00362-021-01270-z。

书稿用
\(Q_n=n\|\bar X\|^2\) 说明 Gaussian oracle 的 weighted-chi-square law，并在
normal-reference 小节称其为 one-sample “starting statistic”。这对解释已知
\(\Sigma\) 的 oracle law 是正确的，但不能作为样本 covariance 未知时的可行 API。
对 null-centered rows \(Z_i=X_i-\mu_0\) 和分母为 \(n-1\) 的样本 covariance \(S\)，
primary 的实际统计量是

\[
 T=n\|\bar Z\|^2-\operatorname{tr}(S)
  =\frac{2}{n-1}\sum_{i<j}Z_i^\top Z_j.
\]

第二个等号显示它是 centered ordered-cross-product/U-statistic，而不是未中心化的
\(n\|\bar Z\|^2\)。两者的均值、前三阶 cumulants 和 matched chi-square 参数不同。

**实现决策。** `zhang_zhou_guo_one_sample_test()` 及 paired/same-unit linear-
hypothesis wrappers 使用上述 centered statistic；oracle quadratic 只作为诊断返回。
trace-power unbiased estimators、三 cumulant matching 和右尾 chi-square calibration
均针对 \(T\)，测试显式锁定 `oracle - trace(S) = U-statistic`，不得把 oracle quantity
重新接到主 p 值。

## CH2-ZZZ-NRSI-01：ZZZ 被混成 raw-\(L_2\) 与 F-type normal reference

**来源。** `chapters/ch2_location.tex:2622-2640`；Zhang、Zhu 和 Zhang (2023)，
doi:10.1080/02664763.2020.1834516，尤其 primary Eqs. (24)--(25)。

书稿先写 raw \(\|\bar X_1-\bar X_2\|^2\)，随后称该文使用 F-type reference，并总结为
“retain the classical \(L_2\) statistic”。这混合了三条不同的 normal-reference 工作线。
该 citation 所对应的 scale-invariant primary 实际定义

\[
 \widehat\Omega_n=\frac{n_2}{n}S_1+\frac{n_1}{n}S_2,
 \qquad \widehat D_n=\operatorname{diag}(\widehat\Omega_n),
\]
\[
 T_{n,p}=\frac{n_1n_2}{np}
 (\bar X_1-\bar X_2)^\top\widehat D_n^{-1}
 (\bar X_1-\bar X_2).
\]

它既不是 raw-\(L_2\)，也不是 F ratio；可行参考为 Welch--Satterthwaite
\(\chi^2_{\widehat d}/\widehat d\)，其中 \(\widehat d=p^2/\widehat q\)，
\(\widehat q\) 由 Eqs. (24)--(25) 的两组 bias-corrected squared traces 和 cross trace
组成。论文的 \(c_{n,p}\) threshold 只调整 df，不改变或 bias-subtract \(T_{n,p}\)。

**实现决策。** `zhang_zhu_zhang_two_sample_test()` 使用 crossed covariance weights、
对角标准化、Eq. (24)--(25) 和右尾 scaled-chi-square calibration；同时返回 adjusted 与
unadjusted df/p/critical。`df_adjustment = "paper"` 只应用论文 threshold rule。软件明确
声明 `raw.L2.normal.reference = FALSE` 和 `F.type.reference = FALSE`，不从名称相近的
raw-\(L_2\) 或后续 F-type 论文借公式。

## CH2-WX-01：Wang--Xu 是半差 Rademacher 校准，不是标签置换

**来源。** `chapters/ch2_location.tex:2642-2666`；Wang 和 Xu (2022)，
doi:10.1093/biomet/asac014，arXiv:2108.01860。

书稿第 2657--2661 行称该方法把“中心化后的观测重新分配到大小为 \(n_1,n_2\) 的
pseudo-groups”。原文算法没有标签重分组。observed statistic 始终是在两组完整原始
数据上计算的 Chen--Qin statistic。参考分布先在每组内部构造

\[
  \widetilde X_{k,i}
  =\frac{X_{k,2i}-X_{k,2i-1}}{2},
  \qquad m_k=\lfloor n_k/2\rfloor,
\]

然后给每个半差乘独立 Rademacher 符号 \(\varepsilon_{k,i}\)，并用 \(m_1,m_2\)
重新计算 CQ 形式的组内与组间 cross products。奇数样本量时最后一个未配对观测仍参与
observed statistic，但不参与 conditional reference。因而最小实用样本量是每组四行，
使 \(m_k\ge2\)。书稿的 \(\pi\) 更适合改成符号向量 \(E\)，且必须说明 observed 与
reference 使用的不是同一数据矩阵。

完整枚举时 conditional tail 为

\[
  p_{\rm exact}=2^{-(m_1+m_2)}
  \sum_E 1\{T_{CQ}(E;\widetilde X)\ge T_{CQ}(X,Y)\}.
\]

因为 \(T(E)=T(-E)\)，软件可固定第一个符号，只计算一半代表而不改变概率。该枚举只是
精确计算 conditional reference；原文明确说明整个 Behrens--Fisher procedure 不是
有限样本 exact randomization test。Monte Carlo 版本使用

\[
  p_{\rm MC}
  =\frac{1+\sum_{b=1}^{B}1\{T_b^*\ge T_{\rm obs}\}}{B+1},
\]

所以必须保留 `>=` tie rule 与加一修正；最小可报告 p 值为 \(1/(B+1)\)。
“arbitrary covariances”表示不限制平均 covariance 的特征结构，并不删除样本量趋于无穷、
四阶矩和非同分布下单个观测不支配等 regularity conditions。

**实现决策。** `wang_xu_approx_randomization_test()` 将 full-sample observed CQ、
相邻半差 pair map、丢弃的奇数末行、符号维数、评估数、exceedance count、是否加一、
随机种子和 Monte Carlo resolution 全部返回。小符号空间可以完整枚举；较大空间只在
用户选择/允许时做方法本身所要求的 Rademacher Monte Carlo 校准。这里的随机化是检验
定义的一部分，不是论文模拟复现；包仍不提供 size/power 场景。测试使用 \(n_1=n_2=4\)
的全符号逐式枚举、固定随机流和边界契约，不以论文的 Monte Carlo 表格作 release gate。

## CH2-PDQ-01：正半径条件与 finite-\(B\) 软件口径

**来源。** `chapters/ch2_location.tex:2670-3077`；
`FengWang2026PDQ` primary manuscript，重点为模型、inverse-radius 条件和 Rademacher
bootstrap theorem。

书稿第 2685--2689 行把 radial variable 写成 \(\xi_{ki}\ge0\)，但 primary 模型使用
正 radial variable；后续理论还显式需要 inverse-radius moments、small-ball control 和正的
population PDQ quantile。若 \(\Pr(\xi_{ki}=0)>0\)，空间方向、\(G_k\) 中的 inverse radius
以及 coordinatewise quantile scale 都可能退化，不能仅用约定 \(U(0)=0\) 把 primary
条件无声放宽。因此软件边界应表述为 \(\xi_{ki}>0\) a.s.，并把样本零 PDQ scale 或
undefined inverse-radius functional 作为显式失败。

primary 的 bootstrap 结论定义条件 \((1-\beta)\) 分位数并使用严格规则
\(T_n^{\rm PDQ}>c_{1-\beta}^*\)，但没有规定有限 \(B\) 时的 quantile type、ties 或可作为
`htest$p.value` 的 Monte Carlo 公式。为使软件结果可复现，当前约定是

\[
 c_{1-\beta,B}^*=T^*_{(\lceil(1-\beta)B\rceil)},
 \qquad
 p_{\rm MC}=\frac{1+\#\{T_b^*\ge T_n^{\rm PDQ}\}}{B+1}.
\]

前者是 type-1 inverse empirical CDF，并继续使用 paper 的严格 `>`；后者是单独标注的
plus-one 辅助 p 值，ties 进入 `>=` 上尾。有限 \(B\) 时两种 rejection 可以不同，不能
把辅助 p 值反称为论文指定的 finite-\(B\) p 值。

**实现决策。** `feng_wang_pdq_two_sample_test()` 同时返回 primary critical decision 与
辅助 MC p decision、critical rank、exceedance count、seed 和 Monte Carlo resolution；
主文档明确 `paper.finite.B.p.value.specified = FALSE`。固定种子测试逐式复核 Rademacher
统计量与两种 decision。这是方法内生的 bootstrap calibration，不是论文 size/power
模拟复现。

## CH2-ZF-STRONG-01：strong-correlation note 的符号与零 sign 边界

**来源。** `chapters/ch2_location.tex:3084-3199`；Zhao 和 Feng (2026)，
arXiv:2601.08736。

primary Eq. (2.7) 把

\[
 \kappa_4=
 \frac{\mathbb E(U_1^\top U_2)^4}
      {\{\mathbb E(U_1^\top U_2)^2\}^2}
\]

定义为 inner-product fourth-moment kurtosis；Assumption 3.3 却写了一个 trace-fourth
ratio。两者不能在没有说明时互换。Assumption 3.2 的 displayed condition 又使用未在
该处定义的矩阵 \(S\)，不能擅自认作 \(\Sigma_U\)、shape 或 SSCM。书稿重述 theorem
时应保留这一 primary-source ambiguity，而不是把其中一个版本宣称为已解决条件。

另一个边界冲突来自 \(U(0)=0\)。书稿第 3097--3100 行以
\(\|U_i\|=1\) 写成

\[
 S_n=\frac12\left\|\sum_iU_i\right\|^2-\frac n2,
\]

但 primary 同时定义 \(U(0)=0\)，且 bootstrap 的 fitted center 可能恰好命中观测。
始终成立的恒等式是

\[
 S_n=\frac12\left\{
   \left\|\sum_iU_i\right\|^2-\sum_i\|U_i\|^2
 \right\}.
\]

## Chapter 4 implementation errata (audited 2026-08-15)

### CH4-ALPHA-01: regression, robust-sum, and robust-max calibrations

The GRS slopes must come from the joint intercept-plus-factor regression, and
the exact statistic uses the MLE residual covariance divisor. The current book
form loses the matching divisor factor when it writes an unbiased covariance.
The Pesaran--Yamagata sum, Liu--Feng--Ma sign statistic, and
Zhao--Feng--Wang--Wang robust maximum also require their primary feasible
centres, trace estimates, projected residuals, and radial-moment correction;
the shorter book displays do not uniquely supply those quantities.

**Implementation decision.** The six unconditional-alpha APIs use the primary
regression projections and feasible calibrations. Unsupported weighted/INST,
dependent, and generic Lq review paragraphs remain review-only.

### CH4-CONDITIONAL-01: conditional residual and combination definitions

The conditional tests use null-restricted or explicitly supplied nuisance
residuals. The feasible sum traces, CSS split residuals, CSM radial correction,
and primary Fisher or truncated-Cauchy combinations cannot be replaced by the
oracle summaries in the review prose.

**Implementation decision.** Sieve construction, nuisance fits, factor order,
and every tuning value are explicit. The seven APIs reject incompatible
supplied residual/factor objects instead of silently ignoring an argument.

### CH4-CP-01: DMS, spatial-sign, and ERHT change-point pivots

The DMS dense component sums the gamma=1/2 CUSUM squares over the full split
grid; its max pivots and weighted A/D arguments differ from the current book
displays. The spatial-sign L2 pivots divide by the square root of the feasible
trace and use endpoint nuisance fits. ERHT's analytic Cauchy aggregation of
dependent ridge p-values is a practical combination, not an exact joint-limit
calibration without the unknown cross-ridge covariance.

**Implementation decision.** The package exposes corrected DMS/sign/ERHT
components and diagnostics. WBS requires explicit ridge, threshold, and
geometry tuning; functional/temporal prose without an executable statistic
remains review-only.

### CH4-WN-01: feasible white-noise statistics

The Feng--Liu--Ma feasible sum is an off-diagonal inner-product U-statistic,
not the naive Frobenius sum in the book. The spatial-sign variance contains a
single lag-count factor rather than its square. The rank primary establishes
max tests; rank sums and their advertised combinations are future-work or
book-only prose.

**Implementation decision.** Only the primary-backed portmanteau, FLM
SUM/MAX, spatial-sign, and Spearman/Kendall maximum APIs are exposed.

### CH4-COMPLETION-01: formula-complete benchmarks and exact vector U-kernels

**Source.** `chapters/ch4_other_tests.tex:237-384,538-617,1010-1088,2849-2988`;
`ZhaoChenZi2022INSTAlpha`, `LiYang2011ConditionalFactor`,
`AngKristensen2012ConditionalFactor`, and
`WangLiuFeng2026VectorIndep`.

Four previously omitted formula blocks have executable but deliberately narrow
contracts. The general weighted-alpha display is an oracle supplied-score
class, whereas inverse-norm weighting is the primary feasible endpoint. The
Gaussian alpha Cauchy display leaves \(\operatorname{tr}(R^2)\) abstract and
must require it. The low-dimensional conditional Wald display supplies
\(\hat\delta\) and its covariance rather than defining a new nuisance
estimator. Finally, the vector-independence primary gives exact Hoeffding D,
BKR R, and tau-star examples, but those kernels do not fill the separate
mutual-independence studentisation gap.

**Implementation decision.** The five public entries are
`weighted_spatial_sign_alpha_oracle_test()`,
`zhao_chen_zi_inst_alpha_test()`,
`book_gaussian_alpha_cauchy_test()`,
`conditional_factor_wald_test()`, and
`wang_liu_feng_vector_u_independence_test()`. The first records oracle scope;
the second has no radial floor or weight cap; the third requires supplied
`trace_R2`; the fourth uses strict Cholesky; and the fifth exactly enumerates
the order-5/order-6/order-4 symmetrized U-kernels with an explicit workload
guard. Their tests are formula and boundary fixtures with intrinsic
permutations only, not paper simulations.

### CH4-INDEP-01: panel and vector independence methods were conflated

The Long/Feng--Jiang--Liu--Xiong panel test assumes iid time errors and uses a
centred squared-correlation sum plus a maximum. The Wang--Liu--Feng--Ma serial
test instead uses a signed-correlation sum, thresholded temporal covariance,
and Fisher combination. Both panel Gumbel laws in their 4 log N
parameterisation have the 1/sqrt(8 pi) constant, not 1/sqrt(pi). Some vector
summation limits printed as p must be q.

**Implementation decision.** Separate public APIs implement the two panel
models. `wang_liu_feng_vector_independence_test()` exposes the primary-backed
Spearman and Kendall families, while
`wang_liu_feng_vector_u_independence_test()` implements the same primary's
exact Hoeffding D, BKR R and tau-star examples with a deterministic workload
guard and intrinsic permutation variance. The legally inaccessible
mutual-independence studentisation remains blocked; neither vector API is
substituted for it.

All randomization, wild-bootstrap, permutation, or Gaussian-process code in
Chapter 4 is intrinsic to the corresponding callable method's calibration.
It is not a reproduction of paper simulations, size/power studies, tables, or
data analyses.
只有所有 signs 非零时第二项才等于 \(n\)。

**实现决策。** `zhao_feng_strongcorr_sign_test()` 直接计算/等价使用一般 pair-sum
恒等式，零 sign 原样保留；不估计存在符号冲突的 \(\kappa_4\)，也不解释未定义的
\(S\)。observed signs 以 null center 形成，bootstrap signs 以普通 sample spatial median
形成；共同 \(\sqrt\tau\sqrt{\binom n2}\) 因子在比较中抵消，不虚构 \(\tau\) plug-in。

## CH2-ERHT-01：ERHT 只实现正式 feasible 公式

**来源。** `chapters/ch2_location.tex:3207-3712`；Feng、Zhou 和 Wang (2026)，
arXiv:2606.25942。

可核实的 fixed-ridge raw statistic 是

\[
 T_n(\rho)=n(\widehat\theta-\theta_0)^\top
 (\widehat R_n+\rho I)^{-1}(\widehat\theta-\theta_0),
 \qquad \rho>0,
\]

其中 \(\widehat\theta\) 为 sample spatial median，
\(\widehat Y_i=\sqrt p\,U(X_i-\widehat\theta)\) 且
\(\widehat R_n=n^{-1}\sum_i\widehat Y_i\widehat Y_i^\top\)。把
\((w_i,A_n)\) 换成 \((\widehat w_i,\widehat A_n)\) 后，primary 完整定义了
\(\widehat\mu_n(\rho)\) 和 \(\widehat\sigma_{D,n}^2(\rho)\)，正式可行统计量为

\[
 \widehat Z_n(\rho)=
 \frac{T_n(\rho)-n\widehat\mu_n(\rho)}
      {\{n\widehat\sigma_{D,n}^2(\rho)\}^{1/2}}.
\]

论文数值段另提 Bartlett center correction 和 \(a=8n/p\)，但正式文章与公开 arXiv
均没有给出把 \(a\) 映射到 feasible center 的公式。仅凭常数无法判断它是乘法、加法、
df 变换还是作用于哪个中心项，任何实现都会是猜测。

**实现决策。** `elliptical_regularized_hotelling_test()` 和固定-grid Cauchy wrapper
只实现上述正式 companion-matrix center/variance equations；诊断明确记录 Bartlett
correction 未应用及原因。测试以独立 full-matrix/SVD reference 核对 raw statistic、
center、variance 和 \(Z\)，不把论文数值程序中未定义的 modification 当作 release
目标，也不复现其模拟表格。

## CH3-CLASSICAL-01: classical sphericity finite-sample factors

**Source.** `chapters/ch3_matrix.tex:287-405`; Mauchly (1940), John
(1971, 1972), Nagao (1973), and Hallin--Paindaveine (2006).

For unknown mean, the Wishart degrees of freedom are \(m=n-1\). The
Mauchly correction denominator is \(6p\), not \(6(p+1)\); John's
fixed-dimensional pivot is \(mpU/2\), not \(m(p+2)U/2\); and Nagao's
identity statistic is \(m\operatorname{tr}(S-I)^2/2\), so the displayed
book formula is missing a factor \(1/2\). The general
Hallin--Paindaveine statistic is a radial-rank-score weighted sign statistic;
the unweighted sign test is only one special score choice.

**Implementation decision.** `mauchly_sphericity_test()`,
`john_sphericity_test()`, and `nagao_identity_test()` use the primary
factors. `hallin_paindaveine_shape_test()` exposes sign, Wilcoxon,
Spearman, and van der Waerden scores rather than relabelling the SSCM-only
special case as the full rank procedure.

## CH3-GAUSSIAN-01: Gaussian benchmark conventions are source-specific

**Source.** `chapters/ch3_matrix.tex:408-690` and the primary articles
cited in that section.

Bickel--Levina, Rothman--Levina--Zhu, and Cai--Liu use the centred covariance
with divisor \(n\) in their stated estimators. The fourth RLZ rule is
adaptive lasso, not MCP. The corrected high-dimensional sphericity LRT in
this section is the Wang--Yao extension; it should not be attributed to the
different Bai--Jiang--Yao--Zheng given-covariance statistic. In Li--Chen, the
displayed feasible null quantity used in rejection is a standard deviation
even where one display labels it as a variance.

**Implementation decision.** The thresholding APIs default to the primary
divisor and expose any alternate divisor explicitly. The RLZ API implements
hard, soft, SCAD, and adaptive-lasso thresholding. Wang--Yao and Li--Chen
calibrations use their own primary centring and standard-error definitions;
non-positive estimated calibration quantities fail rather than being
absoluted or floored.

## CH3-GAUSSIAN-PREC-01: EC2, Gaussian glasso, and CLIME need distinct contracts

**Source.** `chapters/ch3_matrix.tex:446-534`;
`LiuWangZhao2014EC2`, `YuanLin2007`,
`FriedmanHastieTibshirani2008`, `CaiLiuLuo2011CLIME`,
`LiuLuo2015SCIO`, and `SunZhang2012ScaledLasso`.

The primary convex EC2 branch operates on the sample correlation matrix with
unit diagonal and an explicit minimum-eigenvalue constraint; the shortened
book display looks like a direct covariance program. Graphical lasso uses an
off-diagonal penalty and must certify SPD, descent and the full KKT system.
Primary CLIME certifies the raw column programs and then applies the
smaller-absolute-value symmetrisation; that final symmetrisation is not
automatically SPD and must not be repaired or advertised as such. The short
SCIO and scaled-lasso review does not close their different programs, tuning,
or certificate semantics.

**Implementation decision.** `ec2_covariance()`,
`gaussian_graphical_lasso()`, and `clime_precision()` expose these three
distinct programs through `cpp_ch3gp_ec2_l1`,
`cpp_ch3gp_offdiag_glasso`, and `cpp_ch3gp_clime`. An invalid solver result
has no estimate; its last iterate is diagnostic only. Adaptive/MC+ EC2, SCIO,
and scaled-lasso remain review-only. No ridge, eigenvalue-floor repair,
pseudoinverse, or tolerance relaxation is hidden.

## CH3-SPHERICITY-01: feasible sign bias and adaptive Cauchy definition

**Source.** `chapters/ch3_matrix.tex:691-1059`; Zou--Peng--Feng--Wang,
Feng--Liu, and Zhao--Yang--Zhang--Feng--Wang.

The book invokes an estimated sign-test bias without defining a finite-sample
plug-in. The adaptive display also shows the inner truncated-Cauchy score but
does not consistently distinguish it from the final survival probability;
one indicator is typeset inside a tangent although the primary construction
is symmetric with the indicator outside.

**Implementation decision.** The sign sphericity API exposes the
translation-invariant residual-moment plug-in used by the adaptive primary,
plus literal second-order and normal-limit audit choices. The adaptive API
returns both the Cauchy score and `p.value = 1 - F_C(score)`, with the
indicator outside each tangent. Rank ordered/unordered factors were checked
against the primary and retained.

## CH3-PROP-PREC-01: proportionality scaling and precision certificates

**Source.** `chapters/ch3_matrix.tex:1060-1463`; Cheng et al., Feng et al.,
and Lu--Feng.

The Feng proportionality primary uses a trace-one shape. Rewriting it with
the book's trace-\(p\) convention requires the corresponding powers of
\(p\); the current displays do not apply that conversion consistently.
The precision review also omits that SCLIME is solved columnwise and then
symmetrised by the entry with smaller absolute value, while SGLASSO penalises
the full matrix including its diagonal.

**Implementation decision.** `feng_spatial_rank_proportionality_test()`
uses the primary trace-one scale and exact feasible variance.
`spatial_sign_precision()` accepts only solutions passing explicit
feasibility, stationarity/duality, symmetry, positive-definiteness, and KKT
checks. `threshold_spatial_sign_precision()` uses the primary
greater-than-or-equal absolute threshold. No ridge, pseudoinverse, or
post-hoc constraint relaxation is hidden.

## CH3-TENSOR-01: tensor ordering and threshold sign

**Source.** `chapters/ch3_matrix.tex:1464-1769`; Liu et al.

One primary model display orders the Kronecker factors as
\(\Sigma_1\otimes\cdots\otimes\Sigma_K\), while its column-major
vectorisation, Kolda unfolding, whitening map, and objective all require the
reverse \(K,\ldots,1\) order. A threshold display also omits absolute
values even though the proof allows negative edges and proves sign recovery.
The book's ordering and absolute threshold are the coherent versions.

**Implementation decision.** Tensor APIs use R column-major vectorisation,
Kolda mode unfolding, reverse Kronecker order, off-diagonal penalty
\(p_k\lambda\), and `abs(omega_ij) >= tau`. Mode precisions are
Frobenius-normalised as an identification convention and must pass
SPD/descent/KKT certificates.

## CH3-OLLILA-01: published estimators versus reference-code repairs

**Source.** `chapters/ch3_matrix.tex:1770-1841`; Ollila--Raninen,
Raninen--Ollila, and Raninen--Tyler--Ollila.

The public MATLAB code includes pragmatic operations not specified by the
estimators: moving kurtosis inside its lower bound, deleting zero residuals,
flooring small radii, extrapolating a BASIC lookup table, and adding a small
identity term to pooling. The real-valued BASIC map also contains a factor
\(1/2\), which is easy to miss typographically and is required by the
\(\lambda=1\) endpoint.

**Implementation decision.** The package implements the published moments,
integral map, bounded inversion, and pooling problem directly. It reports
theoretical-boundary equality and fails on undefined zero radii or an
unbracketed inverse; it does not copy undocumented clipping, spline
extrapolation, or identity regularisation.

## CH3-HDHR-01: HDHR is a joint banded fixed point

**Source.** `chapters/ch3_matrix.tex:1842-1968`; Yan--Feng--Zhang,
Algorithm 2.

The primary algorithm jointly updates location and shape. Each iteration
standardises by the current shape, bands that sign SSCM, applies the
congruence update, and trace-normalises. The book instead describes a
fixed-pilot score/raw-SSCM post-processing path and then bands the inverse a
second time; those are not Algorithm 2. The book's HR-centred divisor-\(n\)
scale and stated \(r_n+h^{-\alpha}\) rate are also not the executable
primary construction.

**Implementation decision.** `high_dimensional_hr()` implements the joint
fixed point, performs one banding operation per iteration, and directly
inverts the final trace-\(p\) structured shape. The optional covariance
scale is separately labelled as the primary QDA add-on. Nonconvergence,
zero radii, non-PD banded maps, or invalid pilots return an explicit failure
under non-strict mode and never a repaired estimate.

## CH3-FACTOR-01: factor and tuning quantities remain explicit inputs

**Source.** `chapters/ch3_matrix.tex:1969-2123`; Xu--Ma--Wang--Feng.

The factor upper bound, POET threshold constant, number of factors, and
CLIME/GLASSO tuning parameter are theoretical or user-selected quantities;
the primary does not give a universally observable automatic rule. Values
used in a numerical experiment are not method definitions.

**Implementation decision.** `elliptical_factor_number()` requires an
explicit factor cap; POET APIs require a factor count and either a threshold
or an explicit multiplier of the primary rate; precision fitting requires an
explicit penalty and a certified solver result. The Tyler refinement is the
single primary one-step map, not a newly invented iteration. No simulation
tuning grid or reproduction data is shipped.

## CH5-CLASSIFIER-01: audited classification formula and attribution corrections

**Source.** `chapters/ch5_classification.tex:111-1185`, with the primary
papers cited by the individual package help pages.

The executable Chapter 5 rules require the following distinctions.

- Under a general common elliptical generator and shape, unequal class priors
  do not reduce to the displayed midpoint linear score plus a Gaussian
  log-prior offset. The exact generator likelihood is used; the shortcut is
  allowed only for equal priors or the Gaussian exponential generator.
- The direct Gaussian QDA display is twice the canonical log-likelihood ratio.
  This leaves classes unchanged but not score values. The fixed-dimensional
  LDA covariance contribution also contains both
  `(a' Omega delta)^2` and `(a' Omega a)(delta' Omega delta)`.
- DSDA uses the primary asymmetric response coding and objective scaling.
  The Dantzig equations attributed to Jiang are instead Cai--Zhang SDAR;
  Jiang's method uses separate penalised quadratic and linear losses.
- Any practical LPD/SSLDA ridge changes the fitted constraint operator and is
  therefore exposed and recorded as tuning, never as silent initialisation.
  The SSQDA triple ordered trace equals the ordinary unbiased within-class
  total variation and is evaluated in `O(np)`.
- The primary unequal-prior terms in SDAR/SSQDA are internally inconsistent,
  so the first package interfaces enforce equal priors rather than choosing an
  undocumented factor.
- The cited high-dimensional HR method is GQDA, not the book's ordinary
  Gaussian plug-in QDA. Its threshold `c` is selected by evaluating the
  original inequality at all legal breakpoints, endpoints, and interval
  midpoints; this remains correct for positive, negative, or zero
  log-determinant contrast.

**Implementation decision.** Every classifier stores its score orientation and
scale. Solvers must certify feasibility/KKT and SPD or positive determinant
sign where required; an invalid fit cannot be used for prediction. Operator-
norm covariance consistency alone is not advertised as high-dimensional QDA
log-determinant consistency without uniform eigenvalue and log-determinant
control. Review-only method mentions do not become guessed optimisers.

## Chapter 6 implementation errata (audited 2026-08-15)

### CH6-XCHAPTER-MAP-01: POET and spatial-sign factor number are Chapter 3 reuse

**Source.** `chapters/ch6_pca_factor.tex:385-419,1018-1207`;
`FanLiaoMincheva2013`, `XuMaWangFeng2026EllipticalFactor`, and
`HeKongYuZhang2022FactorNoMoments`.

Chapter 6 explicitly says that Chapter 3 already treated the covariance side
of elliptical factor models, but the old ledger made that reuse look like
missing Chapter 6 code. It also listed only the Kendall factor-number API even
though the earlier sign eigenratio/growth-ratio formula is the distinct
spatial-sign selector implemented in Chapter 3.

**Implementation decision.** `poet_covariance()` is the executable
cross-chapter backend for lines 385--419, and
`elliptical_factor_number()` is the spatial-sign ER/GR mapping for the
sign-based factor-number discussion. `kendall_factor_number()` remains the
separate Kendall MKER/MKTCR implementation. The traceability CSV points these
Chapter 6 rows to their actual Chapter 3 tests, Rd topics and native kernels;
no duplicate wrapper is invented.

### CH6-GSPCA-01: generalized-sign cutoff families and population scope

**Source.** `chapters/ch6_pca_factor.tex:574-661` and the distinct 2019 and
2024 primary methods cited there.

The chapter places the 2024 median/raw-MAD cutoff and the 2019 h-order cutoff
inside one generalized-sign presentation. They are different finite-sample
conventions, not interchangeable constants. The compact population display
also suppresses the general radial variable and the conditions under which a
transformed scatter shares the target eigenspace; it cannot support an
unqualified eigenspace claim for every radial transform.

**Implementation decision.** `generalized_sign_pca()` names and records the
cutoff convention and radial population scope. Median/raw-MAD, h-order and
supplied cutoffs remain separate, with endpoint, tie and zero rules tested.
No simulation tuning grid or paper-specific factor result is used to choose
among them.

### CH6-SSCCA-01: primary metric lasso versus book whitened PMD

**Source.** `chapters/ch6_pca_factor.tex:1208-1342` and the cited
Qian--Liu--Feng primary program.

The book summarizes a whitened constrained-PMD construction, whereas the
primary method is a metric-lasso program with p-scaled blocks, generally
non-unit metric diagonals and BIC selection. Treating the two as the same
algorithm changes both the objective and its tuning interpretation.

**Implementation decision.** `sscca()` implements the primary metric-lasso/BIC
contract and returns KKT and metric diagnostics. The distinct
`sign_whitened_sparse_cca()` implements the book construction with an explicit
ridge. Neither API silently whitens, substitutes a pseudoinverse, or guesses a
paper simulation tuning grid.

## Chapter 7 implementation errata and source boundary (audited 2026-08-15)

### CH7-CHIME-ID-01: component identifiability requires iterative relabeling

**Source.** `chapters/ch7_clustering.tex:145-235`.

The model uses the identifiable representation `omega <= 1/2`, but the printed
iteration does not state how to restore that representation after an M-step.
With `gamma_i = P(Z_i = 2 | X_i)`, exchanging components requires
`omega' = 1 - omega`, `mu1' = mu2`, `mu2' = mu1`,
`gamma_i' = 1 - gamma_i`, and `beta' = -beta`. The linear score and logit
threshold change sign together, so away from an exact decision-boundary tie
the decision regions are unchanged up to the component-name exchange. The
documented exact tie is assigned to canonical component 1.

**Implementation decision.** `chime_clustering()` applies this canonical
exchange after the initial state and every update. Means, responsibilities,
beta, score, classes and history are transformed together, and every stage
records whether a label swap occurred. The supplied lambda path and exact KKT
certificate are otherwise unchanged.

### CH7-SPATIAL-CONTRACT-01: spatial fitting and assignment need distinct contracts

**Source.** `chapters/ch7_clustering.tex:390-568`.

The spatial median fits an unsquared distance objective, while cluster
assignment may compare squared Euclidean or SSCM-metric distances. Squaring a
nonnegative distance preserves its assignment argmin, but it does not change
the center objective or justify reporting a different objective certificate.
The `U(0)` convention, threshold equality, empty-active fallback, empty-cluster
repair, cycle handling, returned final state and selector failure rules are
additional algorithm contracts; together they are not a single automatically
monotone objective.

**Implementation decision.** The spatial APIs expose center and assignment
geometry separately, use `U(0) = 0`, specify threshold equality and ties, and
default to no repair. Any requested farthest-empty action is certified and
cycle/final-state diagnostics are explicit. Output-only coordinate resetting
never feeds the algorithmic state.

### CH7-SEMC-CONTRACT-01: SEMC formulas and software-only controls are separated

**Source.** `chapters/ch7_clustering.tex:671-929`;
`FengZhuang2026SEMC`, arXiv:2605.08995; and the official author repository
`flnankai/GEMcluster` audited at
`10fce04fe690fe274dd5d237cfcd3d5c6a4139f6`.

The former source-blocked conclusion is withdrawn. The primary paper closes
the semiparametric mixture likelihood, posterior, center/shape/generator
updates and equations (2.18)--(2.21) Gap-LSE rule. Practical initialization,
KDE/spline endpoints, clipping and floors, damping, factor/threshold selection,
and a soft posterior dispersion default are identifiable software contracts
rather than unique paper formulas. They must therefore remain explicit
arguments and diagnostics rather than being silently represented as theorem
defaults.

**Implementation decision.** `semc_fit()` exposes all such controls and the
three `shape = "tyler"`, `"poet"`, and `"glasso"` paths. Every
returned valid fit has strict symmetry/SPD and method-specific convergence/KKT
certificates; under non-strict mode an unconverged last iterate is diagnostic
only, and `predict.semc_fit()` rejects it.
`semc_select_k_gap()` defaults to the paper hard-label
`mean(log1p(delta))` dispersion and makes the official-software
posterior-weighted `sum(tau * delta)` alternative explicit.

The six native entries are `cpp_ch7_semc_delta`,
`cpp_ch7_semc_softmax`, `cpp_ch7_semc_weighted_sign_scatter`,
`cpp_ch7_semc_weighted_tyler`, `cpp_ch7_semc_offdiag_glasso`, and
`cpp_ch7_semc_weighted_kde`. Tests lock equations, all shape paths, a frozen
official fixture, RNG isolation, invalid-fit prediction, selector contracts,
and boundary failures; both public entry topics have executable fixed-fixture
examples. No paper simulation, grid, size/power table, or reproduction data is
included.

**Provenance.** The three SEMC implementation files are a complete independent
rewrite; no copyrightable author-repository source expression and no
`huge`/GPL code was copied. The pinned repository declares SPDX `MIT` and
its license metadata names copyright year 2026 and holders Dan Zhuang and Long
Feng. The commit, license, and independent-rewrite statement are retained for
the package-level third-party provenance record.

## TINST-SRC-01：two-sample INST 原始 LaTeX 的组索引

这些问题位于辅助论文源
`相关论文的原始latex/2R/uotwo_cjs_r.tex`，不是最新书稿正文，但会直接影响 tINST
复写：

| 行 | 当前文本 | 正确口径 |
|---:|---|---|
| 251 | denominator `\sum_{j\not=i}^{n}` | `\sum_{j\not=i}^{n_k}` |
| 252 | update 左侧 `\D_i` | `\D_k` |
| 348 | 第二组 \(\widehat{\operatorname{tr}(A_2^2)}\) 的外层 `i=1,...,n_1` | `i=1,...,n_2` |

对应位置还可由第 248--253 行的组别 \(k\) 定义以及第 347--350 行三个 trace
估计量的对称性直接核对。

**实现决策。** 所有 leave-one-out location/diagonal-shape 更新都由组局部的
`n[k]` 驱动；第二组 trace 循环严格使用 `n2`。测试必须包含 `n1 != n2`，否则
第 348 行错误会被平衡样本掩盖。

## SLIDES-SSLDA-01：SSLDA 直接估计目标

**来源。** `slides.tex:1039-1042`。第 1039 行声称
\(w=V_0^{-1}\mu_d\)，但紧接着的约束

\[
\|p\widehat S w-(\widehat\mu_0-\widehat\mu_1)\|_\infty\le\lambda_n
\]

以及该节对 \(V_0\) 的定义均对应 \(w=V_0\mu_d\)，不是
\(V_0^{-1}\mu_d\)。最新书稿第 5 章以 normalized shape
\(\Lambda\) 写作 \(\gamma^\star=\Lambda^{-1}\delta\)，与其自身约束一致；这里的
错误限于 slide 使用的 \(V_0\) 记号。

**实现决策。** 不从 slide 第 1039 行生成目标；SSLDA 的目标向量由优化约束所对应
的 population equation 定义。实现和文档统一使用书稿的
\(\gamma^\star=\Lambda^{-1}\delta\) 记号，并用
[`ZhuangFeng2025SSLDA`](https://arxiv.org/abs/2504.11117) 原文复核符号方向。

## TEX-REF-01：三处 `\ref` 转义损坏

以下位置本应是 `\ref{...}`，但反斜杠被生成工具解释成 carriage return；有的在
同一物理行留下内部 CR，有的被转成换行后的 `ef{...}`：

| 源码位置 | 当前损坏文本 | 应为 |
|---|---|---|
| `ch2_location.tex:3770-3771` | `Section~<CR>ef{sec:ch2-erht}` | `Section~\ref{sec:ch2-erht}` |
| `ch3_matrix.tex:2147-2148` | `Section~<CR>ef{sec:ch3-tensor-elliptical-graph}` | `Section~\ref{sec:ch3-tensor-elliptical-graph}` |
| `ch4_other_tests.tex:3432-3433` | `Section~` 后换行 `ef{sec:ch4-radial-directional}` | `Section~\ref{sec:ch4-radial-directional}` |

**实现决策。** 包文档直接链接规范 section label，不复制损坏字符串。任何自动抽取
脚本都应先扫描：内部 CR、行首 `ef{`、以及未转义的 `ref{`。

## TEX-NEQ-01：第 5 章两处不等号转义丢失

**来源。** `chapters/ch5_classification.tex:213` 和 `:234`。当前均写为

```tex
$\mLambda_1 eq\mLambda_2$
```

根据相邻的 equal-shape LDA 情形和随后 generalized QDA 讨论，两处均应为

```tex
$\mLambda_1 \neq \mLambda_2$
```

该错误也出现在当前 PDF 的可见正文中。

**实现决策。** 分类器分支条件按 `equal_shape = FALSE` / \(\Lambda_1\ne\Lambda_2\)
解释；不得把裸文本 `eq` 解析成等号。测试分别覆盖 common-shape SSLDA 和
unequal-shape SSQDA 路径。

## 采用与关闭规则

每条勘误只有在以下条件全部满足后才能标为“已在软件中关闭”：

1. 对应公开函数或内部 kernel 已采用本文档的实现决策；
2. 至少有一个直接公式测试和一个确定性数值/不变性测试；
3. 帮助页或 vignette 明示样本 covariance denominator、normalization 和退化情形；
4. 测试文件注明本勘误 ID；
5. 与原论文或正式 corrigendum 的数值结果在记录的容差内一致。

书稿后续若修正，只更新本文件的状态和新行号，不删除历史条目；这样 JSS 软件论文
可以准确说明实现为何与 2026-08-14 PDF 的个别公式不同。


## BIB-PRIMARY-01: independent publisher metadata corrects four records

**Independent evidence.** `development/PRIMARY_METHOD_VERIFICATION.csv`,
publisher DOI pages, and official institutional repositories checked on
2026-08-15.

- `BosePalSahaRayNayak2015` has canonical authors Smarajit Bose, Amita Pal,
  Rita SahaRay, and Jitadeepa Nayak. The book record names four different
  people.
- `FisherSunGallagher2010` has third author Colin M. Gallagher, not
  Christopher M. Gallagher.
- `ParkAyyala2013` has canonical authors Junyong Park and Deepak Nag Ayyala,
  not the author forms currently stored in the book bibliography.
- `OllilaRaninen2019Shrinkage` has DOI
  `10.1109/TSP.2019.2908144`; `10.1109/TSP.2019.2906691`, formerly printed
  in package documentation, is not the DOI of this work.

**Implementation decision.** The book source is preserved as an authored
release artifact. The package documentation uses the verified Ollila DOI,
while the traceability ledgers retain source-derived metadata and record the
canonical comparison rather than silently rewriting provenance.

## CH6-GSPCA-PROVENANCE-02: published GSPCA core and unavailable book attribution

**Source.** `chapters/ch6_pca_factor.tex:574-660,942-1017`;
Leyder, Raymaekers, and Verdonck (2024), DOI
`10.1007/s11222-024-10413-9`; Raymaekers and Rousseeuw (2019), DOI
`10.1016/j.jmva.2018.11.010`.

The published sources determine the GSSCM radial transformations and the
GSSCM-to-eigenvector PCA algorithm. The book separately attributes
high-dimensional theory to `WangWangFeng2026GSPCA`, but no local manuscript,
DOI, arXiv record, publisher page, or independently auditable public source
was found. It is therefore not implementation evidence and must not be
labelled `local_source`.

**Implementation decision.** MTR-077 cites the published GSSCM/GSPCA sources.
The unavailable citation remains a book-attributed boundary, and package
documentation does not claim to reproduce its unverified theory.
