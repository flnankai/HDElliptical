# JSS manuscript workspace

This directory contains the submission-ready source for the
*Journal of Statistical Software* manuscript accompanying
`HDElliptical` version 0.1.3, published on CRAN on 2026-09-14. The
canonical package page is
<https://CRAN.R-project.org/package=HDElliptical>, and the package DOI is
<https://doi.org/10.32614/CRAN.package.HDElliptical>. The article covers the
implemented methods from book Chapters 1--7 through representative, executable
workflows and a compact traceability appendix.

## Scope and reproducibility

The audited release exposes 176 interfaces and 140 registered native `.Call`
routines, documented by 223 Rd topics. Its 56 test files contain 688 named
test blocks, and all 6,267 runtime assertions passed at the dated 2026-09-13
release gate. Two fresh generation passes produced the same 226-file manifest.
The manuscript examples use fixed small inputs or fixed local seeds. They do
not reproduce paper-specific Monte Carlo size/power grids, benchmark tables,
or empirical results.

Permutation, multiplier, and empirical-null draws appear only when they are
part of a callable method's own calibration. Such interfaces accept explicit
seeds or fixed draws and isolate the caller's random-number state. SEMC is
implemented through `semc_fit()`, `predict.semc_fit()`, and
`semc_select_k_gap()`: all controls are explicit, the paper hard-label
`mean(log1p(delta))` Gap-LSE rule remains the default, and the independent
rewrite records the pinned MIT-licensed GEMcluster commit as provenance.
Exactly two method families remain source-blocked: Chapter 2 structured-
correlation location combinations and Chapter 4 mutual-independence max-sum
studentization.

## Contents

- `article.Rnw`: authoritative knitr/JSS manuscript source.
- `refs.bib`: bibliography; every entry is cited and every cite key resolves.
- `replication.R`: purl output with `documentation = 2`; it executes all
  displayed workflows in manuscript order.
- `jss.cls`, `jss.bst`, and `jsslogo.jpg`: vendored official JSS style files.
- `../output/benchmarks/HDElliptical-benchmarks.csv` and
  `../output/benchmarks/BENCHMARK-REPORT.md`: deterministic descriptive
  timings for one callable workflow per chapter; these are not paper
  simulations or performance guarantees.
- `../output/pdf/HDElliptical-JSS-article.pdf`: promoted release PDF. Source
  builds are staged separately and promoted only after final visual QA.

The vendored style files were compared byte-for-byte with the official JSS
R/LaTeX template on 2026-08-15. The software is distributed under the
GPL-compatible MIT license stated in `DESCRIPTION` and `LICENSE`.

## Build

Install the exact CRAN release (or the package from the parent directory), then
run from this directory in a clean R session:

```r
install.packages("HDElliptical", repos = "https://cloud.r-project.org")
stopifnot(packageVersion("HDElliptical") == package_version("0.1.3"))
knitr::knit("article.Rnw", output = "article.tex")
knitr::purl(
  "article.Rnw", output = "replication.R", documentation = 2
)
source("replication.R", echo = FALSE, chdir = TRUE)
tools::texi2pdf("article.tex", clean = TRUE)
```

The manuscript itself activates `knitr::render_sweave()` so code input and
output use the JSS environments. A release build must verify that a fresh purl
output is line-for-line identical to the committed `replication.R`, that the
standalone script exits successfully, and that LaTeX/BibTeX reports no
undefined citations or references.

The current official style guide and submission checklist are available at:

- <https://www.jstatsoft.org/style>
- <https://www.jstatsoft.org/about/submissions>
- <https://www.jstatsoft.org/public/journals/1/jss-article-rnw.zip>
