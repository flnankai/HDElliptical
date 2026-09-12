## Maintenance update

This is a maintenance update of HDElliptical 0.1.2, now version 0.1.3,
in response to CRAN's 2026-09-12 notification about the linux-arm64 and
noLD Additional issues.

## Changes made

- On linux-arm64, three failures came from architecture-sensitive numerical
  boundary fixtures. A mathematically degenerate normalized Helmert fixture
  can trigger either the strict trace(Sigma^2) or trace(Sigma^3) guard,
  depending on final-bit rounding; the test now accepts either documented
  no-repair error. The PDQ exact-zero fixture now uses dyadic pooled scaling,
  so its intended zero cross residual is exact with or without fused
  multiply-add.
- On the no-long-double check, two chi-square p-value comparisons differed
  only at approximately 1e-15. Their test tolerance is now 1e-12, which
  remains a stringent comparison while allowing extended-precision and
  libm differences.
- No package implementation, public API, statistical formula, or returned
  result was changed.

## Test environments

- Windows 11 x64 (build 26200), R 4.5.2 (2025-10-31 ucrt), ASCII session,
  GCC 14.3.0, `R CMD check --as-cran --no-manual`

## R CMD check results

The local check of the 0.1.3 source tarball completed with:

0 ERROR | 0 WARNING | 2 NOTEs

One NOTE came from transient connection failures during CRAN incoming remote
checks: the local run could not retrieve a Bioconductor index and could not
connect to three GitHub URLs in README.md (libcurl timeout/connection-reset
errors). The same package URLs are present in version 0.1.2, whose regular CRAN
checks are OK. A prior local run completed without this network NOTE.

The other NOTE is the conservative Windows DLL scan for linked `_exit`, `abort`,
and `exit` symbols. An exact direct-call scan of all R, C, and C++ sources found
no calls to those entry points; the symbols enter through linked
runtime/toolchain libraries. All examples, `donttest` examples, tests, and
vignettes completed successfully. The test suite recorded 6,267 passes with
zero failures, warnings, or skips.

## Downstream dependencies

CRAN lists no reverse dependencies for HDElliptical.
