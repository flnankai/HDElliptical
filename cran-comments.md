## Resubmission

This is a resubmission of HDElliptical 0.1.1, now version 0.1.2.

Thank you for identifying the incomplete authorship and copyright information.
We audited the authorship, provenance, and license metadata of every package
component and made the following changes:

- Authors@R now lists Long Feng with the aut, cre, and cph roles for the
  HDElliptical implementation.
- Authors@R now lists Dan Zhuang with the ctb role for the SEMC methodology
  and the external MIT-licensed GEMcluster reference implementation used
  during validation. No GEMcluster source code is included in HDElliptical.
- inst/COPYRIGHTS now distinguishes the copyright in the external GEMcluster
  repository from the copyright in HDElliptical, and records the immutable
  GEMcluster revision inspected, its MIT license, and its copyright holders,
  Dan Zhuang and Long Feng.
- The SEMC code in HDElliptical is an independent reimplementation from the
  published mathematical specification and documented behavioral contracts.
  No copyrightable source expression was copied from GEMcluster, and no
  GPL-licensed source code was copied or linked.
- During the audit, we found erroneous GPL-3.0-or-later SPDX headers in two
  original package files, src/chapter4_alpha.cpp and
  src/chapter4_alpha_fdr_conditional.cpp. Long Feng confirmed that he holds
  the copyright in these original files and may distribute them under MIT.
  Their headers now identify Long Feng and MIT, consistently with the package
  license; the files contain no GPL-derived code.
- Authors of the statistical publications cited by the package remain
  identified in the documentation references. They are not listed as software
  contributors or copyright holders because no source code from those
  publications was copied or derived.
- The package version was increased from 0.1.1 to 0.1.2 to distinguish this
  resubmission.

These changes make the ownership and provenance of all distributed package
components explicit and consistent with the package's MIT license.

## Test environments

- Windows 11 x64 (build 26200), R 4.6.1 (2026-06-24 ucrt), UTF-8 session,
  GCC 14.3.0, `R CMD check --as-cran --no-manual`

## R CMD check results

The local check of the 0.1.2 source tarball completed with:

0 ERROR | 0 WARNING | 1 NOTE

The NOTE is the conservative Windows DLL scan for linked `_exit`, `abort`, and
`exit` symbols. An exact direct-call scan of all R, C, and C++ sources found no
calls to those entry points; the symbols enter through linked runtime/toolchain
libraries. All examples, `donttest` examples, tests, and vignettes completed
successfully. The test suite recorded 6,267 passes with zero failures,
warnings, or skips.

## Downstream dependencies

This remains the first CRAN submission of HDElliptical, so there are no
downstream dependencies to check.
