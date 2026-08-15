# Implementation and verification policy

`HDElliptical` is an independent implementation of methods described in the
book and its cited papers. This document defines the release gates used while
the package is built chapter by chapter.

## Source fidelity

- The dated book source identifies scope and notation; the original paper is
  authoritative when the book omits a variance estimator, computational step,
  or finite-sample factor.
- Confirmed discrepancies are recorded in `BOOK_ERRATA.md`. The public help
  page states which implementable formula is used.
- Oracle and feasible procedures are separate APIs or explicit modes. A method
  that requires an estimated nuisance quantity must not silently substitute an
  oracle quantity.
- A review-only mention is not marked implemented until its algorithm and
  calibration can be traced to a primary source.

## Code and licensing

- Computationally intensive, self-contained kernels use C++17 through Rcpp and
  RcppArmadillo; R wrappers own user-facing validation and metadata.
- GPL or unlicensed research code may be used as a numerical oracle only. It is
  not copied into this MIT-licensed package. Compatible first-party code is
  still audited before integration.
- Every imported implementation must retain its attribution and compatible
  license notice. Independent rewrites cite the paper rather than source code.

## Practical package scope

- The deliverable is a callable methods package, not a repository of paper
  simulation reproductions. Each paper contributes the statistical method,
  its required nuisance estimation, calibration, diagnostics, documentation,
  and a small usage example.
- Do not add Monte Carlo size/power scripts, paper-specific simulation grids,
  table-reproduction code, or large generated simulation data.
- Small deterministic fixtures and fixed-seed property checks are software
  tests, not simulation studies. They remain required because they verify the
  formula and API without trying to reproduce a paper's empirical section.

## API contract

- Rows are observations and columns are variables throughout.
- Missing or infinite data are rejected unless a method explicitly documents a
  missing-data model.
- Singular covariance matrices, invalid asymptotic variances, and unconverged
  optimization problems are reported; no hidden pseudoinverse, ridge, absolute
  value, or variance floor changes the named method.
- Exact zero directions follow `U(0) = 0`. The default zero tolerance is zero
  so scale equivariance is not lost; a positive absolute tolerance is opt-in.
- High-dimensional tests return raw statistics, calibration quantities, and
  diagnostics in addition to a conventional `htest` interface.

## Numerical requirements

- Quadratic forms use factorizations or triangular solves rather than explicit
  inverses, except where a precision matrix is itself the estimand.
- Direction calculations use scaled Euclidean norms. Pairwise differences and
  shape symmetrization are scaled before arithmetic that can overflow.
- Trace calculations use the smaller primal or dual Gram matrix when this
  avoids allocating a dense `p` by `p` matrix.
- Estimating-equation solvers are declared converged only when their defining
  equation residuals meet tolerance. If a source algorithm specifies only
  fixed-point update stability, expose that status under a distinct name such
  as `iteration.stable`, report the score residual separately, and do not call
  the step-size condition an equation-root certificate.

## Tests required before a method is marked implemented

1. A deterministic, independent R reference for the defining formula.
2. The relevant translation, scale, orthogonal, affine, permutation, or
   group-exchange invariance checks.
3. Degenerate-input and smallest-valid-sample tests.
4. Extreme-scale and near-singular regression tests where applicable.
5. An independent implementation comparison when a legally usable numerical
   oracle already exists and the comparison can remain a compact unit test.

Paper tables and Monte Carlo size/power sections are deliberately outside the
package scope. They are not release gates and are not reproduced in auxiliary
scripts.
