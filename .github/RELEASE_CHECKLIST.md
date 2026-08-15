# Release checklist

This repository is intended to use `HDElliptical` as the GitHub repository
root. Complete this checklist from a clean clone of the release tag.

## Source and provenance

- [ ] Confirm `DESCRIPTION` version, `NEWS.md`, `LICENSE`, `LICENSE.md`, and
  `inst/COPYRIGHTS` agree.
- [ ] Run `development/validate-method-traceability.R` and retain its output.
- [ ] Run `development/validate-bibliography-traceability.R` from both the
  workspace and package roots; require every method-primary row to have
  independently auditable evidence and retain the validator output.
- [ ] Confirm every public API has an Rd alias, a traceability row, a test, and
  a cited source or an explicit book-only/source-blocked label.
- [ ] Confirm `.orig`, `.rej`, conflict markers, compiled objects, and temporary
  PDF files are absent.

## Package gate

- [ ] Run `Rcpp::compileAttributes()` and roxygen twice; require an unchanged
  second-pass generated-file manifest.
- [ ] Run all installed tests and record test files, blocks, passing
  expectations, failures, errors, warnings, and skips.
- [ ] Run all examples, both vignettes, and strict Rd parsing.
- [ ] Build the source tarball and run `R CMD check --as-cran` on the tarball.
- [ ] Run `.github/workflows/R-CMD-check.yaml` on Windows, macOS, and Linux.

## Reproducibility and artifacts

- [ ] Re-run `inst/benchmarks/run-benchmarks.R`; treat results as descriptive
  machine-local evidence, not a paper simulation or performance guarantee.
- [ ] Rebuild the JSS manuscript, regenerate `jss/replication.R` by purl, run
  it against the release package, and visually inspect every PDF page.
- [ ] Build the book PDF/source bundle if either changed.
- [ ] Generate a SHA-256 manifest covering the source tarball, JSS PDF,
  benchmark CSV/report, check log, and book assets.

## GitHub publication

- [ ] Create or confirm the `flnankai/HDElliptical` repository, then make a
  bootstrap push **without** creating the release tag. This first push exists
  only to make the `URL` and `BugReports` fields resolvable.
- [ ] After the repository and issue tracker resolve, update the pre-publication
  404 wording in `README.md`, `NEWS.md`, and `cran-comments.md`; also replace
  the installation comment that says "once the public repository is available."
- [ ] Re-run `R CMD check --as-cran` with incoming remote checks enabled,
  rebuild the source tarball, and regenerate the check summary and release
  manifest from that exact post-bootstrap commit.
- [ ] Create an annotated `v0.1.0` tag only after the remote-enabled check is
  green, and verify that the tag matches the regenerated source manifest.
- [ ] Create a GitHub release from the tag and attach the manifest-listed
  deliverable assets plus `RELEASE-MANIFEST.txt` and its `.sha256` sidecar.
  Do not upload obsolete draft PDFs or development tarballs.
