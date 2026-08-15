# HDElliptical deterministic benchmark report

This report records one representative callable workflow per book chapter. It
is machine-local descriptive evidence, not a reproduction of any paper's
simulation, size/power table, or empirical application, and it is not a
performance guarantee.

- Package: HDElliptical 0.1.0
- R: 4.5.2, `x86_64-w64-mingw32`
- CPU: Intel Core i9-14900HX, 24 cores / 32 logical processors
- Memory: 63.7 GiB visible
- OS: Windows 11 64-bit, version 10.0.26200
- Measurement UTC: 2026-08-15 10:26:33--10:26:38
- Repetitions: 5 measured batches after 1 warm-up
- Minimum batch duration: 0.05 seconds

Reported values are elapsed seconds per call. The harness increases the number
of calls in a batch until the timing exceeds the minimum duration, then divides
by the batch size. Every repeated result fingerprint was identical.

| Chapter | Workflow | n | p | Calls/batch | Min (s) | Median (s) | Max (s) |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | `spatial_median()` | 120 | 40 | 64 | 0.0015625 | 0.00171875 | 0.00203125 |
| 2 | `chen_qin_two_sample_test()` | 145 | 80 | 8 | 0.00625 | 0.0075 | 0.0075 |
| 3 | `poet_covariance()` | 120 | 80 | 32 | 0.0028125 | 0.003125 | 0.0034375 |
| 4 | `classical_cusum_test()` | 800 | 1 | 128 | 0.00046875 | 0.00046875 | 0.000546875 |
| 5 | `classical_lda_classifier()` | 160 | 20 | 256 | 0.0002734375 | 0.0003125 | 0.0003125 |
| 6 | `classical_pca()` | 200 | 100 | 32 | 0.0028125 | 0.003125 | 0.003125 |
| 7 | `lloyd_kmeans()` | 300 | 30 | 128 | 0.000859375 | 0.0009375 | 0.001015625 |

The full numeric output, including deterministic fingerprints, is in
`HDElliptical-benchmarks.csv`. Reproduce it with the installed-package script
under `inst/benchmarks/run-benchmarks.R`.
