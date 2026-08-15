# Deterministic benchmarks

Run the installed-package benchmark with:

```sh
Rscript run-benchmarks.R HDElliptical-benchmarks.csv 5 1 0.05
```

The harness times one representative callable workflow per book chapter and
records a deterministic numeric fingerprint. It uses no package beyond
HDElliptical itself. The timings are descriptive machine-local evidence, not a
paper-specific simulation, size/power experiment, or performance guarantee.

Arguments are the output CSV, the number of measured repetitions, and the
number of warm-up calls, followed by the minimum elapsed time used to calibrate
each measurement batch. The defaults are 5, 1, and 0.05 seconds. The script
reports per-call timings together with the number of calls in each batch,
avoiding zero-resolution timings for fast methods. It fails if repeated calls
produce different numeric fingerprints.
