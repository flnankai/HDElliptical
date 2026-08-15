#!/usr/bin/env Rscript

# Deterministic package-level timing harness.
#
# This script measures representative callable workflows. It is not a
# reproduction of any paper's Monte Carlo design, size/power table, or data
# analysis. Run it against an installed HDElliptical package.

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1L) args[[1L]] else "HDElliptical-benchmarks.csv"
repetitions <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5L
warmups <- if (length(args) >= 3L) as.integer(args[[3L]]) else 1L
minimum_batch_seconds <- if (length(args) >= 4L) {
  as.numeric(args[[4L]])
} else {
  0.05
}
if (is.na(repetitions) || repetitions < 1L ||
    is.na(warmups) || warmups < 0L ||
    length(minimum_batch_seconds) != 1L ||
    !is.finite(minimum_batch_seconds) || minimum_batch_seconds <= 0) {
  stop(paste(
    "Usage: run-benchmarks.R [output.csv] [repetitions>=1] [warmups>=0]",
    "[minimum_batch_seconds>0]"
  ),
       call. = FALSE)
}
if (!requireNamespace("HDElliptical", quietly = TRUE)) {
  stop("Install HDElliptical before running this benchmark.", call. = FALSE)
}

fixture <- function(n, p, phase = 0) {
  i <- seq_len(n)
  j <- seq_len(p)
  sin(outer(i + phase, j + 1, "*") / 37) +
    cos(outer(i + 2 * phase, 2 * j + 1, "+") / 19) +
    (outer(i, j, function(a, b) (a * b + 3 * a + b) %% 23)) / 41
}

numeric_leaves <- function(x) {
  if (is.numeric(x) || is.integer(x) || is.logical(x)) {
    value <- as.numeric(x)
    return(value[is.finite(value)])
  }
  if (is.list(x)) {
    return(unlist(lapply(x, numeric_leaves), use.names = FALSE))
  }
  numeric()
}

fingerprint <- function(x) {
  value <- numeric_leaves(x)
  if (!length(value)) return("no-numeric-output")
  value <- value[seq_len(min(length(value), 256L))]
  weights <- ((seq_along(value) * 104729) %% 1009 + 1) / 1009
  format(sum(signif(value, 12) * weights), digits = 14, scientific = TRUE)
}

x1 <- fixture(120L, 40L)
x2a <- fixture(70L, 80L, 1)
x2b <- fixture(75L, 80L, 2) + 0.015
x3 <- fixture(120L, 80L, 3)
x4 <- sin(seq_len(800L) / 17) + c(rep(0, 400L), rep(0.12, 400L))
y5 <- factor(rep(c("A", "B"), each = 80L), levels = c("A", "B"))
x5 <- fixture(160L, 20L, 5)
x5[y5 == "A", 1:5] <- x5[y5 == "A", 1:5] + 0.35
x6 <- fixture(200L, 100L, 6)
group7 <- rep(seq_len(3L), each = 100L)
x7 <- fixture(300L, 30L, 7)
x7 <- sweep(x7, 1L, c(-1.5, 0, 1.5)[group7], "+")

cases <- list(
  list(
    chapter = 1L, method = "spatial_median", n = 120L, p = 40L,
    run = function() HDElliptical::spatial_median(
      x1, tol = 1e-7, max_iter = 1000L, warn = FALSE
    )
  ),
  list(
    chapter = 2L, method = "chen_qin_two_sample_test",
    n = 145L, p = 80L,
    run = function() HDElliptical::chen_qin_two_sample_test(x2a, x2b)
  ),
  list(
    chapter = 3L, method = "poet_covariance", n = 120L, p = 80L,
    run = function() HDElliptical::poet_covariance(
      x3, factors = 3L, threshold = 0.08
    )
  ),
  list(
    chapter = 4L, method = "classical_cusum_test", n = 800L, p = 1L,
    run = function() HDElliptical::classical_cusum_test(
      x4, variance = "sample", calibration = "asymptotic"
    )
  ),
  list(
    chapter = 5L, method = "classical_lda_classifier",
    n = 160L, p = 20L,
    run = function() HDElliptical::classical_lda_classifier(
      x5, y5, prior = "equal"
    )
  ),
  list(
    chapter = 6L, method = "classical_pca", n = 200L, p = 100L,
    run = function() HDElliptical::classical_pca(x6, components = 5L)
  ),
  list(
    chapter = 7L, method = "lloyd_kmeans", n = 300L, p = 30L,
    run = function() HDElliptical::lloyd_kmeans(
      x7, clusters = 3L, initial = c(1L, 101L, 201L),
      empty_action = "error", solver_max_iter = 100L
    )
  )
)

time_case <- function(specification) {
  if (warmups > 0L) {
    for (iteration in seq_len(warmups)) {
      invisible(specification$run())
    }
  }

  measure_batch <- function(batch_size) {
    result <- NULL
    timing <- system.time({
      for (batch_index in seq_len(batch_size)) {
        result <- specification$run()
      }
    })
    list(
      elapsed = unname(timing[["elapsed"]]),
      result = result
    )
  }

  batch_size <- 1L
  repeat {
    invisible(gc())
    probe <- measure_batch(batch_size)
    if (probe$elapsed >= minimum_batch_seconds || batch_size >= 4096L) {
      break
    }
    batch_size <- min(4096L, batch_size * 2L)
  }

  elapsed <- numeric(repetitions)
  signatures <- character(repetitions)
  for (iteration in seq_len(repetitions)) {
    invisible(gc())
    repeat {
      timing <- measure_batch(batch_size)
      if (timing$elapsed >= minimum_batch_seconds || batch_size >= 4096L) {
        break
      }
      batch_size <- min(4096L, batch_size * 2L)
    }
    elapsed[[iteration]] <- timing$elapsed / batch_size
    signatures[[iteration]] <- fingerprint(timing$result)
  }
  if (length(unique(signatures)) != 1L) {
    stop(sprintf(
      "The deterministic result fingerprint changed for %s.",
      specification$method
    ), call. = FALSE)
  }
  data.frame(
    chapter = specification$chapter,
    method = specification$method,
    n = specification$n,
    p = specification$p,
    repetitions = repetitions,
    warmups = warmups,
    calls_per_measurement = batch_size,
    minimum_batch_seconds = minimum_batch_seconds,
    elapsed_min_seconds = min(elapsed),
    elapsed_median_seconds = stats::median(elapsed),
    elapsed_max_seconds = max(elapsed),
    result_fingerprint = signatures[[1L]],
    package_version = as.character(
      utils::packageVersion("HDElliptical")
    ),
    R_version = as.character(getRversion()),
    platform = R.version$platform,
    measured_at_utc = format(
      Sys.time(), tz = "UTC", usetz = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, lapply(cases, time_case))
output_directory <- dirname(normalizePath(
  output, winslash = "/", mustWork = FALSE
))
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
}
utils::write.csv(result, output, row.names = FALSE)
print(result, row.names = FALSE)
cat("Wrote deterministic benchmark results to", output, "\n")
