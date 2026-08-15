#!/usr/bin/env Rscript

# Validate the machine-readable method ledger without loading HDElliptical.
# The check deliberately accepts a roxygen @export declaration when generated
# NAMESPACE/Rd files have not yet been refreshed in the working tree.

root_candidates <- list(
  list(package = "HDElliptical", workspace = "."),
  list(package = ".", workspace = ".."),
  list(package = "..", workspace = file.path("..", ".."))
)
selected_root <- Filter(function(candidate) {
  file.exists(file.path(candidate$package, "DESCRIPTION")) &&
    file.exists(file.path(
      candidate$package, "development", "METHOD_TRACEABILITY.csv"
    ))
}, root_candidates)
if (!length(selected_root)) {
  stop(
    paste(
      "Cannot locate HDElliptical from the workspace, package, or",
      "development directory."
    ),
    call. = FALSE
  )
}
package_root <- selected_root[[1L]]$package
workspace_root <- selected_root[[1L]]$workspace
ledger_file <- file.path(package_root, "development", "METHOD_TRACEABILITY.csv")
book_root <- file.path(
  workspace_root, "book-source", "book_latest_arxiv_20260814_source"
)
bib_file <- file.path(book_root, "references.bib")

errors <- character()
fail <- function(...) {
  errors <<- c(errors, paste0(...))
}
read_utf8 <- function(path) {
  readLines(path, warn = FALSE, encoding = "UTF-8")
}
split_values <- function(value) {
  if (length(value) == 0L || is.na(value) || !nzchar(trimws(value))) {
    return(character())
  }
  answer <- trimws(strsplit(value, ";", fixed = TRUE)[[1L]])
  answer[nzchar(answer)]
}
escape_regex <- function(value) {
  gsub("([][{}()+*^$|\\?.])", "\\\\\\1", value)
}

required_columns <- c(
  "chapter", "book_range", "method", "status", "primary_keys",
  "public_apis", "native_kernels", "test_files", "rd_topics", "notes"
)
if (!file.exists(ledger_file)) stop("Missing METHOD_TRACEABILITY.csv.", call. = FALSE)
ledger <- tryCatch(
  utils::read.csv(
    ledger_file, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character(), fileEncoding = "UTF-8"
  ),
  error = function(error) stop("Cannot parse METHOD_TRACEABILITY.csv: ",
                               conditionMessage(error), call. = FALSE)
)
if (!identical(names(ledger), required_columns)) {
  fail("CSV columns must be exactly: ", paste(required_columns, collapse = ", "))
}
if (nrow(ledger) == 0L) fail("CSV must contain at least one method row.")
if (anyDuplicated(paste(ledger$chapter, ledger$method, sep = "::"))) {
  duplicates <- unique(paste(ledger$chapter, ledger$method, sep = "::")[
    duplicated(paste(ledger$chapter, ledger$method, sep = "::"))
  ])
  fail("Duplicate chapter/method rows: ", paste(duplicates, collapse = ", "))
}

allowed_status <- c("implemented", "review-only", "source-blocked")
bad_status <- setdiff(unique(ledger$status), allowed_status)
if (length(bad_status)) {
  fail("Unknown status values: ", paste(bad_status, collapse = ", "))
}

r_files <- list.files(file.path(package_root, "R"), pattern = "[.]R$",
                      full.names = TRUE)
r_lines <- lapply(r_files, read_utf8)
names(r_lines) <- r_files
namespace_lines <- if (file.exists(file.path(package_root, "NAMESPACE"))) {
  read_utf8(file.path(package_root, "NAMESPACE"))
} else character()
namespace_text <- paste(namespace_lines, collapse = "\n")

source_exported <- function(api) {
  definition <- paste0("^", escape_regex(api), "[[:space:]]*<-[[:space:]]*function\\b")
  for (lines in r_lines) {
    at <- grep(definition, lines, perl = TRUE)
    for (index in at) {
      cursor <- index - 1L
      block <- character()
      while (cursor >= 1L &&
             (grepl("^#'", lines[[cursor]]) || !nzchar(trimws(lines[[cursor]])))) {
        block <- c(lines[[cursor]], block)
        cursor <- cursor - 1L
      }
      if (any(grepl("^#'[[:space:]]+@export(?:[[:space:]]|$)", block,
                    perl = TRUE))) return(TRUE)
    }
  }
  FALSE
}
namespace_exported <- function(api) {
  direct <- paste0("export\\(", escape_regex(api), "\\)")
  if (grepl(direct, namespace_text, perl = TRUE)) return(TRUE)
  parts <- strsplit(api, ".", fixed = TRUE)[[1L]]
  if (length(parts) >= 2L) {
    generic <- parts[[1L]]
    class <- paste(parts[-1L], collapse = ".")
    s3 <- paste0("S3method\\(", escape_regex(generic), ",",
                 escape_regex(class), "\\)")
    if (grepl(s3, namespace_text, perl = TRUE)) return(TRUE)
  }
  FALSE
}

cpp_files <- list.files(file.path(package_root, "src"), pattern = "[.](cpp|cc|cxx)$",
                        full.names = TRUE)
cpp_lines <- lapply(cpp_files, read_utf8)
names(cpp_lines) <- cpp_files
registered_text <- paste(
  unlist(cpp_lines, use.names = FALSE),
  if (file.exists(file.path(package_root, "R", "RcppExports.R")))
    read_utf8(file.path(package_root, "R", "RcppExports.R")) else character(),
  collapse = "\n"
)
kernel_exported <- function(kernel) {
  registered <- grepl(
    paste0("_HDElliptical_", escape_regex(kernel), "\\b"),
    registered_text, perl = TRUE
  )
  if (registered) return(TRUE)
  definition <- paste0("\\b", escape_regex(kernel), "[[:space:]]*\\(")
  for (lines in cpp_lines) {
    at <- grep(definition, lines, perl = TRUE)
    for (index in at) {
      start <- max(1L, index - 4L)
      if (any(grepl("\\[\\[Rcpp::export\\]\\]", lines[start:index],
                    perl = TRUE))) return(TRUE)
    }
  }
  FALSE
}

bib_lines <- if (file.exists(bib_file)) read_utf8(bib_file) else character()
bib_matches <- regexec("^@[[:alpha:]]+\\{([^,]+),", bib_lines, perl = TRUE)
bib_keys <- vapply(regmatches(bib_lines, bib_matches), function(value) {
  if (length(value) >= 2L) value[[2L]] else NA_character_
}, character(1L))
bib_keys <- unique(stats::na.omit(bib_keys))

validate_book_range <- function(value, row_label) {
  entries <- split_values(value)
  if (!length(entries)) {
    fail(row_label, ": missing book_range")
    return(invisible(NULL))
  }
  for (entry in entries) {
    if (entry %in% c("PRIMARY_ONLY", "SOURCE_BLOCKED")) next
    match <- regexec("^([^:]+[.]tex):([0-9]+)-([0-9]+)$", entry, perl = TRUE)
    pieces <- regmatches(entry, match)[[1L]]
    if (length(pieces) != 4L) {
      fail(row_label, ": invalid book range token `", entry, "`")
      next
    }
    chapter_file <- file.path(book_root, "chapters", pieces[[2L]])
    if (!file.exists(chapter_file)) {
      fail(row_label, ": missing manuscript file for `", entry, "`")
      next
    }
    first <- as.integer(pieces[[3L]])
    last <- as.integer(pieces[[4L]])
    number_lines <- length(read_utf8(chapter_file))
    if (is.na(first) || is.na(last) || first < 1L || last < first ||
        last > number_lines) {
      fail(row_label, ": out-of-range manuscript token `", entry,
           "` (file has ", number_lines, " lines)")
    }
  }
}

for (row in seq_len(nrow(ledger))) {
  label <- paste0("row ", row + 1L, " [Ch", ledger$chapter[[row]], ": ",
                  ledger$method[[row]], "]")
  validate_book_range(ledger$book_range[[row]], label)
  citations <- split_values(ledger$primary_keys[[row]])
  if (!length(citations)) {
    fail(label, ": missing primary_keys")
  } else {
    for (citation in citations) {
      allowed_marker <- citation %in% c("BOOK_ONLY", "SOURCE_BLOCKED")
      allowed_url <- grepl("^https?://", citation)
      if (!allowed_marker && !allowed_url && !citation %in% bib_keys) {
        fail(label, ": citation key not in book references: `", citation, "`")
      }
    }
  }
  if (identical(ledger$status[[row]], "implemented")) {
    apis <- split_values(ledger$public_apis[[row]])
    kernels <- split_values(ledger$native_kernels[[row]])
    tests <- split_values(ledger$test_files[[row]])
    topics <- split_values(ledger$rd_topics[[row]])
    if (!length(apis)) fail(label, ": implemented row has no public API")
    if (!length(kernels)) fail(label, ": implemented row has no kernel/PURE_R marker")
    if (!length(tests)) fail(label, ": implemented row has no test file")
    if (!length(topics)) fail(label, ": implemented row has no Rd topic")
    for (api in apis) {
      if (!namespace_exported(api) && !source_exported(api)) {
        fail(label, ": API is neither in NAMESPACE nor roxygen-exported: `",
             api, "`")
      }
    }
    for (kernel in setdiff(kernels, "PURE_R")) {
      if (grepl("[*]", kernel)) {
        fail(label, ": wildcard native token is forbidden: `", kernel, "`")
      } else if (!kernel_exported(kernel)) {
        fail(label, ": native kernel is not registered or Rcpp-exported: `",
             kernel, "`")
      }
    }
    for (test_file in tests) {
      path <- file.path(package_root, test_file)
      if (!file.exists(path)) fail(label, ": missing test file `", test_file, "`")
    }
    for (topic in topics) {
      rd_file <- file.path(package_root, "man", paste0(topic, ".Rd"))
      if (!file.exists(rd_file) && !source_exported(topic)) {
        fail(label, ": missing Rd topic and exported roxygen source `", topic, "`")
      }
    }
  }
}

if (length(errors)) {
  cat("METHOD_TRACEABILITY validation failed:\n", file = stderr())
  cat(paste0("- ", unique(errors), "\n"), file = stderr())
  quit(save = "no", status = 1L)
}

cat(sprintf(
  "METHOD_TRACEABILITY OK: %d rows, %d implemented, %d review-only, %d source-blocked.\n",
  nrow(ledger), sum(ledger$status == "implemented"),
  sum(ledger$status == "review-only"),
  sum(ledger$status == "source-blocked")
))


