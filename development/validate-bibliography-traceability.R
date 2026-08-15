#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

current_ctype <- Sys.getlocale("LC_CTYPE")
if (!grepl("UTF-?8", current_ctype, ignore.case = TRUE)) {
  locale_candidates <- if (.Platform$OS.type == "windows") {
    c(".UTF-8", "English_United States.utf8",
      "Chinese (Simplified)_China.utf8")
  } else {
    c("C.UTF-8", "en_US.UTF-8")
  }
  for (candidate in locale_candidates) {
    selected <- suppressWarnings(Sys.setlocale("LC_CTYPE", candidate))
    if (nzchar(selected) && grepl("UTF-?8", selected, ignore.case = TRUE)) break
  }
}
if (!grepl("UTF-?8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
  stop("A working UTF-8 LC_CTYPE locale is required.", call. = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

abort <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

check <- function(ok, fmt, ...) {
  if (length(ok) != 1L || is.na(ok) || !ok) abort(fmt, ...)
}

root_candidates <- list(
  list(package = "HDElliptical", workspace = "."),
  list(package = ".", workspace = ".."),
  list(package = "..", workspace = file.path("..", ".."))
)
selected_root <- Filter(function(candidate) {
  file.exists(file.path(candidate$package, "DESCRIPTION")) &&
    file.exists(file.path(
      candidate$package, "development", "METHOD_TRACEABILITY.csv"
    )) &&
    file.exists(file.path(
      candidate$workspace, "book-source",
      "book_latest_arxiv_20260814_source", "references.bib"
    ))
}, root_candidates)
check(length(selected_root) >= 1L,
      "Cannot locate the package and book source from the workspace, package, or development directory.")
package_dir <- selected_root[[1L]]$package
workspace_dir <- selected_root[[1L]]$workspace
development_dir <- file.path(package_dir, "development")
book_dir <- file.path(
  workspace_dir, "book-source", "book_latest_arxiv_20260814_source"
)

paths <- list(
  bibliography = file.path(development_dir, "BIBLIOGRAPHY_TRACEABILITY.csv"),
  primary = file.path(
    development_dir, "PRIMARY_METHOD_VERIFICATION.csv"),
  methods = file.path(development_dir, "METHOD_TRACEABILITY.csv"),
  audit = file.path(development_dir, "LITERATURE_AUDIT.md"),
  bib = file.path(book_dir, "references.bib"),
  bbl = file.path(book_dir, "main.bbl")
)

read_utf8 <- function(path, reject_bom = FALSE) {
  check(file.exists(path), "Required file does not exist: %s", path)
  size <- file.info(path)$size
  raw <- readBin(path, what = "raw", n = size)
  if (reject_bom && length(raw) >= 3L) {
    check(!identical(as.integer(raw[1:3]), c(239L, 187L, 191L)),
          "UTF-8 BOM is not allowed: %s", path)
  }
  text <- rawToChar(raw)
  check(validUTF8(text), "Invalid UTF-8 byte sequence: %s", path)
  check(!is.na(iconv(text, from = "UTF-8", to = "UTF-8", sub = NA)),
        "UTF-8 round-trip failed: %s", path)
  text
}

collapse_space <- function(x) {
  gsub("[[:space:]]+", " ", trimws(x), perl = TRUE)
}

extract_citations <- function(path) {
  lines <- strsplit(read_utf8(path), "\n", fixed = TRUE)[[1L]]
  lines <- sub("(?<!\\\\)%.*$", "", lines, perl = TRUE)
  text <- paste(lines, collapse = "\n")
  pattern <- paste0(
    "\\\\cite[A-Za-z*]*\\s*",
    "(?:\\[[^]]*\\]\\s*){0,2}\\{[^}]+\\}"
  )
  hits <- gregexpr(pattern, text, perl = TRUE)[[1L]]
  calls <- if (identical(hits[1L], -1L)) character() else regmatches(text, list(hits))[[1L]]
  payload <- sub(".*\\{([^}]*)\\}\\s*$", "\\1", calls, perl = TRUE)
  keys <- trimws(unlist(strsplit(payload, ",", fixed = TRUE), use.names = FALSE))
  keys <- keys[nzchar(keys)]
  check(all(grepl("^[A-Za-z0-9_.:+/-]+$", keys)),
        "Malformed citation key(s) in %s: %s", path,
        paste(unique(keys[!grepl("^[A-Za-z0-9_.:+/-]+$", keys)]), collapse = ";"))
  list(keys = keys, unique = sort(unique(keys)), counts = table(keys))
}

section_files <- c(
  `1` = "ch1_foundations.tex",
  `2` = "ch2_location.tex",
  `3` = "ch3_matrix.tex",
  `4` = "ch4_other_tests.tex",
  `5` = "ch5_classification.tex",
  `6` = "ch6_pca_factor.tex",
  `7` = "ch7_clustering.tex",
  `Appendix-P` = "appendix_probability.tex"
)
section_paths <- file.path(book_dir, "chapters", unname(section_files))
names(section_paths) <- names(section_files)
citations <- lapply(section_paths, extract_citations)

expected_chapter_counts <- c(`1` = 12L, `2` = 39L, `3` = 46L, `4` = 44L,
                             `5` = 25L, `6` = 28L, `7` = 29L)
observed_chapter_counts <- vapply(citations[names(expected_chapter_counts)],
                                  function(x) length(x$unique), integer(1L))
check(identical(observed_chapter_counts, expected_chapter_counts),
      "Seven-chapter citation counts changed: observed %s; expected %s.",
      paste(observed_chapter_counts, collapse = ","),
      paste(expected_chapter_counts, collapse = ","))

chapter_union <- sort(unique(unlist(
  lapply(citations[names(expected_chapter_counts)], `[[`, "unique"),
  use.names = FALSE
)))
appendix_keys <- citations[["Appendix-P"]]$unique
manuscript_union <- sort(unique(c(chapter_union, appendix_keys)))
check(length(chapter_union) == 191L,
      "Seven-chapter union must contain 191 keys, found %d.", length(chapter_union))
check(length(appendix_keys) == 4L,
      "Probability appendix must contain 4 keys, found %d.", length(appendix_keys))
check(identical(intersect(chapter_union, appendix_keys), "Tropp2012"),
      "Expected Tropp2012 to be the sole chapter/appendix overlap.")
check(length(manuscript_union) == 194L,
      "Compiled manuscript citation union must contain 194 keys, found %d.",
      length(manuscript_union))

parse_bib_headers <- function(path) {
  lines <- strsplit(read_utf8(path), "\n", fixed = TRUE)[[1L]]
  headers <- grep("^\\s*@[A-Za-z]+\\s*\\{\\s*[^,[:space:]]+\\s*,",
                  lines, value = TRUE, perl = TRUE)
  types <- tolower(sub("^\\s*@([A-Za-z]+).*", "\\1", headers, perl = TRUE))
  keys <- sub("^\\s*@[A-Za-z]+\\s*\\{\\s*([^,[:space:]]+).*", "\\1",
              headers, perl = TRUE)
  check(length(keys) == length(unique(keys)), "Duplicate key(s) in references.bib.")
  setNames(types, keys)
}

parse_bbl_keys <- function(path) {
  text <- read_utf8(path)
  pattern <- paste0(
    "(?s)\\\\bibitem\\s*\\[.*?\\]\\s*\\{\\s*%?\\s*",
    "([^}\\s%]+)\\s*\\}"
  )
  hits <- gregexpr(pattern, text, perl = TRUE)[[1L]]
  check(!identical(hits[1L], -1L), "No bibitem found in main.bbl.")
  items <- regmatches(text, list(hits))[[1L]]
  sub(pattern, "\\1", items, perl = TRUE)
}

bib_types <- parse_bib_headers(paths$bib)
bib_keys <- names(bib_types)
bbl_keys <- parse_bbl_keys(paths$bbl)
check(length(bib_keys) == 200L, "references.bib must contain 200 entries, found %d.",
      length(bib_keys))
check(length(bbl_keys) == 194L && length(unique(bbl_keys)) == 194L,
      "main.bbl must contain 194 unique bibitems, found %d (%d unique).",
      length(bbl_keys), length(unique(bbl_keys)))
check(setequal(bbl_keys, manuscript_union),
      "main.bbl keys differ from the compiled manuscript citation union.")
check(all(manuscript_union %in% bib_keys),
      "Cited key(s) missing from references.bib: %s",
      paste(setdiff(manuscript_union, bib_keys), collapse = ";"))

expected_unused <- sort(c(
  "FengSun2015Note", "FopMurphyScrucca2019",
  "WangFengLiuZhou2021Directional", "WittenTibshirani2011",
  "ZhouPanShen2009", "ZouYinFengWang2014NPMLECP"
))
unused_bib <- sort(setdiff(bib_keys, manuscript_union))
check(identical(unused_bib, expected_unused),
      "The six uncited BibTeX entries changed: %s", paste(unused_bib, collapse = ";"))

check(requireNamespace("jsonlite", quietly = TRUE),
      "Package 'jsonlite' is required for the Pandoc CSL-JSON gate.")
pandoc <- Sys.which("pandoc")
check(nzchar(pandoc), "Pandoc is required for bibliography and CSV gates.")

run_pandoc <- function(args, label) {
  output <- suppressWarnings(system2(pandoc, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status") %||% 0L
  check(status == 0L, "%s failed (exit %d): %s", label, status,
        paste(output, collapse = "\n"))
  invisible(output)
}

csl_path <- tempfile("HDElliptical-bibliography-", fileext = ".json")
on.exit(unlink(csl_path, force = TRUE), add = TRUE)
run_pandoc(c(
  "-f", "biblatex", "-t", "csljson", shQuote(paths$bib),
  "-o", shQuote(csl_path)
), "Pandoc BibTeX-to-CSL conversion")
references <- jsonlite::fromJSON(csl_path, simplifyVector = FALSE)
check(length(references) == 200L,
      "Pandoc CSL-JSON must contain 200 records, found %d.", length(references))
reference_ids <- vapply(references, function(x) x$id %||% "", character(1L))
check(length(unique(reference_ids)) == 200L && setequal(reference_ids, bib_keys),
      "Pandoc CSL identifiers do not exactly match references.bib keys.")
references <- setNames(references, reference_ids)

format_person <- function(person) {
  literal <- person$literal %||% ""
  if (nzchar(literal)) return(collapse_space(literal))
  family <- collapse_space(paste(
    c(person[["non-dropping-particle"]] %||% "", person$family %||% ""),
    collapse = " "
  ))
  given <- collapse_space(person$given %||% "")
  if (nzchar(family) && nzchar(given)) paste0(family, ", ", given) else paste0(family, given)
}

format_authors <- function(reference) {
  authors <- reference$author %||% list()
  collapse_space(paste(vapply(authors, format_person, character(1L)), collapse = "; "))
}

reference_year <- function(reference) {
  parts <- reference$issued[["date-parts"]] %||% list()
  if (!length(parts) || !length(parts[[1L]])) return("")
  as.character(parts[[1L]][[1L]])
}

normalise_identifier <- function(url) {
  if (grepl("^https?://doi\\.org/", url, ignore.case = TRUE)) {
    return(paste0("doi:", sub("^https?://doi\\.org/", "", url,
                              ignore.case = TRUE)))
  }
  if (grepl("^https?://arxiv\\.org/abs/", url, ignore.case = TRUE)) {
    return(paste0("arXiv:", sub("^https?://arxiv\\.org/abs/", "", url,
                                ignore.case = TRUE)))
  }
  paste0("url:", url)
}

identifier_rank <- function(identifier) {
  if (grepl("^doi:", identifier)) 3L else if (grepl("^arXiv:", identifier)) 2L else 1L
}

extract_audit_identifiers <- function(path, valid_keys) {
  lines <- strsplit(read_utf8(path), "\n", fixed = TRUE)[[1L]]
  lines <- lines[grepl("^\\| `", lines) & grepl("https?://", lines)]
  ids <- character()
  sources <- character()
  assign_one <- function(key, identifier) {
    old <- if (key %in% names(ids)) ids[[key]] else ""
    if (!nzchar(old) || identifier_rank(identifier) > identifier_rank(old)) {
      ids[[key]] <<- identifier
      sources[[key]] <<- "LITERATURE_AUDIT"
    }
  }
  for (line in lines) {
    key_tokens <- regmatches(line, gregexpr("`[^`]+`", line, perl = TRUE))[[1L]]
    keys <- gsub("^`|`$", "", key_tokens)
    keys <- keys[keys %in% valid_keys]
    urls <- regmatches(line, gregexpr("https?://[^)]+", line, perl = TRUE))[[1L]]
    if (!length(keys) || !length(urls)) next
    if (length(keys) > 1L && length(keys) == length(urls)) {
      for (i in seq_along(keys)) assign_one(keys[[i]], normalise_identifier(urls[[i]]))
    } else {
      candidates <- vapply(urls, normalise_identifier, character(1L))
      best <- candidates[[which.max(vapply(candidates, identifier_rank, integer(1L)))]]
      for (key in keys) assign_one(key, best)
    }
  }
  list(ids = ids, sources = sources)
}

explicit_arxiv <- function(reference) {
  text <- paste(unlist(reference, recursive = TRUE, use.names = FALSE), collapse = " ")
  hit <- regexec(
    "(?i)arXiv[^0-9]{0,30}([0-9]{4}\\.[0-9]{4,5}(?:v[0-9]+)?)",
    text, perl = TRUE
  )
  parts <- regmatches(text, hit)[[1L]]
  if (length(parts) >= 2L) paste0("arXiv:", parts[[2L]]) else ""
}

method_data <- read.csv(paths$methods, stringsAsFactors = FALSE,
                        check.names = FALSE, fileEncoding = "UTF-8",
                        na.strings = character())
check(nrow(method_data) == 87L, "METHOD_TRACEABILITY.csv must contain 87 rows.")
method_map <- setNames(vector("list", length(manuscript_union)), manuscript_union)
method_keys <- character()
for (i in seq_len(nrow(method_data))) {
  keys <- trimws(strsplit(method_data$primary_keys[[i]], ";", fixed = TRUE)[[1L]])
  keys <- keys[nzchar(keys) & keys != "BOOK_ONLY" & !grepl("^https?://", keys)]
  method_keys <- c(method_keys, keys)
  for (key in keys) {
    check(key %in% manuscript_union,
          "METHOD_TRACEABILITY primary key is not in the 194-key closure: %s", key)
    method_map[[key]] <- c(method_map[[key]], sprintf("MTR-%03d", i))
  }
}
method_keys <- sort(unique(method_keys))

audit_ids <- extract_audit_identifiers(paths$audit, manuscript_union)
stable_ids <- audit_ids$ids
id_sources <- audit_ids$sources
for (key in manuscript_union) {
  candidate <- explicit_arxiv(references[[key]])
  old <- if (key %in% names(stable_ids)) stable_ids[[key]] else ""
  if (nzchar(candidate) && (!nzchar(old) || identifier_rank(candidate) > identifier_rank(old))) {
    stable_ids[[key]] <- candidate
    id_sources[[key]] <- "references.bib"
  }
}
verified_overrides <- c(
  FengWang2026PDQ = "arXiv:2605.03265",
  MaFengWang2026TimeVaryingAlpha = "arXiv:2604.13772",
  ZhaoWang2024ConditionalMaxAlpha = "arXiv:2604.12252"
)
for (key in names(verified_overrides)) {
  stable_ids[[key]] <- verified_overrides[[key]]
  id_sources[[key]] <- "external_verified_2026-08-15"
}

primary_columns <- c(
  "citation_key", "verification_status", "stable_identifier", "official_url",
  "canonical_authors", "canonical_year", "canonical_title",
  "metadata_comparison", "method_trace_rows", "formula_trace", "verified_on"
)
invisible(read_utf8(paths$primary, reject_bom = TRUE))
primary <- read.csv(
  paths$primary, stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8", colClasses = "character", na.strings = character()
)
check(nrow(primary) == 53L,
      "PRIMARY_METHOD_VERIFICATION.csv must contain 53 evidence rows.")
check(identical(names(primary), primary_columns),
      "Primary-method verification columns differ from the required schema.")
check(!anyDuplicated(primary$citation_key),
      "Duplicate citation_key in primary-method verification CSV.")
check(!any(grepl("[\r\n]", as.matrix(primary), perl = TRUE)),
      "Primary-method verification fields must not contain embedded line breaks.")
check(all(primary$verification_status %in% c(
  "independently_verified", "book_attributed_source_unavailable"
)), "Unknown primary-method verification status.")

external_primary_key <- "LeyderRaymaekersVerdonck2024GSPCA"
boundary_key <- "WangWangFeng2026GSPCA"
pre_primary_stable <- names(stable_ids)[nzchar(stable_ids)]
verification_needed <- sort(setdiff(method_keys, pre_primary_stable))
primary_book_keys <- setdiff(
  primary$citation_key, c(external_primary_key, boundary_key)
)
check(length(primary_book_keys) == 51L &&
        all(primary_book_keys %in% method_keys),
      "The 51 book-key primary evidence rows must all map to current method-primary keys.")
check(all(verification_needed %in% primary$citation_key),
      paste0(
        "Every method-primary key not already closed by an independent stable ",
        "identifier must appear in PRIMARY_METHOD_VERIFICATION.csv: ",
        paste(setdiff(verification_needed, primary$citation_key), collapse = ";")
      ))

primary_identifier_pattern <- paste0(
  "^(doi:10\\.[0-9]{4,9}/[^[:space:]]+|",
  "arXiv:[0-9]{4}\\.[0-9]{4,5}(v[0-9]+)?|url:https?://[^[:space:]]+)$"
)
verified_primary <- primary[
  primary$verification_status == "independently_verified", , drop = FALSE
]
check(nrow(verified_primary) == 52L,
      "Exactly 52 independently verified primary evidence rows are required.")
check(all(grepl(primary_identifier_pattern,
                verified_primary$stable_identifier, perl = TRUE)),
      "Every independently verified primary row requires a valid stable identifier.")
check(all(grepl("^https://", verified_primary$official_url)),
      "Every independently verified primary row requires an HTTPS official URL.")
check(all(nzchar(verified_primary$canonical_authors) &
          grepl("^[12][0-9]{3}$", verified_primary$canonical_year) &
          nzchar(verified_primary$canonical_title) &
          nzchar(verified_primary$metadata_comparison) &
          nzchar(verified_primary$method_trace_rows) &
          nzchar(verified_primary$formula_trace) &
          verified_primary$verified_on == "2026-08-15"),
      "Independently verified rows require canonical metadata, method/formula trace, and date.")

verified_book <- verified_primary[
  verified_primary$citation_key %in% manuscript_union, , drop = FALSE
]
check(nrow(verified_book) == 51L,
      "Exactly 51 independently verified rows must belong to the 194-key book closure.")
for (i in seq_len(nrow(verified_book))) {
  key <- verified_book$citation_key[[i]]
  expected_rows <- paste(method_map[[key]] %||% character(), collapse = ";")
  check(identical(verified_book$method_trace_rows[[i]], expected_rows),
        "Primary evidence method-row mapping mismatch for %s.", key)
}
leyder <- primary[
  primary$citation_key == external_primary_key, , drop = FALSE
]
check(nrow(leyder) == 1L &&
        leyder$stable_identifier == "doi:10.1007/s11222-024-10413-9" &&
        leyder$method_trace_rows == "MTR-077",
      "The direct published GSPCA primary must be the Leyder et al. 2024 MTR-077 row.")
boundary <- primary[primary$citation_key == boundary_key, , drop = FALSE]
check(nrow(boundary) == 1L &&
        boundary$verification_status == "book_attributed_source_unavailable" &&
        !nzchar(boundary$stable_identifier) && !nzchar(boundary$official_url) &&
        !nzchar(boundary$method_trace_rows) && nzchar(boundary$formula_trace),
      "The Wang-Wang-Feng GSPCA attribution must remain an explicit unavailable-source boundary.")
mtr77_keys <- trimws(strsplit(
  method_data$primary_keys[[77L]], ";", fixed = TRUE
)[[1L]])
check("RaymaekersRousseeuw2019" %in% mtr77_keys &&
        "https://doi.org/10.1007/s11222-024-10413-9" %in% mtr77_keys &&
        !"WangWangFeng2026GSPCA" %in% mtr77_keys,
      "MTR-077 must cite the published GSSCM/GSPCA sources and exclude the unavailable attribution.")

for (i in seq_len(nrow(verified_book))) {
  key <- verified_book$citation_key[[i]]
  stable_ids[[key]] <- verified_book$stable_identifier[[i]]
  id_sources[[key]] <- "PRIMARY_METHOD_VERIFICATION"
}

notes <- setNames(rep("", length(manuscript_union)), manuscript_union)
notes[["WangWangFeng2026GSPCA"]] <- paste(
  "Book-attributed source unavailable: no local manuscript, DOI, arXiv record,",
  "publisher page, or independently auditable public source was found as of 2026-08-15."
)
notes[["ZhaoWang2024ConditionalMaxAlpha"]] <- paste(
  "Key and references.bib record year are 2024; the public arXiv paper is 2026."
)
notes[["ZhaoWangFeng2025SSPCA"]] <- paste(
  "Duplicate bibliography record for arXiv:2409.13267; BibTeX year is 2024."
)
notes[["ZhaoWangFeng2026SPCA"]] <- paste(
  "Duplicate bibliography record for arXiv:2409.13267; BibTeX year is 2024."
)
notes[["XuMaWangFeng2026EllipticalFactor"]] <- "Citation-key year is 2026; references.bib year is 2025."
notes[["YanFengZhang2025InverseNormMaxsum"]] <- "Citation key says Zhang; CSL authors are Yan, Zhao, and Feng."
notes[["MaFengWang2025DependentAlpha"]] <- "Public-version author/title metadata differs from references.bib; see LITERATURE_AUDIT.md."

for (i in seq_len(nrow(verified_book))) {
  comparison <- verified_book$metadata_comparison[[i]]
  if (!identical(comparison, "book metadata agrees")) {
    key <- verified_book$citation_key[[i]]
    notes[[key]] <- trimws(paste(
      notes[[key]], paste0("Independent verification: ", comparison)
    ))
  }
}

section_for_key <- function(key) {
  paste(names(citations)[vapply(citations, function(x) key %in% x$unique, logical(1L))],
        collapse = ";")
}

occurrences_for_key <- function(key) {
  sum(vapply(citations, function(x) sum(x$keys == key), integer(1L)))
}

build_expected <- function() {
  rows <- lapply(bbl_keys, function(key) {
    reference <- references[[key]]
    stable <- if (key %in% names(stable_ids)) stable_ids[[key]] else ""
    reference_text <- paste(unlist(reference, recursive = TRUE, use.names = FALSE),
                            collapse = " ")
    local_manuscript <- grepl("manuscript|preprint", reference_text,
                              ignore.case = TRUE)
    status <- if (nzchar(stable)) {
      "stable_identifier"
    } else if (identical(key, "WangWangFeng2026GSPCA")) {
      "no_stable_identifier"
    } else if (local_manuscript) {
      "local_source"
    } else {
      "book_metadata_only"
    }
    source <- if (nzchar(stable)) {
      if (key %in% names(id_sources)) {
        id_sources[[key]]
      } else "references.bib"
    } else if (status == "local_source") {
      "local_source"
    } else if (identical(key, "WangWangFeng2026GSPCA")) {
      "book_formula_source"
    } else {
      "references.bib+main.bbl"
    }
    mapped <- method_map[[key]] %||% character()
    data.frame(
      citation_key = key,
      chapters = section_for_key(key),
      citation_occurrences = as.character(occurrences_for_key(key)),
      bib_present = "TRUE",
      bbl_present = "TRUE",
      authors = format_authors(reference),
      year = reference_year(reference),
      title = collapse_space(reference$title %||% ""),
      reference_type = bib_types[[key]],
      stable_identifier = stable,
      identifier_status = status,
      evidence_source = source,
      method_trace_rows = paste(mapped, collapse = ";"),
      relation = if (length(mapped)) "method_primary" else "background_or_review",
      notes = notes[[key]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[] <- lapply(out, as.character)
  out
}

write_csv_utf8 <- function(x, path) {
  write.table(
    x, file = path, sep = ",", row.names = FALSE, col.names = TRUE,
    quote = TRUE, qmethod = "double", na = "", fileEncoding = "UTF-8",
    eol = "\n"
  )
}

expected <- build_expected()
trailing_args <- commandArgs(trailingOnly = TRUE)
if ("--write" %in% trailing_args) {
  write_csv_utf8(expected, paths$bibliography)
  cat(sprintf("WROTE %s (%d rows)\n", paths$bibliography, nrow(expected)))
}

invisible(read_utf8(paths$bibliography, reject_bom = TRUE))
trace <- read.csv(
  paths$bibliography, stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8", colClasses = "character", na.strings = character()
)
check(nrow(trace) == 194L, "BIBLIOGRAPHY_TRACEABILITY.csv must contain 194 rows.")
check(identical(names(trace), names(expected)),
      "Bibliography traceability columns differ from the required schema.")
check(!anyDuplicated(trace$citation_key), "Duplicate citation_key in bibliography CSV.")
check(identical(trace$citation_key, bbl_keys),
      "Bibliography CSV key order must exactly follow main.bbl.")
check(isTRUE(all(trace == expected)),
      "Bibliography CSV differs from source-derived expected values; rerun with --write and review the diff.")
check(all(trace$bib_present == "TRUE" & trace$bbl_present == "TRUE"),
      "Every bibliography trace row must be present in both BibTeX and BBL.")
check(all(nzchar(trace$authors) & nzchar(trace$year) & nzchar(trace$title) &
            nzchar(trace$reference_type)),
      "Every row must contain authors, year, title, and reference_type.")
check(!any(grepl("[\r\n]", as.matrix(trace), perl = TRUE)),
      "CSV fields must not contain embedded line breaks.")

allowed_status <- c("stable_identifier", "book_metadata_only", "local_source",
                    "no_stable_identifier")
check(all(trace$identifier_status %in% allowed_status),
      "Unknown identifier_status in bibliography CSV.")
has_identifier <- nzchar(trace$stable_identifier)
check(all(has_identifier == (trace$identifier_status == "stable_identifier")),
      "stable_identifier presence and identifier_status disagree.")
identifier_pattern <- paste0(
  "^(doi:10\\.[0-9]{4,9}/[^[:space:]]+|",
  "arXiv:[0-9]{4}\\.[0-9]{4,5}(v[0-9]+)?|url:https?://[^[:space:]]+)$"
)
check(all(grepl(identifier_pattern, trace$stable_identifier[has_identifier],
                perl = TRUE)),
      "Malformed stable_identifier value(s): %s",
      paste(trace$stable_identifier[has_identifier &
              !grepl(identifier_pattern, trace$stable_identifier, perl = TRUE)],
            collapse = ";"))

for (key in names(verified_overrides)) {
  row <- trace[trace$citation_key == key, , drop = FALSE]
  check(nrow(row) == 1L && row$stable_identifier == verified_overrides[[key]],
        "Verified identifier mismatch for %s.", key)
}
gspca <- trace[trace$citation_key == "WangWangFeng2026GSPCA", , drop = FALSE]
check(nrow(gspca) == 1L && !nzchar(gspca$stable_identifier) &&
        gspca$identifier_status == "no_stable_identifier" &&
        gspca$evidence_source == "book_formula_source" &&
        gspca$relation == "background_or_review",
      "WangWangFeng2026GSPCA must remain a book-attributed unavailable-source boundary.")
conditional_alpha <- trace[
  trace$citation_key == "ZhaoWang2024ConditionalMaxAlpha", , drop = FALSE
]
check(grepl("public arXiv paper is 2026", conditional_alpha$notes, fixed = TRUE),
      "ZhaoWang2024ConditionalMaxAlpha must retain the key/year mismatch note.")

check(all(method_keys %in% trace$citation_key),
      "Method primary key(s) missing from bibliography CSV: %s",
      paste(setdiff(method_keys, trace$citation_key), collapse = ";"))
for (key in trace$citation_key) {
  expected_rows <- paste(method_map[[key]] %||% character(), collapse = ";")
  observed <- trace$method_trace_rows[trace$citation_key == key]
  check(identical(observed, expected_rows), "Method-row mapping mismatch for %s.", key)
  expected_relation <- if (nzchar(expected_rows)) "method_primary" else "background_or_review"
  check(trace$relation[trace$citation_key == key] == expected_relation,
        "Relation mismatch for %s.", key)
}
method_primary <- trace$relation == "method_primary"
check(all(trace$identifier_status[method_primary] == "stable_identifier"),
      "Every method-primary bibliography row must have an independently auditable stable identifier.")
check(all(!trace$evidence_source[method_primary] %in% c(
  "references.bib+main.bbl", "book_formula_source", "local_source"
)), "Method-primary evidence may not rely only on book/internal metadata.")

roundtrip_path <- tempfile("HDElliptical-bibliography-roundtrip-", fileext = ".csv")
on.exit(unlink(roundtrip_path, force = TRUE), add = TRUE)
write_csv_utf8(trace, roundtrip_path)
roundtrip <- read.csv(roundtrip_path, stringsAsFactors = FALSE, check.names = FALSE,
                      fileEncoding = "UTF-8", colClasses = "character",
                      na.strings = character())
check(identical(trace, roundtrip), "Strict CSV UTF-8 round-trip failed.")

pandoc_csv_path <- tempfile("HDElliptical-bibliography-csv-", fileext = ".json")
on.exit(unlink(pandoc_csv_path, force = TRUE), add = TRUE)
run_pandoc(c(
  "-f", "csv", "-t", "json", shQuote(paths$bibliography),
  "-o", shQuote(pandoc_csv_path)
), "Pandoc CSV parse gate")
check(file.exists(pandoc_csv_path) && file.info(pandoc_csv_path)$size > 0L,
      "Pandoc CSV gate produced no JSON output.")

status_counts <- table(factor(trace$identifier_status, levels = allowed_status))
method_mapped <- sum(trace$relation == "method_primary")
no_stable <- trace$citation_key[trace$identifier_status == "no_stable_identifier"]
cat(sprintf(
  paste0(
    "BIBLIOGRAPHY_TRACEABILITY OK: rows=%d; chapter_union=%d; ",
    "appendix_unique=%d; manuscript_union=%d; bib=%d; bbl=%d; ",
    "stable=%d; book_metadata_only=%d; local_source=%d; ",
    "no_stable_identifier=%d; method_primary=%d; background_or_review=%d; ",
    "unused_bib=%d; pandoc_bib=PASS; pandoc_csv=PASS; utf8_csv=PASS.\n"
  ),
  nrow(trace), length(chapter_union), length(appendix_keys),
  length(manuscript_union), length(bib_keys), length(bbl_keys),
  status_counts[["stable_identifier"]], status_counts[["book_metadata_only"]],
  status_counts[["local_source"]], status_counts[["no_stable_identifier"]],
  method_mapped, nrow(trace) - method_mapped, length(unused_bib)
))
cat("UNCITED_BIB=", paste(unused_bib, collapse = ";"), "\n", sep = "")
cat("NO_STABLE_IDENTIFIER=", paste(no_stable, collapse = ";"), "\n", sep = "")
