#!/usr/bin/env Rscript

# Clean protein FASTA files by retaining only sequences with the 20 standard
# amino acids: ACDEFGHIKLMNPQRSTVWY.

allowed_alphabet <- "ACDEFGHIKLMNPQRSTVWY"
allowed_chars <- strsplit(allowed_alphabet, "", fixed = TRUE)[[1]]
default_wrap <- 60

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
if (is.na(script_path)) {
  script_path <- normalizePath("clean_standard_amino_acids.R", mustWork = FALSE)
}
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
default_output_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
default_input_root <- normalizePath(file.path(default_output_dir, ".."), winslash = "/", mustWork = FALSE)

trailing_args <- commandArgs(trailingOnly = TRUE)

get_option_value <- function(flag, default_value) {
  prefix <- paste0(flag, "=")
  match <- trailing_args[startsWith(trailing_args, prefix)]
  if (length(match) > 0) {
    return(sub(prefix, "", match[1], fixed = TRUE))
  }
  index <- match(flag, trailing_args)
  if (!is.na(index) && index < length(trailing_args)) {
    return(trailing_args[index + 1])
  }
  default_value
}

has_flag <- function(flag) {
  flag %in% trailing_args
}

input_root <- normalizePath(
  get_option_value("--input-root", default_input_root),
  winslash = "/",
  mustWork = FALSE
)
output_dir <- normalizePath(
  get_option_value("--output-dir", default_output_dir),
  winslash = "/",
  mustWork = FALSE
)
wrap_width <- as.integer(get_option_value("--wrap", as.character(default_wrap)))
dry_run <- has_flag("--dry-run")

if (is.na(wrap_width) || wrap_width < 1) {
  stop("--wrap must be a positive integer", call. = FALSE)
}

numeric_prefix <- function(path) {
  name <- basename(path)
  prefix <- sub("^(\\d+)\\..*$", "\\1", name)
  if (identical(prefix, name)) {
    return(10000L)
  }
  as.integer(prefix)
}

protein_folders <- function(root) {
  folders <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  folders <- folders[grepl("^\\d+\\.\\s+", basename(folders))]
  folders[order(vapply(folders, numeric_prefix, integer(1)), tolower(basename(folders)))]
}

find_input_fasta <- function(folder) {
  fasta_files <- list.files(
    folder,
    pattern = "\\.(fa|faa|fasta|fna)$",
    ignore.case = TRUE,
    full.names = TRUE,
    recursive = FALSE
  )
  if (length(fasta_files) == 0) {
    stop(paste("No FASTA file found in", folder), call. = FALSE)
  }
  if (length(fasta_files) > 1) {
    stop(
      paste("More than one FASTA file found in", folder, ":", paste(basename(fasta_files), collapse = ", ")),
      call. = FALSE
    )
  }
  fasta_files[1]
}

read_fasta <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  header_lines <- grep("^>", lines)

  if (length(header_lines) == 0) {
    return(data.frame(header = character(), sequence = character(), stringsAsFactors = FALSE))
  }

  records <- vector("list", length(header_lines))
  for (i in seq_along(header_lines)) {
    start <- header_lines[i]
    end <- if (i < length(header_lines)) header_lines[i + 1] - 1 else length(lines)
    header <- sub("^>", "", lines[start])
    if (start < end) {
      sequence_lines <- lines[(start + 1):end]
    } else {
      sequence_lines <- character()
    }
    sequence_lines <- sequence_lines[!grepl("^>", sequence_lines)]
    sequence <- toupper(paste(gsub("\\s+", "", sequence_lines), collapse = ""))
    records[[i]] <- data.frame(header = header, sequence = sequence, stringsAsFactors = FALSE)
  }

  do.call(rbind, records)
}

accession_from_header <- function(header) {
  parts <- strsplit(header, "\\s+")[[1]]
  if (length(parts) == 0) {
    return("")
  }
  parts[1]
}

invalid_characters <- function(sequence) {
  chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
  invalid <- chars[!(chars %in% allowed_chars)]
  paste(unique(invalid), collapse = "")
}

wrap_sequence <- function(sequence, width) {
  starts <- seq(1, nchar(sequence), by = width)
  substring(sequence, starts, pmin(starts + width - 1, nchar(sequence)))
}

write_fasta <- function(path, records, width) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)

  for (i in seq_len(nrow(records))) {
    writeLines(paste0(">", records$header[i]), connection)
    writeLines(wrap_sequence(records$sequence[i], width), connection)
  }
}

empty_log_frame <- function() {
  data.frame(
    protein_folder = character(),
    input_file = character(),
    cleaned_file = character(),
    record_number = integer(),
    accession = character(),
    header = character(),
    sequence_length = integer(),
    invalid_characters = character(),
    decision = character(),
    stringsAsFactors = FALSE
  )
}

summary_rows <- data.frame(
  protein_folder = character(),
  input_file = character(),
  cleaned_file = character(),
  total_sequences = integer(),
  retained_sequences = integer(),
  removed_sequences = integer(),
  allowed_alphabet = character(),
  invalid_characters_detected = character(),
  removed_accessions = character(),
  stringsAsFactors = FALSE
)
record_rows <- empty_log_frame()
removed_rows <- empty_log_frame()

for (folder in protein_folders(input_root)) {
  input_fasta <- find_input_fasta(folder)
  records <- read_fasta(input_fasta)
  input_name <- basename(input_fasta)
  cleaned_name <- paste0(tools::file_path_sans_ext(input_name), "_cleaned.fasta")
  cleaned_relative <- paste("cleaned_fasta", basename(folder), cleaned_name, sep = "/")

  if (nrow(records) == 0) {
    kept_records <- records
    invalid_by_record <- character()
  } else {
    invalid_by_record <- vapply(records$sequence, invalid_characters, character(1))
    kept_records <- records[invalid_by_record == "", , drop = FALSE]
  }

  accessions <- if (nrow(records) == 0) {
    character()
  } else {
    vapply(records$header, accession_from_header, character(1))
  }
  decisions <- ifelse(invalid_by_record == "", "retained", "removed")
  sequence_lengths <- nchar(records$sequence)

  if (nrow(records) > 0) {
    folder_record_rows <- data.frame(
      protein_folder = basename(folder),
      input_file = input_name,
      cleaned_file = cleaned_relative,
      record_number = seq_len(nrow(records)),
      accession = accessions,
      header = records$header,
      sequence_length = sequence_lengths,
      invalid_characters = invalid_by_record,
      decision = decisions,
      stringsAsFactors = FALSE
    )
    record_rows <- rbind(record_rows, folder_record_rows)
    removed_rows <- rbind(removed_rows, folder_record_rows[folder_record_rows$decision == "removed", , drop = FALSE])
  }

  removed_accessions <- accessions[invalid_by_record != ""]
  detected_invalid <- paste(unique(strsplit(paste(invalid_by_record, collapse = ""), "", fixed = TRUE)[[1]]), collapse = "")

  summary_rows <- rbind(
    summary_rows,
    data.frame(
      protein_folder = basename(folder),
      input_file = input_name,
      cleaned_file = cleaned_relative,
      total_sequences = nrow(records),
      retained_sequences = nrow(kept_records),
      removed_sequences = length(removed_accessions),
      allowed_alphabet = allowed_alphabet,
      invalid_characters_detected = detected_invalid,
      removed_accessions = paste(removed_accessions, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )

  if (!dry_run) {
    write_fasta(file.path(output_dir, cleaned_relative), kept_records, wrap_width)
  }
}

if (dry_run) {
  for (i in seq_len(nrow(summary_rows))) {
    message(
      summary_rows$protein_folder[i],
      ": total=", summary_rows$total_sequences[i],
      ", retained=", summary_rows$retained_sequences[i],
      ", removed=", summary_rows$removed_sequences[i]
    )
  }
  quit(status = 0)
}

reports_dir <- file.path(output_dir, "reports")
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  summary_rows,
  file.path(reports_dir, "fasta_cleaning_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  record_rows,
  file.path(reports_dir, "fasta_cleaning_record_log.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  removed_rows,
  file.path(reports_dir, "fasta_cleaning_removed_sequences.csv"),
  row.names = FALSE,
  quote = TRUE
)

message("Cleaning complete: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
