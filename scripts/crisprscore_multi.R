#!/usr/bin/env Rscript
# 
# Ex: 
#   crisprscore.R <path_to_sgrna_bed_file> <comma_separated_method_numbers> <path_to_output_file> [--chunk-size 100000]
#   crisprscore.R --help
#

args <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat("\n\tUsage: crisprscore.R <path_to_sgrna_bed_file> <comma_separated_method_numbers> <outputfile> <enzyme> <5prime_flank_length> <3prime_flank_length> [--chunk-size <size>] [--debug]\n")
  cat("\tExample: crisprscore.R input.tsv 1,2,5 output_scored.tsv Cas9 25 25\n")
  cat("\tExample: crisprscore.R input.tsv 15,16,17 output_scored.tsv Cas12a 30 15 --chunk-size 10000\n")
  cat("\tExample with debug: crisprscore.R input.tsv 1,2,11 output_scored.tsv Cas9 25 25 --debug\n\n")
  
  cat("\tInput TSV format:\n")
  cat("\tRequired column: 'context' (can be in any position)\n")
  cat("\tAll other columns are preserved in the output\n\n")
  
  cat("\tEnzyme types:\n")
  cat("\tCas9 - Use for SpCas9-based scoring methods (1-14)\n")
  cat("\tCas12a - Use for AsCas12a-based scoring methods (15-17)\n\n")
  
  cat("\tContext format:\n")
  cat("\tFor Cas9: [5' flank] + [20nt spacer] + [3nt PAM] + [3' flank]\n")
  cat("\tFor Cas12a: [5' flank] + [4nt PAM] + [23nt spacer] + [3' flank]\n")
  cat("\tSpecify the lengths of your 5' and 3' flanks as arguments\n")
  cat("\tThe script will automatically trim to the required length for each scoring method\n\n")
  
  cat("\tScoring Methods for Cas9:\n")
  cat("\t1:  RuleSet1 - SpCas9 (Length: 30)\n")
  cat("\t2:  Azimuth - SpCas9 (Length: 30)\n")
  cat("\t3:  DeepHF_WT_U6 - SpCas9 (Length: 23)\n")
  cat("\t4:  DeepHF_WT_T7 - SpCas9 (Length: 23)\n")
  cat("\t5:  DeepHF_ESP_U6 - SpCas9 (Length: 23)\n")
  cat("\t6:  DeepHF_ESP_T7 - SpCas9 (Length: 23)\n")
  cat("\t7:  DeepHF_HF_U6 - SpCas9 (Length: 23)\n")
  cat("\t8:  DeepHF_HF_T7 - SpCas9 (Length: 23)\n")
  cat("\t9:  Lindel - SpCas9 (Length: 65)\n")
  cat("\t10: CRISPRscan - SpCas9 (Length: 35)\n")
  cat("\t11: CRISPRater - SpCas9 (Length: 20, spacer only)\n")
  cat("\t12: DeepSpCas9 - SpCas9 (Length: 30)\n")
  cat("\t13: RuleSet3_Hsu2013 - SpCas9 (Length: 30)\n")
  cat("\t14: RuleSet3_Chen2013 - SpCas9 (Length: 30)\n\n")
  
  cat("\tScoring Methods for Cas12a:\n")
  cat("\t15: DeepCpf1 - AsCas12a (Length: 34, canonical PAM conversion)\n")
  cat("\t16: DeepCpf1_noConvert - AsCas12a (Length: 34, no PAM conversion)\n")
  cat("\t17: EnPAMGB - enAsCas12a (Length: 34)\n\n")
  
  cat("\tNote: CasRx-RF and CRISPRai methods are not currently available\n\n")
  
  cat("\tOptional arguments:\n")
  cat("\t--chunk-size <size>: Process dataframe in chunks of specified size (default: entire file)\n")
  cat("\t--debug: Show detailed trimming information for each method (shows what sequence is sent to each scoring function)\n\n")
}

# Check for help argument or incorrect number of arguments
if (length(args) < 6 || args[1] == "--help" || args[1] == "-h") {
  if (length(args) > 0 && args[1] != "--help" && args[1] != "-h") {
    cat("\n\t\tERROR: Missing arguments\n\n")
  }
  print_help()
  quit(status = 0)
} 

# Extract main arguments
sgrna_file <- args[1]
method_numbers_str <- args[2]
output <- args[3]
enzyme <- tolower(args[4])  # Convert to lowercase for case-insensitive comparison
flank_5prime <- as.integer(args[5])
flank_3prime <- as.integer(args[6])

# Validate enzyme type
if (!enzyme %in% c("cas9", "cas12a")) {
  cat("\n\t\tERROR: Invalid enzyme type. Must be 'Cas9' or 'Cas12a' (case-insensitive).\n\n")
  print_help()
  quit(status = 1)
}

# Validate flank lengths
if (is.na(flank_5prime) || is.na(flank_3prime) || flank_5prime < 0 || flank_3prime < 0) {
  cat("\n\t\tERROR: Invalid flank lengths. Must be non-negative integers.\n\n")
  print_help()
  quit(status = 1)
}

# Parse chunk size if provided
chunk_size <- Inf
if ("--chunk-size" %in% args) {
  chunk_size_index <- match("--chunk-size", args) + 1
  if (chunk_size_index <= length(args)) {
    chunk_size <- as.numeric(args[chunk_size_index])
    if (is.na(chunk_size)) {
      cat("\n\t\tERROR: Invalid chunk size\n\n")
      print_help()
      quit(status = 1)
    }
  } else {
    cat("\n\t\tERROR: Missing value for --chunk-size\n\n")
    print_help()
    quit(status = 1)
  }
}

# Parse debug flag
debug_mode <- "--debug" %in% args
if (debug_mode) {
  cat("DEBUG MODE ENABLED\n\n")
}

# Parse method numbers
method_numbers <- as.integer(unlist(strsplit(method_numbers_str, ",")))

# Validate method numbers
if (any(is.na(method_numbers)) || any(method_numbers < 1) || any(method_numbers > 17)) {
  cat("\n\t\tERROR: Invalid method number(s)\n\n")
  print_help()
  quit(status = 1)
}

# Validate method numbers match enzyme type
cas9_methods <- 1:14
cas12a_methods <- 15:17

if (enzyme == "cas9") {
  invalid_methods <- method_numbers[method_numbers %in% cas12a_methods]
  if (length(invalid_methods) > 0) {
    cat("\n\t\tERROR: Methods", paste(invalid_methods, collapse = ", "), 
        "are for Cas12a but enzyme type is Cas9\n\n")
    cat("\t\tFor Cas9, use methods 1-14\n\n")
    quit(status = 1)
  }
} else if (enzyme == "cas12a") {
  invalid_methods <- method_numbers[method_numbers %in% cas9_methods]
  if (length(invalid_methods) > 0) {
    cat("\n\t\tERROR: Methods", paste(invalid_methods, collapse = ", "), 
        "are for Cas9 but enzyme type is Cas12a\n\n")
    cat("\t\tFor Cas12a, use methods 15-17\n\n")
    quit(status = 1)
  }
}

# Load required libraries
suppressPackageStartupMessages(library(crisprScore, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))

# Define method mapping with all parameter combinations
# Each method specifies: name, function, total length needed, and flanking requirements
method_map <- list(
  list(name = "RuleSet1", func = getRuleSet1Scores, length = 30, 
       flank_5prime = 4, spacer = 20, pam = 3, flank_3prime = 3, params = list()),
  list(name = "Azimuth", func = getAzimuthScores, length = 30, 
       flank_5prime = 4, spacer = 20, pam = 3, flank_3prime = 3, params = list()),
  list(name = "DeepHF_WT_U6", func = getDeepHFScores, length = 23, 
       flank_5prime = 0, spacer = 20, pam = 3, flank_3prime = 0, params = list(enzyme = "WT", promoter = "U6")),
  list(name = "DeepHF_WT_T7", func = getDeepHFScores, length = 23, 
       flank_5prime = 0, spacer = 20, pam = 3, flank_3prime = 0, params = list(enzyme = "WT", promoter = "T7")),
  list(name = "DeepHF_ESP_U6", func = getDeepHFScores, length = 23, 
       flank_5prime = 0, spacer = 20, pam = 3, flank_3prime = 0, params = list(enzyme = "ESP", promoter = "U6")),
  list(name = "DeepHF_ESP_T7", func = getDeepHFScores, length = 23, 
       flank_5prime = 0, spacer = 20, pam = 3, flank_3prime = 0, params = list(enzyme = "ESP", promoter = "T7")),
  list(name = "DeepHF_HF_U6", func = getDeepHFScores, length = 23, 
       flank_5prime = 0, spacer = 20, pam = 3, flank_3prime = 0, params = list(enzyme = "HF", promoter = "U6")),
  list(name = "DeepHF_HF_T7", func = getDeepHFScores, length = 23, 
       flank_5prime = 0, spacer = 20, pam = 3, flank_3prime = 0, params = list(enzyme = "HF", promoter = "T7")),
  list(name = "Lindel", func = getLindelScores, length = 65, 
       flank_5prime = 13, spacer = 20, pam = 3, flank_3prime = 29, params = list()),
  list(name = "CRISPRscan", func = getCRISPRscanScores, length = 35, 
       flank_5prime = 6, spacer = 20, pam = 3, flank_3prime = 6, params = list()),
  list(name = "CRISPRater", func = getCRISPRaterScores, length = 20, 
       flank_5prime = 0, spacer = 20, pam = 0, flank_3prime = 0, spacer_only = TRUE, params = list()),
  list(name = "DeepSpCas9", func = getDeepSpCas9Scores, length = 30, 
       flank_5prime = 4, spacer = 20, pam = 3, flank_3prime = 3, params = list()),
  list(name = "RuleSet3_Hsu2013", func = getRuleSet3Scores, length = 30, 
       flank_5prime = 4, spacer = 20, pam = 3, flank_3prime = 3, params = list(tracrRNA = "Hsu2013")),
  list(name = "RuleSet3_Chen2013", func = getRuleSet3Scores, length = 30, 
       flank_5prime = 4, spacer = 20, pam = 3, flank_3prime = 3, params = list(tracrRNA = "Chen2013")),
  list(name = "DeepCpf1", func = getDeepCpf1Scores, length = 34, 
       flank_5prime = 4, spacer = 23, pam = 4, flank_3prime = 3, pam_position = "5prime", params = list(convertPAM = TRUE)),
  list(name = "DeepCpf1_noConvert", func = getDeepCpf1Scores, length = 34, 
       flank_5prime = 4, spacer = 23, pam = 4, flank_3prime = 3, pam_position = "5prime", params = list(convertPAM = FALSE)),
  list(name = "EnPAMGB", func = getEnPAMGBScores, length = 34, 
       flank_5prime = 4, spacer = 23, pam = 4, flank_3prime = 3, pam_position = "5prime", params = list())
)

# Read the input file
cat("Loading input file:", sgrna_file, "\n")
df <- read.table(sgrna_file, 
                 sep = "\t",
                 header = TRUE,
                 comment.char = "",
                 stringsAsFactors = FALSE,
                 quote = "",
                 check.names = FALSE)  # This preserves original column names

# Store original column names
original_colnames <- colnames(df)

# Check that the context column is present
if (!"context" %in% names(df)) {
  cat("\n\t\tERROR: Input file must have a 'context' column\n\n")
  quit(status = 1)
}

# Report on the columns found
cat("Found", ncol(df), "columns:", paste(names(df), collapse = ", "), "\n")

# Calculate expected context length based on user input and enzyme type
if (enzyme == "cas9") {
  # For SpCas9: 5' flank + 20nt spacer + 3nt PAM + 3' flank
  expected_context_length <- flank_5prime + 20 + 3 + flank_3prime
  cat("Enzyme: Cas9\n")
  cat("Expected context length:", expected_context_length, "nt (", 
      flank_5prime, " + 20 + 3 + ", flank_3prime, ")\n", sep="")
} else {
  # For Cas12a: 5' flank + 4nt PAM + 23nt spacer + 3' flank
  expected_context_length <- flank_5prime + 4 + 23 + flank_3prime
  cat("Enzyme: Cas12a\n")
  cat("Expected context length:", expected_context_length, "nt (", 
      flank_5prime, " + 4 + 23 + ", flank_3prime, ")\n", sep="")
}

cat("User specified flanking lengths: 5'=", flank_5prime, "nt, 3'=", flank_3prime, "nt\n", sep="")

# Validate context lengths in the input file
context_lengths <- nchar(as.character(df$context))
unique_lengths <- unique(context_lengths)

if (length(unique_lengths) > 1) {
  cat("WARNING: Multiple context lengths found:", paste(unique_lengths, collapse = ", "), "\n")
}

# Check if contexts match expected length for the specified enzyme
if (!any(unique_lengths == expected_context_length)) {
  cat("\n\t\tERROR: Context lengths in file (", paste(unique_lengths, collapse = ", "), 
      " nt) don't match expected length\n", sep="")
  if (enzyme == "cas9") {
    cat("\t\tFor Cas9: ", flank_5prime, " + 20 + 3 + ", flank_3prime, " = ", 
        expected_context_length, " nt\n", sep="")
  } else {
    cat("\t\tFor Cas12a: ", flank_5prime, " + 4 + 23 + ", flank_3prime, " = ", 
        expected_context_length, " nt\n", sep="")
  }
  quit(status = 1)
}

# Function to trim context to required length for a specific method
trim_context <- function(context, method_info, user_flank_5prime, user_flank_3prime, debug = FALSE) {
  context_str <- as.character(context)
  context_length <- nchar(context_str)
  
  # Special handling for CRISPRater - only wants the 20nt spacer
  if (!is.null(method_info$spacer_only) && method_info$spacer_only == TRUE) {
    # For CRISPRater: extract just the 20nt spacer, no PAM, no flanking
    spacer_start <- user_flank_5prime + 1
    spacer_end <- spacer_start + 19  # 20nt spacer
    
    if (spacer_end > context_length) {
      return(NA)
    }
    
    spacer_seq <- substr(context_str, spacer_start, spacer_end)
    
    if (debug) {
      cat("    Original context (", context_length, "nt): ", context_str, "\n", sep="")
      cat("    Extracted spacer only (20nt): ", spacer_seq, "\n", sep="")
    }
    
    return(spacer_seq)
  }
  
  # Determine the enzyme type based on PAM position
  pam_position <- ifelse(is.null(method_info$pam_position), "3prime", method_info$pam_position)
  
  if (pam_position == "5prime") {
    # Cas12a-type: 5' flank + PAM + spacer + 3' flank
    user_total_length <- user_flank_5prime + method_info$pam + method_info$spacer + user_flank_3prime
    
    if (context_length != user_total_length) {
      return(NA)  # Context doesn't match expected structure
    }
    
    # Calculate positions in the user's context
    pam_start <- user_flank_5prime + 1
    spacer_start <- pam_start + method_info$pam
    spacer_end <- spacer_start + method_info$spacer - 1
    
    # Extract what we need
    if (method_info$flank_5prime > 0) {
      required_5prime_start <- max(1, pam_start - method_info$flank_5prime)
      flank_5prime_seq <- substr(context_str, required_5prime_start, pam_start - 1)
    } else {
      flank_5prime_seq <- ""
    }
    
    pam_seq <- substr(context_str, pam_start, spacer_start - 1)
    spacer_seq <- substr(context_str, spacer_start, spacer_end)
    
    if (method_info$flank_3prime > 0) {
      required_3prime_end <- min(context_length, spacer_end + method_info$flank_3prime)
      flank_3prime_seq <- substr(context_str, spacer_end + 1, required_3prime_end)
    } else {
      flank_3prime_seq <- ""
    }
    
    trimmed <- paste0(flank_5prime_seq, pam_seq, spacer_seq, flank_3prime_seq)
    
    if (debug) {
      cat("    Original context (", context_length, "nt): ", context_str, "\n", sep="")
      cat("    Structure: [", user_flank_5prime, "nt 5'flank][4nt PAM][23nt spacer][", 
          user_flank_3prime, "nt 3'flank]\n", sep="")
      cat("    Trimmed for method (", nchar(trimmed), "nt): ", trimmed, "\n", sep="")
      cat("      5' flank (", nchar(flank_5prime_seq), "nt): ", flank_5prime_seq, "\n", sep="")
      cat("      PAM (", nchar(pam_seq), "nt): ", pam_seq, "\n", sep="")
      cat("      Spacer (", nchar(spacer_seq), "nt): ", spacer_seq, "\n", sep="")
      cat("      3' flank (", nchar(flank_3prime_seq), "nt): ", flank_3prime_seq, "\n", sep="")
    }
    
  } else {
    # SpCas9-type: 5' flank + spacer + PAM + 3' flank
    user_total_length <- user_flank_5prime + method_info$spacer + method_info$pam + user_flank_3prime
    
    if (context_length != user_total_length) {
      return(NA)  # Context doesn't match expected structure
    }
    
    # Calculate positions in the user's context
    spacer_start <- user_flank_5prime + 1
    spacer_end <- spacer_start + method_info$spacer - 1
    pam_start <- spacer_end + 1
    pam_end <- pam_start + method_info$pam - 1
    
    # Extract what we need based on method requirements
    if (method_info$flank_5prime > 0) {
      required_5prime_start <- max(1, spacer_start - method_info$flank_5prime)
      flank_5prime_seq <- substr(context_str, required_5prime_start, spacer_start - 1)
    } else {
      flank_5prime_seq <- ""
    }
    
    spacer_seq <- substr(context_str, spacer_start, spacer_end)
    
    if (method_info$pam > 0) {
      pam_seq <- substr(context_str, pam_start, pam_end)
    } else {
      pam_seq <- ""
    }
    
    if (method_info$flank_3prime > 0) {
      required_3prime_end <- min(context_length, pam_end + method_info$flank_3prime)
      flank_3prime_seq <- substr(context_str, pam_end + 1, required_3prime_end)
    } else {
      flank_3prime_seq <- ""
    }
    
    trimmed <- paste0(flank_5prime_seq, spacer_seq, pam_seq, flank_3prime_seq)
    
    if (debug) {
      cat("    Original context (", context_length, "nt): ", context_str, "\n", sep="")
      cat("    Structure: [", user_flank_5prime, "nt 5'flank][20nt spacer][3nt PAM][", 
          user_flank_3prime, "nt 3'flank]\n", sep="")
      cat("    Trimmed for method (", nchar(trimmed), "nt): ", trimmed, "\n", sep="")
      cat("      5' flank (", nchar(flank_5prime_seq), "nt): ", flank_5prime_seq, "\n", sep="")
      cat("      Spacer (", nchar(spacer_seq), "nt): ", spacer_seq, "\n", sep="")
      cat("      PAM (", nchar(pam_seq), "nt): ", pam_seq, "\n", sep="")
      cat("      3' flank (", nchar(flank_3prime_seq), "nt): ", flank_3prime_seq, "\n", sep="")
    }
  }
  
  # Verify the trimmed length matches what the method expects
  if (nchar(trimmed) != method_info$length) {
    return(NA)
  }
  
  return(trimmed)
}

# Pre-validate that all selected methods can work with the provided contexts
cat("\nValidating context compatibility with selected methods...\n")
for (method_num in method_numbers) {
  method_info <- method_map[[method_num]]
  
  if (debug_mode) {
    cat("\n  Testing method", method_num, ":", method_info$name, "\n")
  }
  
  # Test trim on first few contexts
  test_contexts <- head(df$context, min(10, nrow(df)))
  trimmed_test <- sapply(test_contexts, function(ctx) {
    trim_context(ctx, method_info, flank_5prime, flank_3prime, debug = debug_mode && which(test_contexts == ctx)[1] == 1)
  })
  
  if (all(is.na(trimmed_test))) {
    cat("\n\t\tERROR: Cannot trim contexts for method '", method_info$name, "'\n", sep="")
    cat("\t\tMethod requires: ", method_info$flank_5prime, "nt 5' flank, ",
        method_info$spacer, "nt spacer, ", method_info$pam, "nt PAM, ",
        method_info$flank_3prime, "nt 3' flank (total: ", method_info$length, "nt)\n", sep="")
    cat("\t\tYour contexts have: ", flank_5prime, "nt 5' flank, ",
        flank_3prime, "nt 3' flank\n", sep="")
    cat("\t\tContext cannot be trimmed to fit this method's requirements.\n")
    quit(status = 1)
  }
  
  cat("  ✓", method_info$name, "- context can be trimmed from", unique_lengths[1], 
      "nt to", method_info$length, "nt\n")
}

# Set chunk size to full dataframe if not specified
if (chunk_size == Inf) {
  chunk_size <- nrow(df)
}

# Function to apply a single scoring method to a chunk
apply_single_score <- function(chunk_df, method_num, debug = FALSE) {
  method_info <- method_map[[method_num]]
  
  if (debug && nrow(chunk_df) > 0) {
    cat("\n  Applying", method_info$name, "to", nrow(chunk_df), "sequences\n")
    cat("  Method requirements: ", method_info$flank_5prime, "nt 5'flank + ",
        method_info$spacer, "nt spacer + ", method_info$pam, "nt PAM + ",
        method_info$flank_3prime, "nt 3'flank = ", method_info$length, "nt total\n", sep="")
    
    # Show first sequence trimming in detail
    cat("  First sequence trimming:\n")
  }
  
  # Trim contexts for this method
  trimmed_contexts <- sapply(seq_along(chunk_df$context), function(i) {
    ctx <- chunk_df$context[i]
    # Only show debug for first sequence to avoid too much output
    show_debug <- debug && i == 1
    trim_context(ctx, method_info, flank_5prime, flank_3prime, debug = show_debug)
  })
  
  # Check if any contexts couldn't be trimmed
  if (any(is.na(trimmed_contexts))) {
    warning(paste("Some contexts could not be trimmed for method", method_info$name, 
                  "- these will have NA scores"))
    # Replace NA with empty string for scoring functions
    trimmed_contexts[is.na(trimmed_contexts)] <- ""
  }
  
  # Apply the scoring function with appropriate parameters
  score_func <- method_info$func
  params <- method_info$params
  
  tryCatch({
    if (length(params) == 0) {
      # No additional parameters
      scores <- score_func(trimmed_contexts)$score
    } else {
      # Call with additional parameters
      scores <- do.call(score_func, c(list(trimmed_contexts), params))$score
    }
    
    # Set NA for contexts that couldn't be trimmed
    scores[trimmed_contexts == ""] <- NA
    
    return(scores)
  }, error = function(e) {
    warning(paste("Error applying", method_info$name, ":", e$message))
    return(rep(NA, nrow(chunk_df)))
  })
}

# Function to apply all selected scoring methods to a chunk
apply_scoring_to_chunk <- function(chunk_df, selected_methods, debug = FALSE) {
  scored_chunk <- chunk_df
  
  for (method_num in selected_methods) {
    method_info <- method_map[[method_num]]
    score_column_name <- paste0(method_info$name, "_score")
    
    if (!debug) {
      cat("  Applying", method_info$name, "to chunk...\n")
    }
    
    scores <- apply_single_score(chunk_df, method_num, debug = debug)
    scored_chunk[[score_column_name]] <- round(scores, 4)
  }
  
  return(scored_chunk)
}

# Process the dataframe
cat("Processing with methods:", paste(method_numbers, collapse = ", "), "\n")
cat("Total rows:", nrow(df), "\n")

# Initialize results dataframe
results_df <- data.frame()

# Process in chunks
total_chunks <- ceiling(nrow(df) / chunk_size)
chunk_count <- 0

for (i in seq(1, nrow(df), by = chunk_size)) {
  chunk_count <- chunk_count + 1
  end_row <- min(i + chunk_size - 1, nrow(df))
  
  cat("Processing chunk", chunk_count, "of", total_chunks, "(rows", i, "to", end_row, ")\n")
  
  chunk <- df[i:end_row, ]
  scored_chunk <- apply_scoring_to_chunk(chunk, method_numbers, debug = debug_mode)
  
  if (nrow(results_df) == 0) {
    results_df <- scored_chunk
  } else {
    results_df <- rbind(results_df, scored_chunk)
  }
}

# Write the output
cat("Writing results to:", output, "\n")

# Restore original column names for non-score columns
num_original_cols <- length(original_colnames)
names(results_df)[1:num_original_cols] <- original_colnames

# Special handling for #chr if it was mangled
if (grepl("chr", names(results_df)[1], ignore.case = TRUE) && !startsWith(names(results_df)[1], "#")) {
  names(results_df)[1] <- "#chr"
}

write.table(results_df, file = output, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Done! Processed", nrow(results_df), "rows with", length(method_numbers), "scoring method(s)\n")
cat("Output saved to:", output, "\n")