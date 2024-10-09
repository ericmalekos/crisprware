#!/usr/bin/env Rscript
# 
# Ex: 
#   crisprscore.R <path_to_sgrna_bed_file> 3 <path_to_output_file> [--chunk-size 100000]
#   crisprscore.R --help
#

args <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat("\n\tUsage: crisprscore.R <path_to_sgrna_bed_file> <method_number> <outputfile> [<additional settings> ... ] \n")
  cat("\tExample: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 2 Azimuth_scored_sgRNAs.bed\n")
  cat("\tMethods:\n")
  cat("\t1:  RuleSet1 - SpCas9 (Length: 30)\n")
  cat("\t2:  Azimuth - SpCas9 (Length: 30)\n")
  cat("\t3:  DeepHF - SpCas9 (Length: 23)\n")
  cat("\t4:  Lindel - SpCas9 (Length: 65)\n")
  cat("\t5:  DeepCpf1 - AsCas12a (Length: 34)\n")
  cat("\t6:  EnPAMGB - enAsCas12a (Length: 34)\n")
  cat("\t7:  CRISPRscan - SpCas9 (Length: 35)\n")
  cat("\t8:  CasRx-RF - CasRx (Length: NA)\n")
  cat("\t9:  CRISPRai - SpCas9 (Length: 22)\n")
  cat("\t10: CRISPRater - SpCas9 (Length: 20)\n")
  cat("\t11: DeepSpCas9 - SpCas9 (Length: 30)\n")
  cat("\t12: RuleSet3 - SpCas9 (Length: 30)\n\n")

  cat("\n\tAdditional optional argument for chunk size:\n")
  cat("\t--chunk-size: Specify the size of chunks for processing the dataframe (optional)\n")
  cat("\tExample with chunk size: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 2 Azimuth_scored_sgRNAs.bed --chunk-size 10000\n\n")
  
  cat("\tAdditional settings for method 3 (DeepHF):\n")
  cat("\tenzyme: Specify the enzyme (options: 'WT', 'ESP', 'HF')\n")
  cat("\tpromoter: Specify the promoter (options: 'U6', 'T7')\n")
  cat("\tExample: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 3 DeepHF_scored_sgRNAs.bed WT U6\n\n")

  cat("\tAdditional setting for method 5 (DeepCpf1):\n")
  cat("\t--no-convertPAM: Specify whether non-canonical PAMs are converted to TTTC [default: TRUE]\n")
  cat("\tExample: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 5 DeepCpf1_scored_sgRNAs.bed --no-convertPAM\n\n")
  
  cat("\tAdditional setting for method 12 (RuleSet3):\n")
  cat("\t--tracrRNA: Specify tracrRNA (options: 'Hsu2013', 'Chen2013')\n\n")
  cat("\tExample: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 12 rs3_scored_sgRNAs.bed Chen2013\n\n")
#  }
}

# Check for help argument or incorrect number of arguments
if (length(args) < 3 || args[1] == "--help") {
  cat("\n\t\tERROR: Missing arguments\n\n")
  print_help()
  quit(status = 0)
} 

# Extract arguments
sgrna_file <- args[1]
method_number <- as.integer(args[2])
output <- args[3]

#print(args)
  
chunk_size <- Inf
if ("--chunk-size" %in% args) {
  chunk_size_index <- match("--chunk-size", args) + 1
  chunk_size <- as.numeric(args[chunk_size_index])
  args <- args[1:(length(args) - 2)]  # Remove chunk size arguments from args
}

#print(args)


if (method_number == 3) {
  if (length(args) != 5) {
    cat("\n\t\tERROR: incorrect number of arguments for method 3\n\n")
    print_help()
    quit(status = 1)
  }
  enzyme <- args[4]
  promoter <- args[5]
  
  # Check if enzyme and promoter are valid
  if (!(enzyme %in% c('WT', 'ESP', 'HF') && promoter %in% c('U6', 'T7'))) {
    cat("\n\t\tERROR: Invalid enzyme or promoter\n\n")
    print_help()
    quit(status = 1)
  }
}


tracrRNA = ""
#print(length(args))
if (method_number == 12) {
  if (length(args) != 4) {
    cat("\n\t\tERROR: incorrect number of arguments for method 12\n\n")
    print_help()
    quit(status = 1)
  }
  tracrRNA = args[4]
  if (!(tracrRNA %in% c('Hsu2013', 'Chen2013'))) {
    cat("\n\t\tERROR: Invalid tracr\n\n")
    print_help()
    quit(status = 1)
  }
}


convertPAM <- TRUE
if (method_number == 5 && "no-convertPAM" %in% args) {
    convertPAM <- FALSE
}

# Check if the method number is valid
if (method_number < 1 || method_number > 12) {
  cat("\n\t\tERROR: Invalid method number\n\n")
  print_help()
  quit(status = 1)
}

# Load the TSV file (assuming the TSV file has a header)


df <- read.table(sgrna_file, 
                 sep = "\t",
                 header = TRUE,
                 comment.char = "")

#colnames(df) <- c("#chr", "start", "stop", "id.sequence.pam.chromosome.position.sense", "context", "strand")

if (chunk_size == Inf) {
  chunk_size <- nrow(df)
}

check_context_length <- function(df, expected_length) {
  if (!is.na(expected_length)) {
    # Check if all context lengths are equal to the expected length
    return(all(nchar(as.character(df$context)) == expected_length))
  }
  TRUE
}

suppressPackageStartupMessages(library(crisprScore, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))

# Function mapping (method number to function name and expected length)
method_map <- list(
  list(func = getRuleSet1Scores, length = 30), #y
  list(func = getAzimuthScores, length = 30), #y
  list(func = getDeepHFScores, length = 23), #N
  list(func = getLindelScores, length = 65), #y
  list(func = getDeepCpf1Scores, length = 34), #N
  list(func = getEnPAMGBScores, length = 34), #y
  list(func = getCRISPRscanScores, length = 35), #y
  list(func = getCasRxRFScores, length = NA), # NA for variable length or unknown
  list(func = getCrispraiScores, length = 22), #NA
  list(func = getCRISPRaterScores, length = 20), #y
  list(func = getDeepSpCas9Scores, length = 30), #y
  list(func = getRuleSet3Scores, length = 30)
)

method_info <- method_map[[method_number]]

method_names <- c("RuleSet1", "Azimuth", "DeepHF", "Lindel", "DeepCpf1",
                  "EnPAMGB", "CRISPRscan", "CasRxRF", "CRISPRai",
                  "CRISPRater", "DeepSpCas9", "RuleSet3")

# # Get the method name based on the method number
selected_method_name <- method_names[method_number]


if (method_number %in% c(8, 9)) {
  cat("\n\tFunctionality for the selected method (", method_number, ") is not available \n\n\n")
  quit(status = 0)
}

if (!check_context_length(df, method_info$length)) {
  cat("\n\t\tERROR: Context length does not match the expected length for the selected method\n\n\n")
  print_help()
  quit(status = 1)
}

cat("Processing", sgrna_file, "with method number", method_number, "\n")


apply_scoring_to_chunk <- function(chunk_df) {
  method_info <- method_map[[method_number]]
  score_func <- method_info$func
  selected_method_name <- method_names[method_number]
  

  
  if (!check_context_length(chunk_df, method_info$length)) {
    stop("Context length does not match the expected length for the selected method")
  }
  
  # Dynamically call the scoring function based on method_number
  # Append the scores as a new column named after the selected method
  score_column_name <- paste0(selected_method_name, "_score")
  
  # Assuming all scoring functions return a dataframe with a 'score' column
  # scored_chunk_df <- chunk_df %>%
  #   mutate(!!score_column_name := score_func(context)$score)
  # 
  if (method_number == 3) {
    scored_chunk_df <- chunk_df %>%
      mutate(!!paste0(selected_method_name, "_score") := score_func(context, enzyme = enzyme, promoter = promoter)$score)
  } else if (method_number == 5) {
    scored_chunk_df <- chunk_df %>%
      mutate(!!paste0(selected_method_name, "_score") := score_func(context, convertPAM = convertPAM)$score)
  } else if (method_number == 12) {
    scored_chunk_df <- chunk_df %>%
      mutate(!!paste0(selected_method_name, "_score") := score_func(context, tracrRNA = tracrRNA)$score)
  } else {
    scored_chunk_df <- chunk_df %>%
      mutate(!!paste0(selected_method_name, "_score") := score_func(context)$score)
  }
  
  return(scored_chunk_df)
}

# Load the dataframe
df <- read.table(sgrna_file, sep = "\t", header = TRUE, comment.char = "")

# Initialize an empty dataframe for the results
results_df <- df[0, ]

#print(chunk_size)

# Process the dataframe in chunks
for (i in seq(1, nrow(df), by = chunk_size)) {
  #print(i)
  chunk <- df[i:min(i + chunk_size - 1, nrow(df)), ]
  scored_chunk <- apply_scoring_to_chunk(chunk)
  results_df <- rbind(results_df, scored_chunk)
}

new_column_names <- c("#chr", "start", "stop", "id,sequence,pam,chromosome,position,sense", "context", "strand")
names(results_df)[1:6] <- new_column_names
# round the score column to 4 decimals
results_df[, ncol(results_df)] <- round(results_df[, ncol(results_df)], 4)
write.table(results_df, file = output, sep = "\t", row.names = FALSE, quote = FALSE)



# TESTING PACKAGE
# library("crisprScore")
# packageVersion("crisprScore")
# 
# spacer  <- "ATCGATGCTGATGCTAGATA" #20bp
# pam     <- "AGG" #3bp
# input   <- paste0(spacer, pam)
# results <- getDeepHFScores(input)
# 
# flank5 <- "ACCT" #4bp
# spacer <- "ATCGATGCTGATGCTAGATA" #20bp
# pam    <- "AGG" #3bp 
# flank3 <- "TTG" #3bp
# input  <- paste0(flank5, spacer, pam, flank3) 
# results <- getAzimuthScores(input)

# flank5 <- "ACC" #3bp
# pam    <- "TTTT" #4bp
# spacer <- "AATCGATGCTGATGCTAGATATT" #23bp
# flank3 <- "AAGT" #4bp
# input  <- paste0(flank5, pam, spacer, flank3) 
# results <- getDeepCpf1Scores(input)
# results
# 
# flank5 <- "ACC" #3bp
# pam    <- "TTTT" #4bp
# spacer <- "AATCGATGCTGATGCTAGATATT" #23bp
# flank3 <- "AAGT" #4bp
# input  <- paste0(flank5, pam, spacer, flank3) 
# results <- getEnPAMGBScores(input)
# results