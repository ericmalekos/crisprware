#!/usr/bin/env Rscript
# 
# Ex: 
#   crisprscore.R <path_to_sgrna_bed_file> 3
#   crisprscore.R --help
#

args <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat("\n\tUsage: crisprscore.R <path_to_sgrna_bed_file> <method_number> <outputfile> ... <additional settings> \n")
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

  cat("\tAdditional settings for method 3 (DeepHF):\n")
  cat("\t--enzyme: Specify the enzyme (options: 'WT', 'ESP', 'HF')\n")
  cat("\t--promoter: Specify the promoter (options: 'U6', 'T7')\n")
  cat("\tExample: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 3 DeepHF_scored_sgRNAs.bed WT U6\n\n")

  cat("\tAdditional setting for method 5 (DeepCpf1):\n")
  cat("\t--no-convertPAM: Specify whether non-canonical PAMs are converted to TTTC [default: TRUE]\n")
  cat("\tExample: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 5 DeepCpf1_scored_sgRNAs.bed no-convertPAM\n\n")
  
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

if (method_number == 12) {
  if (length(args) != 4) {
    cat("\n\t\tERROR: incorrect number of arguments for method 12\n\n")
    print_help()
    quit(status = 1)
  }
  tracrRNA = args[4]
  if (!(tracrRNA %in% c('Hsu2013', 'Chen2013', 'HF'))) {
    cat("\n\t\tERROR: Invalid tracr\n\n")
    print_help()
    quit(status = 1)
  }
}

print(args)

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

colnames(df) <- c("#chr", "start", "stop", "id.sequence.pam.chromosome.position.sense", "context", "strand")

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

#print(head(df,5))

score_func <- method_info$func

if (method_number == 3) {
  df <- df %>%
    mutate(!!paste0(selected_method_name, "_score") := score_func(context, enzyme = enzyme, promoter = promoter)$score)
} else if (method_number == 5) {
  df <- df %>%
    mutate(!!paste0(selected_method_name, "_score") := score_func(context, convertPAM = convertPAM)$score)
} else if (method_number == 12) {
  df <- df %>%
    mutate(!!paste0(selected_method_name, "_score") := score_func(context, tracrRNA = tracrRNA)$score)
} else {
  df <- df %>%
    mutate(!!paste0(selected_method_name, "_score") := score_func(context)$score)
}

#print(head(df,5))

output_path <- "./sgRNAs/scoredSgRNAs/"

# Check if the directory exists, if not create it
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Assuming 'df' is your data frame
output_file <- paste0(output_path, output)

# Save the data frame as a TSV file
write.table(df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)


# TESTING PACKAGE
# library("crisprScore")
# packageVersion("crisprScore")
# 
# spacer  <- "ATCGATGCTGATGCTAGATA" #20bp
# pam     <- "AGG" #3bp
# input   <- paste0(spacer, pam)
# results <- getDeepHFScores(input)
