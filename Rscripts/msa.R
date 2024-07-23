# Load necessary libraries
library(msa)
library(Biostrings)
library(dplyr)
library(purrr)
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to parse command-line arguments
parse_args <- function(args) {
  options <- list(
    input = NULL,
    output = NULL,
    group_index = NULL,
    sequence = NULL,
    seq_index = NULL,
    type = "dna"  # Default type is "dna"
  )
  
  for (i in seq_along(args)) {
    if (args[i] == "--input" || args[i] == "-i") {
      options$input <- args[i + 1]
    } else if (args[i] == "--output" || args[i] == "-o") {
      options$output <- args[i + 1]
    } else if (args[i] == "--group_index" || args[i] == "-g") {
      options$group_index <- args[i + 1]
    } else if (args[i] == "--sequence" || args[i] == "-s") {
      options$sequence <- args[i + 1]
    } else if (args[i] == "--seq_index" || args[i] == "-x") {
      options$seq_index <- args[i + 1]
    } else if (args[i] == "--type" || args[i] == "-t") {
      options$type <- args[i + 1]
    }
  }
  
  if (is.null(options$input) || is.null(options$output) || is.null(options$group_index) ||
      is.null(options$sequence) || is.null(options$seq_index)) {
        print(options)
    stop("Usage: Rscript script.R --input=inputpath --output=outputpath --group_index=group_index_column_name --sequence=sequence_column_name --seq_index=seq_index_column_name [--type=dna/aa]")
  }
  
  return(options)
}

# Parse arguments
opt <- parse_args(args)

# Function to perform MSA and return only msa_result
calMSA <- function(seqs, type){
  if (length(seqs) <=1){
      msa_result = as.character(seqs)
  }
  else if (type == "dna") {
    # Convert sequences to DNAStringSet
    string_set <- DNAStringSet(seqs)
    # Run MSA for DNA sequences
    msa_result <- msa(string_set, type = "dna", order = "input")
  } else if (type == "aa") {
    # Convert sequences to AAStringSet
    string_set <- AAStringSet(seqs)
    # Run MSA for amino acid sequences
    msa_result <- msa(string_set, type = "protein", order = "input")
  } else {
    stop("Unsupported sequence type. Please use 'dna' or 'aa'.")
  }
  
  # Return only msa_result
  as.character(msa_result)
  }

# Read input data
data <- read.csv(opt$input, header = TRUE)

# Perform operations on data
df <- data %>%
  group_by(!!sym(opt$group_index)) %>%
  summarize(sequences = list(data.frame(!!sym(opt$sequence), !!sym(opt$seq_index))))

# Initialize an empty list to store results
result <- list()

# Apply the calMSA function to all sequences and store in result list
result$sequences_aligned <- purrr::map(df$sequences, ~ calMSA(.x[[opt$sequence]], opt$type))

# Extract seq indexes
seq_indexes <- unlist(lapply(df$sequences, function(x) x[[opt$seq_index]]))

# Match and assign aligned sequences
aligned_sequences <- unlist(result$sequences_aligned)

data$aligned_sequence <- aligned_sequences[match(data[[opt$seq_index]], seq_indexes)]

# Write output to CSV file
write.csv(data, file = opt$output, row.names = FALSE)