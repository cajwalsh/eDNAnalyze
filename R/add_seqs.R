#' Add zotu sequences and their lengths to a zotu table
#'
#' @param fasta A character string specifying the path to/name of the zotu fasta file
#' @param df A data frame where each row corresponds to a zotu and includes a column named "zotu" to match names with the fasta file
#' @param new_cols A number specifying where the new columns should be inserted into df
#' @returns The original data frame with columns containing the DNA sequence (df$sequence) and its length (df$seq_length) for all zotus
#' @examples
#' add_seqs(fasta = "~/Desktop/results/zotus.fasta", df = COI_taxa, new_cols = 9)
#' @description Extracts the DNA sequence and its length for each zotu from a fasta file and adds them to two columns of the zotu table data frame
#' @export
add_seqs <- function(fasta, df, new_cols) {

  ## Read the fasta file and format it correctly
  sequences = ape::read.dna(fasta, format = "fasta", as.character = T)
  fasta_sequences = sapply(sequences, function(seq) paste(toupper(seq), collapse = ""))

  ## Get the zotu values from your dataframe
  zotu_values <- df$zotu

  ## Initialize new columns
  seq_list <- list()
  seq_list$sequence = seq_list$seq_length = rep(NA, length(zotu_values))

  ## Iterate over each zotu value
  for(z in 1:length(zotu_values)) {
    ## Find the index of the matching fasta sequence
    index <- match(zotu_values[z], names(fasta_sequences))

    if(!is.na(index)) {
      ## Retrieve the sequence
      sequence <- as.character(fasta_sequences[index])

      ## Update the dataframe with the sequence and its length
      seq_list$sequence[z] <- sequence
      seq_list$seq_length[z] <- nchar(sequence)
    }
  }
  newdf = cbind(df[1:new_cols-1],
                sequence = seq_list$sequence,
                seq_length = seq_list$seq_length,
                df[(new_cols):length(df)])
}
