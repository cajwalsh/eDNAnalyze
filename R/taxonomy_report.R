#' Generate a report of the taxonomic diversity detected within a target group
#'
#' @param df A data frame where each row corresponds to an OTU with columns containing: taxonomic information across 8 levels in descending order
#' @param target_taxon A character string corresponding to the taxon for which a breakdown of nested taxa will be returned
#' @param taxon_rank A character string corresponding to the taxonomic level of the target taxon (one of the eight levels from the columns of df)
#' @param dropped A character string assigned to lower taxonomic levels that couldn't be resolved for an OTU (default = "dropped)
#' @param rm_humans A boolean to specify whether or not human OTUs should be included or removed (default = F)
#' @param file_prefix A character string specifying the prefix of the output file (default = target_taxon)
#' @param outdir A character string specifying the output directory of the text file (default = getwd())
#' @returns A text file listing the names of taxa nested within the target taxon, as well as the number of OTUs assigned to each of the nested taxa
#' @examples
#' eDNA_taxon_report(df = df, target_taxon = "Metazoa", taxon_rank = "kingdom")
#' eDNA_taxon_report(df = coi_edna_data, file_prefix = "COI_eDNA_animal", target_taxon = "Metazoa",\
#'                   dropped = "unknown", rm_humans = T, outdir = "~/Desktop/results")
#' @description Generates a list-style report containing lots of taxonomic information by OTU. The first section provides information about how many taxa were dropped at each level for
#' all OTUs/taxa. The next section does the same thing but only for OTUs in the target taxon. Finally, for each taxonomic level below the rank of the target taxon, a list of all taxa
#' found at this level and the number of OTUs that correspond to each taxon is provided. Unless primates or hominids are target taxa, this function treats all sequences identified as
#' Primates as human (in rm_humans).
#' @export

taxonomy_report <- function(df,
                            target_taxon,
                            taxon_rank,
                            dropped = "dropped",
                            rm_humans = F,
                            file_prefix = target_taxon,
                            outdir = getwd()) {

  ## Remove human sequences if requested
  if(rm_humans==T & target_taxon=="Primates" |
     rm_humans==T & target_taxon=="Hominidae" ) {
    df <- droplevels(df[which(df$genus!="Homo"),])
  } else if(rm_humans==T) {
    df <-  droplevels(df[which(df$order!="Primates"),])
  }

  ## Save console output to object
  output <- capture.output({

    ## Overall taxonomy look
    df <- df[,-which(names(df)=="OTU")]
    if(rm_humans == T) cat("All taxa assignment breakdown (minus human)\n")
    if(rm_humans == F) cat("All taxa assignment breakdown\n")
    for(a in names(df)[1:8]) {
      cat(nrow(df[which(df[,a]==dropped),]), "of", nrow(df),
          "total zOTUs dropped at", a, "level\n")
    }

    ## Target taxonomy look
    target_df <- droplevels(df[which(df[,taxon_rank]==target_taxon),])
    cat(sep = "", "\nTarget subset taxonomic breakdown (within ", target_taxon, ")\n")
    for(b in names(target_df)[1:8]) {
      cat(nrow(target_df[which(target_df[,b]==dropped),]), "of", nrow(target_df),
          "target zOTUs dropped at", b, "level\n")
    }
    ## Produce taxonomy table at every level lower than the target_taxon taxon_rank
    c <- which(names(df)==taxon_rank)+1
    while (c < 9) {
      down1 = names(df)[c]
      down1_table <- as.data.frame(target_df %>% dplyr::count(target_df[,down1], sort = T))
      names(down1_table) <- c(down1, "n")
      cat("\n")
      print(down1_table, row.names = F)
      c <- c + 1
    }
  })
  ## Save console output object to text file
  sink(paste0(outdir, "/", file_prefix, "_taxonomy.txt"))
  cat(output, sep = "\n")
  sink()
}
