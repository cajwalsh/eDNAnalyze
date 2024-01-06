#' Merge multiple eDNAFlow LCA taxonomic assignment tables for more confident taxonomic assignments
#'
#' @param glob_to_pid A character string that specifies the file naming convention for each OTU table with different lca parameters up to the percent identity (pid) value (which needs to be in the file name)
#' @param dir A character string that specifies the directory in which the lca taxonomy (identified in glob_to_pid) files are found (default = getwd())
#' @param lulu A boolean to specify whether you want lulu curated or non-lulu curated files if both can be found in the directory with the same glob_to_pid (default = NA)
#' @returns A data frame with one row for every OTU identified at any pid level with 11 columns including:
#' * taxonomic information across 8 levels in descending order (domain, kingdom, phylum, class, order, family, genus, species)
#' * OTU name (named "OTU")
#' * the number of unique blast hits (numberOfUnq_BlastHits)
#' * the highest percent identity level at which the OTU was classified (pid)
#' @examples
#' lca_progression()
#' lca_progression(glob_to_pid = "taxonomy_", dir = "~/Desktop/results", lulu = T)
#' @description Merges a series of similarly named OTU tables generated from running the eDNAFlow LCA taxonomy assignment process on the same dataset with different minimum percent identity (pid) thresholds to keep only the highest level match per OTU. The main point of this is to provide the approximate pid match for each OTU that can be examined alongside the number of unique blast hits to determine whether there may be a dubious assignment from unclear, incomplete, or incorrect database listings.
#' @export
lca_progression <- function(glob_to_pid = "taxonomy_q[0-9]+_p", dir = getwd(), lulu = NA) {

  ## List the files that contain the user-defined taxonomy table prefixes
  ## and make a big empty one to fill with results from all levels combined
  tables <- list.files(dir, pattern = glob_to_pid)
  if(is.na(lulu)) {
    tables = tables
  } else if(lulu == T) {
    tables = tables[grep("lulu", tables)]
  } else if(lulu == F) {
    tables = tables[-grep("lulu", tables)]
  }
  big_table <- data.frame()

  for(table in tables) {
    ## Define the table name and percent identity from the naming convention
    if(glob_to_pid=="taxonomy_q[0-9]+_p") {
      tname <- gsub("taxonomy_q[0-9]+_|_d.*", "", table)
      pid <- as.numeric(gsub("p", "", tname))
    } else {
      tname <- sub("*.tab|*.txt", "", table)
      pid <- as.numeric(gsub(glob_to_pid, "", tname))
    }
    ## Read the table in, add pid as a column, append them together
    assign(paste(tname), read.table(paste0(dir, "/", table), sep = "\t", header = T), )
    assign(paste(tname), cbind(get(tname),pid = rep(pid, nrow(get(tname)))))
    big_table <- rbind(big_table, get(tname))
  }

  ## Keep only the highest match level for each OTU
  final_table <- big_table %>%
    dplyr::group_by(OTU) %>%
    dplyr::slice_max(pid, with_ties = FALSE)
  non_sample_cols = c("domain", "kingdom", "phylum","class",
                      "order", "family", "genus", "species",
                      "OTU", "numberOfUnq_BlastHits", "pid")
  final_table <- final_table[,c(which(names(final_table) %in% non_sample_cols),
                                which(!names(final_table) %in% non_sample_cols))]
}
