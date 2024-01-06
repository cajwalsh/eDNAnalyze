#' Fill in blank taxonomic assignments
#'
#' @param df A data frame where each row corresponds to an OTU with ONLY 8 columns that contain taxonomic information in descending order (domain, kingdom, phylum, class, order, family, genus, species)
#' @param dropped A character string assigned to lower taxonomic levels that couldn't be resolved for an OTU (default = "dropped)
#' @param unassigned A character string filled into higher order taxa for which the assigned taxon has no assignment (default = NA)
#' @returns The original data frame without blanks
#' @examples
#' fill_blank_taxa(df)
#' fill_blank_taxa(df = COI_taxa, dropped = "unknown", unassigned = "unassigned")
#' @description Fills in blank taxonomic information for OTUs by distinguishing between dropped and unassigned blanks.
#' * Dropped taxonomic levels are defined here as lower level assignments that couldn't be resolved (e.g., an OTU that could only be resolved as a fish (class Actinopteri) with no clear order or lower level assignements would have "dropped" for order-species).
#' * Unassigned taxonomic levels are defined here as higher order taxa to which a lower level taxonomic assignment has no assignment (e.g., an OTU assigned to the family Pomacentridae would be "unassigned" for order as their status at that leve is uncertain).
#' * Writing "NA" (in quotes) for dropped or unassigned will give actual NA values (unassigned = "NA" is the default).
#' @export

fill_blank_taxa = function(df, dropped = "dropped", unassigned = "NA") {
  if(dropped=="NA") dropped = NA
  if(unassigned=="NA") unassigned = NA
  ## Check every row in a data frame for empty or NA taxon values
  for(x in 1:nrow(df)) {
    if("" %in% df[x,] |
       NA %in% df[x,]) {
      ## If all taxonomic levels are empty or NA, replace all of these with dropped
      if(all(df[x,] %in% c("", NA))) {
        df[x,] = dropped
      } else {
        ## If only some are, find the lowest taxonomic assignment level
        lowest_assignment = max(which(!df[x,] %in% c("", "dropped")))
        ## If there is missing information below this level, fill it with dropped
        if(lowest_assignment < 8) {
          df[x,(lowest_assignment+1):8] = dropped
        }
        ## If there is missing information above this level, fill it with NA
        if("" %in% df[x,] |
           NA %in% df[x,]) {
          df[x, c(which(df[x,]==""), which(is.na(df[x,])))] = unassigned
        }
      }
    }
  }
  return(df)
}
