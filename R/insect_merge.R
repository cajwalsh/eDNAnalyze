#' Merge LCA and insect taxonomic assignment data
#' @param lca A data frame where each row corresponds to an OTU with at least 9 columns including:
#' * taxonomic information across 8 levels in descending order (domain, kingdom, phylum, class, order, family, genus, species)
#' * OTU name (named "OTU" or "zotu")
#' * optional: information from blast (unique_hits, taxid, and pid if lca_progression was used)
#' @param insect_path The path to the insect data file (default = "insect_classified.csv)
#' @returns A data frame containing taxonomic assignments from domain to species and information from either blast or insect (whichever assignment was kept) for all OTUs
#' @examples
#' lca_insect_merge(lca = lca, taxdir = "~/Desktop", insect_path = "~/Desktop/results/insect_results.csv")
#' @description Merges BLAST LCA and insect taxonomic assignment data for all OTUs into one data frame. If an OTU was assigned taxonomy by the LCA process (through BLAST), the insect and assignment related paratemers are nullified so the assignment method is clear. PLEASE NOTE THAT THE NCBI TAXONOMY DIRECTORY IS OVER 2GB IN SIZE. Make sure your computer has space for this and refer to this file path each time this function is used to save storage space.
#' @export
lca_insect_merge = function(lca, insect_path = "insect_taxonomy.tsv") {

  ## Load insect results from specified file path
  insect = read.table(insect_path,
                      header = T,
                      sep = "\t")[,c("zotu",
                                     "domain", "kingdom", "phylum", "class",
                                     "order", "family", "genus", "species",
                                     "taxon", "taxid")] # maybe add more columns to be kept here if anyone wants
  
  ## Assign taxonomic information to certain taxa that get missed
  insect[which(insect$taxon=="Rhodophyta"),"phylum"] = "Rhodophyta"
  insect[which(insect$taxon=="PX clade"),"class"] = "Phaeophyceae|Xanthophyceae"

  ## Rename first column to OTU to merge both taxonomy dataframes
  if("zotu" %in% names(insect)) names(insect)[which(names(insect) == "zotu")] = "OTU"
  ids <- merge(lca, insect, by = "OTU", all = T) %>%
    dplyr::mutate(
      taxid = dplyr::coalesce(taxid.x, taxid.y),
      domain = dplyr::coalesce(domain.x, domain.y),
      kingdom = dplyr::coalesce(kingdom.x, kingdom.y),
      phylum = dplyr::coalesce(phylum.x, phylum.y),
      class = dplyr::coalesce(class.x, class.y),
      order = dplyr::coalesce(order.x, order.y),
      family = dplyr::coalesce(family.x, family.y),
      genus = dplyr::coalesce(genus.x, genus.y),
      species = dplyr::coalesce(species.x, species.y)
    ) %>%
    dplyr::select(-ends_with(".x"), -ends_with(".y"))

  ## Remove possibly discordant insect classification information/stats for OTUs identified using BLAST
  ids[!is.na(ids$unique_hits),c("taxon")] = NA

  ## Replace NAs in blast hit column with 0s (as if there had been a blast hit, it would be classified by LCA and not insect)
  ids[is.na(ids$unique_hits),"unique_hits"] = 0

  ## Arrange columns in more logical order
  if("pid" %in% names(ids)) {
  ids = ids[,c("domain", "kingdom", "phylum", "class",
               "order", "family", "genus", "species",
               "OTU", "unique_hits", "pid", "taxon", "taxid")]
  } else {
    ids = ids[,c("domain", "kingdom", "phylum", "class",
                 "order", "family", "genus", "species",
                 "OTU", "unique_hits", "taxon", "taxid")]
  }
}
