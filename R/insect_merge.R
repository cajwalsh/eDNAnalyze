## Downloads the NCBI taxonomy database that is used here to give insect data domain information
get_ncbi <- function(taxdir = "taxdump") {

  ## Download NCBI taxonomy dump if it does not exist already
  if (!fs::file_exists(fs::path(taxdir,"rankedlineage.dmp"))) {
    url <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"

    ## create directory (fails silently if it exists)
    fs::dir_create(taxdir,recurse = TRUE)

    ## download taxonomy dump
    cat(paste0("Downloading NCBI taxonomy files to specified directory: ", taxdir, "\n"))
    resp <- httr::GET(url, httr::write_disk(fs::path(taxdir,"taxdump.zip"), overwrite = T))
    if (resp$status_code == 200) {
      ## unzip it
      unzip(fs::path(taxdir,"taxdump.zip"),exdir=taxdir)
    } else {
      return(NULL)
    }
  }

  ## load ranked lineage
  cat("Loading NCBI ranked lineage file into R")
  readr::read_tsv(
    fs::path(taxdir,"rankedlineage.dmp"),
    col_types = "i_c_c_c_c_c_c_c_c_c_",
    col_names = c("taxid","species","fakespecies","genus","family","order","class","phylum","kingdom","domain")
  ) %>% dplyr::select(taxid,domain,kingdom,phylum,class,order,family,genus,species)
}


#' Merge LCA and insect taxonomic assignment data
#' @param lca_taxa A data frame where each row corresponds to an OTU with 10 or 11 columns including:
#' * taxonomic information across 8 levels in descending order (domain, kingdom, phylum, class, order, family, genus, species)
#' * OTU name (named "OTU")
#' * information from blast (numberOfUnq_BlastHits, and pid if lca_progression was used)
#' @param taxdir The path to the directory to download (or that already contains) the NCBI taxonomy dump (THIS IS OVER 2 GB IN SIZE)
#' @param insect_path The path to the insect data file (default = "insect_classified.csv)
#' @returns A data frame containing taxonomic assignments from domain to species and information from either blast or insect (whichever assignment was kept) for all OTUs
#' @examples
#' lca_insect_merge(lca_taxa = lca[,1:11], taxdir = "~/Desktop", insect_path = "~/Desktop/results/insect_results.csv")
#' @description Merges BLAST LCA and insect taxonomic assignment data for all OTUs into one data frame. If an OTU was assigned taxonomy by the LCA process (through BLAST), the insect and assignment related paratemers are nullified so the assignment method is clear. PLEASE NOTE THAT THE NCBI TAXONOMY DIRECTORY IS OVER 2GB IN SIZE. Make sure your computer has space for this and refer to this file path each time this function is used to save storage space.
#' @export
lca_insect_merge = function(lca_taxa, taxdir, insect_path = "insect_classified.csv") {

  ## Load insect results from specified file path
  insect = read.csv(insect_path, header = T)[,c(1:3, 5:12)] # maybe add more columns to be kept here if anyone wants

  ## Load NCBI taxonomy database
  ncbi = get_ncbi(taxdir)

  ## Merge NCBI taxonomy database with insect results
  insect_ncbi <- insect %>%
    dplyr::left_join(ncbi %>% dplyr::select(taxid,domain), by = c("taxID" = "taxid"))

  ## Assign taxonomic information to certain taxa that get missed
  insect_ncbi[which(insect_ncbi$taxon=="Eukaryota"),"domain"] = "Eukaryota"
  insect_ncbi[which(insect_ncbi$insect_ncbi=="Rhodophyta"),"phylum"] = "Rhodophyta"
  insect_ncbi[which(insect_ncbi$taxon=="PX clade"),"class"] = "Phaeophyceae|Xanthophyceae"

  ## Rename first column to OTU to merge both taxonomy dataframes
  names(insect_ncbi)[1] = "OTU"
  ids <- merge(lca_taxa, insect_ncbi, by = "OTU", all = T) %>%
    dplyr::mutate(
      domain = dplyr::coalesce(domain.x, domain.y),
      kingdom = dplyr::coalesce(kingdom.x, kingdom.y),
      phylum = dplyr::coalesce(phylum.x, phylum.y),
      class = dplyr::coalesce(class.x, class.y),
      order = dplyr::coalesce(order.x, order.y),
      family = dplyr::coalesce(family.x, family.y),
      genus = dplyr::coalesce(genus.x, genus.y),
      species = dplyr::coalesce(species.x, species.y),
    ) %>%
    dplyr::select(-ends_with(".x"), -ends_with(".y"))

  ## Remove possibly discordant insect classification information/stats for OTUs identified using BLAST
  ids[!is.na(ids$numberOfUnq_BlastHits),c("taxon", "score")] = NA

  ## Replace NAs in blast hit column with 0s (as if there had been a blast hit, it would be classified by LCA and not insect)
  ids[is.na(ids$numberOfUnq_BlastHits),"numberOfUnq_BlastHits"] = 0

  ## Arrange columns in more logical order
  if("pid" %in% names(ids)) {
  ids = ids[,c("domain", "kingdom", "phylum", "class",
               "order", "family", "genus", "species",
               "OTU", "numberOfUnq_BlastHits", "pid", "taxon", "score")]
  } else {
    ids = ids[,c("domain", "kingdom", "phylum", "class",
                 "order", "family", "genus", "species",
                 "OTU", "numberOfUnq_BlastHits", "taxon", "score")]
  }
}
