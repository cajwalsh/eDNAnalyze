#' Generate reports of reads and OTUs by sample
#'
#' @param df A data frame where each row corresponds to an OTU with columns containing at least:
#'  * taxonomic information across 8 levels (domain, kingdom, phylum, class, order, family, genus, species);
#'  * read counts per OTU from all samples sequenced (one column per sample) after all taxonomic and other data columns
#' @param first_sample_col A number denoting which column the sample read counts start (and continue to the end of the data frame)
#' @param taxa_mat A three-column matrix of character strings containing information about your taxa of interest (check if non-target taxa with reads you want reported have other sep_ options below):
#' * column 1 - taxon name: like target_taxon in taxonomy_report()
#' * column 2 - taxon rank: like taxon_rank in taxonomy_report()
#' * column 3 - report name: what you want each group to be named (e.g., a common name, defaults to the taxon name in the first column if left blank)
#' * Example:
#'     taxa_mat = matrix(c(
#'       "Metazoa", "kingdom", "Animal",
#'       "Actinopteri", "class", "Fish",
#'       "Viridiplantae", "kingdom", "Plant"),
#'   byrow = T, ncol = 3)
#' @param min_reads_sample A number denoting the minimum quantity of OTU reads in a sample for it to be included in the generated GRAPHS (default = 8)
#' @param sep_human A boolean specifying whether human reads should be separated/tracked in their own category (default = T)
#' @param sep_bacteria A boolean specifying whether bacterial reads should be separated/tracked in their own category (default = F)
#' @param sep_phyta A boolean specifying whether eukaryotic algae/plant reads should be separated/tracked in their own category (default = F)
#' @param sep_eukaryota A boolean specifying whether non-target eukaryotic reads should be separated/tracked in their own category (default = F)
#' @param sep_no_domain A boolean specifying whether reads for which a domain could not be assigned should be separated/tracked in their own category (default = F)
#' @param report_graphs A boolean specifying whether or not to create graphs of the read and OTU counts of specified taxa across samples (default = T)
#' @param custom_palette A character string containing a color palette to be interpreted by ggplot2 for created graphs (default = F)
#' @param report_text A boolean specifying whether or not to create text reports of the read counts of specified taxa across samples (default = T)
#' @param return_list A boolean specifying whether or not to return the sample columns sorted into the groups as a list object is returned in R (default = F)
#' @param file_prefix A character string to be the prefix of the output file (default = NULL)
#' @param outdir A character string to specify the output directory of the output files (default = getwd())
#' @returns
#' * When report_graphs = T, graphs of overall and relative abundance of reads and OTUs per sample (grouped by target taxa specified in the taxa_mat) are returned.
#' * When report_text = T, a text file with the overall and relative abundance of reads per sample is created (grouped by target taxa)
#' * When return_list = T, a list containing two data frames is returned. The first is the sample read count columns with an additional "group" variable. The second object is a pivot_longer() version of this data frame grouped by the "group" variable.
#' @examples
#' sample_report(df = mifish_data, first_sample_col = 9, taxa_mat = matrix(c("Actinopteri", "class", "Fish"), nrow = 1))
#'
#' taxa_groups = sample_report(df = coi_data, first_sample_col = 16,
#'                             taxa_mat = matrix(c("Metazoa", "kingdom", "Animal",
#'                                                 "Arthropoda", "phylum", "Arthropod",
#'                                                 "Rhodophyta", "phylum", "Red algae"),
#'                                               byrow = T, ncol = 3),
#'                              min_reads_sample = 2,
#'                              sep_bacteria = T,
#'                              sep_phyta = T,
#'                              sep_eukaryota = T,
#'                              sep_no_domain = T,
#'                              report_graphs = T,
#'                              custom_palette = ggthemes::colorblind_pal()(8),
#'                              report_text = F,
#'                              return_list = T,
#'                              file_prefix = "coi_samples",
#'                              outdir = "~/Desktop/results")
#' @description Generates graphs, text reports, and data frames of samples sorted by pre-selected groups. This function will use only the 8 columns of taxonomic information (domain to species, in any location
#' across the data frame), as well as sample read counts starting from the first specified column of sample read counts (first_sample_col) through to the end of the data frame (so make sure your data frame only contains
#' columns of sample read counts after the specified first_sample_col).
#' * The purpose of the "sep_..." group of arguments is to more easily sort out common contaminants or non-target sequence data instead of having to
#' put these in your taxa_mat. If you are interested in knowing how many reads were assigned to human, bacteria, land plants + algae, non-specified and non-target eukaryotes, or sequences that could not be to assigned
#' a domain, but one (or all) of these groups were not among your targeted taxa, it is better to use one of the "sep_" arguments. These arguments put reads and OTUs from these categories at the top of each sample's data
#' report, keeping targeted taxa at the bottom where counts and proportions can be more easily found and read. If a subgroup of one of these taxa for which a "sep_" argument exists is a target taxon (e.g. one phylum of
#' algae or genus of bacteria), the sep_ argument as TRUE will exclude target sequences labeling the remaining non-targets as "Other ..." (e.g., algae or bacteria).
#' * If Primates (order) or Hominidae (family) are target taxa (in taxa_mat), the sep_human argument treats all sequences in the genus Homo as human. If neither of these two taxa are in the taxa_mat, all sequences
#' identified only to Primates or Hominidae are also considered human.
#' * Custom color palettes must match the number of groups to be plotted, otherwise the default palette will be used. If you get warned that an incorrect number of colors were provided, check how many were plotted using
#' the default and then run the function again with the correct number of provided colors.
#' @export
sample_report <- function(df,
                         first_sample_col,
                         taxa_mat,
                         min_reads_sample = 8,
                         sep_human = T,
                         sep_bacteria = F,
                         sep_phyta = F,
                         sep_eukaryota = F,
                         sep_no_domain = F,
                         custom_palette = F,
                         report_graphs = T,
                         report_text = T,
                         return_list = F,
                         file_prefix = NULL,
                         outdir = getwd()) {

  if(report_graphs==F & report_text==F & return_list==F) stop("At least one of report_graphs, report_text, or return_list must be set to TRUE.")

  ## Issue warning about looking for particular species if applicable
  if("species" %in% taxa_mat[,2]) cat("Species names are not always precise in BLAST taxonomy output,
                                      I recommend double checking species data using base::grep().\n")

  ## If you want to see phyta, bacteria, (other) eukaryota, or no domain reads, add rows for them in the taxa_mat
  if(sep_phyta==T) taxa_mat = rbind(taxa_mat, c("Phyta", "", ""))
  if(sep_eukaryota==T) taxa_mat = rbind(taxa_mat, c("Eukaryota", "domain", ""))
  if(sep_bacteria==T) taxa_mat = rbind(taxa_mat, c("Bacteria", "domain", ""))
  if(sep_no_domain==T) {
    df[which(df$domain=="dropped"),"domain"] = "No Domain"
    taxa_mat = rbind(taxa_mat, c("No Domain", "domain", "No Domain"))
  }

  ## Add row for human contamination (checking whether other primates are among taxa of interest)
  if(sep_human==T &
     "Primates" %in% taxa_mat[,1] |
     sep_human==T &
     "Hominidae" %in% taxa_mat[,1]) {
    taxa_mat = rbind(taxa_mat, c("Homo", "genus", "Human"))
  } else if(sep_human==T) {
    taxa_mat = rbind(taxa_mat, c("Primates", "order", "Human"))
  } else {
    cat("With sep_humans = FALSE, human reads could be counted in a taxon
        within which humans are nested (e.g. Chordata, Mammalia).\n")
  }

  ## If the plotting name was left blank for a column, fill it with the taxon name
  for(z in 1:nrow(taxa_mat)) {
    if(taxa_mat[z,3]=="") taxa_mat[z,3] <- taxa_mat[z,1]
  }

  ## Assign names hierarchically from lowest to highest taxonomic level
  df_new <- df %>%
    dplyr::mutate(group = dplyr::case_when(
      species %in% taxa_mat[,1] ~ as.character(species),
      genus %in% taxa_mat[,1] ~ as.character(genus),
      family %in% taxa_mat[,1] ~ as.character(family),
      order %in% taxa_mat[,1] ~ as.character(order),
      class %in% taxa_mat[,1] ~ as.character(class),
      phylum %in% taxa_mat[,1] ~ as.character(phylum),
      kingdom %in% taxa_mat[,1] ~ as.character(kingdom),
      domain %in% taxa_mat[,1] ~ as.character(domain),
      TRUE ~ "Other"
    ))

  ## If interested in eukaryotic algae/plants, assign name "Phyta" based on whether
  ## phyla or classes contains botanical suffixes
  if(sep_phyta==T) {
    df_new <- df_new %>%
      dplyr::mutate(group = dplyr::case_when(
        group=="Other" & domain == "Eukaryota" & base::grepl("phyta", phylum) |
          group=="Eukaryota" & base::grepl("phyta", phylum) ~ "Phyta",
        group=="Other" & domain == "Eukaryota" & base::grepl("phyceae", class) |
          group =="Eukaryota" & base::grepl("phyceae", class) ~ "Phyta",
        group != "Other" ~ group,
        TRUE ~ "Other"
      ))
  }

  ## Replace taxonomic names with the other names provided initially
  match_indices <- match(df_new$group, taxa_mat[,1])
  df_new$group <- ifelse(is.na(match_indices), df_new$group, taxa_mat[match_indices, 3])
  lookup = setNames(as.character(taxa_mat[,3]), taxa_mat[,1])

  ## If any taxa of interest are nested in others, clarify the nested one is excluded by naming
  ## the higher-order one with "Other ..." and "Non-human ..." when humans are nested in the
  ## taxon, with some special exceptions for dealing with phyta
  taxa_mat1 = taxa_mat
  if(sep_phyta==T) taxa_mat1 = taxa_mat[-which(taxa_mat[,1]=="Phyta"),]

  for(i in 1:nrow(taxa_mat1)) {
    taxon_i = taxa_mat1[i,1]
    entries_tax_i = which(df_new$group==lookup[taxon_i])

    ## Add "non-human" label to taxa within which humans are nested
    if(sep_human==T &
       taxon_i %in% c("Eukaryota", "Metazoa", "Chordata",
                      "Mammalia", "Primates", "Hominidae") &
       taxa_mat1[i,3] != "Human") {
      df_new[entries_tax_i, "group"] = paste0("Non-human ", lookup[taxon_i])
    }

    ## If the number of OTUs within a group don't match the number belonging to the titular taxon,
    ## there are other taxa nested within it and its name needs to be changed
    if(length(entries_tax_i)!=length(which(df_new[,taxa_mat1[i,2]]==taxon_i))) {

      ## If the "non-human" label hasn't already been added, it can just receive the "other" label
      if(lookup[taxon_i] %in% df_new$group) {
        df_new[which(df_new$group==lookup[taxon_i]), "group"] = paste0("Other ", lookup[taxon_i])
        ## If there are also non-human focal taxa within the current one already labeled non-human,
        ## label it as "other non-human taxon i"
      } else if(length(unique(df_new[which(df_new[,taxa_mat1[i,2]]==taxon_i),"group"]))>2) {
        df_new[which(df_new$group==paste0("Non-human ", lookup[taxon_i])), "group"] = paste0("Other non-human ", lookup[taxon_i])
        ## If the only focal taxon within one already labeled "non-human" is human
        ## (sep_human==T), it can stay labeled "non-human"
      } else {
        df_new[which(df_new$group==paste0("Non-human ", lookup[taxon_i])), "group"] = paste0("Non-human ", lookup[taxon_i])
      }
    }
  }
  ## Since phyta has its own option for being dealt with, isn't monophyletic, and will
  ## never have humans nested within, we have to check it in its own special case here
  if(sep_phyta==T & length(union(grep("phyta", df_new[which(df_new$domain=="Eukaryota"),"phylum"]),
                                        grep("phyceae", df_new[which(df_new$domain=="Eukaryota"),"class"])))!=length(which(df_new$group=="Phyta"))) {
    df_new[which(df_new$group=="Phyta"),"group"] = "Other Phyta"
  }

  ## Remove all information except group and reads per sample
  df_reads <- df_new[,first_sample_col:length(df_new)]

  ## Order the groups for display
  final_groups <- as.character(unique(df_reads$group))
  prefixes <- c("Other ", "Non-human ", "Other non-human ")
  original_groups <- gsub(paste(prefixes, collapse = "|"), "", final_groups)

  ## If a catch-all "other" group is necessary from not defining all sequences, put the "other"
  ## category after the focal ones but before possible phyta, bacteria, and human categories
  if("Other" %in% original_groups) {
    other_spot <- length(final_groups)
    if(sep_phyta==T) other_spot = other_spot-1
    if(sep_eukaryota==T) other_spot = other_spot-1
    if(sep_bacteria==T) other_spot = other_spot-1
    if(sep_no_domain==T) other_spot = other_spot-1
    if(sep_human==T) other_spot = other_spot-1
    new_order <- c(1:(other_spot-1), length(original_groups), other_spot:(length(original_groups)-1))
    ordered_groups <- rev(final_groups[order(match(original_groups, taxa_mat[,3]))][new_order])
  } else {
    ordered_groups <- rev(final_groups[order(match(original_groups, taxa_mat[,3]))])
  }
  ## Order them in the factor
  df_reads$group <- factor(df_reads$group, levels = ordered_groups)

  ## Reformat data to get one entry per group per site
  df_reads_long <- tidyr::pivot_longer(df_reads, cols = c(names(df_reads)[1:length(df_reads)-1]),
                                       names_to = "sample", values_to = "reads")
  df_reads_long <- df_reads_long[which(df_reads_long$reads>=min_reads_sample),]
  df_reads_grouped <- df_reads_long %>% dplyr::group_by(sample, group) %>% dplyr::summarize(.groups = "keep",
                                                                                            reads = sum(reads))
  ## If graphs are desired, set up color palettes and plot them
  if(report_graphs==T) {
    use_colorblind_pal = F
    use_custom_pal = F
    ## If no custom color palette is defined (custom_palette = F by default) and there are 8 or
    ## fewer focal taxa, plot them using the colorblind-friendly palette
    if(class(custom_palette)=="character" & length(final_groups)<9) {
      use_colorblind_pal = T
    } else if(class(custom_palette)=="character") {
      ## If no custom palette was defined and there are >8 focal taxa, warn the user that
      ## the plot can not be made with the default colorblind-friendly palette
      warning("More than 8 groups are to be plotted. The default colorblind-friendly palette could not be used.")
      ## If the user-defined custom palette does not have the same number of colors as focal groups,
      ## warn the user and use the default palette
    } else if(length(custom_palette)!=length(final_groups)) {
      if(length(final_groups)<9) use_colorblind_pal = T
      warning("The length of the palette provided does not match the number of groups,\
              or does not match the required format, plotting with the default color palette.")
    } else {
      ## If the user-defined color palette checks out, advise user to consider whether it is
      ## colorblind-friendly
      cat("Using a custom color palette for plotting. Please consider whether these colors are\
          distinguishable by people who are colorblind.")
      use_custom_pal = T
    }

    ## Make composite ggplots for reads and OTU richness
    total_reads <- ggplot2::ggplot(df_reads_grouped, ggplot2::aes(x = sample, y = reads, fill = group)) +
      ggplot2::geom_col(show.legend = F) +
      {if(use_colorblind_pal) ggthemes::scale_fill_colorblind()} +
      {if(use_custom_pal) ggplot2::scale_fill_manual(values = custom_palette)} +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

    relative_reads <- ggplot2::ggplot(df_reads_grouped, ggplot2::aes(x = sample, y = reads, fill = group)) +
      ggplot2::geom_col(position = "fill") +
      {if(use_colorblind_pal) ggthemes::scale_fill_colorblind()} +
      {if(use_custom_pal) ggplot2::scale_fill_manual(values = custom_palette)} +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

    reads = patchwork::wrap_plots(total_reads, relative_reads, ncol = 1)

    ggplot2::ggsave(paste0(outdir, "/", file_prefix, "_reads.pdf"), reads, device = "pdf",
                    height = 10, width = ifelse(length(df_reads)<50,
                                                10, 0.15*length(df_reads)), units = "in")

    total_otus <- ggplot2::ggplot(df_reads_long, ggplot2::aes(x = sample, fill = group)) +
      ggplot2::geom_bar(show.legend = F) +
      {if(use_colorblind_pal) ggthemes::scale_fill_colorblind()} +
      {if(use_custom_pal) ggplot2::scale_fill_manual(values = custom_palette)} +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

    relative_otus <- ggplot2::ggplot(df_reads_long, ggplot2::aes(x = sample, fill = group)) +
      ggplot2::geom_bar(position = "fill") +
      {if(use_colorblind_pal) ggthemes::scale_fill_colorblind()} +
      {if(use_custom_pal) ggplot2::scale_fill_manual(values = custom_palette)} +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

    otus = patchwork::wrap_plots(total_otus, relative_otus, ncol = 1)

    ggplot2::ggsave(paste0(outdir, "/", file_prefix, "_OTUs.pdf"), otus, device = "pdf",
                    height = 10, width = ifelse(length(df_reads)<50,
                                                10, 0.15*length(df_reads)), units = "in")
  }

  if(report_text==T) {
    ## Create text report of same information for more precision
    output <- capture.output({

      for(e in unique(df_reads_grouped$sample)) {
        sample_groups = dplyr::filter(df_reads_grouped, sample==e)
        total = sum(sample_groups$reads)
        ## Write text file summary
        cat(sep = "", "Out of ", total, " ", e, " reads:\n")

        for(f in nrow(sample_groups):1) {
          cat(sep = "", as.numeric(sample_groups[f,"reads"]),
              " (", round(as.numeric(sample_groups[f,"reads"])/total, digits = 4)*100,
              "%) were assigned to ", as.character(droplevels(sample_groups$group[f])), "\n")
        }
        cat("\n")
      }
    }
    )
    ## Save console output object to text file
    sink(paste0(outdir, "/", file_prefix, "_reads.txt"))
    cat(output, sep = "\n")
    sink()
  }
  ## Return list if requested
  if(return_list==T) return(list(reads = df_reads, grouped = df_reads_grouped))
}
