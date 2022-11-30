#' @title Filtering taxa who is low relative abundance or unclassified
#'
#' @description whether to filter the low relative abundance or unclassified
#' taxa by the threshold. Here, we choose the following criterion:
#' 1. Taxa more than Mean relative abundance across all samples: 0.0001 (1e-04);
#' 2. Taxa more than Minimum relative abundance at least one sample: 0.001 (1e-03).
#'
#' @references Thingholm, Louise B., et al. "Obese individuals with and without
#' type 2 diabetes show different gut microbial functional capacity and
#' composition." Cell host & microbe 26.2 (2019): 252-264.
#'
#' @param object (Required). a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param level (Optional). character. taxonomic level to summarize,
#' default the top level rank of the `ps`. taxonomic level(Kingdom, Phylum,
#' Class, Order, Family, Genus, Species, Strains; default: NULL).
#' @param cutoff_mean (Optional). numeric. Threshold for Mean
#' relative abundance all samples (default, `1e-04`).
#' @param cutoff_one (Optional). numeric. Threshold for Minimum
#' relative abundance at least one sample (default, `1e-03`).
#' @param unclass (Optional). logical. whether to filter the unclassified taxa (default `TRUE`).
#'
#' @return a [`phyloseq::phyloseq-class`] object, where each row represents a
#'   taxa, and each col represents the taxa abundance of each sample.
#'
#' @import phyloseq
#' @importFrom purrr map
#' @importFrom dplyr bind_rows %>% group_split
#' @importFrom tibble as_tibble rownames_to_column column_to_rownames
#'
#' @export
#'
#' @usage filter_abundance(
#'     object,
#'     level = NULL,
#'     cutoff_mean = c(100, 1e-04),
#'     cutoff_one = c(1000, 1e-03),
#'     unclass = TRUE)
#'
#' @examples
#'
#'\dontrun{
#'  data("enterotypes_arumugam")
#'
#'  # absolute abundance
#'  ps <- filter_abundance(
#'    object = enterotypes_arumugam,
#'    level = NULL,
#'    cutoff_mean = 100,
#'    cutoff_one = 1000,
#'    unclass = FALSE)
#'
#'  # relative abundance
#'  enterotypes_arumugam_rb <- normalize(enterotypes_arumugam, method = "TSS")
#'  ps <- filter_abundance(
#'    object = enterotypes_arumugam_rb,
#'    level = "Phylum",
#'    cutoff_mean = 1e-04,
#'    cutoff_one = 1e-03,
#'    unclass = TRUE)
#'
#' }
#'
filter_abundance <- function(
    object,
    level = NULL,
    cutoff_mean,
    cutoff_one,
    unclass = TRUE) {

  # data("enterotypes_arumugam")
  # object = enterotypes_arumugam
  # level = "Phylum"
  # cutoff_mean = 100
  # cutoff_one = 1000
  # unclass = TRUE

  ps_LRA <- LowAbundance_taxa(
      ps = object,
      taxa_level = level,
      cutoff_mean = cutoff_mean,
      cutoff_one = cutoff_one,
      unclass = TRUE)

  ps_LRA_otu_tab <- as(phyloseq::otu_table(ps_LRA), "matrix")

  if (unclass) {
    ps_LRA_otu_tab_filter <- ps_LRA_otu_tab %>% data.frame() %>%
      tibble::rownames_to_column("TaxaID") %>%
      dplyr::filter(!TaxaID %in% grep("Unclassified|LowAbundance", TaxaID, value = T)) %>%
      tibble::column_to_rownames("TaxaID") %>%
      as.matrix()
  } else {
    ps_LRA_otu_tab_filter <- ps_LRA_otu_tab %>% data.frame() %>%
      tibble::rownames_to_column("TaxaID") %>%
      dplyr::filter(!TaxaID %in% grep("LowAbundance", TaxaID, value = T)) %>%
      tibble::column_to_rownames("TaxaID") %>%
      as.matrix()
  }

  if (!is.null(ps_LRA@sam_data)) {
    colnames(ps_LRA_otu_tab_filter) <- rownames(sample_data(ps_LRA) %>%
                                                  data.frame())
  }

  phyloseq::otu_table(ps_LRA) <- phyloseq::otu_table(ps_LRA_otu_tab_filter,
                                                     taxa_are_rows = TRUE)

  return(ps_LRA)
}

#' low abundance taxa whose abundance less than threshold
#' @noRd
LowAbundance_taxa <- function(
    ps,
    taxa_level = NULL,
    cutoff_mean,
    cutoff_one,
    unclass) {

  # ps = object
  # taxa_level = level
  # cutoff_mean = cutoff_mean
  # cutoff_one = cutoff_one
  # unclass = TRUE

  if (!is.null(taxa_level)) {
    ps_taxa <- aggregate_taxa(x = ps, level = taxa_level)
  } else {
    ps_taxa <- ps
  }
  otu_tab <- as(phyloseq::otu_table(ps_taxa), "matrix")

  res <- .Calculate_LowAbundance_taxa(
      ps_taxa,
      threshold_mean = cutoff_mean,
      threshold_one = cutoff_one,
      unclassify = unclass)

  otu_summarized <- phyloseq::otu_table(as.matrix(res$otu),
                                        taxa_are_rows = TRUE)

  if (!is.null(res$tax)) {
    tax_summarized <- phyloseq::tax_table(as.matrix(res$tax))
    rownames(tax_summarized) <- rownames(res$tax)
    colnames(tax_summarized) <- colnames(res$tax)
    rownames(tax_summarized@.Data) <- rownames(otu_summarized)
  } else {
    tax_summarized <- NULL
  }

  if (!is.null(ps@sam_data)) {

    colnames(otu_summarized) <- rownames(sample_data(ps) %>%
                                           data.frame())

    phyloseq_object <- phyloseq(otu_summarized,
                                tax_summarized,
                                sample_data(ps))
  } else {
    phyloseq_object <- phyloseq(otu_summarized,
                                tax_summarized)
  }

  return(phyloseq_object)
}


#' low abundance taxa whose abundance less than threshold
#' @noRd
.Calculate_LowAbundance_taxa <- function(
    ps,
    threshold_mean,
    threshold_one,
    unclassify = FALSE) {

  # ps = ps
  # threshold_mean = cutoff_mean
  # threshold_one = cutoff_one
  # unclassify = unclass

  # OTU table
  otus <- phyloseq::otu_table(ps)
  otus_extend <- slot(otus, ".Data") %>%
    data.frame()
  rownames(otus_extend) <- rownames(otus@.Data)

  ############################################
  # 1. Mean relative abundance across all samples
  taxa_mean_res <- phyloseq::taxa_sums(ps) %>%
    data.frame() %>%
    stats::setNames("SumAbundance") %>%
    tibble::rownames_to_column("TaxaID") %>%
    dplyr::mutate(MeanAbundance = SumAbundance/ncol(otus_extend))

  # 2. Minimum relative abundance at least one sample
  taxa_one_res <- apply(otus_extend, 1, function(x) {
    any(x >= threshold_one)
  }) %>% data.frame() %>%
    stats::setNames("LowAbundanceOrNot") %>%
    tibble::rownames_to_column("TaxaID")

  # 3. combine the results
  taxa_res <- taxa_mean_res %>%
    dplyr::inner_join(taxa_one_res, by = "TaxaID")
  ############################################

  # define the low relative abundance taxa or unclassified taxa
  if (unclassify) {
    taxa_table_define <- taxa_res %>%
      dplyr::mutate(TaxaID2 = ifelse(MeanAbundance <= threshold_mean,
                                     "Others_LowAbundance",
                                     ifelse(LowAbundanceOrNot, TaxaID,
                                            "Others_LowAbundance")))


    taxa_table_define$TaxaID2[grep("unclassified|noname|Other|Unclassified",
                                   taxa_table_define$TaxaID2)] <- "Others_Unclassified"
  } else {
    taxa_table_define <- taxa_res %>%
      dplyr::mutate(TaxaID2 = ifelse(MeanAbundance <= threshold_mean,
                                     "Others_LowAbundance",
                                     ifelse(LowAbundanceOrNot, TaxaID,
                                            "Others_LowAbundance")))
  }

  taxa_table_define_merge <- taxa_table_define %>%
    dplyr::select(TaxaID, TaxaID2) %>%
    dplyr::inner_join(otus_extend %>% tibble::rownames_to_column("TaxaID"),
                      by = "TaxaID") %>%
    dplyr::select(-TaxaID)

  taxa_table_define_merge_aggregate <-  taxa_table_define_merge %>%
    dplyr::group_by(TaxaID2) %>%
    dplyr::summarise(across(everything(), ~ sum(., is.na(.), 0))) %>%
    tibble::column_to_rownames("TaxaID2")

  rownames(otus_extend) <- rownames(otus@.Data)

  if (!is.null(ps@tax_table)) {
    taxas <- tax_table(ps)@.Data %>% data.frame()
    consensus <- taxas[match(rownames(taxa_table_define_merge_aggregate),
                             rownames(taxas)), , F]
    if (unclassify) {
      consensus[nrow(consensus)-1, ] <- rep("Others_LowAbundance", ncol(taxas))
      rownames(consensus)[nrow(consensus)-1] <- "Others_LowAbundance"

      consensus[nrow(consensus), ] <- rep("Others_Unclassified", ncol(taxas))
      rownames(consensus)[nrow(consensus)] <- "Others_Unclassified"
    } else {
      consensus[nrow(consensus), ] <- rep("Others_LowAbundance", ncol(taxas))
      rownames(consensus)[nrow(consensus)] <- "Others_LowAbundance"
    }
    tax_summarized <- consensus
  } else {
    tax_summarized <- NULL
  }

  otu_summarized <- taxa_table_define_merge_aggregate

  res <- list(otu = otu_summarized,
              tax = tax_summarized)

  return(res)
}
