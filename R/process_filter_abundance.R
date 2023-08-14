#' @title Filtering feature who is low relative abundance or unclassified
#'
#' @description whether to filter the low relative abundance or unclassified
#' feature by the threshold. Here, we choose the following criterion:
#' 1. Feature more than Mean absolute or relative abundance across all samples;
#' 2. Feature more than Minimum absolute or relative abundance at least one sample.
#'
#' @references Thingholm, Louise B., et al. "Obese individuals with and without
#' type 2 diabetes show different gut microbial functional capacity and
#' composition." Cell host & microbe 26.2 (2019): 252-264.
#'
#' @author Created by Hua Zou (11/30/2021 Shenzhen China)
#'
#' @param object (Required). a [`phyloseq::phyloseq-class`] or
#' [`SummarizedExperiment::SummarizedExperiment-class`] object.
#' @param level (Optional). character. taxonomic level to summarize,
#' default the top level rank of the `ps`. taxonomic level(Kingdom, Phylum,
#' Class, Order, Family, Genus, Species, Strains; default: NULL).
#' @param cutoff_mean (Optional). numeric. Threshold for Mean
#' absolute (integer) or relative (float) abundance all samples (default, `0`).
#' @param cutoff_one (Optional). numeric. Threshold for Minimum
#' absolute (integer) or relative (float) abundance at least one sample (default, `0`).
#' @param unclass (Optional). logical. whether to filter the unclassified taxa (default `TRUE`).
#'
#' @import phyloseq
#' @importFrom purrr map
#' @importFrom dplyr bind_rows %>% group_split
#' @importFrom tibble as_tibble rownames_to_column column_to_rownames
#'
#' @return a [`phyloseq::phyloseq-class`] or
#' [`SummarizedExperiment::SummarizedExperiment-class`] object,
#' where each row represents a feature and each col represents the
#' feature abundance of each sample.
#'
#' @export
#'
#' @usage filter_abundance(
#'    object,
#'    level = c(NULL, "Kingdom", "Phylum", "Class",
#'            "Order", "Family", "Genus",
#'            "Species", "Strain", "unique"),
#'    cutoff_mean = c(100, 1e-04, 2),
#'    cutoff_one = c(1000, 1e-03, 3),
#'    unclass = TRUE)
#'
#' @examples
#'
#'\dontrun{
#' # phyloseq object
#'  data("Zeybel_2022_gut")
#'  Zeybel_2022_gut_counts <- phyloseq::transform_sample_counts(
#'  Zeybel_2022_gut, function(x) {round(x * 10^7)})
#'
#'  # absolute abundance
#'  ps <- filter_abundance(
#'    object = Zeybel_2022_gut_counts,
#'    level = NULL,
#'    cutoff_mean = 100,
#'    cutoff_one = 1000,
#'    unclass = FALSE)
#'
#'  # relative abundance
#'  ps <- filter_abundance(
#'    object = Zeybel_2022_gut,
#'    level = "Phylum",
#'    cutoff_mean = 1e-04,
#'    cutoff_one = 1e-03,
#'    unclass = TRUE)
#'
#'
#' # SummarizedExperiment object
#' data("Zeybel_2022_protein")
#' filter_abundance(
#'   object = Zeybel_2022_protein,
#'   cutoff_mean = 5,
#'   cutoff_one = 8)
#' }
#'
filter_abundance <- function(
    object,
    level = NULL,
    cutoff_mean = 0,
    cutoff_one = 0,
    unclass = TRUE) {

  # data("Zeybel_2022_gut")
  # object = Zeybel_2022_gut
  # level = "Phylum"
  # cutoff_mean = 1e-04
  # cutoff_one = 1e-03
  # unclass = TRUE

  # data("Zeybel_2022_protein")
  # object = Zeybel_2022_protein
  # level = NULL
  # cutoff_mean = 2
  # cutoff_one = 3
  # unclass = TRUE

  if (all(is.null(cutoff_mean), is.null(cutoff_one))) {
    return(object)
  }

  if (all(!is.null(object), inherits(object, "phyloseq"))) {
    # profile: row->features; col->samples
    stopifnot(inherits(object, "phyloseq"))

    res <- LowAbundance_taxa(
      ps = object,
      taxa_level = level,
      cutoff_mean = cutoff_mean,
      cutoff_one = cutoff_one,
      unclass = TRUE)

    res_otu_tab <- as(phyloseq::otu_table(res), "matrix")

    if (unclass) {
      res_otu_tab_filter <- res_otu_tab %>% data.frame() %>%
        tibble::rownames_to_column("TaxaID") %>%
        dplyr::filter(!TaxaID %in% grep("Unclassified|LowAbundance", TaxaID, value = T)) %>%
        tibble::column_to_rownames("TaxaID") %>%
        as.matrix()
    } else {
      res_otu_tab_filter <- res_otu_tab %>% data.frame() %>%
        tibble::rownames_to_column("TaxaID") %>%
        dplyr::filter(!TaxaID %in% grep("LowAbundance", TaxaID, value = T)) %>%
        tibble::column_to_rownames("TaxaID") %>%
        as.matrix()
    }

    if (!is.null(res@sam_data)) {
      colnames(res_otu_tab_filter) <- rownames(sample_data(res) %>%
                                                    data.frame())
    }

    phyloseq::otu_table(res) <- phyloseq::otu_table(res_otu_tab_filter,
                                                       taxa_are_rows = TRUE)

  } else if (all(!is.null(object), inherits(object, "SummarizedExperiment"))) {
    # profile: row->features; col->samples
    res <- LowAbundance_feature(
      dat = object,
      cutoff_mean = cutoff_mean,
      cutoff_one = cutoff_one)
  }

  return(res)
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

  # OTU table: row->features; col->samples
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

#' low abundance feature whose abundance less than threshold
#' @noRd
LowAbundance_feature <- function(
    dat,
    cutoff_mean,
    cutoff_one) {

  # dat = object
  # cutoff_mean = cutoff_mean
  # cutoff_one = cutoff_one

  # profile: row->features; col->samples
  prof_tab <- as(SummarizedExperiment::assay(dat), "matrix")

  ############################################
  # 1. Mean relative abundance across all samples
  feture_mean <- apply(prof_tab, 1, function(x)sum(x, na.rm = TRUE)) %>%
    data.frame() %>%
    stats::setNames("SumAbundance") %>%
    tibble::rownames_to_column("TaxaID") %>%
    dplyr::mutate(MeanAbundance = SumAbundance/ncol(prof_tab))

  # 2. Minimum relative abundance at least one sample
  feture_one <- apply(prof_tab, 1, function(x) {
    any(x[!is.na(x)] >= cutoff_one)
  }) %>% data.frame() %>%
    stats::setNames("LowAbundanceOrNot") %>%
    tibble::rownames_to_column("TaxaID")

  # 3. combine the results
  feature_res <- feture_mean %>%
    dplyr::inner_join(feture_one, by = "TaxaID")
  ############################################

  feature_define <- feature_res %>%
    dplyr::mutate(TaxaID2 = ifelse(MeanAbundance <= cutoff_mean,
                                   "Others_LowAbundance",
                                   ifelse(LowAbundanceOrNot, TaxaID,
                                          "Others_LowAbundance")))
  feature_merge <- feature_define %>%
    dplyr::select(TaxaID, TaxaID2) %>%
    dplyr::inner_join(prof_tab %>%
                        data.frame() %>%
                        tibble::rownames_to_column("TaxaID"),
                      by = "TaxaID") %>%
    dplyr::select(-TaxaID)

  # feature_aggregate <- feature_merge %>%
  #   dplyr::group_by(TaxaID2) %>%
  #   dplyr::summarise(across(everything(), ~ sum(., is.na(.), 0))) %>%
  #   tibble::column_to_rownames("TaxaID2")

  feature_final <- feature_merge %>%
    dplyr::filter(!TaxaID2 %in% c("Others_LowAbundance")) %>%
    tibble::column_to_rownames("TaxaID2")


  # if (is.null(dat@elementMetadata)) {
  #   res <- import_SE(object = feature_aggregate,
  #                    coldata = dat@colData,
  #                    metadata = dat@metadata)
  # } else {
  #   rowdf <- dat@elementMetadata %>% as.data.frame()
  #   overlap_featureID <- intersect(rownames(feature_aggregate),
  #                                rowdf[, 1])
  #   lowdf <- data.frame(matrix(NA, nrow = 1, ncol = ncol(rowdf))) %>%
  #     stats::setNames(colnames(rowdf))
  #   lowdf[1, 1] <- setdiff(rownames(feature_aggregate),
  #                          rowdf[, 1])
  #
  #   rowdf_cln <- rowdf[rowdf[, 1] %in% overlap_featureID, , F] %>%
  #     rbind(lowdf)
  #   rownames(rowdf_cln) <- rowdf_cln[, 1]
  #
  #   # the latter ordered by the former
  #   data_assay <- feature_aggregate[pmatch(rownames(rowdf_cln),
  #                                          rownames(feature_aggregate)), ]
  #
  #   res <- import_SE(object = data_assay,
  #                    rowdata = rowdf_cln,
  #                    coldata = dat@colData,
  #                    metadata = dat@metadata)
  # }

  res <- dat
  if (nrow(prof_tab) != nrow(feature_final)) {
    res <- base::subset(res, rownames(prof_tab) %in% rownames(feature_final))
  } else {
    SummarizedExperiment::assay(res) <- feature_final
  }

  return(res)
}
