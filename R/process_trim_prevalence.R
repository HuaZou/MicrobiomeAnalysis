#' @title Trimming samples or taxa whose prevalence is less than threshold
#'
#' @description
#' trim samples or taxa in `otu_table` by Prevalence,
#' which means the samples or taxa will be discarded if they could not pass the cutoff.
#'
#' @author Created by Hua Zou (11/30/2021 Shenzhen China)
#'
#' @param object (Required). a [`matrix`] or [`otu_table-class`] or [`phyloseq::phyloseq-class`].
#' @param level (Optional). character. taxonomic level to summarize,
#' default the top level rank of the `ps`. taxonomic level(Kingdom, Phylum,
#' Class, Order, Family, Genus, Species, Strains; default: NULL).
#' @param cutoff (Optional). Numeric. the Prevalence threshold (default: 0.1).
#' @param trim (Optional). Character. trimming to apply, the options include:
#' * "none", return the original data without any actions.
#' * "both", prevalence of taxa and samples more than cutoff.
#' * "feature", prevalence of taxa more than cutoff.
#' * "sample", prevalence of samples more than cutoff.
#'
#' @return
#'  A trimed `object` whose prevalence of taxa or samples more than cutoff.
#'
#' @export
#'
#' @import phyloseq
#' @importFrom stats setNames
#' @importFrom dplyr filter %>%
#'
#' @usage trim_prevalence(
#'     object,
#'     level = NULL,
#'     cutoff = 0.1,
#'     trim = c("none", "both", "feature", "sample"))
#'
#' @examples
#' \donttest{
#'  data("enterotypes_arumugam")
#'  trim_prevalence(object = enterotypes_arumugam,
#'    cutoff = 0.1, trim = "feature")
#' }
#'
trim_prevalence <- function(
    object,
    level = NULL,
    cutoff = 0.1,
    trim = c("none", "both", "feature", "sample")){

  # data("enterotypes_arumugam")
  # object = enterotypes_arumugam
  # level = "Phylum"
  # cutoff = 0.1
  # trim = "both"

  trim <- match.arg(trim, c("none", "both", "feature", "sample"))
  if (inherits(object, "phyloseq")) {

    # taxa level
    if (!is.null(level)) {
      ps <- aggregate_taxa(x = object, level = level)
    } else {
      ps <- object
    }

    prf <- as(phyloseq::otu_table(ps), "matrix")
  } else if (inherits(object, "environment")) {
    prf <- as(object$.Data, "matrix")
  } else {
    prf <- object
  }

  if (trim == "feature") {
    tmp1 <- trim_FeatureOrSample(prf, 1, cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- colnames(prf)
  } else if (trim == "sample") {
    tmp2 <- trim_FeatureOrSample(prf, 2, cutoff)
    remain_features <- rownames(prf)
    remain_samples <- rownames(tmp2)
  } else if(trim == "both") {
    tmp1 <- trim_FeatureOrSample(prf, 1, cutoff)
    tmp2 <- trim_FeatureOrSample(prf, 2, cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- rownames(tmp2)
  } else if(trim == "none") {
    return(object)
  }

  if (length(remain_features) > 1 & length(remain_samples) > 1) {
    prf_remain <- prf[remain_features, remain_samples]
  } else if (length(remain_features) > 1 & length(remain_samples) == 1) {
    prf_remain <- prf[remain_features, remain_samples] %>%
      data.frame()
    colnames(prf_remain) <- remain_samples
  } else if (length(remain_features) == 1 & length(remain_samples) > 1) {
    prf_remain <- prf[remain_features, remain_samples] %>%
      data.frame()
    rownames(prf_remain) <- remain_features
  } else if (length(remain_features) == 1 & length(remain_samples) == 1) {
    prf_remain <- prf[remain_features, remain_samples] %>%
      data.frame()

    colnames(prf_remain) <- remain_samples
    rownames(prf_remain) <- remain_features
  } else {
    stop("No sample or taxa were remained, please reinput your cutoff")
  }

  if (inherits(object, "phyloseq")) {
    phyloseq::otu_table(ps) <- phyloseq::otu_table(prf_remain,
                                                   taxa_are_rows = taxa_are_rows(ps))
    object <- ps
  } else if (inherits(object, "environment")) {
    object <- phyloseq::otu_table(prf_remain, taxa_are_rows = taxa_are_rows(object))
  } else {
    object <- prf_remain
  }

  return(object)
}

# the data is trimmed by threshold
#' @keywords internal
trim_FeatureOrSample <- function(x, nRow, threshold) {

  df_occ <- apply(x, nRow, function(x) {
    length(x[c(which(!is.na(x) & x!=0))]) / length(x)
  }) %>%
    data.frame() %>% stats::setNames("Occ") %>%
    tibble::rownames_to_column("type")
  if(nRow == 1){
    rownames(df_occ) <- rownames(x)
  }else{
    rownames(df_occ) <- colnames(x)
  }
  df_KEEP <- apply(df_occ > threshold, 1, all) %>%
    data.frame() %>% stats::setNames("Status") %>%
    dplyr::filter(Status)

  return(df_KEEP)
}
