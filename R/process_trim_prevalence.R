#' @title Trimming samples or features whose prevalence is less than threshold
#'
#' @description
#' trim samples or features in `profile` by Prevalence,
#' which means the samples or features will be discarded if they could not pass the cutoff.
#'
#' @author Created by Hua Zou (11/30/2021 Shenzhen China)
#'
#' @param object (Required). a [`matrix`], [`otu_table-class`],
#' [`phyloseq::phyloseq-class`] or [`SummarizedExperiment-class`]
#' @param level (Optional). character. taxonomic level to summarize,
#' default the top level rank of the `ps`. taxonomic level(Kingdom, Phylum,
#' Class, Order, Family, Genus, Species, Strains; default: NULL).
#' @param cutoff (Optional). Numeric. the Prevalence threshold (default: 0.1).
#' @param group (Optional). character. filtering features or samples by group
#' (default: NULL).
#' @param trim (Optional). Character. trimming to apply, the options include:
#' * "none", return the original data without any actions.
#' * "both", prevalence of features and samples more than cutoff.
#' * "feature", prevalence of features more than cutoff.
#' * "feature_group", prevalence of features more than cutoff by groups.
#' * "sample", prevalence of samples more than cutoff.
#' @param at_least_one (Optional). Logical. prevalence of at least one group
#' meets cutoff (FALSE means all groups meet cutoff, default: FALSE).
#'
#' @return
#'  A trimed `object` whose prevalence of features or samples more than cutoff.
#'
#' @export
#'
#' @import phyloseq
#' @importFrom stats setNames
#' @importFrom dplyr filter %>%
#' @importFrom SummarizedExperiment colData assay
#'
#' @usage trim_prevalence(
#'     object,
#'     level = c(NULL, "Kingdom", "Phylum", "Class",
#'            "Order", "Family", "Genus",
#'            "Species", "Strain", "unique"),
#'     cutoff = 0.1,
#'     group = NULL,
#'     trim = c("none",
#'       "both",
#'       "feature", "feature_group",
#'       "sample"),
#'     at_least_one = FALSE)
#'
#' @examples
#' \donttest{
#' # phyloseq object
#' data("Zeybel_2022_gut")
#' trim_prevalence(
#'   Zeybel_2022_gut,
#'   level = "Phylum",
#'   cutoff = 0.1,
#'   trim = "feature")
#'
#' # SummarizedExperiment object
#' data("Zeybel_2022_protein")
#' trim_prevalence(
#'   Zeybel_2022_protein,
#'   level = NULL,
#'   cutoff = 0.1,
#'   trim = "feature")
#' }
#'
trim_prevalence <- function(
    object,
    level = NULL,
    cutoff = 0.1,
    group = NULL,
    trim = c("none",
             "both",
             "feature", "feature_group",
             "sample"),
    at_least_one = FALSE) {

  # data(Zeybel_2022_gut)
  # object = Zeybel_2022_gut
  # trim = "sample_group"
  # level = "Genus"
  # cutoff = 0.1
  # group = "LiverFatClass"
  # at_least_one = FALSE

  # data(Zeybel_2022_protein)
  # object = Zeybel_2022_protein
  # trim = "feature"
  # level = NULL
  # cutoff = 0.95
  # group = NULL
  # at_least_one = FALSE

  # trim type
  trim <- match.arg(trim, c("none",
                            "both",
                            "feature", "feature_group",
                            "sample"))
  if (length(grep("group", trim)) > 0) {
    if (is.null(group)) {
      stop("please provide `group` for filtering samples or features by group")
    } else if (any(inherits(object, "phyloseq"),
                   inherits(object, "SummarizedExperiment"))) {
      if (inherits(object, "phyloseq")) {
        metadata <- object@sam_data %>%
          data.frame()
      } else if (inherits(object, "SummarizedExperiment")) {
        metadata <- object@colData %>%
          data.frame()
      }
    } else {
      stop("please using phyloseq or SummarizedExperiment object with
           `group` in metadata for filtering samples or features by group")
    }
  }

  # profile: row->features; col->samples
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
  } else if (inherits(object, "SummarizedExperiment")) {
    # profile: row->features; col->samples
    prf <- SummarizedExperiment::assay(object) %>%
      as.data.frame()
  } else {
    prf <- object
  }

  # trimming
  if (trim == "both") {
    tmp1 <- trim_FeatureOrSample(prf, 1, cutoff)
    tmp2 <- trim_FeatureOrSample(prf, 2, cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- rownames(tmp2)
  } else if (trim == "feature") {
    tmp1 <- trim_FeatureOrSample(prf, 1, cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- colnames(prf)
  } else if(trim == "feature_group") {
    tmp1 <- trim_FeatureByGroup(prf, cutoff,
                                group, metadata, at_least_one)
    remain_features <- rownames(tmp1)
    remain_samples <- colnames(prf)
  } else if(trim == "sample") {
    tmp2 <- trim_FeatureOrSample(prf, 2, cutoff)
    remain_features <- rownames(prf)
    remain_samples <- rownames(tmp2)
  } else if(trim == "none") {
    return(object)
  }

  # remain profile
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

  # subset dataset
  if (inherits(object, "phyloseq")) {
    phyloseq::otu_table(ps) <- phyloseq::otu_table(prf_remain,
                                                   taxa_are_rows = taxa_are_rows(ps))
    object <- ps
  } else if (inherits(object, "SummarizedExperiment")) {

    object <- base::subset(object, rownames(object) %in% rownames(prf_remain))

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
    data.frame() %>% stats::setNames("Occ") #%>%
    #tibble::rownames_to_column("type")
  if (nRow == 1) {
    rownames(df_occ) <- rownames(x)
  } else {
    rownames(df_occ) <- colnames(x)
  }

  df_KEEP <- apply(df_occ >= threshold, 1, all) %>%
    data.frame() %>%
    stats::setNames("Status") %>%
    dplyr::filter(Status)

  return(df_KEEP)
}


# the data is trimmed by threshold and group
#' @keywords internal
trim_FeatureByGroup <- function(x, threshold,
                              group, metadata, at_least_one) {

  # x = prf
  # threshold = cutoff
  # group = group
  # metadata = metadata
  # at_least_one = at_least_one

  # groups for filtering
  colnames(metadata)[which(colnames(metadata) == group)] <- "TempGroupName"
  group_names <- unique(metadata$TempGroupName)
  group_list <- lapply(group_names, function(x){
    index <- which(metadata$TempGroupName == x)
    group_samples <- rownames(metadata)[index]
    return(group_samples)
  })

  # filter by feature
  df_filter <- sapply(1:length(group_list), function(j) {
      df_sub <- x[, colnames(x) %in% group_list[[j]]]
      df_sub_matrix <- as.matrix(df_sub)
      df_occ <- apply(df_sub_matrix, 1, function(x) {
        length(x[c(which(!is.na(x) & x!=0))]) / length(x)
      })
      return(df_occ)
    }) %>%
      data.frame() %>%
      stats::setNames(group_names) #%>%
      #tibble::rownames_to_column("type")

  #rownames(df_filter) <- rownames(x)

  if (at_least_one) {
    df_KEEP <- apply(df_filter, 1, function(x) {any(x >= threshold)}) %>%
      data.frame() %>%
      stats::setNames("Status") %>%
      dplyr::filter(Status)
  } else {
    df_KEEP <- apply(df_filter, 1, function(x) {all(x >= threshold)}) %>%
      data.frame() %>%
      stats::setNames("Status") %>%
      dplyr::filter(Status)
  }

  return(df_KEEP)
}
