# https://github.com/pcastellanoescuder/POMA/blob/master/R/PomaImpute.R
# https://github.com/WandeRum/MVI-evaluation/blob/master/Impute_wrapper.R

#' @title Imputation methods on missing value
#'
#' @description
#' This function offers different methods to impute missing values in data.
#'
#' @references Armitage, E. G., Godzien, J., Alonso‐Herranz, V., López‐Gonzálvez,
#' Á., & Barbas, C. (2015). Missing value imputation strategies for metabolomics
#' data. Electrophoresis, 36(24), 3050-3060.
#'
#' @author Created by Pol Castellano-Escuder;
#' Modified by Hua Zou (12/02/2022 Shenzhen China)
#'
#' @param object (Required). a [`phyloseq::phyloseq-class`] or
#' [`Biobase::ExpressionSet`] object.
#' @param level (Optional). character. Summarization
#' level (from \code{rank_names(pseq)}, default: NULL).
#' @param group_name (Required). character. group for determining missing values.
#' @param ZerosAsNA (Optional). logical. zeros in the data are missing
#' values (default: FALSE).
#' @param RemoveNA (Optional). logical. those features with more than selected
#' cutoff missing values in each group have to be removed (default: TRUE).
#' @param cutoff (Optional). numeric. percentage of missing values allowed in
#' each group. If one of the groups have less missing values than selected
#' cutoff value, these feature will not be removed.
#' @param method (Optional). character. Imputation method. Options are:
#'  * "none": all missing values will be replaced by zero.
#'  * "LOD": specific Limit Of Detection which provides by user.
#'  * "half_min": half minimal values across samples except zero.
#'  * "median": median values across samples except zero.
#'  * "mean": mean values across samples except zero.
#'  * "min": minimal values across samples except zero.
#'  * "knn": k-nearest neighbors samples.
#'  * "rf": nonparametric missing value imputation using Random Forest.
#'  * "global_mean": a normal distribution with a mean that is down-shifted from
#'  the sample mean and a standard deviation that is a
#'  fraction of the standard deviation of the sample distribution.
#'  * "svd": missing values imputation based Singular value decomposition.
#'  * "QRILC": missing values imputation based quantile regression.
#'  (default: "none").
#' @param LOD (Optional). Numeric. limit of detection (default: NULL).
#'
#' @usage impute_abundance(
#'    object,
#'    level = c(NULL, "Kingdom", "Phylum", "Class",
#'            "Order", "Family", "Genus",
#'            "Species", "Strain", "unique"),
#'    group_name,
#'    ZerosAsNA = FALSE,
#'    RemoveNA = TRUE,
#'    cutoff = 20,
#'    method = c("none", "LOD", "half_min", "median",
#'        "mean", "min", "knn", "rf",
#'        "global_mean", "svd", "QRILC"),
#'    LOD = NULL
#'    )
#'
#' @export
#'
#' @return A [`phyloseq::phyloseq-class`] or
#' [`Biobase::ExpressionSet`] object with cleaned data.
#'
#' @importFrom dplyr %>%
#'
#' @examples
#'
#' \dontrun{
#' data("enterotypes_arumugam")
#' impute_abundance(
#'   enterotypes_arumugam,
#'   level = "Phylum",
#'   group_name = "Enterotype",
#'   ZerosAsNA = TRUE,
#'   RemoveNA = TRUE,
#'   cutoff = 20,
#'   method = "knn")
#' }
#'
impute_abundance <- function(
    object,
    level = NULL,
    group_name,
    ZerosAsNA = FALSE,
    RemoveNA = TRUE,
    cutoff = 20,
    method = c("none", "LOD", "half_min", "median",
               "mean", "min", "knn", "rf",
               "global_mean", "svd", "QRILC"),
    LOD = NULL) {

  # data("enterotypes_arumugam")
  # object = enterotypes_arumugam
  # level = "Phylum"
  # group_name = "Enterotype"
  # ZerosAsNA = TRUE
  # RemoveNA = TRUE
  # cutoff = 20
  # method = "QRILC"

  if (base::missing(object)) {
    stop("object argument is empty!")
  }

  if (all(!methods::is(object, "phyloseq"),
          !methods::is(object, "ExpressionSet"))) {
    stop("object is not either a phyloseq or ExpressionSet object.")
  }

  method <- match.arg(
    method, c("none", "LOD", "half_min", "median",
              "mean", "min", "knn", "rf",
              "global_mean", "svd", "QRILC")
  )

  # if (!(method %in% c("none", "LOD", "half_min", "median",
  #                     "mean", "min", "knn", "rf",
  #                     "global_mean", "svd", "QRILC"))) {
  #   stop("Incorrect value for method argument!")
  # }

  if (base::missing(method)) {
    message("method argument is empty! KNN will be used")
  }

  # phyloseq object
  if (all(!is.null(object), inherits(object, "phyloseq"))) {

    ps <- check_sample_names(object = object)

    # taxa level
    if (!is.null(level)) {
      ps <- aggregate_taxa(x = ps, level = level)
    } else {
      ps <- ps
    }
    ## sample table & profile table
    sam_tab <- phyloseq::sample_data(ps) %>%
      data.frame() %>%
      tibble::rownames_to_column("TempRowNames")

    if (phyloseq::taxa_are_rows(ps)) {
      prf_tab <- phyloseq::otu_table(phyloseq::t(ps)) %>%
        data.frame()
    } else {
      prf_tab <- phyloseq::otu_table(ps) %>% data.frame()
    }

  } else if (all(!is.null(object), inherits(object, "ExpressionSet"))) {
    # sample table & profile table
    sam_tab <- Biobase::pData(object) %>%
      data.frame() %>%
      tibble::rownames_to_column("TempRowNames")
    prf_tab <- Biobase::exprs(object) %>%
      data.frame()
  }

  group_index <- which(colnames(sam_tab) == group_name)
  samples_groups <- sam_tab[, group_index]
  to_imp_data <- prf_tab %>% as.matrix()

  if (ZerosAsNA) {
    to_imp_data[to_imp_data == 0] <- NA
    to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))
    colnames(to_imp_data)[2:ncol(to_imp_data)] <- colnames(prf_tab)

  } else {
    to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))
    colnames(to_imp_data)[2:ncol(to_imp_data)] <- colnames(prf_tab)
  }

  percent_na <- sum(is.na(to_imp_data))
  if (percent_na == 0) {
    message("No missing values detected in your data")
    if (method != "none") {
      method <- "none"
    }
  }

  if (isTRUE(RemoveNA)) {
    count_NA <- stats::aggregate(
          . ~ Group,
          data = to_imp_data,
          function(x) { 100 * (sum(is.na(x)) / (sum(is.na(x)) + sum(!is.na(x))) ) },
          na.action = NULL)
    count_NA <- count_NA %>%
      dplyr::select(-Group)
    correct_names <- names(count_NA)
    supress <- unlist(as.data.frame(lapply(count_NA, function(x) all(x > cutoff))))
    names(supress) <- correct_names
    correct_names <- names(supress[supress == "FALSE"])
    depurdata <- to_imp_data[, 2:ncol(to_imp_data)][!supress]
    depurdata <- sapply(depurdata, function(x) as.numeric(as.character(x)))
  } else {
    depurdata <- to_imp_data[, 2:ncol(to_imp_data)]
    depurdata <- sapply(depurdata, function(x) as.numeric(as.character(x)))
    correct_names <- colnames(prf_tab)
  }

  # Row->feature;Col->sample
  if (method == "none") {
    depurdata[is.na(depurdata)] <- 0
  } else if (method == "LOD") {
    if (is.null(LOD)) {
      message("No LOD provided, regard one-tenth mininal value as LOD")
      depurdata_withoutNA <- depurdata[!is.na(depurdata)]
      LOD <- min(depurdata_withoutNA[depurdata_withoutNA != 0]) / 10
    }
    depurdata[is.na(depurdata)] <- LOD
    depurdata[depurdata == 0] <- LOD
  } else if (method == "half_min") {
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), min(x, na.rm = TRUE)/2, x) else x})
  } else if (method == "median") {
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), median(x, na.rm = TRUE), x) else x})
  } else if (method == "mean") {
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), mean(x, na.rm = TRUE), x) else x})
  } else if (method == "min") {
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), min(x, na.rm = TRUE), x) else x})
  } else if (method == "knn") {
    depurdata <- t(depurdata)
    datai <- impute::impute.knn(depurdata)
    depurdata <- t(datai$data)
  } else if (method == "rf") {
    # depurdata <- data.frame(group = samples_groups, depurdata)
    # depurdata <- randomForest::rfImpute(group ~ ., depurdata)
    # depurdata <- depurdata %>%
    #   dplyr::select(-group)
    fit <- missForest::missForest(t(depurdata))
    depurdata <- fit$ximp %>%
      t()
  } else if (method == "global_mean") {
    depurdata <- .GlobalMean(object = t(depurdata)) %>%
      t()
  } else if (method == "svd") {
    depurdata <- .SVD_wrapper(depurdata)
  } else if (method == "QRILC") {
    fit <- log(t(depurdata)) %>%
      imputeLCMD::impute.QRILC() #%>%
      #extract2(1) %>% exp
    depurdata <- t(fit[[1]])
  }

  colnames(depurdata) <- correct_names
  rownames(depurdata) <- rownames(prf_tab)

  if (methods::is(object, "phyloseq")) {
    res <- ps
    phyloseq::otu_table(res) <- phyloseq::otu_table(t(depurdata),
                                                   taxa_are_rows = TRUE)
  }

  if (methods::is(object, "ExpressionSet")) {
    res <- object
    Biobase::exprs(res) <- depurdata
  }

  return(res)
}


#' Global Standard imputation
#'
#' Global Standard imputation is used the Global mean and standard as the
#' cutoff to impute the missing value. We calculate the mean  for
#' all non missing values of that variable then replace missing
#' value with mean.
#'
#' @keywords internal
#' @noRd
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
.GlobalMean <- function(
    object,
    width = 0.3,
    downshift = 1.8) {

  # Row->Features; Column->Samples
  if (length(object[is.na(object)]) > 1) {
    for (i in 1:ncol(object)) {
      temp <- object[, i]
      temp_sd <- width * sd(temp, na.rm = TRUE)
      temp_mean <- mean(temp, na.rm = TRUE) - downshift * sd(temp, na.rm = TRUE)
      n_missing <- sum(is.na(temp))
      object[is.na(temp), i] <- rnorm(n_missing, mean = temp_mean, sd = temp_sd)
      object_imputed <- object
    }
  } else {
    object_imputed <- object
    message("No NAs were found\n")
  }
  return(object_imputed)
}

#' @keywords internal
#' @noRd
.SVD_wrapper <- function(
    object,
    K = 5) {

  data_sc_res <- .scale_recover(object, method = "scale")
  data_sc <- data_sc_res[[1]]
  data_sc_param <- data_sc_res[[2]]

  result <- data_sc %>%
    imputeLCMD::impute.wrapper.SVD(., K = K) %>%
    .scale_recover(., method = "recover", param_df = data_sc_param) #%>%
    #extract2(1)

  res <- result[[1]]
  return(res)
}

# Scale and recover -------------
#' @keywords internal
#' @noRd
.scale_recover <- function(
    object,
    method = "scale",
    param_df = NULL) {

  results <- list()
  data_res <- object

  if (!is.null(param_df)) {
    if (method == "scale") {
      data_res[] <- scale(object, center=param_df$mean, scale=param_df$std)
    } else if (method=="recover") {
      data_res[] <- t(t(object) * param_df$std + param_df$mean)
    }
  } else {
    if (method == "scale") {
      param_df <- data.frame(mean=apply(object, 2, function(x) mean(x, na.rm=T)),
                             std=apply(object, 2, function(x) sd(x, na.rm=T)))
      data_res[] <- scale(object, center=param_df$mean, scale=param_df$std)
    } else {stop("no param_df found for recover...")}
  }
  results[[1]] <- data_res
  results[[2]] <- param_df

  return(results)
}
