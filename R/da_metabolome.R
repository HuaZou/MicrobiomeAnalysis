#' @title Perform differential analysis on metabolomic data
#'
#' @description
#' Differential expression analysis contains foldchange analysis, VIP and t-test
#' on metabolomic data.
#'
#' @details
#' To identify the potential significant metabolites is important
#' in metabolomics. Here, we use three methods to obtain the results:
#'
#' 1. FoldChange: the raw intensity values;
#' 2. VIP (Variable important in projection) by PLS-DA:
#' the normalized intensity values;
#' 3. T-test: the normalized intensity values.
#'
#' Combining all the results and providing more choice for users to determine
#' the differential metabolites (Recommend: Log2FoldChange and AdjustedPvalue).
#'
#' @author Created by Hua Zou (12/2/2023 Shenzhen China)
#'
#' @param object_raw (Optional). a [`SummarizedExperiment::SummarizedExperiment`]
#' object with raw intensity values (default: NULL).
#' @param object_norm (Optional). a [`SummarizedExperiment::SummarizedExperiment`]
#' object with normalized intensity values (default: NULL).
#' @param variable (Required). character. grouping variable for test.
#' @param variable_name (Required). two characters. variable' names.
#' @param DA_method (Optional). character. method for differential analysis
#' (default: "all"). options include:
#'   * "all": three methods.
#'   * "fc": foldchange analysis.
#'   * "vip": VIP analysis.
#'   * "t": t-test analysis.
#' @param method_VIP (Optional). character. method for VIP (default: "PLS"),
#' options include:
#'   * "PLS": Partial Least Squares Discriminant Analysis.
#'   * "OPLS": Orthogonal Partial Least Square Discriminant Analysis.
#' @param cutoff_prev (Optional). Numeric.
#' the Prevalence threshold (default: 0.1).
#' @param p_adjust (Optional). character. method to adjust p-values by.
#' Options include "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none". See [`stats::p.adjust()`]
#' for more details (default: "BH").
#'
#' @return
#'   A data frame of the differential results.
#'
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames
#'
#' @usage run_metabolomeDA(
#'     object_raw = NULL,
#'     object_norm = NULL,
#'     variable,
#'     variable_name,
#'     DA_method = c("all", "fc", "vip", "t"),
#'     method_VIP = c("PLS", "OPLS"),
#'     cutoff_prev = 0.1,
#'     p_adjust = c("none", "fdr", "bonferroni", "holm",
#'                  "hochberg", "hommel", "BH", "BY"))
#'
#' @export
#'
#' @examples
#'
#' data("Zeybel_2022_protein")
#' Zeybel_2022_protein_imput <- impute_abundance(
#'       Zeybel_2022_protein,
#'       group = "LiverFatClass",
#'       method = "knn")
#' Zeybel_2022_protein_norm <- scale_variables(
#'       Zeybel_2022_protein_imput,
#'       method == "zscore")
#' DA_results <- run_metabolomeDA(
#'   object_raw = Zeybel_2022_protein,
#'   object_norm = Zeybel_2022_protein_norm,
#'   variable = "LiverFatClass",
#'   variable_name = c("None", "Severe"))
#'
#' \dontrun{
#'
#' }
#'
run_metabolomeDA <- function(
    object_raw = NULL,
    object_norm = NULL,
    variable,
    variable_name,
    DA_method = c("all", "fc", "vip", "t"),
    method_VIP = "PLS",
    cutoff_prev = 0.1,
    p_adjust = "BH") {

  # object_raw = Zeybel_2022_protein
  # Zeybel_2022_protein_imput <- impute_abundance(Zeybel_2022_protein,
  #                                               group = "LiverFatClass",
  #                                               method = "knn")
  # object_norm = scale_variables(Zeybel_2022_protein_imput, method = "zscore")
  # variable = "LiverFatClass"
  # variable_name = c("None", "Severe")
  # DA_method = "all"
  # method_VIP = "PLS"
  # cutoff_prev = 0.1
  # p_adjust = "BH"

  DA_method <- match.arg(
    DA_method, c("all", "fc", "vip", "t"))

  method_VIP <- match.arg(
    method_VIP, c("PLS", "OPLS"))

  p_adjust <- match.arg(
    p_adjust, c("none", "fdr", "bonferroni", "holm",
                "hochberg", "hommel", "BH", "BY"))

  if (DA_method == "all") {

    # fold change analysis
    tryCatch(
      expr = {
        fc_res <- DA_FoldChange(
          x = object_raw,
          group = variable,
          group_names = variable_name,
          occ_cutoff = cutoff_prev)
      },
      error = function(e){
        message('DA_FoldChange Caught an error!')
        print(e)
      }
    )

    # VIP analysis
    tryCatch(
      expr = {
        VIP_res <- DA_VIP(
          x = object_norm,
          group = variable,
          group_names = variable_name,
          method = method_VIP)
      },
      error = function(e){
        message('DA_VIP Caught an error!')
        print(e)
      }
    )

    # t-test analysis
    tryCatch(
      expr = {
        t_res <- DA_ttest(
          x = object_norm,
          group = variable,
          group_names = variable_name,
          padjut = p_adjust)
      },
      error = function(e){
        message('DA_ttest Caught an error!')
        print(e)
      }
    )

    res <- DA_mergedResults(
      fc_result = fc_res,
      vip_result = VIP_res,
      test_result = t_res,
      group_names = variable_name)
  } else if (DA_method == "fc") {

    # fold change analysis
    tryCatch(
      expr = {
        fc_res <- DA_FoldChange(
          x = object_raw,
          group = variable,
          group_names = variable_name,
          occ_cutoff = cutoff_prev)
      },
      error = function(e){
        message('DA_FoldChange Caught an error!')
        print(e)
      }
    )
    res <- fc_res
  } else if (DA_method == "vip") {

    # VIP analysis
    tryCatch(
      expr = {
        VIP_res <- DA_VIP(
          x = object_norm,
          group = variable,
          group_names = variable_name,
          method = method_VIP)
      },
      error = function(e){
        message('DA_VIP Caught an error!')
        print(e)
      }
    )

    res <- VIP_res
  } else if (DA_method == "t") {

    # t-test analysis
    tryCatch(
      expr = {
        t_res <- DA_ttest(
          x = object_norm,
          group = variable,
          group_names = variable_name,
          padjut = p_adjust)
      },
      error = function(e){
        message('DA_ttest Caught an error!')
        print(e)
      }
    )

    res <- t_res
  }

  return(res)
}


#' fold change between two groups.
#' @keywords internal
#' @noRd
DA_FoldChange <- function(
    x,
    group,
    group_names,
    occ_cutoff) {

  # dataseat
  metadata <- SummarizedExperiment::colData(x) %>%
    as.data.frame()
  profile <- SummarizedExperiment::assay(x) %>%
    as.data.frame()

  tryCatch(
    expr = {
      feature <- SummarizedExperiment::rowData(x) %>%
        as.data.frame()
    },
    error = function(e){
      message('SummarizedExperiment::rowData Caught an error!')
      print(e)
    }
  )

  colnames(metadata)[which(colnames(metadata) == group)] <- "CompVar"
  phenotype <- metadata %>%
    dplyr::filter(CompVar %in% group_names) %>%
    dplyr::mutate(CompVar = as.character(CompVar)) %>%
    dplyr::mutate(CompVar = factor(CompVar, levels = group_names))

  if (length(unique(phenotype$CompVar)) < 2) {
    stop("FoldChange analysis requires only two groups.")
  }

  sid <- intersect(rownames(phenotype), colnames(profile))
  phen <- phenotype[pmatch(sid, rownames(phenotype)), , ]
  prof <- profile %>%
    dplyr::select(dplyr::all_of(sid))

  if (!all(colnames(prof) == rownames(phen))) {
    stop("Wrong Order")
  }

  .FeatureOrSample_trim <- function(x, nRow, threshold) {

    df_occ <- apply(x, nRow, function(x) {
      length(x[c(which(!is.na(x) & x!=0))]) / length(x)
    }) %>%
      data.frame() %>% stats::setNames("Occ") %>%
      tibble::rownames_to_column("type")
    if (nRow == 1) {
      rownames(df_occ) <- rownames(x)
    } else {
      rownames(df_occ) <- colnames(x)
    }
    df_KEEP <- apply(df_occ > threshold, 1, all) %>%
      data.frame() %>% stats::setNames("Status") %>%
      dplyr::filter(Status)

    res <- x %>%
      tibble::rownames_to_column("featureid") %>%
      dplyr::filter(featureid %in% rownames(df_KEEP)) %>%
      tibble::column_to_rownames("featureid")

    return(res)
  }
  tryCatch(
    expr = {
      prof_cln <- .FeatureOrSample_trim(prof, 1, occ_cutoff)
    },
    error = function(e){
      message('test_fun Caught an error!')
      print(e)
    }
  )

  fc_res <- apply(prof_cln, 1, function(x1, y1) {

    dat <- data.frame(value = as.numeric(x1), group = y1)
    mn <- tapply(dat$value, dat$group, mean, na.rm = TRUE) %>%
      data.frame() %>%
      stats::setNames("value") %>%
      tibble::rownames_to_column("Group")

    mn1 <- with(mn, mn[Group %in% group_names[1], "value"])
    mn2 <- with(mn, mn[Group %in% group_names[2], "value"])
    mnall <- mean(dat$value, na.rm = TRUE)

    if (all(mn1 != 0, mn2 != 0)) {
      fc <- mn1 / mn2
    } else {
      fc <- NA
    }
    logfc <- log2(fc)

    res <- c(fc, logfc, mnall, mn1, mn2)
    return(res)
  }, phen$CompVar) %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Feature")

  colnames(fc_res) <- c("FeatureID", "FoldChange",
                        "Log2FoldChange",
                        "Mean Abundance\n(All)",
                        paste0("Mean Abundance\n", group_names))

  # Number of Group
  dat_status <- table(phen$CompVar)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  fc_res$Block <- paste(paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                        "vs",
                        paste(dat_status_number[2], dat_status_name[2], sep = "_"))

  if (ncol(feature) != 0) {
    res <- fc_res %>%
      dplyr::select(FeatureID, Block, everything()) %>%
      dplyr::inner_join(feature %>%
                          tibble::rownames_to_column("FeatureID"),
                        by = "FeatureID")
  } else {
    res <- fc_res %>%
      dplyr::select(FeatureID, Block, everything())
  }

  return(res)
}

#' VIP between two groups.
#' @keywords internal
#' @noRd
DA_VIP <- function(
    x,
    group,
    group_names,
    method) {

  # dataseat
  metadata <- SummarizedExperiment::colData(x) %>%
      as.data.frame()
  profile <- SummarizedExperiment::assay(x) %>%
      as.data.frame()
  tryCatch(
    expr = {
      feature <- SummarizedExperiment::rowData(x) %>%
        as.data.frame()
    },
    error = function(e){
      message('SummarizedExperiment::rowData Caught an error!')
      print(e)
    }
  )

  colnames(metadata)[which(colnames(metadata) == group)] <- "CompVar"
  phenotype <- metadata %>%
    dplyr::filter(CompVar %in% group_names) %>%
    dplyr::mutate(CompVar = as.character(CompVar)) %>%
    dplyr::mutate(CompVar = factor(CompVar, levels = group_names))

  if (length(unique(phenotype$CompVar)) < 2) {
    stop("VIP analysis requires only two groups.")
  }

  sid <- intersect(rownames(phenotype), colnames(profile))
  phen <- phenotype[pmatch(sid, rownames(phenotype)), , ]
  prof <- profile %>%
    dplyr::select(dplyr::all_of(sid))

  if (!all(colnames(prof) == rownames(phen))) {
    stop("Wrong Order")
  }

  prof_cln <- prof

  dataMatrix <- prof_cln %>% t() # row->sampleID; col->features
  sampleMetadata <- phen # row->sampleID; col->features
  variableMetadata <- feature # row->features; col->character

  comparsionVn <- sampleMetadata[, "CompVar"]

  pvaVn <- apply(dataMatrix, 2,
                 function(feaVn) cor.test(as.numeric(comparsionVn), feaVn)[["p.value"]])

  if (method == "OPLS") {
    vipVn <- ropls::getVipVn(
          ropls::opls(dataMatrix,
                      comparsionVn,
                      predI = 1,
                      orthoI = NA,
                      fig.pdfC = "none"))
  } else if (method == "PLS") {
    vipVn <- ropls::getVipVn(
          ropls::opls(dataMatrix,
                      comparsionVn,
                      predI = 1,
                      fig.pdfC = "none"))
  }
  quantVn <- qnorm(1 - pvaVn / 2)
  rmsQuantN <- sqrt(mean(quantVn^2))

  res_temp <- data.frame(
      FeatureID = names(vipVn),
      VIP = vipVn,
      CorPvalue = pvaVn) %>%
    dplyr::arrange(desc(VIP))


  # Number of Group
  dat_status <- table(phen$CompVar)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  res_temp$Block <- paste(paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                       "vs",
                       paste(dat_status_number[2], dat_status_name[2], sep = "_"))

  if (ncol(feature) != 0) {
    res <- res_temp %>%
      dplyr::select(FeatureID, Block, everything()) %>%
      dplyr::inner_join(feature %>%
                          tibble::rownames_to_column("FeatureID"),
                        by = "FeatureID")
  } else {
    res <- res_temp %>%
      dplyr::select(FeatureID, Block, everything())
  }

  return(res)
}

#' t-test between two groups.
#' @keywords internal
#' @noRd
DA_ttest <- function(
    x,
    group,
    group_names,
    padjut) {

  # dataseat
  metadata <- SummarizedExperiment::colData(x) %>%
    as.data.frame()
  profile <- SummarizedExperiment::assay(x) %>%
    as.data.frame()
  tryCatch(
    expr = {
      feature <- SummarizedExperiment::rowData(x) %>%
        as.data.frame()
    },
    error = function(e){
      message('SummarizedExperiment::rowData Caught an error!')
      print(e)
    }
  )

  colnames(metadata)[which(colnames(metadata) == group)] <- "CompVar"
  phenotype <- metadata %>%
    dplyr::filter(CompVar %in% group_names) %>%
    dplyr::mutate(CompVar = as.character(CompVar)) %>%
    dplyr::mutate(CompVar = factor(CompVar, levels = group_names))

  if (length(unique(phenotype$CompVar)) < 2) {
    stop("t-test analysis requires only two groups.")
  }

  sid <- intersect(rownames(phenotype), colnames(profile))
  phen <- phenotype[pmatch(sid, rownames(phenotype)), , ]
  prof <- profile %>%
    dplyr::select(dplyr::all_of(sid))

  if (!all(colnames(prof) == rownames(phen))) {
    stop("Wrong Order")
  }

  prof_cln <- prof
  t_res <- apply(prof_cln, 1, function(x1, y1) {
    dat <- data.frame(value = as.numeric(x1), group = y1)

    rest <- t.test(data = dat, value ~ group)

    res <- c(rest$statistic, rest$p.value)
    return(res)
  }, phen$CompVar) %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("Feature")

  colnames(t_res) <- c("FeatureID", "Statistic", "Pvalue")
  t_res$AdjustedPvalue <- p.adjust(as.numeric(t_res$Pvalue), method = padjut)

  # Number of Group
  dat_status <- table(phen$CompVar)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  t_res$Block <- paste(paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                       "vs",
                       paste(dat_status_number[2], dat_status_name[2], sep = "_"))

  if (ncol(feature) != 0) {
    res <- t_res %>%
      dplyr::select(FeatureID, Block, everything()) %>%
      dplyr::inner_join(feature %>%
                          tibble::rownames_to_column("FeatureID"),
                        by = "FeatureID")
  } else {
    res <- t_res %>%
      dplyr::select(FeatureID, Block, everything())
  }

  return(res)
}

#' Combination of results.
#' @keywords internal
#' @noRd
DA_mergedResults <- function(
    fc_result,
    vip_result,
    test_result,
    group_names) {

  overlap_cols <- intersect(
    intersect(colnames(fc_result),
              colnames(vip_result)),
    colnames(test_result))

  overlap_cols <- overlap_cols[overlap_cols != "FeatureID"]

  mdat <- fc_result %>%
    dplyr::mutate(Block2 = paste(group_names, collapse = " vs ")) %>%
    dplyr::mutate(FeatureID = make.names(FeatureID)) %>%
    dplyr::inner_join(vip_result %>%
                        dplyr::select(-dplyr::all_of(overlap_cols)) %>%
                        dplyr::mutate(FeatureID = make.names(FeatureID)),
                      by = "FeatureID") %>%
    dplyr::inner_join(test_result %>%
                        dplyr::select(-dplyr::all_of(overlap_cols)) %>%
                        dplyr::mutate(FeatureID = make.names(FeatureID)),
                      by = "FeatureID")

  res <- mdat %>%
    dplyr::select(FeatureID, Block2, Block,
                  FoldChange, Log2FoldChange,
                  VIP, CorPvalue,
                  Statistic, Pvalue, AdjustedPvalue,
                  everything()) %>%
    dplyr::arrange(AdjustedPvalue, Log2FoldChange)

  return(res)
}
