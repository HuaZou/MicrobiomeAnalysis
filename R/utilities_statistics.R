#' @title Calculate the Median abundance of features per group (two groups)
#'
#' @param profile matrix; (Required).
#' @param metadata matrix; (Required).
#' @param groups character; (Required).
#'
#' @importFrom dplyr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @import stats
#'
#' @aliases calculate_median_abundance
#'
#' @export
#'
calculate_median_abundance <- function(profile, metadata, groups){

  if (!all(rownames(metadata) == colnames(profile))) {
    stop("Order of sampleID between colData and proData is wrong please check your inputdata")
  }

  res <- apply(profile, 1, function(x, y){
    dat <- data.frame(value=as.numeric(x), group=y)
    mn <- tapply(dat$value, dat$group, median) %>%
      data.frame() %>% setNames("value") %>%
      tibble::rownames_to_column("Group")
    mn1 <- with(mn, mn[Group%in%groups[1], "value"])
    mn2 <- with(mn, mn[Group%in%groups[2], "value"])
    mnall <- median(dat$value)

    if (all(mn1 != 0, mn2 != 0)) {
      Log2median <- log2(mn1/mn2)
    } else {
      Log2median <- NA
    }

    res <- c(Log2median, mnall, mn1, mn2)
    return(res)
  }, metadata$Compvar) %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("TaxaID")

  colnames(res) <- c("TaxaID",
                     paste0("Log2FoldChange (Median)\n", paste(groups, collapse = "_vs_")),
                     "Median Abundance\n(All)",
                     paste0("Median Abundance\n", groups))

  return(res)
}


#' @title Calculate the Mean abundance of features per group (two groups)
#'
#' @param profile (Required). matrix.
#' @param metadata (Required). matrix.
#' @param groups (Required). character.
#'
#' @importFrom dplyr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @import stats
#'
#' @aliases calculate_mean_abundance
#'
#' @export
#'
calculate_mean_abundance <- function(profile, metadata, groups){

  if (!all(rownames(metadata) == colnames(profile))) {
    stop("Order of sampleID between colData and proData is wrong please check your inputdata")
  }

  res <- apply(profile, 1, function(x, y){
    dat <- data.frame(value=as.numeric(x), group=y)
    mn <- tapply(dat$value, dat$group, mean) %>%
      data.frame() %>% setNames("value") %>%
      tibble::rownames_to_column("Group")
    mn1 <- with(mn, mn[Group%in%groups[1], "value"])
    mn2 <- with(mn, mn[Group%in%groups[2], "value"])
    mnall <- mean(dat$value)

    if (all(mn1 != 0, mn2 != 0)) {
      Log2mean <- log2(mn1/mn2)
    } else {
      Log2mean <- NA
    }

    res <- c(Log2mean, mnall, mn1, mn2)
    return(res)
  }, metadata$Compvar) %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("TaxaID")

  colnames(res) <- c("TaxaID",
                     paste0("Log2FoldChange (Mean)\n", paste(groups, collapse = "_vs_")),
                     "Mean Abundance\n(All)",
                     paste0("Mean Abundance\n", groups))

  return(res)
}


#' @title Calculate the Rank Mean abundance of features per group(two groups)
#'
#' @param profile (Required). matrix.
#' @param metadata (Required). matrix.
#' @param groups (Required). character.
#'
#' @importFrom dplyr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @import stats
#'
#' @aliases calculate_MeanRank_abundance
#'
#' @export
#'
calculate_MeanRank_abundance <- function(profile, metadata, groups) {

  if (!all(rownames(metadata) == colnames(profile))) {
    stop("Order of sampleID between colData and proData is wrong please check your inputdata")
  }

  res <- apply(profile, 1, function(x, y) {
    dat <- data.frame(value=as.numeric(x), group=y)
    rk <- rank(dat$value)
    rnk <- round(tapply(rk, dat$group, mean), 2)  %>%
      data.frame() %>% setNames("value")  %>%
      tibble::rownames_to_column("Group")
    rnk1 <- with(rnk, rnk[Group%in%groups[1], "value"])
    rnk2 <- with(rnk, rnk[Group%in%groups[2], "value"])
    Log2rnk <- log2(rnk1/rnk2)

    res <- c(Log2rnk, rnk1, rnk2)
    return(res)
  }, metadata$Compvar) %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("TaxaID")

  colnames(res) <- c("TaxaID",
                     paste0("Log2FoldChange (Rank)\n", paste(groups, collapse = "_vs_")),
                     paste0("Mean Rank Abundance\n", groups))

  return(res)
}


#' @title Calculate the Mean abundance of features per group (two groups)
#'
#' @param profile (Required). matrix.
#' @param metadata (Required). matrix.
#' @param groups (Required). character.
#'
#' @importFrom dplyr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @import stats
#'
#' @aliases calculate_geometricmean_abundance
#'
#' @export
#'
calculate_geometricmean_abundance <- function(profile, metadata, groups) {

  if (!all(rownames(metadata) == colnames(profile))) {
    stop("Order of sampleID between colData and proData is wrong please check your inputdata")
  }

  res <- apply(profile, 1, function(x, y) {
    dat <- data.frame(value=as.numeric(x), group=y)
    dat$value_scale <- scale(dat$value, center = TRUE, scale = TRUE)
    mn_GM <- tapply(dat$value_scale, dat$group, compositions::geometricmean) %>%
      data.frame() %>% stats::setNames("value") %>%
      rownames_to_column("Group")
    mn_GM1 <- with(mn_GM, mn_GM[Group%in%groups[1], "value"])
    mn_GM2 <- with(mn_GM, mn_GM[Group%in%groups[2], "value"])
    mn_GM_all <- compositions::geometricmean(dat$value_scale)

    if (all(mn_GM1 != 0, !is.nan(mn_GM1), mn_GM2 != 0, !is.nan(mn_GM2))) {
      Log2FC_GM <- log2(mn_GM1/mn_GM2)
    } else {
      Log2FC_GM <- NA
    }

    res <- c(Log2FC_GM, mn_GM_all, mn_GM1, mn_GM2)
    return(res)
  }, metadata$Compvar) %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("TaxaID")

  colnames(res) <- c("TaxaID",
                     paste0("Log2FoldChange (geometricmean)\n", paste(groups, collapse = "_vs_")),
                     "GeometricMean Abundance\n(All)",
                     paste0("GeometricMean Abundance\n", groups))

  return(res)
}


#' @title Calculate the Occurrence of features per group (two groups)
#'
#' @param profile (Required). matrix.
#' @param metadata (Required). matrix.
#' @param groups (Required). character.
#'
#' @importFrom dplyr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @import stats
#'
#' @aliases calculate_occurrence_taxa
#'
#' @export
#'
calculate_occurrence_taxa <- function(profile, metadata, groups) {

  if (!all(rownames(metadata) == colnames(profile))) {
    stop("Order of sampleID between colData and proData is wrong please check your inputdata")
  }

  res <- apply(profile, 1, function(x, y){
    dat <- data.frame(value=as.numeric(x), group=y)
    occ <- tapply(dat$value, dat$group, function(x){
      length(x[c(which(!is.na(x)&x!=0))])/length(x)
    }) %>%
      data.frame() %>% setNames("occ") %>%
      tibble::rownames_to_column("Group")
    occ1 <- occ[occ$Group%in%groups[1], "occ"]
    occ2 <- occ[occ$Group%in%groups[2], "occ"]
    occall <- length(dat$value[c(which(!is.na(dat$value)&dat$value!=0))]) / length(dat$value)

    # 100%
    occ1_percentage <- round(occ1, 4) * 100
    occ2_percentage <- round(occ2, 4) * 100
    occall_percentage <- round(occall, 4) * 100

    res <- c(occall_percentage, occ1_percentage, occ2_percentage)
    return(res)
  }, metadata$Compvar) %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("TaxaID")

  colnames(res) <- c("TaxaID",
                     "Occurrence (100%)\n(All)",
                     paste0("Occurrence (100%)\n", groups))

  return(res)
}


#' @title 95% confidential interval Odds Ratio
#'
#' @description
#' Odds ratio shows the ratio between the probability of Group1 and Group2,
#' indicating the risk between groups. `run_OddRatio` is used glm with
#' scaled values to calculate coefficients of groups to obtain the
#' 95% confidential interval.
#'
#' @author Created by Hua Zou (12/13/2021 Guangzhou China)
#'
#' @param datx (Required). Data.frame. Metadata with SampleID and
#' Group information.
#' @param daty (Required). Data.frame. Profile with SampleID,
#' ordered by datx (Row->Features; Column->SampleID).
#' @param GroupName (Required). Character. the contrast group.
#'
#' @return
#'
#' 95% confidential interval Odds Ratio per feature
#'
#' @importFrom dplyr %>% select all_of
#' @importFrom tibble rownames_to_column
#' @importFrom stats glm
#'
#' @usage run_OddRatio(datx, daty, GroupName)
#'
#' @aliases run_OddRatio
#'
#' @export
#'
run_OddRatio <- function(datx, daty, GroupName){

  if (length(GroupName) != 2) {
    stop("Levels of Group don't equal 2")
  }

  # glm result for odd ratios 95%CI
  mdat <- dplyr::inner_join(datx %>% tibble::rownames_to_column("SampleID") %>%
                              dplyr::select(dplyr::all_of(c("SampleID",
                                                            "Compvar"))),
                            daty %>% t() %>% data.frame() %>%
                              tibble::rownames_to_column("SampleID"),
                            by = "SampleID") %>%
    tibble::column_to_rownames("SampleID")

  # Judge whether the data are not 0 per groups
  occ_fun <- function(x){
    length(x[c(which(!is.na(x)&x != 0))]) / length(x)
  }
  mdat_zero_occ <- mdat %>% dplyr::group_by(Compvar) %>%
    dplyr::summarise(across(everything(), ~ occ_fun(.))) %>%
    tibble::column_to_rownames("Compvar") %>%
    t() %>% data.frame()
  mdat_zero_occ$KEEP <- ifelse(mdat_zero_occ[, 1] > 0 &
                                 mdat_zero_occ[, 2] > 0, TRUE, FALSE)
  nfeature <- rownames(mdat_zero_occ)[mdat_zero_occ$KEEP]
  mdat_remain <- mdat[, colnames(mdat)%in%c("Compvar", nfeature)]

  # Judge whether the data are identical (such as all equal=1)
  idential_res <- apply(mdat_remain[, -1], 2, function(x){
    length( unique(x)[unique(x) != 0] ) != 1
  }) %>%
    data.frame()
  nfeature_idential <- rownames(idential_res)[idential_res$.]
  mdat_remain_final <- mdat_remain[, colnames(mdat_remain) %in%
                                     c("Compvar", nfeature_idential)]

  dat_phe <- mdat_remain_final %>% dplyr::select(Compvar) %>%
    dplyr::mutate(Compvar=ifelse(Compvar == GroupName[2], 1, 0))
  dat_prf <- mdat_remain_final %>% dplyr::select(-Compvar)

  glmFun <- function(GroupN, MarkerN){

    MarkerN[MarkerN==0] <- min(MarkerN[MarkerN!=0])
    dat_glm <- data.frame(group=GroupN,
                          marker=scale(MarkerN, center=TRUE, scale=TRUE)) %>%
      na.omit()
    model <- summary(stats::glm(group ~ marker, data = dat_glm,
                                family = binomial(link = "logit")))
    res <- signif(exp(model$coefficients["marker", 1]) +
                  qnorm(c(0.025, 0.5, 0.975)) * model$coefficients["marker", 1], 2)

    return(res)
  }

  glm_res <- t(apply(dat_prf, 2, function(x, group){
    res <- glmFun(group, as.numeric(x))
    return(res)
  }, dat_phe$Compvar))

  Odd <- glm_res %>% data.frame() %>%
    setNames(c("upper", "expected","lower")) %>%
    dplyr::mutate("Odds Ratio (95% CI)" = paste0(expected,
                                                 " (", lower, ";", upper, ")"))
  Odd$TaxaID <- rownames(glm_res)

  res_ratain <- Odd[, c(5, 4)]

  # drop features
  drop_features <- dplyr::setdiff(rownames(daty), nfeature_idential)
  if (length(drop_features) > 0) {
    res_drop <- data.frame(TaxaID = drop_features,
                           Value = NA) %>%
      stats::setNames(c("TaxaID", "Odds Ratio (95% CI)"))

    res <- rbind(res_ratain, res_drop)
  } else {
    res <- res_ratain
  }

  return(res)
}


#' @title Confidence Interval
#'
#' @description
#' Calculates the confidence interval of a vector of data.
#'
#' @references https://github.com/cran/Rmisc/
#'
#' @author Created by Hua Zou (5/19/2022 Shenzhen China)
#'
#' @param x (Required). vector. a vector of data.
#' @param ci (optional). numeric. confidence interval to be calculated.
#'
#' @return
#' \item{upper}{Upper bound of interval.}
#' \item{mean}{Mean of data.}
#' \item{lower}{Lower bound of interval.}
#'
#' @export
#'
#' @examples
#' run_CI(rnorm(100))
#'
run_CI <- function(x, ci = 0.95) {

  mean_value <- mean(x)
  sd_value <- sd(x)
  n <- length(x)
  error <- stats::qt(ci + (1 - ci) / 2, df = n - 1) * sd_value/ sqrt(n)
  res <- c(upper=mean_value + error,
           mean=mean_value,
           lower=mean_value - error)

  return(res)
}


#' @title Group Confidence Interval
#'
#' @description
#' Calculates the confidence interval of grouped data
#'
#' @references https://github.com/cran/Rmisc/
#'
#' @author Created by Hua Zou (5/19/2022 Shenzhen China)
#'
#' @param x (Required). an `aggregate` compatible formula.
#' @param data (Required). a data frame (or list) from which the variables in formula should be taken.
#' @param ci (optional). numeric. confidence interval to be calculated.
#'
#' @return
#' A data frame consisting of one column for each grouping factor plus
#' three columns for the upper bound, mean and lower bound of the confidence interval
#' for each level of the grouping factor
#'
#' @export
#'
#' @examples
#' run_group.CI(Sepal.Length~Species, iris, 0.95)
#'
run_group.CI <- function(x, data, ci = .95) {

  return(.group.UCL(x,
                    data,
                    FUN=run_CI,
                    ci=ci))
}

.group.UCL <- function(x, data, FUN, ...) {

  d <- stats::aggregate(x, data, FUN=FUN, ...)
  y <- colnames(d)[ncol(d)]
  n <- as.data.frame(d[, y])
  colnames(n) <- sapply(list("upper", "mean", "lower"),
                        function(l){
                          paste(y, l, sep=".")
                        })
  d[ncol(d)] <- NULL

  return(cbind(d, n))
}


#' @title Standard Error
#'
#' @description
#' Calculates the standard error interval of a vector of data
#'
#' @references https://github.com/cran/Rmisc/
#'
#' @author Created by Hua Zou (5/19/2022 Shenzhen China)
#'
#' @param x (Required). vector. a vector of data.
#'
#' @return
#' \item{upper}{Upper bound of interval.}
#' \item{mean}{Mean of data.}
#' \item{lower}{Lower bound of interval.}
#'
#' @export
#'
#' @examples
#' run_STDERR(rnorm(100))
#'
run_STDERR <- function(x) {

  mean_value <- mean(x)
  sd_value <- sd(x)
  n <- length(x)
  error <- (sd_value / sqrt(n))
  res <- c(upper=mean_value + error,
           mean=mean_value,
           lower=mean_value - error)

  return(res)
}


#' @title Group Standard Error Interval
#'
#' @description
#' Calculates the standard error interval of grouped data.
#'
#' @references https://github.com/cran/Rmisc/
#'
#' @author Created by Hua Zou (5/19/2022 Shenzhen China)
#'
#' @param x (Required). an `aggregate` compatible formula.
#' @param data (Required). a data frame (or list) from which the variables in formula should be taken.
#'
#' @return A data frame consisting of one column for each grouping factor plus
#' three columns for the upper bound, mean and lower bound of the standard error
#' interval for each level of the grouping factor.
#'
#' @export
#'
#' @examples
#' run_group.STDERR(Sepal.Length~Species, iris)
#'
run_group.STDERR <- function(x, data) {

  return(.group.UCL(x,
                    data,
                    FUN=run_STDERR))

}


#' @title Summarizes data for mean median sd etc
#'
#' @description
#' Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#'
#' @references https://github.com/cran/Rmisc/
#'
#' @author Created by Hua Zou (5/19/2022 Shenzhen China)
#'
#' @param data (Required). a data frame.
#' @param measurevar (Required). character. the name of a column that contains the variable to be summariezed.
#' @param groupvars (Required). vector. a vector containing names of columns that contain grouping variables.
#' @param na.rm (Required). logical. a boolean that indicates whether to ignore NA's (default: FALSE).
#' @param conf.interval (Required). numeric. the percent range of the confidence interval (default: 0.95).
#' @param .drop (Required). logical. should combinations of variables that do not appear in the input data be
#' preserved (FALSE) or dropped (default: TRUE).
#'
#' @return a data frame with count, mean, standard deviation,
#' standard error of the mean, and confidence interval (default 95%).
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' library(dose)
#' run_summarySE(ToothGrowth, measurevar="len", groupvars=c("supp", "dose"))
#' }
#'
run_summarySE <- function(data,
                          measurevar,
                          groupvars=NULL,
                          na.rm=FALSE,
                          conf.interval=0.95,
                          .drop=TRUE) {

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) {
      sum(!is.na(x))
    } else {
      length(x)
    }
  }

  # This is does the summary; it's not easy to understand...
  # datac <- plyr::ddply(data, groupvars, .drop=.drop,
  #                .fun= function(xx, col, na.rm) {
  #                  c( N    = length2(xx[,col], na.rm=na.rm),
  #                     mean = mean   (xx[,col], na.rm=na.rm),
  #                     sd   = sd     (xx[,col], na.rm=na.rm)
  #                  )
  #                },
  #                measurevar,
  #                na.rm
  # )

  colnames(data)[which(colnames(data) == measurevar)] <- "measure"
  datac <- data %>% dplyr::group_by(.dots = groupvars, .drop=.drop) %>%
    dplyr::summarise(N = length2(measure, na.rm=na.rm),
                     mean = mean(measure, na.rm=na.rm),
                     median = median(measure, na.rm=na.rm),
                     sd = sd(measure, na.rm=na.rm),
                     .groups = "drop")

  # Rename the "mean" column
  colnames(datac)[which(colnames(datac) == "mean")] <- measurevar

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}


#' @title Insert Row into a Matrix
#'
#' @description
#' Insert Row into a Matrix
#'
#' @references https://github.com/cran/miscTools/
#'
#' @author Created by Hua Zou (5/19/2022 Shenzhen China)
#'
#' @param m (Required). a matrix.
#' @param r (Required). numeric. row number where the new row should be inserted.
#' @param v (optional). numeric. values for the new row.
#' @param rName (optional). character. the name of the new row.
#'
#' @return a matrix with one more row than the provided matrix m.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' m <- matrix(1:4, 2)
#' insertRow(m, 2, 5:6)
#' }
#'
insertRow <- function(m, r, v = NA, rName = "") {

  if (!inherits(m, "matrix")) {
    stop("argument 'm' must be a matrix")
  }
  if (r == as.integer(r)) {
    r <- as.integer(r)
  }
  else {
    stop("argument 'r' must be an integer")
  }
  if (length(r) != 1) {
    stop("argument 'r' must be a scalar")
  }
  if (r < 1) {
    stop("argument 'r' must be positive")
  }
  if (r > nrow(m) + 1) {
    stop("argument 'r' must not be larger than the number of rows",
         " of matrix 'm' plus one")
  }
  if (!is.character(rName)) {
    stop("argument 'rName' must be a character string")
  }
  if (length(rName) != 1) {
    stop("argument 'rName' must be a be a single character string")
  }
  nr <- nrow(m)
  nc <- ncol(m)
  rNames <- rownames(m)
  if (is.null(rNames) & rName != "") {
    rNames <- rep("", nr)
  }
  if (r == 1) {
    m2 <- rbind(matrix(v, ncol = nc), m)
    if (!is.null(rNames)) {
      rownames(m2) <- c(rName, rNames)
    }
  }
  else if (r == nr + 1) {
    m2 <- rbind(m, matrix(v, ncol = nc))
    if (!is.null(rNames)) {
      rownames(m2) <- c(rNames, rName)
    }
  }
  else {
    m2 <- rbind(m[1:(r - 1), , drop = FALSE], matrix(v, ncol = nc),
                m[r:nr, , drop = FALSE])
    if (!is.null(rNames)) {
      rownames(m2) <- c(rNames[1:(r - 1)], rName, rNames[r:nr])
    }
  }
  return(m2)
}


#' @title Insert Row into a Matrix
#'
#' @description
#' Insert Row into a Matrix
#'
#' @references https://github.com/cran/miscTools/
#'
#' @author Created by Hua Zou (5/19/2022 Shenzhen China)
#'
#' @param m (Required). a matrix.
#' @param c (Required). numeric. column number where the new column should be inserted.
#' @param v (optional). numeric. values for the new column.
#' @param cName (optional). character. the name of the new column.
#'
#' @return a matrix with one more column than the provided matrix m.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' m <- matrix(1:4, 2)
#' insertCol(m, 2, 5:6)
#' }
#'
insertCol <- function (m, c, v = NA, cName = ""){
  if (!inherits(m, "matrix")) {
    stop("argument 'm' must be a matrix")
  }
  if (c == as.integer(c)) {
    c <- as.integer(c)
  }
  else {
    stop("argument 'c' must be an integer")
  }
  if (length(c) != 1) {
    stop("argument 'c' must be a scalar")
  }
  if (c < 1) {
    stop("argument 'c' must be positive")
  }
  if (c > ncol(m) + 1) {
    stop("argument 'c' must not be larger than the number of columns",
         " of matrix 'm' plus one")
  }
  if (!is.character(cName)) {
    stop("argument 'cName' must be a character string")
  }
  if (length(cName) != 1) {
    stop("argument 'cName' must be a be a single character string")
  }
  nr <- nrow(m)
  nc <- ncol(m)
  cNames <- colnames(m)
  if (is.null(cNames) & cName != "") {
    cNames <- rep("", nc)
  }
  if (c == 1) {
    m2 <- cbind(matrix(v, nrow = nr), m)
    if (!is.null(cNames)) {
      colnames(m2) <- c(cName, cNames)
    }
  }
  else if (c == nc + 1) {
    m2 <- cbind(m, matrix(v, nrow = nr))
    if (!is.null(cNames)) {
      colnames(m2) <- c(cNames, cName)
    }
  }
  else {
    m2 <- cbind(m[, 1:(c - 1), drop = FALSE], matrix(v, nrow = nr),
                m[, c:nc, drop = FALSE])
    if (!is.null(cNames)) {
      colnames(m2) <- c(cNames[1:(c - 1)], cName, cNames[c:nc])
    }
  }
  return(m2)
}
