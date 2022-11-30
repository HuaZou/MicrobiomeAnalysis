#' @title Permutational multivariate analysis of variance (PERMANOVA)
#'
#' @description
#' The `run_PERMANOVA` function is used to identify the association between
#' the community and environmental variables.
#'
#' @details
#' The `run_PERMANOVA` function is used to identify the association between
#' the community and environmental variables, applying the distance in profile
#' and calculating the F statistic between community and variable by permutation
#' test to determine the significance. It can be applied to both
#' [`phyloseq::phyloseq-class`] and [`Biobase::ExpressionSet`] object.
#'
#' @references Anderson, Marti J. "Permutational multivariate analysis of
#' variance." Department of Statistics, University of Auckland,
#' Auckland 26 (2005): 32-46.
#'
#' @author Created by Hua Zou (5/14/2022 Shenzhen China)
#'
#' @param object (Required). a [`phyloseq::phyloseq-class`] or
#' [`Biobase::ExpressionSet`] object.
#' @param level (Optional). character. Summarization
#' level (from \code{rank_names(pseq)}, default: NULL).
#' @param variables (Optional). vector. variables for test (default: all).
#' @param method (Optional). character. Provide one of the currently supported
#' options. See `distanceMethodList` for a detailed list of the supported options
#' and links to accompanying documentation. Options include:
#'  * "unifrac": unweighted UniFrac distance.
#'  * "wunifrac": weighted-UniFrac distance.
#'  * "GUniFrac": The variance-adjusted weighted UniFrac distances (default: alpha=0.5).
#'  * "bray": bray crutis distance.
#'  * "dpcoa": sample-wise distance used in Double Principle Coordinate Analysis.
#'  * "jsd": Jensen-Shannon Divergence.
#'  Alternatively, you can provide a character string that defines a custom
#'  distance method, if it has the form described in `designdist` (default: "bray").
#' @param mode (Optional). character. Since there are missing values
#' in some columns, providing two test model:
#'  * "one": test on one variable per time.
#'  * "all": test on all variables once.
#' (default: "one").
#' @param seedNum (Optional). numeric. specify seeds for reproduction (default: 123).
#' @param alpha (Optional). numeric. the parameter for "GUniFrac" controlling
#' weight on abundant lineages (default: 0.5).
#'
#' @return
#'   A data.frame of PERMANOVA result.
#'
#' @importFrom vegan adonis vegdist
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom Biobase pData exprs
#' @importFrom stats setNames p.adjust
#'
#' @usage run_PERMANOVA(
#'    object,
#'    level = NULL,
#'    variables = "all',
#'    method = "bray",
#'    mode = c("one", "all"),
#'    seedNum = 123,
#'    alpha = 0.5)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data("enterotypes_arumugam")
#' run_PERMANOVA(enterotypes_arumugam,
#'   method = "bray",
#'   mode = "one")
#'
#' data("caporaso")
#' run_PERMANOVA(caporaso,
#'   method = "bray",
#'   mode = "all",
#'   variables = c("Month", "SampleType"))
#' }
#'
run_PERMANOVA <- function(
    object,
    level = NULL,
    variables = "all",
    method = "bray",
    mode = "one",
    seedNum = 123,
    alpha = 0.5) {

  # data("enterotypes_arumugam")
  # object = enterotypes_arumugam
  # level = "Phylum"
  # variables = "all"
  # mode = "all"
  # method = "bray"
  # seedNum = 123
  # alpha = 0.5

  # phyloseq object
  if (all(!is.null(object), inherits(object, "phyloseq"))) {
    ps <- check_sample_names(object = object)

    # taxa level
    if (!is.null(level)) {
      ps <- aggregate_taxa(x = ps, level = level)
    } else {
      ps <- ps
    }

    if (!is.null(ps@phy_tree) & (method %in%
                                 c("unifrac", "wunifrac", "GUniFrac"))) {
      method <- match.arg(
        method,
        c("unifrac", "wunifrac", "GUniFrac")
      )
    } else if (method %in% c("unifrac", "wunifrac", "GUniFrac")) {
      message("It enforces to use Bray-Curtis because no phy_tree")
      method <- "bray"
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

  # distance
  disMatrix <- run_distance(object = object, method = method, alpha = alpha)

  # variables for test
  if (all(length(variables) == 1, variables == "all")) {
    sam_tab <- sam_tab
  } else {
    sam_tab <- sam_tab %>%
      dplyr::select(dplyr::all_of(c("TempRowNames", variables)))
  }

  # set seed
  if (!is.null(seedNum)) {
    set.seed(seedNum)
  }

  if (mode == "one") {
    res <- .one_permanova(x = sam_tab,
                          y = disMatrix,
                          z = prf_tab,
                          m = method)
  } else {
    res <- .all_permanova(x = sam_tab,
                          y = prf_tab,
                          m = method)
  }

  return(res)
}

.one_permanova <- function(x, y, z, m) {

  # x = sam_tab
  # y = disMatrix
  # z = prf_tab
  # m = method

  dat_x <- x %>%
    dplyr::select(-TempRowNames)

  res <- data.frame()
  for (i in 1:ncol(dat_x)) {
    dat <- data.frame(value = dat_x[, i], z)
    na_index <- which(is.na(dat$value))

    if (!length(na_index) == 0) {
      # missing values
      datphe <- dat$value[-na_index]

      if (length(datphe) == 0 | length(unique(datphe)) == 1) {

        res <- c(length(datphe), rep(NA, 6))
      }
      if (length(unique(datphe)) < 6) {
        datphe <- as.factor(datphe)
      }

      dat_cln <- dat[-na_index, ]
      datprf <- dat_cln[, -1, F]
      dis <- vegan::vegdist(datprf, method = m)

    } else {

      datphe <- dat$value

      if (length(datphe) == 0 | length(unique(datphe)) == 1) {
        res <- c(length(datphe), rep(NA, 6))
      }

      if (length(unique(datphe)) < 6) {
        datphe <- as.factor(datphe)
      }

      dis <- y
    }

    # https://github.com/joey711/phyloseq/issues/1457
    ad <- vegan::adonis(unname(dis) ~ datphe, permutations = 999)
    tmp <- as.data.frame(ad$aov.tab) %>% dplyr::slice(1)
    out <- c(length(datphe), as.numeric(tmp[, c(1:6)])) %>%
      data.frame() %>%
      t()

    res <- rbind(res, out)
  }

  colnames(res) <- c("SumsOfSample", "Df", "SumsOfSqs",
                     "MeanSqs", "F.Model", "R2", "Pr(>F)")
  res$AdjustedPvalue <- stats::p.adjust(res$`Pr(>F)`, method = "BH")
  rownames(res) <- colnames(dat_x)

  return(res)
}

.all_permanova <- function(x, y, m) {

  # x = sam_tab
  # y = prf_tab
  # m = method

  # remove missing values of all columns
  phen <- na.omit(x)

  if (nrow(phen) != nrow(x)) {
    # profile
    prof <- y[rownames(y) %in% rownames(phen), , ]
  } else {
    prof <- y
  }

  # remove TempRowNames
  dat_x <- phen %>%
    dplyr::select(-TempRowNames)

  # formula
  fma <- formula(paste("prof",
                       paste0(colnames(dat_x), collapse = "+"),
                       sep = "~"))
  print(fma)

  fit <- vegan::adonis2(formula = fma,
                        data = dat_x,
                        permutations = 999,
                        method = m,
                        by = "margin")

  print(fit)

  temp_res <- data.frame(fit)[1:ncol(dat_x), , F]
  colnames(temp_res) <- c("Df", "SumsOfSqs",
                     "R2", "F.Model", "Pr(>F)")
  temp_res$SumsOfSample <- nrow(dat_x)
  temp_res$AdjustedPvalue <- stats::p.adjust(temp_res$`Pr(>F)`, method = "BH")

  res <- temp_res %>%
    dplyr::select(SumsOfSample, everything())

  return(res)
}