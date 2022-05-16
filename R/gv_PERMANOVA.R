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
#' @param variables (Optional). vector. variables for test.
#' @param method (Optional). character. Provide one of the currently supported
#' options. See `distanceMethodList` for a detailed list of the supported options
#' and links to accompanying documentation. Options include:
#'  * "unifrac" : unweighted UniFrac distance.
#'  * "wunifrac": weighted-UniFrac distance.
#'  * "GUniFrac": The variance-adjusted weighted UniFrac distances (default: alpha=0.5).
#'  * "bray": bray crutis distance.
#'  * "dpcoa": sample-wise distance used in Double Principle Coordinate Analysis.
#'  * "jsd": Jensen-Shannon Divergence.
#'  Alternatively, you can provide a character string that defines a custom
#'  distance method, if it has the form described in `designdist` (default: "bray").
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
#' @usage run_PERMANOVA(object,
#'                      variables,
#'                      method,
#'                      seedNum,
#'                      alpha)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data("enterotypes_arumugam")
#' run_PERMANOVA(enterotypes_arumugam, method = "bray")
#' }
#'
run_PERMANOVA <- function(
              object,
              variables = NULL,
              method = "bray",
              seedNum = 123,
              alpha = 0.5) {

  # phyloseq object
  if (all(!is.null(object), inherits(object, "phyloseq"))) {
    ps <- check_sample_names(object = object)
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
    sam_tab <- phyloseq::sample_data(ps) %>% data.frame() %>%
      tibble::rownames_to_column("SampleID")
    if (phyloseq::taxa_are_rows(ps)) {
      prf_tab <- phyloseq::otu_table(phyloseq::t(ps)) %>%
        data.frame()
    } else {
      prf_tab <- phyloseq::otu_table(ps) %>% data.frame()
    }
  } else if (all(!is.null(object), inherits(object, "ExpressionSet"))) {
    # sample table & profile table
    sam_tab <- Biobase::pData(object) %>% data.frame() %>%
      tibble::rownames_to_column("SampleID")
    prf_tab <- Biobase::exprs(object) %>% data.frame()
  }
  # distance
  disMatrix <- run_distance(object = object, method = method, alpha = alpha)
  # variables for test
  if (!is.null(variables)) {
    sam_tab <- sam_tab %>%
      dplyr::select(dplyr::all_of(c("SampleID", variables)))
  }

  res <- apply(sam_tab, 2, function(x, data_distance, profile){
    dat <- data.frame(value = x, profile)
    na_index <- which(is.na(dat$value))

    if (!length(na_index) == 0) {
      datphe <- dat$value[-na_index]

      if (length(datphe) == 0 | length(unique(datphe)) == 1) {
        res <- c(length(datphe), rep(NA, 6))
      }
      if (length(unique(datphe)) < 6) {
        datphe <- as.factor(datphe)
      }

      dat_cln <- dat[-na_index, ]
      datprf <- dat_cln[, -1, F]
      dis <- vegan::vegdist(datprf, method = method)

    } else {
      datphe <- dat$value
      if (length(datphe) == 0 | length(unique(datphe)) == 1) {
        res <- c(length(datphe), rep(NA, 6))
      }
      if (length(unique(datphe)) < 6) {
        datphe <- as.factor(datphe)
      }
      dis <- data_distance
    }

    # set seed
    if (!is.null(seedNum)) {
      set.seed(seedNum)
    }

    # https://github.com/joey711/phyloseq/issues/1457
    ad <- vegan::adonis(unname(dis) ~ datphe, permutations = 999)
    tmp <- as.data.frame(ad$aov.tab) %>% dplyr::slice(1)
    out <- c(length(datphe), as.numeric(tmp[, c(1:6)]))
    return(out)
  }, disMatrix, prf_tab) %>%
    base::t() %>% data.frame()

  colnames(res) <- c("SumsOfSample", "Df", "SumsOfSqs",
                     "MeanSqs", "F.Model", "R2", "Pr(>F)")
  res$AdjustedPvalue <- stats::p.adjust(res$`Pr(>F)`, method = "BH")

  return(res)
}
