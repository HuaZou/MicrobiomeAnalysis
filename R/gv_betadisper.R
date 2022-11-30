#' @title Multivariate homogeneity of groups dispersions (variances)
#'
#' @description
#' The `run_betadisper` function is a multivariate analogue of Levene's test
#' for homogeneity of variances.
#'
#' @details
#' One measure of multivariate dispersion (variance) for a group of samples is
#' to calculate the average distance of group members to the group centroid or
#' spatial median (both referred to as 'centroid' from now on unless stated
#' otherwise) in multivariate space. To test if the dispersions (variances)
#' of one or more groups are different, the distances of group members to
#' the group centroid are subject to ANOVA. This is a multivariate analogue
#' of Levene's test for homogeneity of variances if the distances between
#' group members and group centroids is the Euclidean distance.
#' See `?vegan::betadisper` for more details. It can be applied to both
#' [`phyloseq::phyloseq-class`] and [`Biobase::ExpressionSet`] object.
#'
#' @references Anderson, Marti J., Kari E. Ellingsen, and Brian H. McArdle.
#' "Multivariate dispersion as a measure of beta diversity."
#' Ecology letters 9.6 (2006): 683-693.
#'
#' @author Created by Hua Zou (5/15/2022 Shenzhen China)
#'
#' @param object (Required). a [`phyloseq::phyloseq-class`] or
#' [`Biobase::ExpressionSet`] object.
#' @param level (Optional). character. Summarization
#' level (from \code{rank_names(pseq)}, default: NULL).
#' @param variable (Required). character. grouping variable for test.
#' @param type (Optional). character. the type of analysis to perform. Use the
#' spatial "median" or the group "centroid" (default: "median").
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
#'   A betadisper model.
#'
#' @importFrom vegan vegdist betadisper permutest
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom Biobase pData exprs
#' @importFrom stats setNames p.adjust
#'
#' @usage run_betadisper(
#'    object,
#'    level = c(NULL, "Kingdom", "Phylum", "Class",
#'            "Order", "Family", "Genus",
#'            "Species", "Strain", "unique"),
#'    variable,
#'    type = c("median","centroid"),
#'    method = c("unifrac", "wunifrac", "GUniFrac", "bray", "dpcoa", "jsd"),
#'    seedNum = 123,
#'    alpha = 0.5)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data("enterotypes_arumugam")
#' run_betadisper(enterotypes_arumugam,
#'   variable = "Enterotype",
#'   method = "bray")
#' }
#'
run_betadisper <- function(
    object,
    level = NULL,
    variable,
    type = "median",
    method = "bray",
    seedNum = 123,
    alpha = 0.5) {

  # data("enterotypes_arumugam")
  # object = enterotypes_arumugam
  # level = "Phylum"
  # variable = "Enterotype"
  # type = "median"
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
    sam_tab <- phyloseq::sample_data(ps) %>% data.frame() %>%
      tibble::rownames_to_column("TempRowNames")
    if (phyloseq::taxa_are_rows(ps)) {
      prf_tab <- phyloseq::otu_table(phyloseq::t(ps)) %>%
        data.frame()
    } else {
      prf_tab <- phyloseq::otu_table(ps) %>% data.frame()
    }
  } else if (all(!is.null(object), inherits(object, "ExpressionSet"))) {
    # sample table & profile table
    sam_tab <- Biobase::pData(object) %>% data.frame() %>%
      tibble::rownames_to_column("TempRowNames")
    prf_tab <- Biobase::exprs(object) %>% data.frame()
  }

  # distance
  disMatrix <- run_distance(object = object, method = method, alpha = alpha)

  # grouping variable for test
  if (variable %in% colnames(sam_tab)) {
    colnames(sam_tab)[which(colnames(sam_tab) == variable)] <- "GroupVar"
  } else {
    stop(variable, " no exists in columns of Metadata please check your input")
  }

  betadisper_core <- function(x, data_distance, profile){

    dat <- data.frame(value = x, profile)
    na_index <- which(is.na(dat$value))

    if (!length(na_index) == 0) {
      datphe <- dat$value[-na_index]

      if (length(datphe) == 0 | length(unique(datphe)) == 1) {
        res <- NA
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
        res <- NA
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

    # Checking the homogeneity condition
    mod <- vegan::betadisper(dis, datphe, type = type)
    betadisper_res <- vegan::permutest(mod, pairwise = TRUE,
                                       permutations = permute::how(nperm = 999))
    return(betadisper_res)
  }

  res <- betadisper_core(x = sam_tab$GroupVar,
                         data_distance = disMatrix,
                         profile = prf_tab)

  return(res)
}
