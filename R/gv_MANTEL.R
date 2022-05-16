#' @title Mantel Test (MANTEL)
#'
#' @description
#' The `run_MANTEL` function is used to identify the association between
#' the community and environmental variables.
#'
#' @details
#' The `run_MANTEL` function is used to test the correlation between two
#' distance matrices, calculating the Z statistic between community distance
#' and variables distance to determine the significance. While the Mantel test
#' is used to compare between two distance (dissimilarity) matrices, such as
#' X and Y, the partial Mantel test is used to estimate the correlation
#' between these two matrices, while controlling for the effect of
#' a control matrix Z. It can be applied to both
#' [`phyloseq::phyloseq-class`] and [`Biobase::ExpressionSet`] object.
#'
#' @references Mantel, Nathan. "The detection of disease clustering and
#' a generalized regression approach." Cancer research 27.2 Part 1
#' (1967): 209-220.
#'
#' @author Created by Hua Zou (5/15/2022 Shenzhen China)
#'
#' @param object (Required). a [`phyloseq::phyloseq-class`] or
#' [`Biobase::ExpressionSet`] object.
#' @param y_variables (Required). vector. variables for a control matrix Y.
#' @param z_variables (Optional). vector. variables for a control matrix Z (default: NULL).
#' @param norm (Optional). logical. whether to norm y and z matrix into unit by
#' `scale` (default: TRUE).
#' @param method (Required). character. Mantel test method, "mantel" or
#' "mantel.partial" or "mantel.randtest", "mantel.rtest" (default: "mantel").
#' @param method_cor (Optional). character. Correlation method, as
#' accepted by cor: "pearson", "spearman" or "kendall" (default: "spearman").
#' @param method_dist (Required). character. methods for three(X, Y, Z matrix)
#' distance matrices. Provide one of the currently supported
#' options. See `distanceMethodList` for a detailed list of the supported options
#' and links to accompanying documentation. Options include:
#'  * "unifrac" : unweighted UniFrac distance.
#'  * "wunifrac": weighted-UniFrac distance.
#'  * "GUniFrac": The variance-adjusted weighted UniFrac distances (default: alpha=0.5).
#'  * "bray": bray crutis distance.
#'  * "dpcoa": sample-wise distance used in Double Principle Coordinate Analysis.
#'  * "jsd": Jensen-Shannon Divergence.
#'  Alternatively, you can provide a character string that defines a custom
#'  distance method, if it has the form described in `designdist`
#'  (default: c("bray", "euclidean", "jaccard")).
#' @param seedNum (Optional). numeric. specify seeds for reproduction (default: 123).
#' @param alpha (Optional). numeric. the parameter for "GUniFrac" controlling
#' weight on abundant lineages (default: 0.5).
#'
#' @return
#'   A Mantel test model.
#'
#' @importFrom vegan vegdist mantel mantel.partial
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom Biobase pData exprs
#' @importFrom stats setNames p.adjust
#'
#' @usage run_MANTEL(object,
#'                   y_variables,
#'                   z_variables,
#'                   norm,
#'                   method,
#'                   method_cor,
#'                   method_dist,
#'                   seedNum,
#'                   alpha)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data("enterotypes_arumugam")
#' run_MANTEL(enterotypes_arumugam,
#'           y_variables = c("Enterotype", "Clinical.Status"),
#'           z_variables = c("Nationality", "Gender"),
#'           norm = FALSE,
#'           method = "mantel.partial",
#'           method_dist = c("bray", "euclidean", "jaccard"))
#' }
#'
run_MANTEL <- function(
              object,
              y_variables,
              z_variables = NULL,
              norm = TRUE,
              method = c("mantel", "mantel.partial", "mantel.randtest", "mantel.rtest"),
              method_cor = "spearman",
              method_dist = c("bray", "euclidean", "jaccard"),
              seedNum = 123,
              alpha = 0.5) {

  # Mantel test method
  method <- match.arg(
    method,
    c("mantel", "mantel.partial", "mantel.randtest", "mantel.rtest")
  )
  if (!is.null(z_variables)) {
    if (method == "mantel") {
      message("Force to choose mantel.partial")
    }
    method <- "mantel.partial"
  }
  # Correlation method
  method_cor <- match.arg(method_cor, c("pearson", "spearman","kendall"))
  # distance method
  if (!is.null(z_variables)) {
    if (is.na(method_dist[3])) {
      stop("Please provide the method to calculate Z matrix distance")
    }
  }

  # phyloseq object
  if (all(!is.null(object), inherits(object, "phyloseq"))) {
    ps <- check_sample_names(object = object)
    if (!is.null(ps@phy_tree) & (method_dist[1] %in%
                                 c("unifrac", "wunifrac", "GUniFrac"))) {
      method_dist[1] <- match.arg(
        method_dist[1],
        c("unifrac", "wunifrac", "GUniFrac")
      )
    } else if (method_dist[1] %in% c("unifrac", "wunifrac", "GUniFrac")) {
      message("It enforces to use Bray-Curtis because no phy_tree")
      method_dist[1] <- "bray"
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
  # x distance
  x_dis <- run_distance(object = object, method = method_dist[1], alpha = alpha)
  # y distance
  y_dis <- .cal_variable_distance(dat = sam_tab,
                                  vars = y_variables,
                                  approach = method_dist[2],
                                  norm = norm)

  # z distance
  if (!is.null(z_variables)) {
    z_dis <- .cal_variable_distance(dat = sam_tab,
                                    vars = z_variables,
                                    approach = method_dist[3],
                                    norm = norm)
  }

  # set seed
  if (!is.null(seedNum)) {
    set.seed(seedNum)
  }

  if (method == "mantel") {
    res <- vegan::mantel(
              xdis = x_dis,
              ydis = y_dis,
              method = method_cor,
              permutations = 999)
  } else if (method == "mantel.partial") {
    res <- vegan::mantel.partial(
              xdis = x_dis,
              ydis = y_dis,
              zdis = z_dis,
              method = method_cor,
              permutations = 999)
  } else if (method == "mantel.randtest") {
    res <- ade4::mantel.randtest(
              m1 = x_dis,
              m2 = y_dis,
              nrepet = 999)
  } else if (method == "mantel.rtest") {
    res <- ade4::mantel.rtest(
              m1 = x_dis,
              m2 = y_dis,
              nrepet = 999)
  }

  return(res)
}

#' calculate the distance of Y and Z matrix
#' @noRd
.cal_variable_distance <- function(dat,
                                   vars,
                                   approach,
                                   norm) {

  var_df <- dat %>% dplyr::select(dplyr::all_of(vars))
  # only numeric data for calculating distance
  var_df_numeric <- apply(var_df, 2, function(x){
    if(is.character(x)) {as.numeric(as.factor(x))} else {x}
  })
  rownames(var_df_numeric) <- rownames(var_df)
  # environmental variables were all measured using different metrics
  # that are not comparable to each other
  if (norm) {
    var_scale <- scale(var_df_numeric, center = TRUE, scale = TRUE)
  } else {
    var_scale <- var_df_numeric
  }

  if (any(data.frame(var_scale) < 0)) {
    if (approach %in% c("bray", "jaccard")) {
      message(approach, " is not suitable for Negative values and Replace it by euclidean")
      approach <- "euclidean"
    }
  }

  if (all(approach == "distHaversine", ncol(var_df_numeric) == 2)) {
    # geographic data frame - haversine distance (longitude and latitude)
    dis <- geosphere::distm(var_df_numeric, fun = geosphere::distHaversine)
  } else {
    dis <- vegan::vegdist(as.matrix(var_scale), method = approach)
  }

  return(dis)
}
