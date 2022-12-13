#' @title Multi-response Permutation Procedures (MRPP)
#'
#' @description
#' The `run_MRPP` function is used to identify the association between
#' the community and environmental variables by average distance.
#'
#' @details
#' The `run_MRPP` function is used to identify the association between
#' the community and environmental variables by average distance, applying the
#' distance in profile and calculating the delta statistic between community and
#' variable by permutation test to determine the significance. It can be applied
#' to both [`phyloseq::phyloseq-class`] and [`SummarizedExperiment::SummarizedExperiment`] object.
#'
#' @references Mielke Jr, Paul W. "The application of multivariate permutation
#' methods based on distance functions in the earth sciences."
#' Earth-Science Reviews 31.1 (1991): 55-71.
#'
#' @author Created by Hua Zou (5/15/2022 Shenzhen China)
#'
#' @param object (Required). a [`phyloseq::phyloseq-class`] or
#' [`SummarizedExperiment::SummarizedExperiment`] object.
#' @param level (Optional). character. Summarization
#' level (from \code{rank_names(pseq)}, default: NULL).
#' @param variable (Required). character. grouping variable for test.
#' @param method (Optional). character. Provide one of the currently supported
#' options. See `distanceMethodList` for a detailed list of the supported options
#' and links to accompanying documentation. Options include:
#'  * "bray": bray crutis distance.
#'  * "unifrac" : unweighted UniFrac distance.
#'  * "wunifrac": weighted-UniFrac distance.
#'  * "GUniFrac": The variance-adjusted weighted UniFrac distances (default: alpha=0.5).
#'  * "dpcoa": sample-wise distance used in Double Principle Coordinate Analysis.
#'  * "jsd": Jensen-Shannon Divergence.
#'  Alternatively, you can provide a character string that defines a custom
#'  distance method, if it has the form described in `designdist` (default: "bray").
#' @param seedNum (Optional). numeric. specify seeds for reproduction (default: 123).
#' @param alpha (Optional). numeric. the parameter for "GUniFrac" controlling
#' weight on abundant lineages (default: 0.5).
#'
#' @return
#'   A MRPP model.
#'
#' @importFrom vegan vegdist mrpp
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom SummarizedExperiment colData assay
#' @importFrom stats setNames p.adjust
#'
#' @usage run_MRPP(
#'    object,
#'    level = c(NULL, "Kingdom", "Phylum", "Class",
#'            "Order", "Family", "Genus",
#'            "Species", "Strain", "unique"),
#'    variable,
#'    method = c("unifrac", "wunifrac", "GUniFrac", "bray", "dpcoa", "jsd"),
#'    seedNum = 123,
#'    alpha = 0.5)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # phyloseq object
#' data("Zeybel_2022_gut")
#' run_MRPP(Zeybel_2022_gut,
#'           variable = "LiverFatClass",
#'           method = "bray")
#'
#' # SummarizedExperiment object
#' data("Zeybel_2022_protein")
#' run_MRPP(Zeybel_2022_protein,
#'           variable = "LiverFatClass",
#'           method = "bray")
#' }
#'
run_MRPP <- function(
    object,
    level = NULL,
    variable,
    method = c("bray", "unifrac", "wunifrac",
               "GUniFrac", "dpcoa", "jsd"),
    seedNum = 123,
    alpha = 0.5) {

  # data("Zeybel_2022_gut")
  # object = Zeybel_2022_gut
  # level = "Phylum"
  # variable = "LiverFatClass"
  # method = "bray"
  # seedNum = 123
  # alpha = 0.5

  # data("Zeybel_2022_protein")
  # object = Zeybel_2022_protein
  # level = NULL
  # variable = "LiverFatClass"
  # method = "bray"
  # seedNum = 123
  # alpha = 0.5

  method <- match.arg(
    method, c("bray", "unifrac", "wunifrac",
              "GUniFrac", "dpcoa", "jsd")
  )

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
  } else if (all(!is.null(object), inherits(object, "SummarizedExperiment"))) {
    # sample table & profile table
    sam_tab <- SummarizedExperiment::colData(object) %>%
      data.frame() %>%
      tibble::rownames_to_column("TempRowNames")
    prf_tab <- SummarizedExperiment::assay(object) %>%
      data.frame() %>%
      t()
  }
  # distance
  disMatrix <- run_distance(object = object, method = method, alpha = alpha)
  # grouping variable for test
  if (variable %in% colnames(sam_tab)) {
    colnames(sam_tab)[which(colnames(sam_tab) == variable)] <- "GroupVar"
  } else {
    stop(variable, " no exists in columns of Metadata please check your input")
  }

  MRPP_core <- function(x, data_distance, profile){

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

    MRPP_res <- vegan::mrpp(dis, datphe, permutations = 999)
    return(MRPP_res)
  }

  res <- MRPP_core(x=sam_tab$GroupVar,
                   data_distance=disMatrix,
                   profile=prf_tab)

  return(res)
}
