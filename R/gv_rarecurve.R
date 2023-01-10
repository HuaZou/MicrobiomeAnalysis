#' @title Calculating result of rarefaction curve
#'
#' @description
#' The function is to result of rarefaction curve
#'
#' @author Created by Hua Zou (5/14/2022 Shenzhen China)
#'
#' @param object (Required). a [`phyloseq::phyloseq-class`] object.
#' @param level (Optional). character. Summarization
#' level (from \code{rank_names(pseq)}, default: NULL).
#' @param index_name (Required). character.
#' alpha-diversity indices (default: "Observed").
#' Values must be among those supported:
#' c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
#' "InvSimpson", "Fisher", "Evenness", "readsNums").
#' "Chao1", "ACE", "Fisher" and "readsNums" only supported by
#' counts matrix (integers).
#' @param chunks (Optional). integer. the number of index in a sample (default: 200).
#' @param rng_seed (Optional). integer. random seed (default: 123).
#'
#' @return
#'   an object, which could be applied for rarefaction curve plotting.
#'
#' @import dplyr
#' @importFrom stats setNames as.dist
#'
#' @usage run_rarecurve(
#'    object,
#'    level = c(NULL, "Kingdom", "Phylum", "Class",
#'            "Order", "Family", "Genus",
#'            "Species", "Strain", "unique"),
#'    index_name = c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
#'                  "InvSimpson", "Fisher", "Evenness", "readsNums"),
#'    chunks = 200,
#'    rng_seed = 123)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # phyloseq object
#' data("enterotypes_arumugam")
#' run_rarecurve(
#'    object = enterotypes_arumugam,
#'    level = "Phylum",
#'    index_name = "Observed")
#' }
#'
run_rarecurve <- function(
      object,
      level = NULL,
      index_name = "Observed",
      chunks = 200,
      rng_seed = 123) {

  # data("enterotypes_arumugam")
  # object = enterotypes_arumugam
  # level = "Phylum"
  # index_name = "Observed"
  # chunks = 200
  # rng_seed = NULL

  index_name <- match.arg(
    index_name, c("Observed", "Chao1", "ACE", "Shannon",
                  "Simpson", "InvSimpson", "Fisher", "Evenness")
  )
  # phyloseq object
  ps <- check_sample_names(object = object)

  # taxa level
  if (!is.null(level)) {
    ps <- aggregate_taxa(x = ps, level = level)
  } else {
    ps <- ps
  }

  otu_tab <- ps@otu_table %>%
    data.frame()

  # rarefaction
  sdepth <- sum(otu_tab)
  step <- base::trunc(sdepth/chunks)
  if (step == 0) {
    message("The amounts of some samples is less than chunks = ", chunks, "!")
    step = 2
  }

  n <- seq(0, sdepth, by=step)[-1]
  n <- c(n, sdepth)
  out <- lapply(n, function(x) {
    tmp <- withr::with_seed(rng_seed, get_alphaindex(ps = ps, mindepth=x, indices = index_name))
    tmp$readsNums <- x
    return(tmp)
    })

  res <- do.call("rbind", out)
  res[is.na(res)] <- 0

  return(res)
}
