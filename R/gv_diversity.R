#' @title Global view on Alpha diversity on microbiota data
#'
#' @description
#' The function is to use `phyloseq::estimate_richness()` to calculate
#' the alpha diversity.
#' Performs a number of standard alpha diversity estimates, and returns
#' results as data.frame.
#' **NOTE**: You must use untrimmed datasets for meaningful results,
#' as these estimates (and even the "observed" richness) are highly dependent
#' on the number of singletons. You can always trim the data later on if needed,
#' just not before using this function.
#'
#' @references
#' https://scienceparkstudygroup.github.io/microbiome-lesson/04-alpha-diversity/index.html
#'
#' @author Created by Hua Zou (12/02/2022 Shenzhen China)
#'
#' @param ps (Required). a [`phyloseq::phyloseq-class`] object.
#' @param level (Optional). character. taxonomic level to summarize,
#' default the top level rank of the `ps`. taxonomic level(Kingdom, Phylum,
#' Class, Order, Family, Genus, Species, Strains; default: NULL).
#' @param indexes (Optional). character, meaning that all available
#' alpha-diversity indexes will be included. Alternatively, you can specify
#' one or more indexes as a character vector of measure names (default: all).
#' Values must be among those supported:
#' c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
#' "InvSimpson", "Fisher", "Evenness").
#' "Chao1", "ACE" and "Fisher" only supported by
#' counts matrix (integers).
#'
#' @return
#'   A data.frame of alpha index with metadata.
#'
#' @import phyloseq
#' @import dplyr
#'
#' @usage get_alphaindex(
#'     ps,
#'     level = c(NULL, "Kingdom", "Phylum", "Class",
#'            "Order", "Family", "Genus",
#'            "Species", "Strain", "unique"),
#'     indexes = c("all",
#'       "Observed", "Chao1", "ACE", "Shannon",
#'       "Simpson", "InvSimpson", "Fisher", "Evenness")
#'    )
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' # relative abundance
#' data("enterotypes_arumugam")
#' get_alphaindex(
#'   ps = enterotypes_arumugam,
#'   level = "Genus",
#'   indexes = c("Shannon", "Observed"))
#'
#' # absolute abundance
#' data("caporaso")
#' get_alphaindex(
#'   ps = caporaso,
#'   level = "Genus",
#'   indexes = c("Shannon", "Chao1"))
#' }
#'
#'
get_alphaindex <- function(
    ps,
    level = NULL,
    indexes = c("Observed", "Chao1", "ACE", "Shannon",
                "Simpson", "InvSimpson", "Fisher", "Evenness")) {

  # ps = enterotypes_arumugam
  # level = "Genus"
  # indexes = "Shannon"

  # ps = caporaso
  # level = "Genus"
  # indexes = "Shannon"

  stopifnot(inherits(ps, "phyloseq"))
  ps <- preprocess_ps(ps)

  # taxa level
  if (!is.null(level)) {
    ps <- aggregate_taxa(x = ps, level = level)
  } else {
    ps <- ps
  }

  # alpha diversity indexes
  measures <- c("Observed", "Chao1", "ACE", "Shannon",
                "Simpson", "InvSimpson", "Fisher", "Evenness")
  if (all(length(indexes) == 1, indexes == "all")) {
    indexes <- as.character(measures)
  }
  if (any(indexes %in% measures)) {
    indexes[indexes %in% measures] <- measures[measures %in% indexes]
  }
  if (!any(indexes %in% measures)) {
    stop("None of the `indexes` you provided are supported. Try default `NULL` instead.")
  }

  # otu table
  otu_tab <- as(phyloseq::otu_table(ps), "matrix") %>%
    t()

  # shannon Simpson & Evenness (Pielou) (Row sample * Col feature)
  Shannon <- vegan::diversity(otu_tab)
  Simpson <- vegan::diversity(otu_tab, index = "simpson")
  InvSimpson <- vegan::diversity(otu_tab, index = "invsimpson")
  Evenness <- Shannon/log(vegan::specnumber(otu_tab))

  if (identical(all.equal(otu_tab, round(otu_tab)), TRUE)) {
    spn <- vegan::estimateR(otu_tab)
    Fisher <- vegan::fisher.alpha(otu_tab)
    Observed <- spn[1,]
    Chao1 <- spn[2, ]
    ACE <- spn[4, ]
    alpha <- data.frame(Fisher = Fisher,
                        Observed = Observed,
                        Chao1 = Chao1,
                        ACE = ACE,
                        Shannon = Shannon,
                        Simpson = Simpson,
                        InvSimpson = InvSimpson,
                        Evenness = Evenness)
  } else {
    message("Chao1, ACE and Fisher could not be calculated on integer values.")

    Observed <- apply(otu_tab, 1, function(x) {sum(x > 0)})
    # Fisher <- NULL
    # Chao1 <- NULL
    # ACE <- NULL
    alpha <- data.frame(Observed = Observed,
                        Shannon = Shannon,
                        Simpson = Simpson,
                        InvSimpson = InvSimpson,
                        Evenness = Evenness)
  }

  res <- phyloseq::sample_data(ps) %>%
    data.frame() %>%
    tibble::rownames_to_column("TempRowNames") %>%
    dplyr::inner_join(alpha %>%
                        dplyr::select(all_of(indexes)) %>%
                        tibble::rownames_to_column("TempRowNames"),
                      by = "TempRowNames") %>%
    tibble::column_to_rownames("TempRowNames")

  return(res)
}
