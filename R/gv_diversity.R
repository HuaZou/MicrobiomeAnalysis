#' @title Calculating index of Alpha diversity on microbiota data
#'
#' @description
#' The function is to calculate alpha diversity estimates.
#'
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
#' @param indices (Optional). character, meaning that all available
#' alpha-diversity indices will be included. Alternatively, you can specify
#' one or more indices as a character vector of measure names (default: all).
#' Values must be among those supported:
#' c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
#' "InvSimpson", "Fisher", "Evenness", "TaxaNumber").
#' "Chao1", "ACE" and "Fisher" only supported by
#' counts matrix (integers).
#' @param mindepth (Optional). numeric,
#' Subsample size for rarefying community (default: all).
#' @param force (Optional). logical,
#' whether to rarefy samples (default: FASLE).
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
#'     indices = c("all",
#'       "Observed", "Chao1", "ACE", "Shannon",
#'       "Simpson", "InvSimpson", "Fisher",
#'       "Evenness", "TaxaNumber"),
#'     mindepth = NULL,
#'     force = FALSE
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
#'   indices = c("Shannon", "Observed"))
#'
#' # absolute abundance
#' data("caporaso")
#' get_alphaindex(
#'   ps = caporaso,
#'   level = "Genus",
#'   indices = c("Shannon", "Chao1"))
#' }
#'
#'
get_alphaindex <- function(
    ps,
    level = NULL,
    indices = c("Observed", "Chao1", "ACE", "Shannon",
                "Simpson", "InvSimpson", "Fisher",
                "Evenness", "TaxaNumber"),
    mindepth = NULL,
    force = FALSE) {

  # ps = caporaso
  # level = "Genus"
  # indices = "Shannon"
  # mindepth = NULL
  # force = FALSE

  stopifnot(inherits(ps, "phyloseq"))
  ps <- preprocess_ps(ps)

  # taxa level
  if (!is.null(level)) {
    ps <- aggregate_taxa(x = ps, level = level)
  } else {
    ps <- ps
  }

  # rarefy
  sample_LibrarySize <- colSums(ps@otu_table %>% data.frame())
  if (missing(mindepth) || is.null(mindepth)){
    mindepth <- min(sample_LibrarySize)
  }
  if (sample_LibrarySize %>% var != 0 && force){
    ps <- norm_rarefy(object = ps, size = mindepth)
  }

  # alpha diversity indices
  measures <- c("Observed", "Chao1", "ACE", "Shannon",
                "Simpson", "InvSimpson", "Fisher",
                "Evenness", "TaxaNumber")
  if (all(length(indices) == 1, indices == "all")) {
    indices <- as.character(measures)
  }
  if (any(indices %in% measures)) {
    indices[indices %in% measures] <- measures[measures %in% indices]
  }
  if (!any(indices %in% measures)) {
    stop("None of the `indices` you provided are supported. Try default `all` instead.")
  }

  # otu table
  otu_tab <- as(phyloseq::otu_table(ps), "matrix") %>%
    t()

  # shannon Simpson & Evenness (Pielou) (Row sample * Col feature)
  Shannon <- vegan::diversity(otu_tab)
  Simpson <- vegan::diversity(otu_tab, index = "simpson")
  InvSimpson <- vegan::diversity(otu_tab, index = "invsimpson")
  Evenness <- Shannon/log(vegan::specnumber(otu_tab))
  TaxaNumber <- vegan::specnumber(otu_tab)

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
                        Evenness = Evenness,
                        TaxaNumber = TaxaNumber)
  } else {
    #message("Chao1, ACE and Fisher could not be calculated without integer values.")
    Observed <- apply(otu_tab, 1, function(x) {sum(x > 0)})
    # Fisher <- NULL
    # Chao1 <- NULL
    # ACE <- NULL
    alpha <- data.frame(Observed = Observed,
                        Shannon = Shannon,
                        Simpson = Simpson,
                        InvSimpson = InvSimpson,
                        Evenness = Evenness,
                        TaxaNumber = TaxaNumber)
  }

  # overlap of index
  index_overlap <- dplyr::intersect(colnames(alpha), indices)

  res <- phyloseq::sample_data(ps) %>%
    data.frame() %>%
    tibble::rownames_to_column("TempRowNames") %>%
    dplyr::inner_join(alpha %>%
                        dplyr::select(all_of(index_overlap)) %>%
                        tibble::rownames_to_column("TempRowNames"),
                      by = "TempRowNames") %>%
    tibble::column_to_rownames("TempRowNames")

  return(res)
}
