#' check whether all names of taxonomic ranks include in available_ranks
#'
#' @param ps a [`phyloseq::phyloseq-class`] object; (Required).
#'
#' @export
#'
check_rank_names <- function(ps) {

  if (!is.null(ps@tax_table)) {
    ps_ranks <- rank_names(ps)
    if (!all(ps_ranks %in% available_ranks)) {
      stop(
        "ranks of taxonimic profile must be one of ",
        paste(available_ranks, collapse = ", ")
      )
    }
  } else {
    message("There is no taxa table in phyloseq object")
  }

  invisible(ps)
}


#' check whether samples' names start with numeric and then paste "S_"
#'
#' @param object (Optional). a [`phyloseq::phyloseq-class`] or
#' [`Biobase::ExpressionSet`] or [`phyloseq::otu_table-class`] or
#' [`phyloseq::otu_table-class`] or matrix/metadata.
#'
#' @export
#'
check_sample_names <- function(object) {

  # phyloseq
  if (inherits(object, "phyloseq")) {
    if(length(grep("^\\d+", phyloseq::sample_names(object))) > 0) {
      message("samples' names are start with numeric")
      Names_index <- grep("^\\d+", phyloseq::sample_names(object))
      phyloseq::sample_names(object)[Names_index] <-
        paste0("S_", phyloseq::sample_names(object)[Names_index])
    }
  }

  # ExpressionSet
  if (inherits(object, "ExpressionSet")) {
    if(length(grep("^\\d+", Biobase::sampleNames(object))) > 0) {
      message("samples' names are start with numeric")
      Names_index <- grep("^\\d+", Biobase::sampleNames(object))
      Biobase::sampleNames(object)[Names_index] <-
        paste0("S_", Biobase::sampleNames(object)[Names_index])
    }
  }

  # otu table
  if (inherits(object, "otu_table")) {
    if (phyloseq::taxa_are_rows(object)) {
      if(length(grep("^\\d+", colnames(object))) > 0) {
        Names_index <- grep("^\\d+", colnames(object))
        colnames(object)[Names_index] <- paste0("S_", colnames(object)[Names_index])
      }
    } else {
      if(length(grep("^\\d+", rownames(object))) > 0) {
        Names_index <- grep("^\\d+", rownames(object))
        rownames(object)[Names_index] <- paste0("S_", rownames(object)[Names_index])
      }
    }
  }

  # matrix metadata
  if(length(grep("^\\d+", rownames(object))) > 0) {
    Names_index <- grep("^\\d+", rownames(object))
    rownames(object)[Names_index] <- paste0("S_", rownames(object)[Names_index])
  }
  if(length(grep("^\\d+", colnames(object))) > 0) {
    Names_index <- grep("^\\d+", colnames(object))
    colnames(object)[Names_index] <- paste0("S_", colnames(object)[Names_index])
  }

  return(object)
}


#' check whether the object is phyloseq-class object
#'
#' @param ps (Required). a [`phyloseq::phyloseq-class`] object.
#' @param fill_na_taxa (Optional). whether there is taxa table (default: TRUE).
#'
#' @export
#'
check_phyloseq <- function (ps, fill_na_taxa = TRUE) {

  # Sanity checks for a phyloseq object. Required with some methods.
  if (!taxa_are_rows(ps)) {
    x@otu_table <- otu_table(t(otu_table(ps)), taxa_are_rows = TRUE)
  }

  if (fill_na_taxa || is.character(fill_na_taxa)) {
    M <- as.matrix(tax_table(ps))
    if (!taxa_are_rows(ps)) {
      M <- t(M)
    }

    if (!is.character(fill_na_taxa)) {
      fill_na_taxa <- "Unknown"
    }

    M[is.na(M)] <- fill_na_taxa
    ps@tax_table <- tax_table(M)
  }

  return(ps)
}

