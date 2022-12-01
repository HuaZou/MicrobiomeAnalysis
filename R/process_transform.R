#' @title Transform the taxa abundances in `otu_table` sample by sample
#'
#' @description
#' Transform the taxa abundances in `otu_table` sample by sample, which means
#' the counts of each sample will be transformed individually.
#'
#' @author Created by Yang Cao; modified by Hua Zou (11/30/2022 Shenzhen China)
#'
#' @param object  [`otu_table-class`], [`phyloseq-class`], or
#'   [`microbiomeMarker-class`].
#' @param transform transformation to apply, the options inclulde:
#' * "identity", return the original data without any transformation.
#' * "log10", the transformation is `log10(object)`, and if the data contains
#'   zeros the transformation is `log10(1 + object)`.
#' * "log10p", the transformation is `log10(1 + object)`.
#' @param level (Optional). character. Summarization
#' level (from \code{rank_names(pseq)}, default: NULL).
#'
#' @importFrom phyloseq t otu_table<-
#'
#' @usage transform_abundances(
#'    object,
#'    transform = c("identity", "log10", "log10p"),
#'    level = c(NULL, "Kingdom", "Phylum", "Class",
#'            "Order", "Family", "Genus",
#'            "Species", "Strain", "unique"))
#'
#' @return A object matches the class of argument `object` with the transformed
#'   `otu_table`.
#'
#' @export
#'
#' @seealso [`abundances()`]
#'
#' @examples
#' \dontrun{
#' data(enterotypes_arumugam)
#' x1 <- transform_abundances(object = enterotypes_arumugam)
#' head(otu_table(x1), 10)
#'
#' x2 <- transform_abundances(enterotypes_arumugam, transform = "log10")
#' head(otu_table(x2), 10)
#'
#' x3 <- transform_abundances(enterotypes_arumugam, transform = "log10p")
#' head(otu_table(x3), 10)
#' }
#'
transform_abundances <- function(
    object,
    transform = c("identity", "log10", "log10p"),
    level = NULL) {

  transform <- match.arg(transform, c("identity", "log10", "log10p"))

  if (any(inherits(object, "environment"), inherits(object, "phyloseq"))) {

    # taxa level
    if (!is.null(level)) {
      object <- aggregate_taxa(x = object, level = level)
    } else {
      object <- object
    }

    otu <- as(otu_table(object), "matrix")
  } else {
    otu <- as.matrix(object)
  }

  if (transform == "identity") {
    abd <- otu
  } else if (transform == "log10") {
    abd <- transform_log10(otu)
  } else {
    abd <- transform_log10p(otu)
  }

  if (any(inherits(object, "environment"), inherits(object, "phyloseq"))) {
    otu_table(object) <- otu_table(abd, taxa_are_rows = taxa_are_rows(object))
  } else {
    object <- abd
  }

  return(object)
}

# the data is transformed using log10(1 + x) if the data contains zeroes
transform_log10 <- function(x) {

  x_norm <- x
  if (min(x) == 0) {
    warning("OTU table contains zeroes. Using log10(1 + x) instead.")
    for(i in 1:nrow(x_norm)) {
      for (j in 1:ncol(x_norm)) {
        if (x_norm[i, j] == 0) {
          x_norm[i, j] <- log10(x_norm[i, j] + 1)
        } else {
          x_norm[i, j] <- log10(x_norm[i, j])
        }
      }
    }
  } else {
    x_norm <- log10(x)
  }

  return(x_norm)
}

# the data is transformed using log10(1 + x)
transform_log10p <- function(x) {

  x_norm <- x
  for(i in 1:nrow(x_norm)) {
    for (j in 1:ncol(x_norm)) {
      if (x_norm[i, j] == 0) {
        x_norm[i, j] <- log10(x_norm[i, j] + 1)
      } else {
        x_norm[i, j] <- log10(x_norm[i, j])
      }
    }
  }
  return(x_norm)
}