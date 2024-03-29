# This function is inspired from microbiome::abundances

#' Extract taxa abundances
#'
#' Extract taxa abundances from phyloseq objects.
#'
#' @param object [`otu_table-class`], [`phyloseq-class`], or
#'   [`microbiomeMarker-class`].
#' @param transform transformation to apply, the options inclulde:
#' * "identity", return the original data without any transformation.
#' * "log10", the transformation is `log10(object)`, and if the data contains
#'   zeros the transformation is `log10(1 + object)`.
#' * "log10p", the transformation is `log10(1 + object)`.
#' * "log2", the transformation is `log2(object)`, and if the data contains
#'   zeros the transformation is `log2(1 + object)`.
#' * "log2p", the transformation is `log2(1 + object)`.
#' * "SquareRoot", the transformation is `Square Root`.
#' * "CubicRoot", the transformation is `Cubic Root`.
#' * "logit", the transformation is `Zero-inflated Logit Transformation`
#' (Does not work well for microbiome data).
#' @param norm logical, indicating whether or not to return the normalized
#'   taxa abundances.
#' @return abundance matrix with taxa in rows and samples in columns.
#' @seealso [`otu_table-class`], [`phyloseq-class`],
#' [`microbiomeMarker-class`],[`transform_abundances`]
#' @docType methods
#' @rdname abundances-methods
#' @export
#' @examples
#' data(caporaso)
#' abd <- abundances(caporaso)
setGeneric(
    "abundances",
    function(object,
    transform = c("identity", "log10", "log10p",
                  "log2", "log2p",
                  "SquareRoot", "CubicRoot", "logit"),
    norm = FALSE) {
        standardGeneric("abundances")
    }
)

# otu_table object
#' @aliases abundances, otu_table-method
#' @rdname abundances-methods
setMethod(
    abundances, "otu_table",
    function(object,
        transform = c("identity", "log10", "log10p",
                      "log2", "log2p",
                      "SquareRoot", "CubicRoot", "logit"),
        norm = FALSE) {

        transform <- match.arg(transform, c("identity", "log10", "log10p",
                                            "log2", "log2p",
                                            "SquareRoot", "CubicRoot", "logit"))
        obj_transed <- transform_abundances(object, transform = transform)
        abd <- as(otu_table(obj_transed), "matrix")

        if (norm) {
            nf <- get_norm_factors(object)
            if (!is.null(nf)) {
                abd <- sweep(abd, 2, nf, "/")
            }
        }

        return(abd)
    }
)

# phyloseq object
#' @aliases abundances,phyloseq-method
#' @rdname abundances-methods
setMethod(
    abundances, "phyloseq",
    function(object,
        transform = c("identity", "log10", "log10p",
                      "log2", "log2p",
                      "SquareRoot", "CubicRoot", "logit"),
        norm = FALSE) {

        transform <- match.arg(transform, c("identity", "log10", "log10p",
                                            "log2", "log2p",
                                            "SquareRoot", "CubicRoot", "logit"))
        otu <- otu_table(object)
        if (norm) {
            nf <- get_norm_factors(object)
            if (!is.null(nf)) {
                attr(otu, "norm_factor") <- nf
            }
        }
        otu <- abundances(otu, transform = transform, norm = norm)

        return(otu)
    }
)

# microbiomeMarker object
#' @aliases abundances,microbiomeMarker-method
#' @rdname abundances-methods
setMethod(
    abundances, "microbiomeMarker",
    function(object,
        transform = c("identity", "log10", "log10p",
                      "log2", "log2p",
                      "SquareRoot", "CubicRoot", "logit")) {

        transform <- match.arg(transform, c("identity", "log10", "log10p",
                                            "log2", "log2p",
                                            "SquareRoot", "CubicRoot", "logit"))
        otu <- otu_table(object)
        otu <- abundances(otu, transform = transform, norm = FALSE)

        return(otu)
    }
)

# assay object
#' @aliases abundances, assay-method
#' @rdname abundances-methods
setMethod(
  abundances, "assay",
  function(object,
           transform = c("identity", "log10", "log10p",
                         "log2", "log2p",
                         "SquareRoot", "CubicRoot", "logit"),
           norm = FALSE) {

    transform <- match.arg(transform, c("identity", "log10", "log10p",
                                        "log2", "log2p",
                                        "SquareRoot", "CubicRoot", "logit"))
    obj_transed <- transform_abundances(object, transform = transform)
    abd <- as(SummarizedExperiment::assay(obj_transed), "matrix")

    if (norm) {
      nf <- get_norm_factors(object)
      if (!is.null(nf)) {
        abd <- sweep(abd, 2, nf, "/")
      }
    }

    return(abd)
  }
)

# SummarizedExperiment object
#' @aliases abundances,SummarizedExperiment-method
#' @rdname abundances-methods
setMethod(
  abundances, "SummarizedExperiment",
  function(object,
           transform = c("identity", "log10", "log10p",
                         "log2", "log2p",
                         "SquareRoot", "CubicRoot", "logit"),
           norm = FALSE) {

    transform <- match.arg(transform, c("identity", "log10", "log10p",
                                        "log2", "log2p",
                                        "SquareRoot", "CubicRoot", "logit"))
    otu <- SummarizedExperiment::assay(object)
    if (norm) {
      nf <- get_norm_factors(object)
      if (!is.null(nf)) {
        attr(otu, "norm_factor") <- nf
      }
    }
    otu <- abundances(otu, transform = transform, norm = norm)

    return(otu)
  }
)
