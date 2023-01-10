#' @title Convert the output of dada2 into phyloseq object
#'
#' @description
#' Convert inputs into SummarizedExperiment object.
#'
#' @details
#' The profile of inputs is a feature table of RNA-seq or metabolites or others:
#' A matrix with rows corresponding to features and
#' columns to samples, in which the value of each entry is the number of times
#' that features was observed in that sample.
#'
#' @author Created by Hua Zou (1/8/2023 Shenzhen China)
#'
#' @param object (Required). matrix-like, features table with expressed values
#' (row: featureID; columns: samples).
#' @param rowdata (Optional). DataFrame,
#' A DataFrame object describing the rows (default: NULL).
#' @param rowranges (Optional).
#' A GRanges or GRangesList object describing the ranges of interest. (default: NULL).
#' @param coldata (Optional). DataFrame,
#' An optional DataFrame describing the samples (default: NULL).
#' @param metadata (Optional). List, An optional list of arbitrary content
#' describing the overall experiment (default: NULL).
#'
#' @return [`SummarizedExperiment::SummarizedExperiment-class`] object hold expressed profile,
#' information related row and column or experimental design.
#'
#' @import SummarizedExperiment
#'
#' @usage import_SE(
#'     object,
#'     rowdata = NULL,
#'     rowranges = GRangesList(),
#'     coldata = NULL,
#'     metadata = list())
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("Zeybel_2022_protein")
#' assay <- SummarizedExperiment::assay(Zeybel_2022_protein) %>%
#'   data.frame()
#' rowData <- SummarizedExperiment::rowData(Zeybel_2022_protein) %>%
#'   data.frame()
#' colData <- SummarizedExperiment::colData(Zeybel_2022_protein) %>%
#'   data.frame()
#' metadata <- list(lab="hua", type="protein")
#'
#' assay <- assay[1:10, 1:10]
#' import_SE(
#'     object = assay,
#'     rowdata = rowData,
#'     coldata = colData,
#'     metadata = metadata)
#' }
#'
import_SE <- function(
    object,
    rowdata = NULL,
    rowranges = NULL,
    coldata = NULL,
    metadata = list()) {

  # data(Zeybel_2022_protein)
  # object = SummarizedExperiment::assay(Zeybel_2022_protein) %>%
  #   data.frame()
  # rowdata = SummarizedExperiment::rowData(Zeybel_2022_protein) %>%
  #   data.frame()
  # rowranges = GenomicRanges::GRangesList()
  # coldata = SummarizedExperiment::colData(Zeybel_2022_protein) %>%
  #   data.frame()
  # metadata = list(lab="hua")


  # overlap of samples
  if (all(!is.null(object), !is.null(coldata))) {
    overlap_sample <- intersect(colnames(object), rownames(coldata))
    if (length(overlap_sample) == 0) {
      stop("No overlap of samples between assay and colData, please check your data")
    }
    object <- object %>%
      dplyr::select(dplyr::all_of(overlap_sample))
    #coldata <- coldata[rownames(coldata) %in% overlap_sample, , F]
    coldata <- coldata[colnames(object), , F]
  }

  # overlap of features
  if (all(!is.null(object), !is.null(rowdata))) {
    overlap_feature <- intersect(rownames(object), rownames(rowdata))
    if (length(overlap_feature) == 0) {
      stop("No overlap of features between assay and rowData, please check your data")
    }
    object <- object[rownames(object) %in% overlap_feature, , F]
    #rowdata <- rowdata[rownames(rowdata) %in% overlap_feature, , F]
    rowdata <- rowdata[rownames(object), , F]
  }

  # overlap of features
  if (all(!is.null(object), length(rowranges) != 0)) {
    overlap_feature <- intersect(rownames(object), rownames(rowranges))
    if (length(overlap_feature) == 0) {
      stop("No overlap of features between assay and rowRanges, please check your data")
    }
    object <- object[rownames(object), , F]
    #rowranges <- rowranges[rownames(rowranges) %in% overlap_feature, , F]
    rowranges <- rowranges[pmatch(rownames(rowranges), rownames(object)), , F]
  }

  if (all(!is.null(rowdata), length(rowranges) != 0)) {
    stop("Only on of rowData and rowRanges to be used as input, please check your data")
  }

  # SummarizedExperiment object
  if (!is.null(rowdata)) {
    res <- SummarizedExperiment::SummarizedExperiment(
      assays = object,
      rowData = rowdata,
      colData = coldata,
      metadata = metadata)
  } else if (length(rowranges) != 0) {
    res <- SummarizedExperiment::SummarizedExperiment(
      assays = object,
      rowRanges = rowranges,
      colData = coldata,
      metadata = metadata)
  } else {
    res <- SummarizedExperiment::SummarizedExperiment(
      assays = object,
      #rowData = rowdata,
      #rowRanges = rowranges,
      colData = coldata,
      metadata = metadata)
  }

  return(res)
}
