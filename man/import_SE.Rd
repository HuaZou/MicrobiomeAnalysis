% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_import_SE.R
\name{import_SE}
\alias{import_SE}
\title{Convert inputs into SummarizedExperiment object}
\usage{
import_SE(
    object,
    rowdata = NULL,
    rowranges = GRangesList(),
    coldata = NULL,
    metadata = list())
}
\arguments{
\item{object}{(Required). matrix-like, features table with expressed values
(row: featureID; columns: samples).}

\item{rowdata}{(Optional). DataFrame,
A DataFrame object describing the rows (default: NULL).}

\item{rowranges}{(Optional).
A GRanges or GRangesList object describing the ranges of interest. (default: NULL).}

\item{coldata}{(Optional). DataFrame,
An optional DataFrame describing the samples (default: NULL).}

\item{metadata}{(Optional). List, An optional list of arbitrary content
describing the overall experiment (default: NULL).}
}
\value{
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment::SummarizedExperiment}} object hold expressed profile,
information related row and column or experimental design.
}
\description{
Convert inputs into SummarizedExperiment object.
}
\details{
The profile of inputs is a feature table of RNA-seq or metabolites or others:
A matrix with rows corresponding to features and
columns to samples, in which the value of each entry is the number of times
that features was observed in that sample.
}
\examples{
\dontrun{
data("Zeybel_2022_protein")
assay <- SummarizedExperiment::assay(Zeybel_2022_protein) \%>\%
  data.frame()
rowData <- SummarizedExperiment::rowData(Zeybel_2022_protein) \%>\%
  data.frame()
colData <- SummarizedExperiment::colData(Zeybel_2022_protein) \%>\%
  data.frame()
metadata <- list(lab="hua", type="protein")

assay <- assay[1:10, 1:10]
import_SE(
    object = assay,
    rowdata = rowData,
    coldata = colData,
    metadata = metadata)
}

}
\author{
Created by Hua Zou (1/8/2023 Shenzhen China)
}
