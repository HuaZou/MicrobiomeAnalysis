% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{marker_table}
\alias{marker_table}
\alias{marker_table,data.frame-method}
\alias{marker_table,microbiomeMarker-method}
\title{Build or access the marker_table}
\usage{
marker_table(object)

\S4method{marker_table}{data.frame}(object)

\S4method{marker_table}{microbiomeMarker}(object)
}
\arguments{
\item{object}{an object among the set of classes defined by the
microbiomeMarker package that contain \code{marker_table}}
}
\value{
a \code{\linkS4class{marker_table}} object.
}
\description{
This is the recommended function for both building and accessing microbiome
marker table (\code{\linkS4class{marker_table}}).
}
\examples{

data(enterotypes_arumugam)
mm <- run_limma_voom(
    enterotypes_arumugam,
    "Enterotype",
    contrast = c("Enterotype 3", "Enterotype 2"),
    pvalue_cutoff = 0.05,
    p_adjust = "fdr"
)
marker_table(mm)
}
