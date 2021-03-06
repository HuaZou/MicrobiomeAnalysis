% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gv_distance.R
\name{run_distance}
\alias{run_distance}
\title{Calculate the distance among samples}
\usage{
run_distance(
                 object,
                 method,
                 alpha)
}
\arguments{
\item{object}{(Required). a \code{\link[phyloseq:phyloseq-class]{phyloseq::phyloseq}} or
\code{\link[Biobase:class.ExpressionSet]{Biobase::ExpressionSet}} object.}

\item{method}{(Required). character. Provide one of the currently supported
options. See \code{distanceMethodList} for a detailed list of the supported options
and links to accompanying documentation. Options include:
\itemize{
\item "unifrac" : unweighted UniFrac distance.
\item "wunifrac": weighted-UniFrac distance.
\item "GUniFrac": The variance-adjusted weighted UniFrac distances (default: alpha=0.5).
\item "bray": bray crutis distance.
\item "dpcoa": sample-wise distance used in Double Principle Coordinate Analysis.
\item "jsd": Jensen-Shannon Divergence.
Alternatively, you can provide a character string that defines a custom
distance method, if it has the form described in \code{designdist}.
}}

\item{alpha}{(Optional). numeric. the parameter for "GUniFrac" controlling
weight on abundant lineages (default: 0.5).}
}
\value{
distance object, which could be applied for ordination analysis.
}
\description{
The function is to calculate the distance among samples
}
\examples{

\dontrun{
data("enterotypes_arumugam")
res <- run_distance(enterotypes_arumugam,
                    method = "bray")
}

}
\author{
Created by Hua Zou (5/14/2022 Shenzhen China)
}
