% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/da_metagenomeSeq.R
\name{phyloseq2metagenomeSeq}
\alias{phyloseq2metagenomeSeq}
\alias{otu_table2metagenomeSeq}
\title{Convert phyloseq data to MetagenomeSeq \code{MRexperiment} object}
\usage{
phyloseq2metagenomeSeq(ps, ...)

otu_table2metagenomeSeq(ps, ...)
}
\arguments{
\item{ps}{\code{\link[phyloseq:phyloseq-class]{phyloseq::phyloseq}} object for
\code{phyloseq2metagenomeSeq()}, or \code{\link[phyloseq:otu_table-class]{phyloseq::otu_table}} object
for \code{otu_table2metagenomeseq()}.}

\item{...}{optional, additional named arguments passed  to
\code{\link[metagenomeSeq:newMRexperiment]{metagenomeSeq::newMRexperiment()}}. Most users will not need to pass
any additional arguments here.}
}
\value{
A \code{\link[metagenomeSeq:MRexperiment-class]{metagenomeSeq::MRexperiment}} object.
}
\description{
The phyloseq data is converted to the relevant
\code{\link[metagenomeSeq:MRexperiment-class]{metagenomeSeq::MRexperiment}} object, which can then be tested in
the zero-inflated mixture model framework in the metagenomeSeq package.
}
\examples{
data(caporaso)
phyloseq2metagenomeSeq(caporaso)
}
\seealso{
\code{\link[metagenomeSeq:fitTimeSeries]{metagenomeSeq::fitTimeSeries()}},
\code{\link[metagenomeSeq:fitLogNormal]{metagenomeSeq::fitLogNormal()}},\code{\link[metagenomeSeq:fitZig]{metagenomeSeq::fitZig()}},
\code{\link[metagenomeSeq:MRtable]{metagenomeSeq::MRtable()}},\code{\link[metagenomeSeq:MRfulltable]{metagenomeSeq::MRfulltable()}}
}
