% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/da_lefse.R
\name{run_lefse}
\alias{run_lefse}
\title{Liner discriminant analysis (LDA) effect size (LEFSe) analysis}
\usage{
run_lefse(
  ps,
  group,
  subgroup = NULL,
  taxa_rank = "all",
  transform = c("identity", "log10", "log10p", "SquareRoot", "CubicRoot", "logit"),
  norm = "CPM",
  norm_para = list(),
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = c("0", "1", "2"),
  sample_min = 10,
  only_same_subgrp = FALSE,
  curv = FALSE
)
}
\arguments{
\item{ps}{a \code{\link[phyloseq]{phyloseq-class}} object}

\item{group}{character, the column name to set the group}

\item{subgroup}{character, the column name to set the subgroup}

\item{taxa_rank}{character to specify taxonomic rank to perform
differential analysis on. Should be one of
\code{phyloseq::rank_names(phyloseq)}, or "all" means to summarize the taxa by
the top taxa ranks (\code{summarize_taxa(ps, level = rank_names(ps)[1])}), or
"none" means perform differential analysis on the original taxa
(\code{taxa_names(phyloseq)}, e.g., OTU or ASV).}

\item{transform}{character, the methods used to transform the microbial
abundance. See \code{\link[=transform_abundances]{transform_abundances()}} for more details. The
options include:
\itemize{
\item "identity", return the original data without any transformation
(default).
\item "log10", the transformation is \code{log10(object)}, and if the data contains
zeros the transformation is \code{log10(1 + object)}.
\item "log10p", the transformation is \code{log10(1 + object)}.
\item "SquareRoot", the transformation is \verb{Square Root}.
\item "CubicRoot", the transformation is \verb{Cubic Root}.
\item "logit", the transformation is \verb{Zero-inflated Logit Transformation}
(Does not work well for microbiome data).
}}

\item{norm}{the methods used to normalize the microbial abundance data. See
\code{\link[=normalize]{normalize()}} for more details.
Options include:
\itemize{
\item "none": do not normalize.
\item "rarefy": random subsampling counts to the smallest library size in the
data set.
\item "TSS": total sum scaling, also referred to as "relative abundance", the
abundances were normalized by dividing the corresponding sample library
size.
\item "TMM": trimmed mean of m-values. First, a sample is chosen as reference.
The scaling factor is then derived using a weighted trimmed mean over the
differences of the log-transformed gene-count fold-change between the
sample and the reference.
\item "RLE", relative log expression, RLE uses a pseudo-reference calculated
using the geometric mean of the gene-specific abundances over all
samples. The scaling factors are then calculated as the median of the
gene counts ratios between the samples and the reference.
\item "CSS": cumulative sum scaling, calculates scaling factors as the
cumulative sum of gene abundances up to a data-derived threshold.
\item "CLR": centered log-ratio normalization.
\item "CPM": pre-sample normalization of the sum of the values to 1e+06.
}}

\item{norm_para}{named \code{list}. other arguments passed to specific
normalization methods.  Most users will not need to pass any additional
arguments here.}

\item{kw_cutoff}{numeric, p value cutoff of kw test, default 0.05}

\item{lda_cutoff}{numeric, lda score cutoff, default 2}

\item{bootstrap_n}{integer, the number of bootstrap iteration for LDA,
default 30}

\item{bootstrap_fraction}{numeric, the subsampling fraction value for each
bootstrap iteration, default \code{2/3}}

\item{wilcoxon_cutoff}{numeric, p value cutoff of wilcoxon test, default 0.05}

\item{multigrp_strat}{logical, for multiple group tasks, whether the test is
performed in a one-against one (more strict) or in a one-against all
setting, default \code{FALSE}.}

\item{strict}{multiple testing options, 0 for no correction (default), 1 for
independent comparisons, 2 for independent comparison.}

\item{sample_min}{integer, minimum number of samples per subclass for
performing wilcoxon test, default 10}

\item{only_same_subgrp}{logical, whether perform the wilcoxon test only
among the subgroups with the same name, default \code{FALSE}}

\item{curv}{logical, whether perform the wilcoxon test using the
Curtis's approach, defalt \code{FALSE}}
}
\value{
a \linkS4class{microbiomeMarker} object, in which the \code{slot} of
\code{marker_table}
contains four variables:
\itemize{
\item \code{feature}, significantly different features.
\item \code{enrich_group}, the class of the differential features enriched.
\item \code{lda}, logarithmic LDA score (effect size)
\item \code{pvalue}, p value of kw test.
}
}
\description{
Perform Metagenomic LEFSe analysis based on phyloseq object.
}
\examples{
data(Zeybel_2022_gut)
Zeybel_2022_gut_small <- phyloseq::subset_taxa(
    Zeybel_2022_gut,
    Phylum == "p__Firmicutes"
)
mm_lefse <- run_lefse(
    Zeybel_2022_gut_small,
    wilcoxon_cutoff = 0.01,
    group = "LiverFatClass",
    kw_cutoff = 0.01,
    multigrp_strat = TRUE,
    lda_cutoff = 4
)
}
\references{
Segata, Nicola, et al. Metagenomic biomarker discovery and
explanation. Genome biology 12.6 (2011): R60.
}
\seealso{
\link{normalize}
}
\author{
Yang Cao
}
