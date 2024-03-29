% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/da_multiple_groups.R
\name{run_test_multiple_groups}
\alias{run_test_multiple_groups}
\title{Statistical test for multiple groups}
\usage{
run_test_multiple_groups(
  ps,
  group,
  taxa_rank = "all",
  transform = c("identity", "log10", "log10p", "SquareRoot", "CubicRoot", "logit"),
  norm = "TSS",
  norm_para = list(),
  method = c("anova", "kruskal"),
  p_adjust = c("none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY"),
  pvalue_cutoff = 0.05,
  effect_size_cutoff = NULL
)
}
\arguments{
\item{ps}{a \code{\link[phyloseq:phyloseq-class]{phyloseq::phyloseq}} object}

\item{group}{character, the variable to set the group}

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

\item{norm_para}{arguments passed to specific normalization methods}

\item{method}{test method, must be one of "anova" or "kruskal"}

\item{p_adjust}{method for multiple test correction, default \code{none},
for more details see \link[stats:p.adjust]{stats::p.adjust}.}

\item{pvalue_cutoff}{numeric, p value cutoff, default 0.05.}

\item{effect_size_cutoff}{numeric, cutoff of effect size default \code{NULL}
which means no effect size filter. The eta squared is used to measure the
effect size for anova/kruskal test.}
}
\value{
a \code{\linkS4class{microbiomeMarker}} object.
}
\description{
Differential expression analysis for multiple groups.
}
\examples{
data(enterotypes_arumugam)
ps <- phyloseq::subset_samples(
    enterotypes_arumugam,
    Enterotype \%in\% c("Enterotype 3", "Enterotype 2", "Enterotype 1")
)
mm_anova <- run_test_multiple_groups(
    ps,
    group = "Enterotype",
    method = "anova"
)
}
\seealso{
\code{\link[=run_posthoc_test]{run_posthoc_test()}},\code{\link[=run_test_two_groups]{run_test_two_groups()}},\code{\link[=run_simple_stat]{run_simple_stat()}}
}
