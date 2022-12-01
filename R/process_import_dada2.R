# This function is modified from import-dada2.R in microbiomeMarker
# https://github.com/yiluheihei/microbiomeMarker/tree/master/R/import-dada2.R

#' @title Convert the output of dada2 into phyloseq object
#'
#' @description
#' Convert the output of dada2 into phyloseq object.
#'
#' @details
#' The output of the dada2 pipeline is a feature table of amplicon sequence
#' variants (an ASV table): A matrix with rows corresponding to samples and
#' columns to ASVs, in which the value of each entry is the number of times
#' that ASV was observed in that sample. This table is analogous to the
#' traditional OTU table. Conveniently, taxa names are saved as ASV1, ASV2,
#' ..., in the returned phyloseq object.
#'
#' @author Created by Yang Cao; modified by Hua Zou (5/14/2022 Shenzhen China)
#'
#' @param seq_tab (Required). matrix-like, ASV table, the output of
#' `dada2::removeBimeraDenovo`.
#' @param tax_tab (Optional). matrix, taxonomy table, the output of
#' `dada2::assignTaxonomy` or `dada2::addSpecies`.
#' @param sam_tab (Optional). data.frame or [`phyloseq::sample_data-class`],
#' sample data.
#' @param phy_tree (Optional). [`ape::phylo`] class or character represents
#' the path of the tree file.
#' @param keep_taxa_rows (Optional). logical. whether keep taxa in rows or not
#' in the `otu_table` of the returned `phyloseq` object, default `TRUE`.
#'
#' @return [`phyloseq::phyloseq-class`] object hold the taxonomy info,
#'   sample metadata, number of reads per ASV.
#'
#' @importFrom phyloseq sample_data read_tree
#' @importMethodsFrom phyloseq t
#'
#' @export
#'
#' @examples
#' \dontrun{
#' seq_tab <- readRDS(
#'   system.file("extdata", "dada2_seqtab.rds",
#'               package = "MicrobiomeAnalysis"))
#' tax_tab <- readRDS(
#'   system.file("extdata", "dada2_taxtab.rds",
#'               package = "MicrobiomeAnalysis"))
#' sam_tab <- read.table(
#'   system.file("extdata", "dada2_samdata.txt",
#'               package = "MicrobiomeAnalysis"),
#'   sep = "\t", header = TRUE, row.names = 1)
#' ps <- import_dada2(
#'    seq_tab = seq_tab,
#'    tax_tab = tax_tab,
#'    sam_tab = sam_tab)
#' ps
#' }
#'
import_dada2 <- function(
    seq_tab,
    tax_tab = NULL,
    sam_tab = NULL,
    phy_tree = NULL,
    keep_taxa_rows = TRUE) {

  # seq_tab <- readRDS(
  #   system.file("extdata", "dada2_seqtab.rds",
  #               package = "MicrobiomeAnalysis"))
  # tax_tab <- readRDS(
  #   system.file("extdata", "dada2_taxtab.rds",
  #               package = "MicrobiomeAnalysis"))
  # sam_tab <- read.table(
  #   system.file("extdata", "dada2_samdata.txt",
  #               package = "MicrobiomeAnalysis"),
  #   sep = "\t", header = TRUE, row.names = 1)
  #
  # phy_tree = NULL
  # keep_taxa_rows = TRUE

  # refseq
  refseq <- colnames(seq_tab)

  # set refseq and taxa names to ASV_1, ASV_2,...
  refseq_nm <- paste0("ASV", seq_along(refseq))
  colnames(seq_tab) <- refseq_nm
  names(refseq) <- refseq_nm

  if (!is.null(tax_tab)) {
    if (!identical(refseq_nm, row.names(tax_tab))) {
      tax_tab_temp <- .import_dada2_taxa(dada2_taxa = tax_tab)
      tax_tab <- tax_tab_temp[match(refseq,
                                    row.names(tax_tab_temp)), ,
                drop = FALSE]
    }
    row.names(tax_tab) <- refseq_nm
    tax_tab <- phyloseq::tax_table(as.matrix(tax_tab))
  }

  # refseq to XStringSet
  refseq <- Biostrings::DNAStringSet(refseq)

  if (!is.null(sam_tab)) {
    sam_tab <- phyloseq::sample_data(sam_tab)
  }

  if (!is.null(phy_tree) && inherits(phy_tree, "character")) {
    phy_tree <- phyloseq::read_tree(phy_tree)
  }

  asv_tab <- phyloseq::otu_table(seq_tab, taxa_are_rows = FALSE)

  ps <- phyloseq::phyloseq(
          otu_tab = asv_tab,
          tax_tab = tax_tab,
          sam_tab = sam_tab,
          phy_tree = phy_tree,
          refseq = refseq)

  if (keep_taxa_rows) {
    ps <- phyloseq::t(ps)
  }

  return(ps)
}

#' @title replace the NA with unclassified and give taxa the prefix
#' @keywords internal
#' @noRd
#'
.import_dada2_taxa <- function(dada2_taxa) {

  # tax_tab <- readRDS(
  #   system.file("extdata", "dada2_taxtab.rds",
  #               package = "MicrobiomeAnalysis"))
  #
  # dada2_taxa = tax_tab

  if (is.data.frame(dada2_taxa)) {
    taxa_tab <- dada2_taxa %>%
      as.matrix()
  } else {
    taxa_tab <- dada2_taxa
  }

  tax_prefix <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

  # i-> taxa; j->OTU
  res <- taxa_tab
  for (i in 1:nrow(res)) {
    for (j in 1:ncol(res)) {
      tax_name <- as.character(make.names(res[i, j]))
      if (length(which(tax_name == "NA.")) == 0) {
        if (j == 7) {
          upper_tax_name <- as.character(make.names(res[i, j-1]))
          res[i, j] <- paste0(upper_tax_name, "_", tax_name)
        } else {
          res[i, j] <- tax_name
        }
      } else if (length(which(tax_name == "NA.")) != 0) {
        upper_tax_name <- as.character(make.names(res[i, j-1]))

        if (length(grep("_unclassified", upper_tax_name)) != 0) {
          res[i, j] <- upper_tax_name
        } else {
          res[i, j] <- paste0(upper_tax_name, "_unclassified")
        }
      }

      # replace dot with "_"
      res[i, j] <- gsub("\\.", "_", res[i, j])

    }
  }

  for (K in 1:ncol(res)) {
    res[, K] <- paste0(tax_prefix[K], res[, K])
  }

  return(res)
}

