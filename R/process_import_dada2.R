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
#' @usage
#' import_dada2(
#'     seq_tab,
#'     tax_tab = NULL,
#'     sam_tab = NULL,
#'     phy_tree = NULL,
#'     keep_taxa_rows = TRUE)
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
#' ps <- import_dada2(seq_tab = seq_tab, tax_tab = tax_tab, sam_tab = sam_tab)
#' ps
#' }
#'
import_dada2 <- function(
    seq_tab,
    tax_tab = NULL,
    sam_tab = NULL,
    phy_tree = NULL,
    keep_taxa_rows = TRUE) {

    # refseq
    refseq <- colnames(seq_tab)
    # set refseq and taxa names to ASV_1, ASV_2,...
    refseq_nm <- paste0("ASV", seq_along(refseq))
    colnames(seq_tab) <- refseq_nm
    names(refseq) <- refseq_nm

    if (!is.null(tax_tab)) {
        if (!identical(refseq_nm, row.names(tax_tab))) {
            tax_tab <- tax_tab[match(refseq, row.names(tax_tab)), ,
                drop = FALSE]
        }
        row.names(tax_tab) <- refseq_nm
        tax_tab <- tax_table(as.matrix(tax_tab))
    }

    # refseq to XStringSet
    refseq <- Biostrings::DNAStringSet(refseq)

    if (!is.null(sam_tab)) {
        sam_tab <- sample_data(sam_tab)
    }

    if (!is.null(phy_tree) && inherits(phy_tree, "character")) {
        phy_tree <- read_tree(phy_tree)
    }

    asv_tab <- otu_table(seq_tab, taxa_are_rows = FALSE)
    ps <- phyloseq(asv_tab, tax_tab, sam_tab, phy_tree, refseq)

    if (keep_taxa_rows) {
        ps <- phyloseq::t(ps)
    }

    return(ps)
}


