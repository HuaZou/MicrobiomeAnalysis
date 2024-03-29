library(MicrobiomeAnalystR)
library(phyloseq)
library(magrittr)

# build package
# roxygen2::roxygenise(getwd())
# devtools::check(document = FALSE)
# devtools::build(binary = FALSE, manual = TRUE, quiet = FALSE)

# build pkgdown
# pkgdown::build_site() # Run to build the website

# Human Moving Picture from MicrobiomeAnalyst server ----------------------
download.file(
  "https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/data/treebiom.zip",
  "data-raw/caporaso.zip"
)
unzip("data-raw/caporaso.zip", exdir = "data-raw/")
file.rename("data-raw/treebiom/", "data-raw/caporaso/")

ps <- import_biom(
  "data-raw/caporaso/otu_table_mc2_w_tax_no_pynast_failures.biom",
  treefilename = "data-raw/caporaso/rep_set.tre",
)

colnames(tax_table(ps)) <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
)

sampledata <- read.delim("data-raw/caporaso/map.txt", row.names = 1) %>%
  sample_data()
caporaso <- merge_phyloseq(ps, sampledata)

usethis::use_data(caporaso, overwrite = TRUE)
unlink("data-raw/cap*", recursive = TRUE)


# cid data from github.com/ying14/yingtools2 (times series data)-----------
download.file(
  "https://github.com/ying14/yingtools2/raw/master/data/cid.phy.rda",
  "data-raw/cid.phy.rda"
)
load("data-raw/cid.phy.rda")
cid_ying <- cid.phy
tax_table(cid_ying) <- tax_table(cid.phy)[, -7]
usethis::use_data(cid_ying, overwrite = TRUE)
unlink("data-raw/cid*")

# pediatric ibd (cross sectional data)-------------------------------------

# https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/data/ibd_data.zip
download.file(
  "https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/data/ibd_data.zip",
  "data-raw/pediatric_idb.zip"
)
unzip("data-raw/pediatric_idb.zip", exdir = "data-raw/")
asv_abundance <- readr::read_tsv("data-raw/ibd_data/IBD_data/ibd_asv_table.txt") %>%
  tibble::column_to_rownames("#NAME")
asv_table <- readr::read_tsv("data-raw/ibd_data/IBD_data/ibd_taxa.txt") %>%
  tibble::column_to_rownames("#TAXONOMY")
sample_table <- readr::read_csv("data-raw/ibd_data/IBD_data/ibd_meta.csv") %>%
  tibble::column_to_rownames("#NAME")
pediatric_ibd <- phyloseq(
  otu_table(asv_abundance, taxa_are_rows = TRUE),
  tax_table(as.matrix(asv_table)),
  sample_data(sample_table)
)
tree <- read_tree(treefile = "data-raw/ibd_data/IBD_data/ibd_tree.tre")
phy_tree(pediatric_ibd) <- tree
# usethis::use_data(pediatric_ibd, overwrite = TRUE)
# unlink("data-raw/ibd_data", recursive = TRUE)
# unlink("data-raw/pediatric_idb.zip")


# oxygen availability -----------------------------------------------------
# a small subset of the HMP 16S dataset for finding biomarkers characterizing
# different level of oxygen availability in different bodysites

# oxygen_dat <- readr::read_tsv(
#   "https://raw.githubusercontent.com/biobakery/biobakery/master/demos/biobakery_demos/data/lefse/input/hmp_small_aerobiosis.txt",
#   col_names = FALSE
# )

# sample_meta <- data.frame(
#   oxygen_availability	= c(oxygen_dat[1, ][-1], recursive = TRUE),
#   body_site = c(oxygen_dat[2, ][-1], recursive = TRUE),
#   subject_id = c(oxygen_dat[3, ][-1], recursive = TRUE)
# )
#   # tibble::rownames_to_column() %>%
#   # tidyr::pivot_longer(-rowname) %>%
#   #tidyr::pivot_wider(names_from = "rowname", values_from = "value") %>%
#   #tibble::column_to_rownames("name")
# tax_dat <- oxygen_dat$X1[-(1:3)]
#
# sample_abd <- dplyr::slice(oxygen_dat, -(1:3)) %>%
#   dplyr::select(-1) %>%
#   purrr::map_df(as.numeric)
# row.names(sample_abd) <- tax_dat
#
# tax_mat <- as.matrix(tax_dat)
# row.names(tax_mat) <- tax_dat
# colnames(tax_mat) <-  "Summarize"
#
# oxygen <- phyloseq(
#   otu_table(sample_abd, taxa_are_rows = TRUE),
#   tax_table(tax_mat),
#   sample_data(sample_meta)
# )
#
download.file("https://raw.githubusercontent.com/biobakery/biobakery/master/demos/biobakery_demos/data/lefse/input/hmp_small_aerobiosis.txt", "data-raw/oxygen.txt")
oxygen <- import_biobakery_lefse_in(
  "data-raw/oxygen.txt",
  ranks_prefix = c("k", "p", "c", "o", "f", "g"),
  meta_rows = 1:3,
)
# unlink("data-raw/oxygen.txt")
# usethis::use_data(oxygen, overwrite = TRUE)

# data from lefse galaxy --------------------------------------------------
# Fecal microbiota in a mouse model of spontaneous colitis. The dataset contains
# 30 abundance profiles (obtained processing the 16S reads with RDP) belonging
# to 10 rag2 (control) and 20 truc (case) mice
# spontaneous_colitis <- readr::read_tsv(
#   "https://raw.githubusercontent.com/biobakery/galaxy_lefse/master/test-data/lefse_input",
#   col_names = FALSE
# )
# class <- spontaneous_colitis[1, ]
# taxas <- spontaneous_colitis[, 1]
#
# sample_meta <- data.frame(
#   class = unlist(class[-1]),
#   stringsAsFactors = FALSE
# )
# tax_dat <- as.matrix(taxas[-1, ])
# row.names(tax_dat) <- tax_dat
# colnames(tax_dat) <- "Summarize"
# tax_abd <- spontaneous_colitis[-1, -1] %>%
#   purrr::map_df(as.numeric)
# row.names(tax_abd) <- tax_dat[,1]
#
# spontaneous_colitis <- phyloseq(
#   otu_table(tax_abd, taxa_are_rows = TRUE),
#   tax_table(tax_dat),
#   sample_data(sample_meta)
# )
download.file("https://raw.githubusercontent.com/biobakery/galaxy_lefse/master/test-data/lefse_input", "data-raw/lefse_in")
spontaneous_colitis <- import_biobakery_lefse_in(
  "data-raw/lefse_in",
  ranks_prefix = c("k", "p", "c", "o", "f", "g"),
  meta_rows = 1,
)
unlink("data-raw/lefse_in")
# usethis::use_data(spontaneous_colitis, overwrite = TRUE)

# Enterotypes data from Arumugam's paper from stamp -----------------------

enterotypes_arumugam <- readr::read_tsv("https://github.com/yiluheihei/STAMP/raw/master/examples/EnterotypesArumugam/Enterotypes.profile.spf")

enterotypes_arumugam_meta <- readr::read_tsv("https://github.com/yiluheihei/STAMP/raw/master/examples/EnterotypesArumugam/Enterotypes.metadata.tsv") %>% as.data.frame()
row.names(enterotypes_arumugam_meta) <- enterotypes_arumugam_meta$`Sample Id`

enterotype_abd <- dplyr::select(enterotypes_arumugam, -Phyla, -Genera)
enterotype_tax <- dplyr::select(enterotypes_arumugam, Phylum = Phyla, Genus = Genera)

enterotypes_arumugam <- phyloseq(
  otu_table(enterotype_abd, taxa_are_rows = TRUE),
  tax_table(as.matrix(enterotype_tax)),
  sample_data(enterotypes_arumugam_meta)
)

usethis::use_data(enterotypes_arumugam, overwrite = TRUE)


# kostic crc --------------------------------------------------------------
# data from https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html

# A publicly available data from a study on colorectal cancer:
# Genomic analysis identifies association of Fusobacterium with colorectal
# carcinoma. Kostic, A. D., Gevers, D., Pedamallu, C. S., Michaud, M., Duke,
# ., Earl, A. M., et al. (2012). Genome research, 22(2), 292-298.
#
filepath <- system.file(
  "extdata",
  "study_1457_split_library_seqs_and_mapping.zip",
  package = "phyloseq"
)
kostic <- phyloseq::microbio_me_qiime(filepath)

# remove the 5 samples that had no DIAGNOSIS attribute assigned
kostic <- subset_samples(kostic, DIAGNOSIS != "None")
# remove samples with less than 500 reads (counts)
kostic_crc <- prune_samples(sample_sums(kostic) > 500, kostic)
usethis::use_data(kostic_crc, overwrite = TRUE)


# ecam from ancom paper ---------------------------------------------------
ecam_meta <- readr::read_tsv("https://raw.githubusercontent.com/FrederickHuangLin/ANCOM/master/data/ecam-sample-metadata.tsv")
# remove var types: #q2:types
ecam_meta <- ecam_meta[-1, ]
ecam_feature_table <- readr::read_tsv("https://raw.githubusercontent.com/FrederickHuangLin/ANCOM/master/data/ecam-table-taxa.tsv", skip = 1)
taxa <- ecam_feature_table$`feature-id`
feature_table <- dplyr::select(ecam_feature_table, -`feature-id`)

# taxa table
taxa_table <- lapply(taxa, strsplit, split = ";", fixed = TRUE)
taxa_table <- lapply(taxa_table, unlist)
taxa_table <- do.call(rbind, taxa_table)
taxa_table <- data.frame(taxa_table)
names(taxa_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# make sure all unknown taxa as prefix__
# set the NA|unknown|unclassified to __, and then add prefix
prefixes <- c("k", "p", "c", "o", "f", "g")
taxa_table <- purrr::map2_df(
  taxa_table, prefixes,
  ~ ifelse(.x == "__", paste0(.y, .x), .x)
)
ecam_meta <- phyloseq::sample_data(ecam_meta)
row.names(ecam_meta) <- names(feature_table)


ecam <- phyloseq::phyloseq(
  phyloseq::otu_table(as(feature_table, "matrix"), taxa_are_rows = TRUE),
  phyloseq::tax_table(as(taxa_table, "matrix")),
  phyloseq::sample_data(ecam_meta)
)
usethis::use_data(ecam, overwrite = TRUE)


##### Zeybel et al. - 2022
# ../../bookdown/Microbiota_notes/dataset/Zeybel-2022/result/
#cross-section data
Zeybel_2022_gut <- readRDS("Zeybel_2022_gut_MGS_ps.RDS")
Zeybel_2022_oral <- readRDS("Zeybel_2022_oral_MGS_ps.RDS")
Zeybel_2022_metabolite <- readRDS("Zeybel_2022_fecal_metabolite_se.RDS")
Zeybel_2022_protein <- readRDS("Zeybel_2022_plasma_protein_se.RDS")

usethis::use_data(Zeybel_2022_gut, overwrite = TRUE)
usethis::use_data(Zeybel_2022_oral, overwrite = TRUE)
usethis::use_data(Zeybel_2022_metabolite, overwrite = TRUE)
usethis::use_data(Zeybel_2022_protein, overwrite = TRUE)

# longitudinal data
Zeybel_2022_gut_paired <- readRDS("Zeybel_2022_gut_MGS_ps_Paired.RDS")
Zeybel_2022_oral_paired <- readRDS("Zeybel_2022_oral_MGS_ps_Paired.RDS")
Zeybel_2022_metabolite_paired <- readRDS("Zeybel_2022_fecal_metabolite_se_Paired.RDS")
Zeybel_2022_protein_paired <- readRDS("Zeybel_2022_plasma_protein_se_Paired.RDS")

usethis::use_data(Zeybel_2022_gut_paired, overwrite = TRUE)
usethis::use_data(Zeybel_2022_oral_paired, overwrite = TRUE)
usethis::use_data(Zeybel_2022_metabolite_paired, overwrite = TRUE)
usethis::use_data(Zeybel_2022_protein_paired, overwrite = TRUE)
