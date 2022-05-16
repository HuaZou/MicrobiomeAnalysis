#' 16S rRNA data from "Moving pictures of the human microbiome"
#'
#' 16S read counts and phylogenetic tree file of 34 Illumina samples derived
#' from Moving Pictures of the Human Microbiome (Caporaso et al.) Group label:
#' gut, left palm, right palm, and tongue - indicating different sampled body
#' sites.
#'
#' @format a [phyloseq::phyloseq] object
#' @references
#' Caporaso, et al. Moving pictures of the human microbiome. Genome Biol 12,
#' R50 (2011).
#'
#' \url{https://doi.org/10.1186/gb-2011-12-5-r50}
#'
#' @source Data was downloaded from https://www.microbiomeanalyst.ca
#'
#'
#' @name data-caporaso
#' @aliases caporaso
#' @docType data
#' @author Yang Cao
NA

#' 16S rRNA data of 94 patients from CID 2012
#'
#' Data from a cohort of 94 Bone Marrow Transplant patients previously published
#' on in CID
#'
#' @format a [phyloseq::phyloseq] object
#' @references
#' Ying, et al. Intestinal Domination and the Risk of Bacteremia in Patients
#' Undergoing Allogeneic Hematopoietic Stem Cell Transplantation,
#' Clinical Infectious Diseases, Volume 55, Issue 7, 1 October 2012,
#' Pages 905–914,
#'
#' \url{https://academic.oup.com/cid/article/55/7/905/428203}
#'
#' @source \url{https://github.com/ying14/yingtools2/tree/master/data}
#' @name data-cid_ying
#' @aliases cid_ying
#' @docType data
#' @author Yang Cao
NA


#' Enterotypes data of 39 samples
#'
#' The data contains 22 European metagenomes from Danish, French, Italian,
#' and Spanish individuals, and 13 Japanese and 4 American.
#'
#' @format a [`phyloseq::phyloseq-class`] object
#' @references
#' Arumugam, Manimozhiyan, et al. Enterotypes of the human gut microbiome.
#' nature 473.7346 (2011): 174-180.
#' @name data-enterotypes_arumugam
#' @aliases enterotypes_arumugam
#' @docType data
#' @author Yang Cao
NA


#' Data from a study on colorectal cancer (kostic 2012)
#'
#' The data from a study on colorectal cancer. Samples that had no `DIAGNOSIS`
#' attribute assigned and with less than 500 reads (counts) were removed, and
#' 191 samples remains (91 healthy and 86 Tumors).
#'
#' @format a [`phyloseq::phyloseq-class`] object
#'
#' @references
#' Kostic et al. Genomic analysis identifies association of Fusobacterium with
#' colorectal carcinoma. Genome research, 2012, 22(2), 292-298.
#' @name data-kostic_crc
#' @aliases kostic_crc
#' @docType data
#' @author Yang Cao
NA


#' Data from Early Childhood Antibiotics and the Microbiome (ECAM) study
#'
#' The data from a subset of the Early Childhood Antibiotics and the
#' Microbiome (ECAM) study, which tracked the microbiome composition and
#' development of 43 infants in the United States from birth to 2 years of age,
#' identifying microbiome associations with antibiotic exposure, delivery mode,
#' and diet.
#'
#' @format a [`phyloseq::phyloseq-class`] object
#'
#' @references
#' Bokulich, Nicholas A., et al. "Antibiotics, birth mode, and diet shape
#' microbiome maturation during early life." Science translational medicine
#' 8.343 (2016): 343ra82-343ra82.
#'
#' \url{https://github.com/FrederickHuangLin/ANCOM/tree/master/data}
#' @name data-ecam
#' @aliases ecam
#' @docType data
NA
