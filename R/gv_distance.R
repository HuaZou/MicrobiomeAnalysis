#' @title Calculate the distance among samples
#'
#' @description
#' The function is to calculate the distance among samples
#'
#' @author Created by Hua Zou (5/14/2022 Shenzhen China)
#'
#' @param object (Required). a [`phyloseq::phyloseq-class`] or
#' [`Biobase::ExpressionSet`] object.
#' @param level (Optional). character. Summarization
#' level (from \code{rank_names(pseq)}, default: NULL).
#' @param method (Required). character. Provide one of the currently supported
#' options. See `distanceMethodList` for a detailed list of the supported options
#' and links to accompanying documentation. Options include:
#'  * "unifrac" : unweighted UniFrac distance.
#'  * "wunifrac": weighted-UniFrac distance.
#'  * "GUniFrac": The variance-adjusted weighted UniFrac distances (default: alpha=0.5).
#'  * "bray": bray crutis distance.
#'  * "dpcoa": sample-wise distance used in Double Principle Coordinate Analysis.
#'  * "jsd": Jensen-Shannon Divergence.
#'  Alternatively, you can provide a character string that defines a custom
#'  distance method, if it has the form described in `designdist`.
#' @param alpha (Optional). numeric. the parameter for "GUniFrac" controlling
#' weight on abundant lineages (default: 0.5).
#'
#' @return
#'   distance object, which could be applied for ordination analysis.
#'
#' @import phyloseq
#' @import dplyr
#' @import vegan
#' @importFrom stats setNames as.dist
#'
#' @usage run_distance(
#'    object,
#'    level = c(NULL, "Kingdom", "Phylum", "Class",
#'            "Order", "Family", "Genus",
#'            "Species", "Strain", "unique"),
#'    method = c("unifrac", "wunifrac", "GUniFrac", "bray", "dpcoa", "jsd"),
#'    alpha = 0.5)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data("enterotypes_arumugam")
#' res <- run_distance(enterotypes_arumugam,
#'                     method = "bray")
#' }
#'
run_distance <- function(
      object,
      level = NULL,
      method = "bray",
      alpha = 0.5) {


  # data("enterotypes_arumugam")
  # object = enterotypes_arumugam
  # level = "Phylum"
  # method = "bray"
  # alpha = 0.5

  # phyloseq object
  if (all(!is.null(object), inherits(object, "phyloseq"))) {
    ps <- check_sample_names(object = object)

    # taxa level
    if (!is.null(level)) {
      ps <- aggregate_taxa(x = ps, level = level)
    } else {
      ps <- ps
    }

    if (!is.null(ps@phy_tree) &
        (method %in% c("unifrac", "wunifrac", "GUniFrac"))) {
      method <- match.arg(
        method,
        c("unifrac", "wunifrac", "GUniFrac")
      )
      if (method == "GUniFrac") {
        otu_tab <- phyloseq::otu_table(ps)
        tree_tab <- phyloseq::phy_tree(ps)
      }
    }
    # distance
    if (method == "GUniFrac") {
      res_temp <- .get_GUniFrac(otu.tab = otu_tab,
                               tree = tree_tab,
                               alpha = alpha)
      disMatrix <- stats::as.dist(res_temp$unifracs[, , paste0("d_", alpha)])
    } else {
      disMatrix <- phyloseq::distance(physeq = ps, method = method)
    }
  } else if (all(!is.null(object), inherits(object, "ExpressionSet"))) {
    prf_tab <- Biobase::exprs(object) %>% data.frame()
    if (any(prf_tab < 0)) {
      if (method %in% c("bray", "jaccard")) {
        message(method, " is not suitable for Negative values and Replace it by euclidean")
        method <- "euclidean"
      }
    }

    disMatrix <- vegan::vegdist(prf_tab, method = method)
  }

  return(disMatrix)
}

#' @title Generalized UniFrac distances for comparing microbial communities
#'
#' @description
#' A generalized version of commonly used UniFrac distances.
#' The unweighted and weighted UniFrac, and variance-adjusted weighted UniFrac distances are also implemented.
#'
#' @references https://github.com/cran/GUniFrac/blob/master/R/GUniFrac.R
#'
#' @author Created by Jun Chen; modified by Hua Zou (5/14/2022 Shenzhen China)
#'
#' @param otu.tab (Required). a matrix, an OTU count table, row - n sample, column - q OTU.
#' @param tree (Required). tree. a rooted phylogenetic tree of R class “phylo”.
#' @param alpha (optional). numeric. the parameter controlling weight on abundant lineages.
#'
#' @return a three dimensional array containing all the UniFrac distance matrices.
#'
#' unifracs: three dimensional array containing the generalized
#'           UniFrac distances, unweighted UniFrac distance and
#'           variance adjusted UniFrac distances.
#'
#' @usage .get_GUniFrac(otu.tab,
#'                  tree,
#'                  alpha)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data("cid_ying")
#' otu.tab <- phyloseq::otu_table(cid_ying)
#' tree <- phyloseq::phy_tree(cid_ying)
#'
#' unifracs <- .get_GUniFrac(otu.tab, tree, alpha=c(0, 0.5, 1))$unifracs
#'
#' dw <- unifracs[, , "d_1"]    # Weighted UniFrac
#' du <- unifracs[, , "d_UW"]   # Unweighted UniFrac
#' dv <- unifracs[, , "d_VAW"]  # Variance adjusted weighted UniFrac
#' d0 <- unifracs[, , "d_0"]    # GUniFrac with alpha 0
#' d5 <- unifracs[, , "d_0.5"]  # GUniFrac with alpha 0.5
#'
#' }
#'
.get_GUniFrac <- function(
      otu.tab,
      tree,
      alpha = c(0, 0.5, 1)) {

  if (!ape::is.rooted(tree)) {
    stop("Rooted phylogenetic tree required!")
  }

  if (phyloseq::taxa_are_rows(otu.tab)) {
    otu.tab <- phyloseq::t(otu.tab)
  }

  # Convert into proportions
  otu.tab <- as.matrix(otu.tab)
  row.sum <- rowSums(otu.tab)
  otu.tab <- otu.tab / row.sum
  n <- nrow(otu.tab)

  # Construct the returning array
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste("comm", 1:n, sep="_")
  }
  # d_UW: unweighted UniFrac, d_VAW: weighted UniFrac
  dimname3 <- c(paste("d", alpha, sep="_"), "d_UW", "d_VAW")
  unifracs <- array(NA, c(n, n, length(alpha) + 2),
                    dimnames=list(rownames(otu.tab), rownames(otu.tab), dimname3))
  for (i in 1:(length(alpha)+2)) {
    for (j in 1:n) {
      unifracs[j, j, i] <- 0
    }
  }

  # Check OTU name consistency
  if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
    stop("The OTU table contains unknown OTUs! OTU names
						in the OTU table and the tree should match!" )
  }

  # Get the subtree if tree contains more OTUs
  absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
  if (length(absent) != 0) {
    tree <- ape::drop.tip(tree, absent)
    warning("The tree has more OTU than the OTU table!")
  }

  # Reorder the otu.tab matrix if the OTU orders are different
  tip.label <- tree$tip.label
  otu.tab <- otu.tab[, tip.label]

  ntip <- length(tip.label)
  nbr <- nrow(tree$edge)
  edge <- tree$edge
  edge2 <- edge[, 2]
  br.len <- tree$edge.length

  #  Accumulate OTU proportions up the tree
  cum <- matrix(0, nbr, n)							# Branch abundance matrix
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)
    cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
    node <- edge[tip.loc, 1]						# Assume the direction of edge
    node.loc <- which(edge2 == node)
    while (length(node.loc)) {
      cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }

  # Calculate various UniFrac distances
  cum.ct <- round(t(t(cum) * row.sum)) 	# For VAW
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      cum1 <- cum[, i]
      cum2 <- cum[, j]
      ind <- (cum1 + cum2) != 0
      cum1 <- cum1[ind]
      cum2 <- cum2[ind]
      br.len2 <- br.len[ind]
      mi <- cum.ct[ind, i] + cum.ct[ind, j]
      mt <- row.sum[i] + row.sum[j]
      diff <- abs(cum1 - cum2) / (cum1 + cum2)

      # Generalized UniFrac distance
      for(k in 1:length(alpha)){
        w <- br.len2 * (cum1 + cum2)^alpha[k]
        unifracs[i, j, k] <- unifracs[j, i, k] <- sum(diff * w) / sum(w)
      }

      #	Variance Adjusted UniFrac Distance
      ind2 <- (mt != mi)
      w <- br.len2 * (cum1 + cum2) / sqrt(mi * (mt - mi))
      unifracs[i, j, (k + 2)] <- unifracs[j, i, (k + 2)] <-
        sum(diff[ind2] * w[ind2]) / sum(w[ind2])

      #	Unweighted UniFrac Distance
      cum1 <- (cum1 != 0)
      cum2 <- (cum2 != 0)
      unifracs[i, j, (k + 1)] <- unifracs[j, i, (k + 1)] <-
        sum(abs(cum1 - cum2) / (cum1 + cum2) * br.len2) / sum(br.len2)
    }
  }
  return(list(unifracs=unifracs))
}
