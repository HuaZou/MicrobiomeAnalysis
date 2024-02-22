#' @title Ordination for microbiota data
#'
#' @description
#' The primary goal of ordination was considered “exploratory” (Gauch 1982a, b),
#' with the introduction of canonical correspondence analysis (CCA), ordination
#' has gone beyond mere “exploratory” analysis (ter Braak 1985) and become
#' hypothesis testing as well.
#'
#' @details
#' The primary aim of ordination is to represent multiple samples (subjects) in
#' a reduced number of orthogonal (i.e., independent) axes,
#' where the total number of axes is less than or equal to the number of samples
#'
#' @references
#' Xia, Y., Sun, J., & Chen, D. G. (2018). Statistical analysis of microbiome
#' data with R (Vol. 847). Singapore: Springer.
#'
#' @author Created by Hua Zou (8/9/2023 Shenzhen China)
#'
#' @param object (Required). a [`phyloseq::phyloseq-class`] or
#' [`SummarizedExperiment::SummarizedExperiment`] object.
#' @param level (Optional). character. Summarization
#' level (from \code{rank_names(pseq)}, default: NULL).
#' @param variable (Required). character. grouping variable for test.
#' @param transform character, the methods used to transform the microbial
#'   abundance. See [`transform_abundances()`] for more details. The
#'   options include:
#'   * "identity", return the original data without any transformation
#'     (default).
#'   * "log10", the transformation is `log10(object)`, and if the data contains
#'     zeros the transformation is `log10(1 + object)`.
#'   * "log10p", the transformation is `log10(1 + object)`.
#'   * "SquareRoot", the transformation is `Square Root`.
#'   * "CubicRoot", the transformation is `Cubic Root`.
#'   * "logit", the transformation is `Zero-inflated Logit Transformation`
#' (Does not work well for microbiome data).
#' @param norm the methods used to normalize the microbial abundance data. See
#'   [`normalize()`] for more details.
#'   Options include:
#'   * "none": do not normalize.
#'   * "rarefy": random subsampling counts to the smallest library size in the
#'     data set.
#'   * "TMM": trimmed mean of m-values. First, a sample is chosen as reference.
#'     The scaling factor is then derived using a weighted trimmed mean over the
#'     differences of the log-transformed gene-count fold-change between the
#'     sample and the reference.
#'   * "RLE", relative log expression, RLE uses a pseudo-reference calculated
#'     using the geometric mean of the gene-specific abundances over all
#'     samples. The scaling factors are then calculated as the median of the
#'     gene counts ratios between the samples and the reference.
#'   * "CSS": cumulative sum scaling, calculates scaling factors as the
#'     cumulative sum of gene abundances up to a data-derived threshold.
#' @param method (Optional). character. Ordination method (default: "PCoA"),
#' options include:
#'   * "PCA": Principal Component Analysis.
#'   * "PCoA": Principal Coordinate Analysis.
#'   * "tSNE": t-distributed stochastic neighbor embedding.
#'   * "UMAP": Uniform Manifold Approximation and Projection.
#'   * "NMDS": Non-metric Multidimensional Scaling.
#' @param distance (Optional). character. Provide one of the currently supported
#' options. See `vegan::vegdist` for a detailed list of the supported options
#' and links to accompanying documentation (default: "bray"). Options include:
#'  * "bray": bray crutis distance.
#'  * "unifrac" : unweighted UniFrac distance.
#'  * "wunifrac": weighted-UniFrac distance.
#'  * "GUniFrac": The variance-adjusted weighted UniFrac distances (default: alpha=0.5).
#'  * "dpcoa": sample-wise distance used in Double Principle Coordinate Analysis.
#'  * "jsd": Jensen-Shannon Divergence.
#'  Alternatively, you can provide a character string that defines a custom
#'  distance method, if it has the form described in `designdist`.
#' @param para (Optional). list. the additional parameters for methods.
#'   * "Perplexity": numeric; Perplexity parameter (should not be bigger than
#'   3 perplexity < nrow(X) - 1.
#'   * "Y_vars": Constraining matrix, typically of environmental variables.
#'   * "Z_vars": Conditioning matrix, the effect of which is removed ("partial out")
#'   before next step.
#'   * "scale": Scale features to unit variance (like correlations).
#'   * "center": Scale features to unit variance (like correlations).
#' @param ... (Optional). additional parameters.
#'
#' @return
#'   A list of the ordination's results.
#'
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom vegan adonis vegdist
#' @importFrom stats setNames
#'
#' @usage run_ord(
#'     object,
#'     level = NULL,
#'     variable,
#'     transform = c("identity", "log10", "log10p",
#'                   "SquareRoot", "CubicRoot", "logit"),
#'     norm = c("none", "rarefy", "TSS", "TMM",
#'              "RLE", "CSS", "CLR", "CPM"),
#'     method = c("PCA", "PCoA", "tSNE", "UMAP", "NMDS",
#'                "CA", "RDA", "CCA", "CAP"),
#'     distance = c("bray", "unifrac", "wunifrac",
#'                  "GUniFrac", "dpcoa", "jsd"),
#'     para = list(Perplexity = NULL,
#'                 Y_vars = NULL,
#'                 Z_vars = NULL,
#'                 scale = TRUE,
#'                 center = TRUE,),
#'     ...)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' # phyloseq object
#' data("Zeybel_2022_gut")
#' ps_zeybel <- summarize_taxa(Zeybel_2022_gut, level = "Genus")
#' ord_result <- run_ord(
#'   object = ps_zeybel,
#'   variable = "LiverFatClass",
#'   method = "PCoA")
#'
#' # SummarizedExperiment object
#' data("Zeybel_2022_protein")
#' Zeybel_2022_protein_imp <- impute_abundance(
#'   Zeybel_2022_protein,
#'   group = "LiverFatClass",
#'   ZerosAsNA = TRUE,
#'   RemoveNA = TRUE,
#'   cutoff = 20,
#'   method = "knn")
#' ord_result <- run_ord(
#'   object = Zeybel_2022_protein_imp,
#'   variable = "LiverFatClass",
#'   method = "PCA")
#'
#' }
#'
run_ord <- function(
    object,
    level = NULL,
    variable,
    transform = "identity",
    norm = "none",
    method = "PCoA",
    distance = "bray",
    para = list(Perplexity = NULL,
                Y_vars = NULL,
                Z_vars = NULL,
                scale = TRUE,
                center = TRUE),
    ...) {

  # data("Zeybel_2022_gut")
  # ps = summarize_taxa(Zeybel_2022_gut, level = "Genus")
  # level = NULL
  # variable = "LiverFatClass"
  # transform = "identity"
  # norm = "none"
  # method = "PCoA"
  # distance = "bray"
  # para = list(Perplexity = NULL,
  #             Y_vars = NULL,
  #             Z_vars = NULL,
  #             scale = FALSE)

  # data("Zeybel_2022_protein")
  # object = Zeybel_2022_protein
  # level = NULL
  # variable = "LiverFatClass"
  # transform = "identity"
  # norm = "none"
  # method = "PCA"
  # distance = "bray"
  # para = list(scale = TRUE,
  #             center = TRUE)

  # transform
  transform <- match.arg(
    transform, c("identity", "log10", "log10p",
                 "SquareRoot", "CubicRoot", "logit"))
  # normalization
  norm <- match.arg(
    norm, c("none", "rarefy", "TSS", "TMM",
            "RLE", "CSS", "CLR", "CPM")
  )

  # ordination approaches
  method <- match.arg(
    method,
    c("PCA", "PCoA", "tSNE", "UMAP",
      "NMDS", "CA", "RDA", "CCA", "CAP")
  )

  if (!is.null(object)) {
    # stopifnot(inherits(object, "phyloseq"))
    if (inherits(object, "phyloseq")) {
      # add tryCatch (2023/12/2 update)
      tryCatch(
        expr = {
          ps <- preprocess_ps(object)
        },
        error = function(e){
          message('preprocess_ps Caught an error!')
          print(e)
        }
      )
      if (!is.null(level)) {
        ps <- summarize_taxa(ps, level = level)
      } else {
        ps <- ps
      }
      # preprocess phyloseq object
      ps <- transform_abundances(ps, transform = transform)
      # normalize the data
      ps_normed <- normalize(object, method = norm)
      # otu table: row -> taxa; column -> SampleID
      if (!(otu_table(ps_normed)@taxa_are_rows)) {
        otu_tab <- phyloseq::otu_table(ps_normed) %>%
          data.frame() %>%
          t() %>% data.frame()
      } else {
        otu_tab <- phyloseq::otu_table(ps_normed) %>%
          data.frame()
      }
      # sample data
      sam_tab <- phyloseq::sample_data(ps_normed) %>%
        data.frame()
    } else if (inherits(object, "SummarizedExperiment")) {
      ps_normed <- transform_abundances(object, transform = transform)
      otu_tab <- SummarizedExperiment::assay(ps_normed) %>%
        data.frame()
      sam_tab <- SummarizedExperiment::colData(ps_normed) %>%
        data.frame()

    }
  }

  # check whether group is valid
  if (!variable %in% colnames(sam_tab)) {
    stop(
      variable, " are not contained in the `sample_data` of `ps`",
      call. = FALSE
    )
  }

  # whether there is negative value in the otu table
  if (any(otu_tab < 0)) {
    distance <- "euclidean"
  }

  # filter taxa with a constant/zero value
  cutoff_value <- 1/nrow(otu_tab)
  if (cutoff_value > 0.1) {
    cutoff_value <- 0.01
  }
  ps_trim <- trim_prevalence(
    object = ps_normed,
    cutoff = cutoff_value,
    trim = "feature")

  # beta dispersion
  bd_fit <- run_betadisper(
    object = ps_trim,
    variable = variable,
    method = distance)

  bd_pva <- bd_fit$tab$`Pr(>F)`[1]

  # PERMANOVA
  if (bd_pva < 0.05) {
    print("Pvalue of beta dispersion less than 0.05")
  }

  bd_sig_lab <- paste(cut(bd_pva,
                           breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                           label = c("***", "**", "*", "ns")))

  bd_res <- bquote(atop(atop("beta dispersion"),
                          atop("p-value="~.(bd_pva)~.(bd_sig_lab), phantom())))

  # permanova
  ## add tryCatch (2023/12/2 update)
  tryCatch(
    expr = {
      per_res <- .calculatePERMANOVA(
        object = ps_trim,
        distance = distance,
        group = variable)
    },
    error = function(e){
      message('.calculatePERMANOVA Caught an error!')
      print(e)
    }
  )

  # metadata and profile
  if (inherits(object, "phyloseq")) {
    meta_ps <- phyloseq::sample_data(ps_trim) %>%
      data.frame()
    prof_ps <- phyloseq::otu_table(ps_trim) %>%
      data.frame()
  } else if (inherits(object, "SummarizedExperiment")) {
    meta_ps <- SummarizedExperiment::colData(ps_trim) %>%
      data.frame()
    prof_ps <- SummarizedExperiment::assay(ps_trim) %>%
      data.frame()
  }

  # ordination methods
  ## add tryCatch (2023/12/2 update)
  tryCatch(
    expr = {
      res <- .calculateOrdination(
        metadata = meta_ps,
        profile = prof_ps,
        method = method,
        distance = distance,
        group = variable,
        permanova = per_res,
        betadisperion = bd_res,
        parameters = para)
    },
    error = function(e){
      message('.calculateOrdination Caught an error!')
      print(e)
    }
  )

  return(res)
}


#' calculate results of permanova
#' @noRd
#' @keywords internal
#'
.calculatePERMANOVA <- function(object, distance, group) {

  set.seed(123)

  per <- run_PERMANOVA(
    object = object,
    method = distance,
    variables = group)

  #per_temp <- per[rownames(per)%in%"Compvar", ,]

  pvalue <- signif(per[, 7], 3)
  r2 <- signif(per[, 6], 3)

  signi_label <- paste(cut(pvalue,
                           breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                           label = c("***", "**", "*", "ns")))

  res_PERMANOVA <- bquote(atop(atop("PERMANOVA", R^2==~.(r2)),
                               atop("p-value="~.(pvalue)~.(signi_label), phantom())))

  return(res_PERMANOVA)
}


#' calculate results of Ordination
#' @noRd
#' @keywords internal
#'
.calculateOrdination <- function(
    metadata,
    profile,
    method,
    distance,
    group,
    permanova,
    betadisperion,
    parameters) {

  # metadata = meta_ps
  # profile = prof_ps
  # method = method
  # distance = distance
  # group = variable
  # permanova = per_res
  # betadisperion = bd_res
  # parameters = para

  # whether the parameters list is null
  if (method %in% c("tSNE")) {
    if(is.null(parameters$Perplexity)) {
      warning("please provide the additional parameters for ordination")
      Perplexity <- 2
    } else {
      Perplexity <- parameters$Perplexity
    }
  }

  # ordination methods
  ord_list <- switch(
    method,
    "PCA" = PCA_fun(x = profile, y = metadata,
                    center_ = parameters$center,
                    scale_ = parameters$scale),
    "PCoA" = PCoA_fun(x = profile, y = metadata,
                     n = distance),
    "tSNE" = tSNE_fun(x = profile, y = metadata,
                      n = Perplexity,
                      center_ = parameters$center,
                      scale_ = parameters$scale,
                      m = method),
    "UMAP" = UMAP_fun(x = profile, y = metadata,
                      center_ = parameters$center,
                      scale_ = parameters$scale,
                      m = method),
    "NMDS" = NMDS_fun(x = profile, y = metadata,
                      n = distance,
                      m = method),
    "CA" = CA_fun(x = profile, y = metadata),
    "RDA" = RDA_fun(x = profile, y = metadata,
                    n = parameters))


  # res <- list(fit = datfit,
  #             dat = score,
  #             explains = explains,
  #             eigvalue = eig,
  #             PERMANOVA = permanova,
  #             BETADISPER = betadisperion)

  res <- list(fit = ord_list[[1]],
              dat = ord_list[[2]],
              explains = ord_list[[3]],
              eigvalue = ord_list[[4]],
              PERMANOVA = permanova,
              BETADISPER = betadisperion)

  return(res)
}

#' @title PCA
#' @noRd
#' @keywords internal
#'
PCA_fun <- function(x, y, center_, scale_) {

  datfit <- stats::prcomp(scale(t(x), center = center_, scale = scale_))
  eig <- get_eigValue(datfit)
  explains <- paste0(paste0("PCA", seq(2)), "(", paste0(round(eig[1:2, 2], 2), "%"), ")")
  # plotdata
  score <- dplyr::inner_join(datfit$x %>% data.frame() %>%
                               dplyr::select(c(1:4)) %>%
                               stats::setNames(paste0("Axis", seq(4))) %>%
                               tibble::rownames_to_column("TempRowNames"),
                             y %>% tibble::rownames_to_column("TempRowNames"),
                             by = "TempRowNames")

  return(list(datfit, score, explains, eig))
}

#' @title PCoA
#' @noRd
#' @keywords internal
#'
PCoA_fun <- function(x, y, n) {

  ps_pcoa <- phyloseq::phyloseq(
    phyloseq::otu_table(x, taxa_are_rows = TRUE),
    phyloseq::sample_data(y))

  dat_dis <- run_distance(object = ps_pcoa, method = n)
  datfit <- ape::pcoa(dat_dis)
  eig <- datfit$values[, "Eigenvalues"]
  eig_var <- eig[1:2]
  eig_var_explain <- round(eig_var / sum(eig), 4) * 100
  # explains variable
  explains <- paste0(paste0("PCoA", seq(2)), " (", paste0(eig_var_explain, "%"), ")")
  # plotdata
  score <- dplyr::inner_join(datfit$vectors[, c(1:4)] %>% data.frame() %>%
                               stats::setNames(paste0("Axis", seq(4))) %>%
                               tibble::rownames_to_column("TempRowNames"),
                             y %>% tibble::rownames_to_column("TempRowNames"),
                             by = "TempRowNames")

  return(list(datfit, score, explains, eig))
}

#' @title Rtsne
#' @noRd
#' @keywords internal
#'
tSNE_fun <- function(x, y, n, center_, scale_, m) {

  datfit <- Rtsne::Rtsne(scale(t(x), center = center_, scale = scale_),
                         dims = 3,
                         perplexity = n,
                         verbose = FALSE,
                         max_iter = 500,
                         eta = 200)
  explains <- paste0(m, c(1, 2))
  eig <- NULL
  point <- datfit$Y %>% data.frame() %>%
    dplyr::select(c(1:3)) %>%
    stats::setNames(paste0("Axis", seq(3)))
  rownames(point) <- colnames(x)
  score <- dplyr::inner_join(point %>% tibble::rownames_to_column("TempRowNames"),
                             y %>% tibble::rownames_to_column("TempRowNames"),
                             by = "TempRowNames")

  return(list(datfit, score, explains, eig))
}

#' @title umap
#' @noRd
#' @keywords internal
#'
UMAP_fun <- function(x, y, center_, scale_, m) {

  datfit <- umap::umap(scale(t(x), center = center_, scale = scale_))
  explains <- paste0(m, c(1, 2))
  eig <- NULL
  point <- datfit$layout %>% data.frame() %>%
    dplyr::select(c(1:2)) %>%
    stats::setNames(c("Axis1", "Axis2"))
  rownames(point) <- colnames(x)
  score <- dplyr::inner_join(point %>% tibble::rownames_to_column("TempRowNames"),
                             y %>% tibble::rownames_to_column("TempRowNames"),
                             by = "TempRowNames")

  return(list(datfit, score, explains, eig))
}

#' @title metaMDS
#' @noRd
#' @keywords internal
#'
NMDS_fun <- function(x, y, n, m) {

  datfit <- vegan::metaMDS(t(x), dist = n)
  explains <- paste0(m, c(1, 2))
  eig <- NULL
  point <- datfit$points %>% data.frame() %>%
    dplyr::select(c(1:2)) %>%
    stats::setNames(c("Axis1", "Axis2"))
  rownames(point) <- colnames(x)
  score <- dplyr::inner_join(point %>% tibble::rownames_to_column("TempRowNames"),
                             y %>% tibble::rownames_to_column("TempRowNames"),
                             by = "TempRowNames")

  return(list(datfit, score, explains, eig))
}

#' @title CA
#' @noRd
#' @keywords internal
#'
CA_fun <- function(x, y) {

  datfit <- CA(x, graph = FALSE) # FactoMineR R package
  eig_var_explain <- c(signif(datfit$eig[1, 2], 4), signif(datfit$eig[2, 2], 4))
  eig <- eig_var_explain
  explains <- paste0(paste0("CA", seq(2)), " (", paste0(eig_var_explain, "%"), ")")
  # default: samples points
  point <- datfit$col$coord %>% data.frame() %>%
    dplyr::select(c(1:2)) %>%
    stats::setNames(c("Axis1", "Axis2"))
  rownames(point) <- colnames(x)
  score <- dplyr::inner_join(point %>% tibble::rownames_to_column("TempRowNames"),
                             y %>% tibble::rownames_to_column("TempRowNames"),
                             by = "TempRowNames")

  return(list(datfit, score, explains, eig))
}

#' @title RDA
#' @noRd
#' @keywords internal
#'
RDA_fun <- function(x, y, n) {

  if (!is.null(n$Y_vars)) {
    Y_var <- y[, n$Y_vars]
  } else {
    Y_var <- NULL
  }

  if (!is.null(n$Z_vars)) {
    Z_var <- y[, n$Z_vars]
  } else {
    Z_var <- NULL
  }

  datfit <- vegan::rda(X = t(x),
                       Y = Y_var,
                       Z = Z_var,
                       scale = n$scale)
  if (all(is.null(n$Y_vars), is.null(n$Z_vars))) {
    eig_var_explain <- as.numeric(c(signif(datfit$CA$eig[1], 4),
                                    signif(datfit$CA$eig[2], 4)))
    eig <- eig_var_explain
    explains <- paste0(paste0("CA", seq(2)), " (", paste0(eig_var_explain, "%"), ")")
    # default: samples points
    point <- datfit$CA$u %>% data.frame() %>%
      dplyr::select(c(1:2)) %>%
      stats::setNames(c("Axis1", "Axis2"))
  } else {
    eig_var_explain <- as.numeric(c(signif(datfit$CCA$eig[1], 4),
                                    signif(datfit$CCA$eig[2], 4)))
    eig <- eig_var_explain
    explains <- paste0(paste0("RDA", seq(2)), " (", paste0(eig_var_explain, "%"), ")")
    # default: samples points
    point <- datfit$CCA$u %>% data.frame() %>%
      dplyr::select(c(1:2)) %>%
      stats::setNames(c("Axis1", "Axis2"))
  }

  rownames(point) <- colnames(x)
  score <- dplyr::inner_join(point %>% tibble::rownames_to_column("TempRowNames"),
                             y %>% tibble::rownames_to_column("TempRowNames"),
                             by = "TempRowNames")

  return(list(datfit, score, explains, eig))
}


# https://github.com/husson/FactoMineR/blob/master/R/CA.R
#' @title Correspondence Analysis
#'
#' @noRd
#' @keywords internal
#'
CA <- function (X,
                ncp = 5,
                row.sup = NULL,
                col.sup = NULL,
                quanti.sup = NULL,
                quali.sup = NULL,
                graph = FALSE,
                axes = c(1, 2),
                row.w = NULL,
                excl = NULL) {

  fct.eta2 <- function(vec,x,weights) {
    VB <- function(xx) {
      return(sum((colSums((tt * xx) * weights)^2) / ni))
    }
    tt <- tab_disjonctif(vec)
    ni <- colSums(tt * weights)
    unlist(lapply(as.data.frame(x), VB)) / colSums(x * x * weights)
  }

  if (is.table(X)) {
    X <- matrix(as.vector(X), nrow(X), dimnames=dimnames(X))
  }

  if (is.null(rownames(X))) {
    rownames(X) <- 1:nrow(X)
  }

  if (is.null(colnames(X))) {
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "V")
  }

  if (!is.null(row.sup) & !is.numeric(row.sup)) {
    row.sup <- which(rownames(X)%in%row.sup)
  }

  if (!is.null(col.sup) & !is.numeric(col.sup)) {
    col.sup <- which(colnames(X)%in%col.sup)
  }

  if (!is.null(quali.sup) & !is.numeric(quali.sup)) {
    quali.sup <- which(colnames(X)%in%quali.sup)
  }

  if (!is.null(quanti.sup) & !is.numeric(quanti.sup)) {
    quanti.sup <- which(colnames(X)%in%quanti.sup)
  }

  X <- as.data.frame(X)
  is.quali <- which(!unlist(lapply(X, is.numeric)))
  X[, is.quali] <- lapply(X[, is.quali, drop=FALSE], as.factor)

  for (i in is.quali) {
    X[, i] <- as.factor(X[, i])
  }

  X <- droplevels(X)
  Xtot <- X
  if (any(!unlist(lapply(X, is.numeric)))) {
    auxi <- NULL
    for (j in (1:ncol(X))[!((1:ncol(X))%in%quali.sup)]) {
      if (!is.numeric(X[, j])) {
        auxi <- c(auxi, colnames(X)[j])
        }
    }

    if (!is.null(auxi)) {
      stop(paste("\nThe following variables are not quantitative: ", auxi))
    }
  }

  if (!inherits(X, "data.frame")) {
    stop("X is not a data.frame")
  }

  if (!is.null(row.sup)) {
    X <- as.data.frame(X[-row.sup, ])
  }

  if ((!is.null(col.sup))||(!is.null(quanti.sup))||(!is.null(quali.sup))) {
    X <- as.data.frame(X[, -c(col.sup, quanti.sup, quali.sup)])
  }

  if (any(apply(X, 1, sum) == 0)) {
    warning(paste0("The rows ", paste(rownames(X)[which(apply(X, 1, sum) == 0)], collapse = ", "),
                   " sum at 0. They were suppressed from the analysis"))
    X <- X[-which(apply(X, 1, sum) == 0), , drop = FALSE]
  }
  if (any(apply(X, 2, sum) == 0)){
    warning(paste0("The columns ", paste(colnames(X)[which(apply(X, 2, sum) == 0)], collapse = ", "),
                   " sum at 0. They were suppressed from the analysis"))
    X <- X[, -which(apply(X, 2, sum) == 0), drop = FALSE]
  }
  ### 3 lignes rajoutees
  if (is.null(row.w)) {
    row.w <- rep(1, nrow(X))
  }

  row.w.init <- row.w

  if (length(row.w) != nrow(X)) {
    stop("length of vector row.w should be the number of active rows")
  }

  total <- sum(X * row.w)
  F <- as.matrix(X) * (row.w/total)
  marge.col <- colSums(F)
  marge.row <- rowSums(F)
  ncp <- min(ncp, (nrow(X) - 1), (ncol(X) - 1))
  Tc <- t(t(F/marge.row)/marge.col) - 1
  if(!is.null(excl)) {
    marge.col[excl] <- 1e-15
  }
  tmp <- svd_triplet(Tc, row.w = marge.row, col.w = marge.col, ncp = ncp)
  if(!is.null(excl)) marge.col[excl] <- 0
  eig <- tmp$vs^2
  vp <- matrix(NA, length(eig), 3)
  rownames(vp) <- paste("dim", 1:length(eig))
  colnames(vp) <- c("eigenvalue", "percentage of variance", "cumulative percentage of variance")
  vp[, "eigenvalue"] <- eig
  vp[, "percentage of variance"] <- (eig/sum(eig)) * 100
  vp[, "cumulative percentage of variance"] <- cumsum(vp[, "percentage of variance"])
  V <- tmp$V
  U <- tmp$U
  eig <- eig[1:ncol(U)]
  coord.col <- t(t(V) * sqrt(eig))
  coord.row <- t(t(U) * sqrt(eig))
  dist2.col <- colSums(Tc^2 * marge.row)
  contrib.col <- t(t(coord.col^2 * marge.col)/eig)
  cos2.col <- coord.col^2/dist2.col
  colnames(coord.col) <- colnames(contrib.col) <- colnames(cos2.col) <- paste("Dim", 1:length(eig))
  rownames(coord.col) <- rownames(contrib.col) <- rownames(cos2.col) <- attributes(X)$names
  dist2.row <- rowSums(t(t(Tc^2) * marge.col))
  contrib.row <- t(t(coord.row^2 * marge.row)/eig)
  cos2.row <- coord.row^2/dist2.row
  colnames(coord.row) <- colnames(contrib.row) <- colnames(cos2.row) <- paste("Dim", 1:length(eig))
  rownames(coord.row) <- rownames(contrib.row) <- rownames(cos2.row) <- attributes(X)$row.names
  inertia.row <- marge.row*dist2.row
  inertia.col <- marge.col*dist2.col
  names(inertia.col) <- attributes(coord.col)$row.names
  names(inertia.row) <- attributes(coord.row)$row.names

  res.call <- list(X = X,
                   marge.col = marge.col,
                   marge.row = marge.row,
                   ncp = ncp,
                   row.w = row.w,
                   excl = excl,
                   call = match.call(),
                   Xtot = Xtot,
                   N = sum(row.w*rowSums(X)))
  res.col <- list(coord = as.matrix(coord.col[, 1:ncp]),
                  contrib = as.matrix(contrib.col[, 1:ncp] * 100),
                  cos2 = as.matrix(cos2.col[, 1:ncp]),
                  inertia = inertia.col)
  res.row <- list(coord = coord.row[, 1:ncp],
                  contrib = contrib.row[, 1:ncp] * 100,
                  cos2 = cos2.row[, 1:ncp],
                  inertia = inertia.row)
  res <- list(eig = vp[1:min(nrow(X) - 1, ncol(X) - 1), , drop = FALSE],
              call = res.call,
              row = res.row,
              col = res.col,
              svd = tmp)
  if (!is.null(row.sup)){
    X.row.sup <- as.data.frame(Xtot[row.sup, ])
    if ((!is.null(col.sup))||(!is.null(quanti.sup))||(!is.null(quali.sup))) {
      X.row.sup <- as.data.frame(X.row.sup[, -c(col.sup, quanti.sup, quali.sup)])
  }
    somme.row <- rowSums(X.row.sup)
    X.row.sup <- X.row.sup/somme.row
    coord.row.sup <- crossprod(t(as.matrix(X.row.sup)), V)
    # modif
    dist2.row <- rowSums(t((t(X.row.sup) - marge.col)^2/marge.col))
    cos2.row.sup <- coord.row.sup^2/dist2.row
    coord.row.sup <- coord.row.sup[, 1:ncp, drop = FALSE]
    cos2.row.sup <- cos2.row.sup[, 1:ncp, drop = FALSE]
    colnames(coord.row.sup) <- colnames(cos2.row.sup) <- paste("Dim", 1:ncp)
    rownames(coord.row.sup) <- rownames(cos2.row.sup) <- rownames(X.row.sup)
    res.row.sup <- list(coord = coord.row.sup, cos2 = cos2.row.sup)
    res$row.sup <- res.row.sup
    res$call$row.sup <- row.sup
  }
  if (!is.null(col.sup)) {
    X.col.sup <- as.data.frame(Xtot[, col.sup])
    if (!is.null(row.sup)) {
      X.col.sup <- as.data.frame(X.col.sup[-row.sup, ])
    }
    ## 1 ligne rajoutee
    X.col.sup <- X.col.sup * row.w
    colnames(X.col.sup) <- colnames(Xtot)[col.sup]
    somme.col <- colSums(X.col.sup)
    X.col.sup <- t(t(X.col.sup)/somme.col)
    coord.col.sup <- crossprod(as.matrix(X.col.sup), U)

    dist2.col <- colSums((X.col.sup-marge.row)^2/marge.row)
    coord.col.sup <- as.matrix(coord.col.sup[, 1:ncp, drop = FALSE])
    cos2.col.sup <- coord.col.sup^2/dist2.col
    cos2.col.sup <- cos2.col.sup[, 1:ncp, drop = FALSE]
    colnames(coord.col.sup) <- colnames(cos2.col.sup) <- paste("Dim", 1:ncp)
    rownames(coord.col.sup) <- rownames(cos2.col.sup) <- colnames(X.col.sup)
    res.col.sup <- list(coord = coord.col.sup, cos2 = cos2.col.sup)
    res$col.sup <- res.col.sup
    res$call$col.sup <- col.sup
  }
  ## Ajout variable quanti supp.
  if (!is.null(quanti.sup)) {
    coord.quanti.sup <- matrix(NA, length(quanti.sup), ncp)
    if (is.null(row.sup)) {
      coord.quanti.sup <- cov.wt(cbind.data.frame(res$row$coord,Xtot[, quanti.sup, drop = FALSE]),
                                 cor = TRUE,
                                 wt = marge.row,method="ML")$cor[-(1:ncp), 1:ncp, drop = FALSE]
    } else {coord.quanti.sup <- cov.wt(cbind.data.frame(res$row$coord, Xtot[-row.sup,quanti.sup, drop = FALSE]),
                                       wt=  marge.row,
                                       cor = TRUE,
                                       method = "ML")$cor[-(1:ncp), 1:ncp, drop = FALSE]
    }
    dimnames(coord.quanti.sup) <- list(colnames(Xtot)[quanti.sup], paste("Dim", 1:ncp, sep = "."))
    res$quanti.sup$coord <- coord.quanti.sup
    res$quanti.sup$cos2 <- coord.quanti.sup^2
    res$call$quanti.sup <- quanti.sup
  }
  ## Ajout variable quali supp.
  if (!is.null(quali.sup)) {

    if (!is.null(row.sup)) {
      X.del <- as.data.frame(Xtot[-row.sup, -c(col.sup, quanti.sup, quali.sup)])
    } else {
      X.del <- Xtot[, -c(col.sup, quanti.sup, quali.sup)]
    }
    X.quali.sup <- NULL
    Xtot2 <- Xtot
    if (!is.null(row.sup)) {
      Xtot2 <- Xtot[-row.sup, ]
    }
    for (j in 1:length(quali.sup)) {
      Xtot2[, quali.sup[j]] <- droplevels(Xtot2[, quali.sup[j]] , reorder = FALSE)
      X.quali.sup <- rbind(X.quali.sup, matrix(unlist(by(X.del,
                                                         Xtot2[, quali.sup[j]], colSums)),
                                               ncol = ncol(X.del),
                                               byrow = T))
    }
    somme.quali <- rowSums(X.quali.sup)
    X.quali.sup <- X.quali.sup/somme.quali
    coord.quali.sup <- crossprod(t(as.matrix(X.quali.sup)),V)
    dist2.quali <- rowSums(t((t(X.quali.sup) - marge.col)^2/marge.col))
    cos2.quali.sup <- coord.quali.sup^2/dist2.quali
    coord.quali.sup <- coord.quali.sup[, 1:ncp, drop = FALSE]
    cos2.quali.sup <- cos2.quali.sup[, 1:ncp, drop = FALSE]
    rownames(coord.quali.sup) <- rownames(cos2.quali.sup) <- paste(rep(colnames(Xtot2)[quali.sup],
                                                                       lapply(Xtot2[, quali.sup, drop = FALSE],
                                                                              nlevels)),
                                                                   unlist(lapply(Xtot2[, quali.sup, drop = FALSE],
                                                                                 levels)), sep = ".")
    colnames(coord.quali.sup) <- colnames(cos2.quali.sup) <- paste("Dim", 1:ncp)
    res$quali.sup <- list(coord = coord.quali.sup, cos2 = cos2.quali.sup)

    Zqs <- tab_disjonctif(Xtot2[, quali.sup])
    Nj <- colSums(Zqs * row.w)
    Nj <- colSums(Zqs * marge.row)*total
    if (total>1) coef <- sqrt(Nj * ((total - 1)/(total - Nj)))
    else coef <- sqrt(Nj)
    res$quali.sup$v.test <- res$quali.sup$coord*coef

    eta2 <- matrix(NA, length(quali.sup), ncp)
    eta2 <- sapply(as.data.frame(Xtot2[, quali.sup, drop = FALSE]),
                   fct.eta2,
                   res$row$coord,
                   weights = marge.row)
    eta2 <- t(as.matrix(eta2, ncol=ncp))
    colnames(eta2) <- paste("Dim", 1:ncp)
    rownames(eta2) <- colnames(Xtot)[quali.sup]

    res$quali.sup$eta2 <- eta2
    res$call$quali.sup <- quali.sup
  }

  class(res) <- c("CA", "list")
  if (graph & (ncp > 1)) {
    print(plot(res, axes=axes))
    if (!is.null(quanti.sup)) {
      print(plot(res, choix="quanti.sup", axes = axes, new.plot = TRUE))
    }
  }
  return(res)
}


# FactoMineR::svd.triplet
#' @title Singular Value Decomposition of a Matrix
#'
#' @noRd
#' @keywords internal
#'
svd_triplet <- function (X,
                         row.w = NULL,
                         col.w = NULL,
                         ncp = Inf) {
  tryCatch.W.E <- function(expr) {
    W <- NULL
    w.handler <- function(w) {
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler), warning = W)
  }
  if (is.null(row.w))
    row.w <- rep(1/nrow(X), nrow(X))
  if (is.null(col.w))
    col.w <- rep(1, ncol(X))
  ncp <- min(ncp, nrow(X) - 1, ncol(X))
  row.w <- row.w/sum(row.w)
  X <- t(t(X) * sqrt(col.w)) * sqrt(row.w)
  if (ncol(X) < nrow(X)) {
    svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$val
    if (names(svd.usuelle)[[1]] == "message") {
      svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp,
                                      nv = ncp))$val
      if (names(svd.usuelle)[[1]] == "d") {
        aux <- svd.usuelle$u
        svd.usuelle$u <- svd.usuelle$v
        svd.usuelle$v <- aux
      }
      else {
        bb <- eigen(crossprod(X, X), symmetric = TRUE)
        svd.usuelle <- vector(mode = "list", length = 3)
        svd.usuelle$d[svd.usuelle$d < 0] = 0
        svd.usuelle$d <- sqrt(svd.usuelle$d)
        svd.usuelle$v <- bb$vec[, 1:ncp]
        svd.usuelle$u <- t(t(crossprod(t(X), svd.usuelle$v))/svd.usuelle$d[1:ncp])
      }
    }
    U <- svd.usuelle$u
    V <- svd.usuelle$v
    if (ncp > 1) {
      mult <- sign(as.vector(crossprod(rep(1, nrow(V)),
                                       as.matrix(V))))
      mult[mult == 0] <- 1
      U <- t(t(U) * mult)
      V <- t(t(V) * mult)
    }
    U <- U/sqrt(row.w)
    V <- V/sqrt(col.w)
  }
  else {
    svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp, nv = ncp))$val
    if (names(svd.usuelle)[[1]] == "message") {
      svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$val
      if (names(svd.usuelle)[[1]] == "d") {
        aux <- svd.usuelle$u
        svd.usuelle$u <- svd.usuelle$v
        svd.usuelle$v <- aux
      }
      else {
        bb <- eigen(crossprod(t(X), t(X)), symmetric = TRUE)
        svd.usuelle <- vector(mode = "list", length = 3)
        svd.usuelle$d[svd.usuelle$d < 0] = 0
        svd.usuelle$d <- sqrt(svd.usuelle$d)
        svd.usuelle$v <- bb$vec[, 1:ncp]
        svd.usuelle$u <- t(t(crossprod(X, svd.usuelle$v))/svd.usuelle$d[1:ncp])
      }
    }
    U <- svd.usuelle$v
    V <- svd.usuelle$u
    mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
    mult[mult == 0] <- 1
    V <- t(t(V) * mult)/sqrt(col.w)
    U <- t(t(U) * mult)/sqrt(row.w)
  }
  vs <- svd.usuelle$d[1:min(ncol(X), nrow(X) - 1)]
  num <- which(vs[1:ncp] < 1e-15)
  if (length(num) == 1) {
    U[, num] <- U[, num, drop = FALSE] * vs[num]
    V[, num] <- V[, num, drop = FALSE] * vs[num]
  }
  if (length(num) > 1) {
    U[, num] <- t(t(U[, num]) * vs[num])
    V[, num] <- t(t(V[, num]) * vs[num])
  }
  res <- list(vs = vs, U = U, V = V)
  return(res)
}


# FactoMineR::tab.disjonctif
#' @title Make a disjonctif table
#'
#' @noRd
#'
tab_disjonctif <- function (tab) {
  tab <- as.data.frame(tab)
  modalite.disjonctif <- function(i) {
    moda <- as.factor(tab[, i])
    n <- length(moda)
    x <- matrix(0L, n, nlevels(moda))
    x[(1:n) + n * (unclass(moda) - 1L)] <- 1L
    return(x)
  }
  if (ncol(tab) == 1) {
    res <- modalite.disjonctif(1)
    dimnames(res) <- list(attributes(tab)$row.names, levels(tab[,
                                                                1]))
  }
  else {
    variable <- rep(attributes(tab)$names, lapply(tab, nlevels))
    listModa <- unlist(lapply(tab, levels))
    wlistModa <- which((listModa) %in% c("y", "n", "Y",
                                         "N"))
    if (!is.null(wlistModa))
      listModa[wlistModa] <- paste(variable[wlistModa],
                                   listModa[wlistModa], sep = ".")
    numlistModa <- which(unlist(lapply(listModa, is.numeric)))
    if (!is.null(numlistModa))
      listModa[numlistModa] <- paste(variable[numlistModa],
                                     listModa[numlistModa], sep = ".")
    res <- lapply(1:ncol(tab), modalite.disjonctif)
    res <- as.matrix(data.frame(res, check.names = FALSE))
    dimnames(res) <- list(attributes(tab)$row.names, listModa)
  }
  return(res)
}
