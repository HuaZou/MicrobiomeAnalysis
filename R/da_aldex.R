#' @title Perform differential analysis using ALDEx2
#'
#' @description
#' The `run_aldex` function is used to identify the significant features.
#'
#' @details
#' The `run_aldex` function is used to identify the significant features.
#' It can be applied to both [`phyloseq::phyloseq-class`]
#' and [`Biobase::ExpressionSet`] object.
#'
#' @references Fernandes, A.D., Reid, J.N., Macklaim, J.M. et al. Unifying the
#'   analysis of high-throughput sequencing datasets: characterizing RNA-seq,
#'   16S rRNA gene sequencing and selective growth experiments by compositional
#'   data analysis. Microbiome 2, 15 (2014).
#'
#' @author Created by Yang Cao; modified by Hua Zou (5/17/2022 Shenzhen China)
#'
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param group character, the variable to set the group
#' @param taxa_rank character to specify taxonomic rank to perform
#'   differential analysis on. Should be one of
#'   `phyloseq::rank_names(phyloseq)`, or "all" means to summarize the taxa by
#'   the top taxa ranks (`summarize_taxa(ps, level = rank_names(ps)[1])`), or
#'   "none" means perform differential analysis on the original taxa
#'   (`taxa_names(phyloseq)`, e.g., OTU or ASV).
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
#'   * "TSS": total sum scaling, also referred to as "relative abundance", the
#'     abundances were normalized by dividing the corresponding sample library
#'     size.
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
#'   * "CLR": centered log-ratio normalization.
#'   * "CPM": pre-sample normalization of the sum of the values to 1e+06.
#' @param norm_para arguments passed to specific normalization methods
#' @param method test method, options include: "t.test" and "wilcox.test"
#'   for two groups comparison,  "kruskal" and "glm_anova" for multiple groups
#'   comparison.
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param pvalue_cutoff cutoff of p value, default 0.05.
#' @param mc_samples integer, the number of Monte Carlo samples to use for
#'   underlying distributions estimation, 128 is usually sufficient.
#' @param denom character string, specifiy which features used to as the
#'   denominator for the geometric mean calculation. Options are:
#'   * "all", with all features.
#'   * "iqlr", accounts for data with systematic variation and centers the
#'    features on the set features that have variance that is between the lower
#'    and upper quartile of variance.
#'   *  "zero", a more extreme case where there are many non-zero features in
#'     one condition but many zeros in another. In this case the geometric mean
#'     of each group is calculated using the set of per-group non-zero features.
#'   * "lvha", with house keeping features.
#' @param paired logical, whether to perform paired tests, only worked for
#'   method "t.test" and "wilcox.test".
#'
#' @return a [`microbiomeMarker-class`] object.
#'
#' @export
#'
#' @seealso [`ALDEx2::aldex()`]
#'
#' @examples
#' \dontrun{
#' data(caporaso)
#' ps <- phyloseq::subset_samples(
#'     caporaso,
#'     SampleType %in% c("gut", "tongue", "right palm")
#' )
#' run_aldex(ps, group = "SampleType", method = "kruskal")
#' }
#'
run_aldex <- function(
    ps,
    group,
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p",
                  "SquareRoot", "CubicRoot", "logit"),
    norm = "none",
    norm_para = list(),
    method = c("t.test", "wilcox.test",
               "kruskal", "glm_anova"),
    p_adjust = c("none", "fdr", "bonferroni", "holm",
                 "hochberg", "hommel", "BH", "BY"),
    pvalue_cutoff = 0.05,
    mc_samples = 128,
    denom = c("all", "iqlr", "zero", "lvha"),
    paired = FALSE) {

  # ps = phyloseq::subset_samples(
  #   caporaso,
  #   SampleType %in% c("gut", "tongue", "right palm"))
  # group = "SampleType"
  # taxa_rank = "Genus"
  # transform = "identity"
  # norm = "none"
  # norm_para = list()
  # method = "kruskal"
  # p_adjust = "none"
  # pvalue_cutoff = 0.05
  # mc_samples = 128
  # denom = "all"
  # paired = FALSE

  stopifnot(inherits(ps, "phyloseq"))
  ps <- check_rank_names(ps) %>%
    check_taxa_rank(taxa_rank)

  denom <- match.arg(denom, c("all", "iqlr", "zero", "lvha"))
  p_adjust <- match.arg(
    p_adjust,
    c(
      "none", "fdr", "bonferroni", "holm",
      "hochberg", "hommel", "BH", "BY"
    )
  )

  # trans method as argument test in ALDEx2::aldex
  method <- match.arg(
    method,
    c("t.test", "wilcox.test", "kruskal", "glm_anova")
  )
  if (method %in% c("t.test", "wilcox.test")) {
    test <- "t"
  } else {
    test <- "kw"
  }

  # check whether group is valid, write a function
  sample_meta <- sample_data(ps)
  meta_nms <- names(sample_meta)
  if (!group %in% meta_nms) {
    stop(
      group, " are not contained in the `sample_data` of `ps`",
      call. = FALSE
    )
  }

  transform <- match.arg(
    transform, c("identity", "log10", "log10p",
                 "SquareRoot", "CubicRoot", "logit")
  )

  # preprocess phyloseq object
  ps <- preprocess_ps(ps)
  ps <- transform_abundances(ps, transform = transform)

  # normalize the data
  norm_para <- c(norm_para, method = norm, object = list(ps))
  ps_normed <- do.call(normalize, norm_para)
  ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
  groups <- sample_meta[[group]]
  abd <- abundances(ps_summarized, norm = TRUE)

  test_fun <- ifelse(test == "t", aldex_t, aldex_kw)
  test_para <- list(
    reads = abd,
    conditions = groups,
    method = method,
    mc_samples = mc_samples,
    denom = denom,
    p_adjust = p_adjust
  )
  if (test == "t") {
    test_para <- c(test_para, paired = paired)
  }

# test_fun(reads = abd,
#          conditions = groups,
#          method = method,
#          mc_samples = mc_samples,
#          denom = denom,
#          p_adjust = p_adjust)

  test_out <- tryCatch(
    do.call(test_fun, test_para),
    error = function(e) e
  )

  # check whether counts are integers
  if (inherits(test_out, "error") &&
      conditionMessage(test_out) == "not all reads are integers") {
    warning(
      "Not all reads are integers, the reads are ceiled to integers.\n",
      "   Raw reads is recommended from the ALDEx2 paper.",
      call. = FALSE
    )
    test_para$reads <- ceiling(abd)
    test_out <- do.call(test_fun, test_para)
  }

  sig_feature <- dplyr::filter(test_out, .data$padj <= pvalue_cutoff)
  marker <- return_marker(sig_feature, test_out)

  feature <- test_out$feature
  tax <- matrix(feature) %>%
    tax_table()
  row.names(tax) <- row.names(abd)

  mm <- microbiomeMarker(
    marker_table = marker,
    norm_method = get_norm_method(norm),
    diff_method = paste0("ALDEx2_", method),
    sam_data = sample_data(ps_summarized),
    otu_table = otu_table(abd, taxa_are_rows = TRUE),
    tax_table = tax
  )

  mm
}

# aldex t test, wilcox test
# In the original version of ALDEx2, each p value is corrected using the
# Benjamini-Hochberg method. Here, we add a new argument `p_adjust` to
# make aldex support for other correction methods.
aldex_t <- function(reads,
                    conditions,
                    mc_samples,
                    method = c("t.test", "wilcox.test"),
                    denom = c("all", "iqlr", "zero", "lvha"),
                    p_adjust = c(
                      "none", "fdr", "bonferroni", "holm",
                      "hochberg", "hommel", "BH", "BY"
                    ),
                    paired = FALSE) {
  method <- match.arg(method, c("t.test", "wilcox.test"))
  demon <- match.arg(denom, c("all", "iqlr", "zero", "lvha"))
  p_adjust <- match.arg(
    p_adjust,
    c(
      "none", "fdr", "bonferroni", "holm",
      "hochberg", "hommel", "BH", "BY"
    )
  )
  conditions <- as.factor(conditions)

  if (!inherits(reads, "aldex.clr")) {
    reads_clr <- ALDEx2::aldex.clr(
      reads = reads,
      conds = conditions,
      mc.samples = mc_samples,
      denom = denom
    )
    feature <- row.names(reads)
  } else {
    reads_clr <- reads
    feature <- row.names(reads@reads)
  }

  mc_instance <- reads_clr@analysisData
  mc_instance_ldf <- convert_instance(mc_instance, mc_samples)

  if (method == "t.test") {
    pvalue <- purrr::map_dfc(
      mc_instance_ldf,
      t_fast,
      group = conditions, paired = paired
    )
  } else {
    pvalue <- purrr::map_dfc(
      mc_instance_ldf,
      wilcox_fast,
      group = conditions, paired = paired
    )
  }

  padj <- purrr::map_dfc(pvalue, p.adjust, method = p_adjust)
  # expect value
  e_pvalue <- rowMeans(pvalue)
  e_padj <- rowMeans(padj)

  # effect size
  ef <- ALDEx2::aldex.effect(
    reads_clr,
    include.sample.summary = FALSE,
    verbose = FALSE
  )
  # enrich group
  cds <- gsub("rab.win.", "", names(ef)[2:3])
  ef <- ef$effect
  enrich_group <- ifelse(ef > 0, cds[1], cds[2])

  res <- data.frame(
    feature = feature,
    enrich_group = enrich_group,
    ef_aldex = ef,
    pvalue = e_pvalue,
    padj = e_padj
  )

  res
}

# aldex kruskal-wallis test and glm anova statistics
#' @importFrom stats kruskal.test glm drop1
aldex_kw <- function(reads,
                     conditions,
                     method = c("kruskal", "glm_anova"),
                     mc_samples = 128,
                     denom = c("all", "iqlr", "zero", "lvha"),
                     p_adjust = c(
                       "none", "fdr", "bonferroni", "holm",
                       "hochberg", "hommel", "BH", "BY"
                     )) {
  method <- match.arg(method, c("kruskal", "glm_anova"))
  demon <- match.arg(denom, c("all", "iqlr", "zero", "lvha"))
  p_adjust <- match.arg(
    p_adjust,
    c(
      "none", "fdr", "bonferroni", "holm",
      "hochberg", "hommel", "BH", "BY"
    )
  )
  conditions <- as.factor(conditions)

  if (!inherits(reads, "aldex.clr")) {
    reads_clr <- ALDEx2::aldex.clr(
      reads = reads,
      conds = conditions,
      mc.samples = mc_samples,
      denom = denom
    )
    feature <- row.names(reads)
  } else {
    reads_clr <- reads
    feature <- row.names(reads@reads)
  }

  mc_instance <- reads_clr@analysisData
  # convert mc_instance to a list of data frame, each element represents a mc
  # sample for all samples.
  mc_instance_ldf <- convert_instance(mc_instance, mc_samples)

  if (method == "kruskal") {
    pvalue <- purrr::map_dfc(
      mc_instance_ldf,
      function(x) {
        apply(
          x, 1,
          function(y) {
            stats::kruskal.test(y, g = factor(conditions))[[3]]
          }
        )
      }
    )
  } else {
    pvalue <- purrr::map_dfc(
      mc_instance_ldf,
      function(x) {
        apply(
          x, 1,
          function(y) {
            stats::glm(as.numeric(y) ~ factor(conditions)) %>%
              stats::drop1(test = "Chis") %>%
              purrr::pluck(5, 2)
          }
        )
      }
    )
  }

  padj <- purrr::map_dfc(pvalue, p.adjust, method = p_adjust)
  e_pvalue <- rowMeans(pvalue)
  e_padj <- rowMeans(padj)

  # f statistic
  ef_F_statistic <- purrr::map_dfc(
    mc_instance_ldf,
    function(x) {
      apply(
        x, 1,
        function(y) {
          summary(aov(y ~ factor(conditions)))[[1]]$`F value`[1]
        }
      )
    }
  ) %>%
    rowMeans()

  enrich_group <- get_aldex_kwglm_enrich_group(mc_instance_ldf, conditions)

  res <- data.frame(
    feature = feature,
    enrich_group = enrich_group,
    ef_F_statistic = ef_F_statistic,
    pvalue = e_pvalue,
    padj = e_padj
  )

  return(res)
}

# enriched group for kw and glm anova
get_aldex_kwglm_enrich_group <- function(mc_instance_ldf, conditions) {
  instance_split <- purrr::map(
    mc_instance_ldf,
    ~ split(data.frame(t(.x)), conditions)
  )
  instance_mean <- purrr::map(
    instance_split,
    ~ purrr::map_dfc(.x, colMeans)
  )
  instance_mean <- Reduce("+", instance_mean)
  max_idx <- apply(instance_mean, 1, which.max)
  enrich_group <- names(instance_mean)[max_idx]

  return(enrich_group)
}

# Each element of mc instances of a clr object represents all instances of a
# sample, this function convert mc instances to list data frames where each
# element represents a mc instance for all samples
convert_instance <- function(mc_instance, mc_samples) {
  mc_instance_ldf <- purrr::map(
    seq.int(mc_samples),
    function(x) {
      res <- purrr::map_dfc(mc_instance, function(y) y[, x])
      names(res) <- names(mc_instance)
      res
    }
  )

  return(mc_instance_ldf)
}


# fast test function modified from ALDEx2
#' @importFrom stats pt
t_fast <- function(x, group, paired = FALSE) {
  grp1 <- group == unique(group)[1]
  grp2 <- group == unique(group)[2]
  n1 <- sum(grp1)
  n2 <- sum(grp2)

  if (paired) {
    # Order pairs for the mt.teststat function
    if (n1 != n2) stop("Cannot pair uneven groups.")
    idx1 <- which(grp1)
    idx2 <- which(grp2)
    paired_order <- unlist(
      lapply(
        seq_along(idx1),
        function(i) c(idx1[i], idx2[i])
      )
    )

    t <- multtest::mt.teststat(
      x[, paired_order],
      as.numeric(grp1)[paired_order],
      test = "pairt",
      nonpara = "n"
    )
    df <- length(idx1) - 1
    res <- pt(abs(t), df = df, lower.tail = FALSE) * 2
  } else {
    t <- multtest::mt.teststat(x,
                               as.numeric(grp1),
                               test = "t",
                               nonpara = "n"
    )
    s1 <- apply(x[, grp1], 1, sd)
    s2 <- apply(x[, grp2], 1, sd)
    df <- ((s1^2 / n1 + s2^2 / n2)^2) / ((s1^2 / n1)^2 / (n1 - 1) +
                                           (s2^2 / n2)^2 / (n2 - 1))
    res <- pt(abs(t), df = df, lower.tail = FALSE) * 2
  }

  return(res)
}

# wilcox.fast function replaces wilcox.test
#  * runs much faster
#  * uses exact distribution for ties!
#    * this differs from ?wilcox.test
#  * optional paired test
#    * equivalent to wilcox.test(..., correct = FALSE)
#  * uses multtest
#' @importFrom stats psignrank pnorm pwilcox wilcox.test
wilcox_fast <- function(x, group, paired = FALSE) {
  grp1 <- group == unique(group)[1]
  grp2 <- group == unique(group)[2]
  n1 <- sum(grp1)
  n2 <- sum(grp2)

  # Check for ties in i-th Monte-Carlo instance
  xt <- t(x)
  if (paired) {
    any_ties <- any(
      apply(
        xt[grp1, ] - xt[grp2, ], 2,
        function(y) length(unique(y))
      ) != ncol(x) / 2
    )
  } else {
    any_ties <- any(
      apply(
        xt, 2,
        function(y) length(unique(y))
      ) != ncol(x)
    )
  }

  # Ties trigger slower, safer wilcox.test function
  if (any_ties) {
    res <- apply(
      xt, 2,
      function(i) {
        wilcox.test(
          i[grp1], i[grp2],
          paired = paired, correct = FALSE
        )$p.value
      }
    )
  }

  if (paired) {
    if (n1 != n2) stop("Cannot pair uneven groups.")
    x_diff <- xt[grp1, ] - xt[grp2, ]
    v <- apply(x_diff, 2, function(y) sum(rank(abs(y))[y > 0]))
    topscore <- (n1 * (n1 + 1)) / 2
    v_lower <- ifelse(v > topscore / 2, topscore - v, v)
    if (sum(grp1) < 50) {
      # as per wilcox test, use exact -- ASSUMES NO TIES!!
      v_p <- psignrank(v_lower, n1) * 2
      # psignrank returns non-zero for W = mean
      res <- ifelse(v_p > 1, 1, v_p)
    } else { # Use normal approximation
      v_std <- (topscore / 2 - v_lower) /
        sqrt(n1 * (n1 + 1) * (2 * n1 + 1) / 24)
      res <- pnorm(v_std, lower.tail = FALSE) * 2
    }
  } else {
    w_std <- multtest::mt.teststat(x, as.numeric(grp1), test = "wilcoxon")
    if (sum(grp1) < 50 && sum(grp2) < 50) {
      # as per wilcox test, use exact -- ASSUMES NO TIES!!
      w_var <- sqrt((n1 * n2) * (n1 + n2 + 1) / 12)
      w <- abs(w_std) * w_var + (n1 * n2) / 2
      w_p <- pwilcox(w - 1, n1, n2, lower.tail = FALSE) * 2
      # pwilcox returns non-zero for W = mean
      res <- ifelse(w_p > 1, 1, w_p)
    } else { # Use normal approximation
      res <- pnorm(abs(w_std), lower.tail = FALSE) * 2
    }
  }

  return(res)
}
