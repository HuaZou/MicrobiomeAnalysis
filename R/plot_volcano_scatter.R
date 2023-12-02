#' @title plot the differential analysis results by volcano plot
#'
#' @description
#' The function is to use ggplot2 to draw volcano.
#'
#' @author  Created by Hua Zou (12/3/2022 Shenzhen China)
#'
#' @param da_res (Required). a data frame from differential analysis.
#' @param group_names (Optional). character. the two groups' names (default: NULL).
#' @param x_index (Required). character. the variable for x-axis.
#' @param x_index_cutoff (Optional). numeric. the cutoff of x-axis (default: 1).
#' @param rev_x_index (Optional). logical. whether to reverse the x index.
#' The effectsize, logFC or coefficients are controversial in `aldex2`, `ancom`,
#' `edgeR`, `limma_voom`, `masslin2` and `mbzinb`. For example,
#' the negative value indicates group1 and the positive indicates group2
#' when the order of group names is group1, group2. In fact, we expect that
#' the negative value indicates group2 and the positive indicates group1.
#'  (default: FALSE).
#' @param y_index (Required). character. the variable for y-axis.
#' @param y_index_cutoff (Optional). numeric. the cutoff of y-axis (default: 0.05).
#' @param group_colors (Optional). character. the color for plotting.
#' (default: upregulated->"red", downregulated->"blue", nonsignificant->"grey").
#' @param topN (Optional). integer. the number of top significant taxa to
#' plot (default: NULL).
#' @param feature_name (Optional). character. the interesting significant
#' feature to show(default: NULL).
#' @param p_size (Optional). numeric. point size of the scatterplot (default: 1).
#' @param line_size (Optional). numeric. line size (default: 0.7).
#' @param geom_text_size (Optional). numeric. geom: text size (default: 6).
#' @param geom_label_repel_size (Optional). numeric. main geom: label
#' repel size (default: 4).
#' @param theme_text_size (Optional). numeric. main theme: text size (default: 10).
#' @param theme_title_size (Optional). numeric. main theme: title size (default: 12).
#' @param theme_legend_size (Optional). numeric. main legend: text size (default: 8).
#' @param theme_strip_size (Optional). numeric. main strip: text size (default: 14).
#' @param add_enrich_arrow (Optional). logical. whether to add enriched
#' arrow. (default: FALSE).
#' @param x_gap_times (Optional). numeric. times of x axis gap (default: 1.2).
#' @param segment_size (Optional). numeric. arrow segment size (default: 8).
#' @param segment_text_size (Optional). numeric. size of text in arrow
#' segment (default: 4).
#'
#' @import ggplot2
#' @importFrom dplyr %>% select filter all_of inner_join
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @usage plot_volcano(
#'     da_res,
#'     group_names = NULL,
#'     x_index,
#'     x_index_cutoff = 1,
#'     rev_x_index = FALSE,
#'     y_index,
#'     y_index_cutoff = 0.05,
#'     group_colors = c("red", "grey", "blue"),
#'     topN = 5,
#'     feature_name = NULL,
#'     p_size = 1,
#'     line_size = 0.7,
#'     geom_text_size = 5,
#'     geom_label_repel_size = 4,
#'     theme_text_size = 10,
#'     theme_title_size = 12,
#'     theme_legend_size = 8,
#'     theme_strip_size = 14,
#'     add_enrich_arrow = FALSE,
#'     x_gap_times = 1.2,
#'     segment_size = 8,
#'     segment_text_size = 4)
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' data("Zeybel_2022_protein")
#' Zeybel_2022_protein_imput <- impute_abundance(
#'       Zeybel_2022_protein,
#'       group = "LiverFatClass",
#'       method = "knn")
#' Zeybel_2022_protein_norm <- scale_variables(
#'       Zeybel_2022_protein_imput,
#'       method == "zscore")
#' DA_results <- run_metabolomeDA(
#'   object_raw = Zeybel_2022_protein,
#'   object_norm = Zeybel_2022_protein_norm,
#'   variable = "LiverFatClass",
#'   variable_name = c("None", "Severe"))
#'
#' plot_volcano(
#'    da_res = DA_results,
#'    group_names = c("None", "Severe"),
#'    x_index = "Log2FoldChange",
#'    x_index_cutoff = 0.1,
#'    rev_x_index = TRUE,
#'    y_index = "AdjustedPvalue",
#'    y_index_cutoff = 0.5,
#'    group_color = c("red", "grey", "blue"),
#'    topN = 5)

#' }
#'
plot_volcano <- function(
    da_res,
    group_names = NULL,
    x_index,
    x_index_cutoff = 1,
    rev_x_index = FALSE,
    y_index,
    y_index_cutoff = 0.05,
    group_colors = c("red", "grey", "blue"),
    topN = 5,
    feature_name = NULL,
    p_size = 1,
    line_size = 0.7,
    geom_text_size = 5,
    geom_label_repel_size = 4,
    theme_text_size = 10,
    theme_title_size = 12,
    theme_legend_size = 8,
    theme_strip_size = 14,
    add_enrich_arrow = FALSE,
    x_gap_times = 1.2,
    segment_size = 8,
    segment_text_size = 4) {

  # data("Zeybel_2022_protein")
  # Zeybel_2022_protein_imput <- impute_abundance(
  #   Zeybel_2022_protein,
  #   group = "LiverFatClass",
  #   method = "knn")
  # Zeybel_2022_protein_norm <- scale_variables(
  #   Zeybel_2022_protein_imput,
  #   method == "zscore")
  # DA_results <- run_metabolomeDA(
  #   object_raw = Zeybel_2022_protein,
  #   object_norm = Zeybel_2022_protein_norm,
  #   variable = "LiverFatClass",
  #   variable_name = c("None", "Severe"))
  #
  # da_res = DA_results
  # group_names = c("None", "Severe")
  # x_index = "Log2FoldChange"
  # x_index_cutoff = 0
  # rev_x_index = TRUE
  # y_index = "AdjustedPvalue"
  # y_index_cutoff = 0.5
  # group_colors = c("red", "grey", "blue")
  # topN = 5
  # feature_name = NULL
  # p_size = 1
  # line_size = 0.7
  # geom_text_size = 5
  # geom_label_repel_size = 4
  # theme_text_size = 10
  # theme_title_size = 12
  # theme_legend_size = 8
  # theme_strip_size = 14
  # add_enrich_arrow = TRUE
  # x_gap_times = 1.2
  # segment_size = 8
  # segment_text_size = 4


  # feature ID into taxaID
  colnames(da_res)[1] <- "TaxaID"

  # group
  if (any(colnames(da_res) %in% "Block")) {
    da_res_group_names <- gsub("\\d+_", "", unlist(strsplit(da_res$Block[1], " vs ")))
    if (is.null(group_names)) {
      group_names <- da_res_group_names
    } else {
      if (!all(group_names == da_res_group_names)) {
        message("group names are in wrong order, and reoder them")
        group_names <- da_res_group_names
      } else {
        group_names <- group_names
      }
    }
  } else {
    if (is.null(group_names)) {
      group_names <- da_res_group_names
      stop("please provide groups for comparison")
    } else {
      group_names <- group_names
    }
  }

  if (!x_index %in% colnames(da_res)) {
    stop("No x_index matched the DA results' column, please check out your inputdata")
  }
  if (!y_index %in% colnames(da_res)) {
    stop("No y_index matched the DA results' column, please check out your inputdata")
  }
  colnames(da_res)[which(colnames(da_res) == x_index)] <- "Xindex"
  colnames(da_res)[which(colnames(da_res) == y_index)] <- "Yindex"

  # ancom
  if (y_index %in% c("(W)q-values < alpha", "W_ratio")) {
    da_res <- da_res %>%
      dplyr::filter(!is.na(Yindex))
  }

  # whether to reverse the Xindex (aldex2, ancom, edgeR, limma_voom, masslin2 and mbzinb)
  if (rev_x_index) {
    da_res$Xindex <- -(da_res$Xindex)
  }

  # enrichment by new x_index_cutoff and y_index_cutoff
  da_res[which(da_res$Xindex > x_index_cutoff & da_res$Yindex < y_index_cutoff),
         "EnrichedDir"] <- group_names[1]
  da_res[which(da_res$Xindex < -x_index_cutoff & da_res$Yindex < y_index_cutoff),
         "EnrichedDir"] <- group_names[2]
  da_res[which(abs(da_res$Xindex) <= x_index_cutoff | da_res$Yindex >= y_index_cutoff),
         "EnrichedDir"] <- "Nonsignif"

  # ancom
  if (y_index %in% c("(W)q-values < alpha", "W_ratio")) {
    # enrichment by new x_index_cutoff and y_index_cutoff
    da_res[which(da_res$Xindex > x_index_cutoff & da_res$Yindex > y_index_cutoff),
           "EnrichedDir"] <- group_names[1]
    da_res[which(da_res$Xindex < -x_index_cutoff & da_res$Yindex > y_index_cutoff),
           "EnrichedDir"] <- group_names[2]
    da_res[which(abs(da_res$Xindex) <= x_index_cutoff | da_res$Yindex <= y_index_cutoff),
           "EnrichedDir"] <- "Nonsignif"
  }

  # print(table(dat$color))
  df_status <- table(da_res$EnrichedDir) %>% data.frame() %>%
    stats::setNames(c("Group", "Number"))

  grp1_number <- with(df_status, df_status[Group %in% group_names[1], "Number"])
  grp2_number <- with(df_status, df_status[Group %in% group_names[2], "Number"])
  nsf_number <- with(df_status, df_status[Group %in% "Nonsignif", "Number"])

  if (all(length(grp1_number) > 0, length(nsf_number) > 0, length(grp2_number) > 0)) {
    legend_label <- c(paste0(group_names[1], " (", grp1_number, ")"),
                      paste0("Nonsignif", " (", nsf_number, ")"),
                      paste0(group_names[2], " (", grp2_number, ")"))
  } else if (all(!length(grp1_number) > 0, length(nsf_number) > 0, length(grp2_number) > 0)) {
    legend_label <- c(paste0("Nonsignif", " (", nsf_number, ")"),
                      paste0(group_names[2], " (", grp2_number, ")"))
  } else if (all(length(grp1_number) > 0, length(nsf_number) > 0, !length(grp2_number) > 0)) {
    legend_label <- c(paste0(group_names[1], " (", grp1_number, ")"),
                      paste0("Nonsignif", " (", nsf_number, ")"))
  } else if (all(length(grp1_number) > 0, !length(nsf_number) > 0, length(grp2_number) > 0)) {
    legend_label <- c(paste0(group_names[1], " (", grp1_number, ")"),
                      paste0(group_names[2], " (", grp2_number, ")"))
  } else if (all(!length(grp1_number) > 0, length(nsf_number) > 0, !length(grp2_number) > 0)) {
    legend_label <- paste0("Nonsignif", " (", nsf_number, ")")
  } else if (all(length(grp1_number) > 0, !length(nsf_number) > 0, !length(grp2_number) > 0)) {
    legend_label <- paste0(group_names[1], " (", grp1_number, ")")
  } else if (all(!length(grp1_number) > 0, !length(nsf_number) > 0, length(grp2_number) > 0)) {
    legend_label <- paste0(group_names[2], " (", grp2_number, ")")
  }

  da_res_signif <- da_res %>%
    dplyr::arrange(Yindex, Xindex) %>%
    dplyr::filter(Yindex < y_index_cutoff) %>%
    dplyr::filter(abs(Xindex) > x_index_cutoff)

  # ancom
  if (y_index %in% c("(W)q-values < alpha", "W_ratio")) {
    da_res_signif <- da_res %>%
      dplyr::arrange(Yindex, Xindex) %>%
      dplyr::filter(Yindex > y_index_cutoff) %>%
      dplyr::filter(abs(Xindex) > x_index_cutoff) %>%
      dplyr::arrange(desc(Yindex))
  }

  if (nrow(da_res_signif) != 0) {
    if (!is.null(topN)) {
      if (nrow(da_res_signif) < topN) {
        message("warning: topN no more than total significant taxa")
        dat_top_signif <- da_res_signif
      } else if (topN > 50) {
        message("warning: Don't recommand to show more than 50 significant taxa in the volcanplot")
        dat_top_signif <- da_res_signif %>%
          dplyr::slice(1:topN)
      } else {
        dat_top_signif <- da_res_signif %>%
          dplyr::slice(1:topN)
      }

    } else {
      dat_top_signif <- NULL
    }

  } else {
    dat_top_signif <- NULL
    #stop("No significant taxa has been chosed, please reset your cutoff")
  }

  if (!is.null(feature_name)) {
    dat_taxa_show <- da_res %>%
      dplyr::filter(TaxaID %in% feature_name)
  } else {
    dat_taxa_show <- NULL
  }

  if (!is.null(dat_taxa_show)) {
    if (is.null(dat_top_signif)) {
      dat_final_tax <- dat_taxa_show
    } else {
      intersect_taxa <- dplyr::intersect(dat_taxa_show$TaxaID,
                                         dat_top_signif$TaxaID)

      if (length(intersect_taxa) == 0) {
        dat_final_tax <- rbind(dat_top_signif,
                               dat_taxa_show)
      } else {
        dat_final_tax <- dat_top_signif
      }
    }
  } else {
    if (is.null(dat_top_signif)) {
      dat_final_tax <- NULL
    } else {
      dat_final_tax <- dat_top_signif
    }
  }

  if (!is.null(group_colors)) {
    plot_group_color <- group_colors
    names(plot_group_color) <- c(group_names[1], "Nonsignif", group_names[2])
  } else {
    plot_group_color <- c("red", "grey", "blue")
    names(plot_group_color) <- c(group_names[1], "Nonsignif", group_names[2])
  }
  # 2023/5/17 update for legend label
  if (all(length(grp1_number) > 0, length(nsf_number) > 0, length(grp2_number) > 0)) {
    plot_group_color <- plot_group_color
  } else if (all(!length(grp1_number) > 0, length(nsf_number) > 0, length(grp2_number) > 0)) {
    plot_group_color <- plot_group_color[c(2, 3)]
  } else if (all(length(grp1_number) > 0, length(nsf_number) > 0, !length(grp2_number) > 0)) {
    plot_group_color <- plot_group_color[c(1, 2)]
  } else if (all(length(grp1_number) > 0, !length(nsf_number) > 0, length(grp2_number) > 0)) {
    plot_group_color <- plot_group_color[c(1, 3)]
  } else if (all(!length(grp1_number) > 0, length(nsf_number) > 0, !length(grp2_number) > 0)) {
    plot_group_color <- plot_group_color[2]
  } else if (all(length(grp1_number) > 0, !length(nsf_number) > 0, !length(grp2_number) > 0)) {
    plot_group_color <- plot_group_color[1]
  } else if (all(!length(grp1_number) > 0, !length(nsf_number) > 0, length(grp2_number) > 0)) {
    plot_group_color <- plot_group_color[3]
  }

  if (x_index == "logFC"){
    xlabel <- paste0("Log2FoldChange (", paste(group_names, collapse = "/"), ")")
  } else {
    xlabel <- paste0(x_index, " (", paste(group_names, collapse = "/"), ")")
  }

  if (y_index == "AdjustedPvalue") {
    ylable <- expression(-log[10]("Adjusted p-value"))
  } else if (y_index == "Pvalue") {
    ylable <- expression(-log[10]("P-value"))
  } else {
    ylable <- y_index
  }

  # ancom
  if (y_index %in% c("(W)q-values < alpha", "W_ratio")) {
    xlabel <- paste0("CLR mean difference (", paste(group_names, collapse = "/"), ")")
    if (length(grep("alpha", y_index)) == 1) {
      ylable <- "W statistic" # (W)q-values < alpha
    } else if (length(grep("W_ratio", y_index)) == 1) {
      ylable <- "Ratio of W statistic" # W_ratio
    } else {
      stop("No column matched the DA result")
    }
  }

  da_res$EnrichedDir <- factor(da_res$EnrichedDir, levels = names(plot_group_color))

  if (y_index == "AdjustedPvalue") {
    p <- ggplot(da_res, aes(x = Xindex, y = -log10(Yindex), color = EnrichedDir)) +
      geom_point(size = p_size, alpha = 1, stroke = 1)
  } else if (y_index == "Pvalue") {
    p <- ggplot(da_res, aes(x = Xindex, y = -log10(Yindex), color = EnrichedDir)) +
      geom_point(size = p_size, alpha = 1, stroke = 1)
  } else {
    p <- ggplot(da_res, aes(x = Xindex, y = Yindex, color = EnrichedDir)) +
      geom_point(size = p_size, alpha = 1, stroke = 1)
  }

  p1 <- p +
    scale_color_manual(name = NULL,
                       values = plot_group_color,
                       labels = legend_label) +
    xlab(xlabel) +
    ylab(ylable) +
    geom_vline(xintercept = x_index_cutoff, alpha = .8, linetype = 2, linewidth = line_size) +
    geom_vline(xintercept = -x_index_cutoff, alpha = .8, linetype = 2, linewidth = line_size) +
    annotate("text",
             x = x_index_cutoff, y = 0,
             label = x_index_cutoff,
             size = geom_text_size, color = "red") +
    annotate("text", x = -x_index_cutoff, y = 0,
             label = -x_index_cutoff,
             size = geom_text_size, color = "red") +
    scale_x_continuous(position = "top") +
    # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + # 2023/5/18 update
    guides(color = guide_legend(override.aes = list(size = geom_text_size - 3)))

  if (y_index == "AdjustedPvalue") {
    p2 <- p1 +
      scale_y_continuous(trans = "log1p") +
      geom_hline(yintercept = -log10(y_index_cutoff), alpha = .8, linetype = 2, linewidth = line_size) +
      annotate("text",
               x = min(da_res$Xindex),
               y = -log10(y_index_cutoff),
               label = y_index_cutoff,
               size = geom_text_size, color = "red")
  } else if (y_index == "Pvalue") {
    p2 <- p1 +
      scale_y_continuous(trans = "log1p") +
      geom_hline(yintercept = -log10(y_index_cutoff), alpha = .8, linetype = 2, linewidth = line_size) +
      annotate("text",
               x = min(da_res$Xindex),
               y = -log10(y_index_cutoff),
               label = y_index_cutoff,
               size = geom_text_size, color = "red")
  } else {
    p2 <- p1 +
      geom_hline(yintercept = y_index_cutoff, alpha = .8, linetype = 2, linewidth = line_size) +
      annotate("text",
               x = min(da_res$Xindex),
               y = y_index_cutoff,
               label = y_index_cutoff,
               size = geom_text_size, color = "red")
  }

  if (!is.null(dat_final_tax)) {
    p3 <- p2 +
      ggrepel::geom_text_repel(data = dat_final_tax,
                               aes(label = TaxaID),
                               size = geom_label_repel_size,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 80),
                               segment.linetype = 1,
                               segment.curvature = -1e-20,
                               box.padding = unit(0.35, "lines"),
                               point.padding = unit(0.3, "lines"),
                               arrow = arrow(length = unit(0.005, "npc")),
                               bg.r = 0.15)
  } else {
    p3 <- p2
  }

  p4 <- p3 +
    theme_bw() +
    theme(axis.title = element_text(size = theme_title_size, face = "bold", color = "black"),
          axis.text = element_text(size = theme_text_size, color = "black"),
          text = element_text(size = theme_text_size - 2, color = "black"),
          #panel.grid = element_blank(),
          legend.position = c(.15, .1),
          legend.key.height = unit(0.6, "cm"),
          # transparent legend (2023/11/30 update)
          legend.key = element_rect(color = NA, fill = NA),
          legend.text = element_text(size = theme_legend_size, face = "bold", color = "black"),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent", color = NA),
          strip.text = element_text(size = theme_strip_size, face = "bold"))

  if (add_enrich_arrow) {
    # add arrow segment for enriched direction (2023/11/30 update)
    x_gap <- (max(da_res$Xindex) - min(da_res$Xindex)) / 3
    x_left_end <- min(da_res$Xindex)
    x_left_start <- min(da_res$Xindex) + x_gap * x_gap_times
    x_right_end <- max(da_res$Xindex)
    x_right_start <- max(da_res$Xindex) - x_gap * x_gap_times
    y_h <- mean(da_res$Yindex)

    x_left_mean <- (x_left_end + x_left_start) / 2
    x_right_mean <- (x_right_end + x_right_start) / 2

    enriched_group1 <- paste0("Enriched in ", group_names[1])
    enriched_group2 <- paste0("Enriched in ", group_names[2])

    xlabel_pl <- ggplot(da_res, aes(x = Xindex, y = Yindex)) +
      annotate("segment", x = x_right_start, y = y_h, xend = x_right_end, yend = y_h,
               size = segment_size, linejoin = "mitre",
               arrow = arrow(type = "closed", length = unit(0.02, "npc")),
               color = group_colors[1]) +
      annotate("text", x = x_right_mean, y = y_h, label = enriched_group1, color = "white",
               size = segment_text_size, fontface = "bold") +
      annotate("segment", x = x_left_start, y = y_h, xend = x_left_end, yend = y_h,
               size = segment_size, linejoin = "mitre",
               arrow = arrow(type = "closed", length = unit(0.02, "npc")),
               color = group_colors[3]) +
      annotate("text", x = x_left_mean, y = y_h, label = enriched_group2, color = "white",
               size = segment_text_size, fontface = "bold") +
      theme_void() +
      theme(axis.text = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())

    p5 <- cowplot::insert_xaxis_grob(
      p4, xlabel_pl,
      grid::unit(0.2, "null"),
      position = "bottom")

    pl_final <- cowplot::ggdraw(p5)
  } else {
    pl_final <- p4 + scale_x_continuous(position = "bottom")
  }

  return(pl_final)
}
