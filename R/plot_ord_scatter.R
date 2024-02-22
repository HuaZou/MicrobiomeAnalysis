#' @title Visualization of Ordination results with scatterplot
#'
#' @description
#' Show the result of Ordination using scatterplot.
#'
#' @author  Created by Hua Zou (8/9/2022 Shenzhen China)
#'
#' @param reslist (Required). list. Results of Ordination.
#' @param variable (Optional). character. the variable for x-axis.
#' @param variable_name (Optional). character. variable' names (default: NULL).
#' @param variable_color (Optional). character. the color for plotting (default: NULL).
#' @param var_shape (Optional). character. the variable for shape' column (default: NULL).
#' @param var_shape_name (Optional). character. the shape' names (default: NULL).
#' @param var_shape_value (Optional). character. the shape values (default: NULL).
#' @param display_test (Optional). logical. whether to show pvalue of PERMANOVA &
#' beta dispersion  (default: FALSE).
#' @param sample_label (Optional). logical. whether to show sample names (default: FALSE).
#' @param ellipse_type (Optional). character. how to show scatter plot, including
#' "none", "ellipse", "ellipse_CI" , "ellipse_groups" and "ellipse_line"
#' (default: "none").
#' @param sideboxplot (Optional). logical. whether to show side boxplot of axis
#'  (default: FALSE).
#' @param point_size (Optional). numeric. point size of the scatterplot (default: 3).
#' @param line_size (Optional). numeric. line size (default: 0.3).
#' @param ellipse_line_size (Optional). numeric. ellipse line size (default: 0.5).
#' @param geom_text_size (Optional). numeric. geom: text size (default: 5).
#' @param theme_text_size (Optional). numeric. main theme: text size (default: 8).
#' @param theme_strip_size (Optional). numeric. main theme: strip size (default: 6).
#' @param legend_position (Optional). numeric. main legend: position (default: c(0, 0)).
#' @param legend_justification (Optional). numeric. main legend: justification (default: c(-0.02, -0.02)).
#' @param legend_title_size (Optional). numeric. main legend: title size (default: 7).
#' @param legend_text_size (Optional). numeric. main legend: text size (default: 6).
#' @param test_annotate_size (Optional). numeric. PERMANOVA text size (default: 4).
#' @param geom_label_repel_size (Optional). numeric. main geom: label repel size (default: 2).
#' @param ... additional parameters.
#'
#' @return
#' A scatterplot with boxplot to show the results of ordination
#'
#' @importFrom tibble column_to_rownames column_to_rownames
#' @import ggplot2
#'
#' @export
#'
#' @usage plot_ord(
#'     reslist,
#'     variable,
#'     variable_name = NULL,
#'     variable_color = NULL,
#'     var_shape = NULL,
#'     var_shape_name = NULL,
#'     var_shape_value = NULL,
#'     display_test = FALSE,
#'     sample_label = FALSE,
#'     ellipse_type = c("none", "ellipse", "ellipse_CI" ,
#'                      "ellipse_groups", "ellipse_line"),
#'     sideboxplot = FALSE,
#'     point_size = 3,
#'     line_size = 0.3,
#'     ellipse_line_size = 0.5,
#'     geom_text_size = 5,
#'     theme_text_size = 8,
#'     theme_strip_size = 6,
#'     legend_position = c(0, 0),
#'     legend_justification = c(-0.02, -0.02),
#'     legend_title_size = 7,
#'     legend_text_size = 6,
#'     test_annotate_size = 4,
#'     geom_label_repel_size = 2,
#'     ...)
#' @examples
#'
#' \donttest{
#'
#' # phyloseq object
#' data("Zeybel_2022_gut")
#' ps_zeybel <- summarize_taxa(Zeybel_2022_gut, level = "Genus")
#' ord_result <- run_ord(
#'   object = ps_zeybel,
#'   variable = "LiverFatClass",
#'   method = "PCoA")
#'
#'  pl <- plot_ord(
#'   reslist = ord_result,
#'   variable = "LiverFatClass",
#'   ellipse_type = "ellipse",
#'   sideboxplot = TRUE,
#'   sample_label = TRUE
#'   )
#'  pl
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
#'  pl <- plot_ord(
#'   reslist = ord_result,
#'   variable = "LiverFatClass",
#'   ellipse_type = "ellipse",
#'   sideboxplot = TRUE,
#'   sample_label = TRUE
#'   )
#'  pl
#'
#' }
#'
plot_ord <- function(
    reslist,
    variable,
    variable_name = NULL,
    variable_color = NULL,
    var_shape = NULL,
    var_shape_name = NULL,
    var_shape_value = NULL,
    display_test = FALSE,
    sample_label = FALSE,
    ellipse_type = "none",
    sideboxplot = FALSE,
    point_size = 3,
    line_size = 0.3,
    ellipse_line_size = 0.5,
    geom_text_size = 5,
    theme_text_size = 8,
    theme_strip_size = 6,
    legend_position = c(0, 0),
    legend_justification = c(-0.02, -0.02),
    legend_title_size = 7,
    legend_text_size = 6,
    test_annotate_size = 4,
    geom_label_repel_size = 2,
    ...) {


  # data("Zeybel_2022_gut")
  # ps_zeybel <- summarize_taxa(Zeybel_2022_gut, level = "Genus")
  # ord_result <- run_ord(
  #   object = ps_zeybel,
  #   variable = "LiverFatClass",
  #   method = "PCoA")
  # reslist = ord_result
  # variable = "LiverFatClass"
  # variable_name = NULL
  # variable_color = NULL
  # var_shape = NULL
  # var_shape_name = NULL
  # var_shape_value = NULL
  # display_test = FALSE
  # sample_label = FALSE
  # ellipse_type = "ellipse_groups"
  # sideboxplot = FALSE
  # point_size = 3
  # line_size = 0.3
  # ellipse_line_size = 0.5
  # geom_text_size = 5
  # theme_text_size = 8
  # theme_strip_size = 6
  # legend_position = c(0, 0)
  # legend_justification = c(-0.02, -0.02)
  # legend_title_size = 7
  # legend_text_size = 6
  # test_annotate_size = 4
  # geom_label_repel_size = 2

  if ("dat" %in% names(reslist)) {
    plotdata <- reslist$dat
  } else {
    stop("There is no plotdata, please check your inputdata")
  }

  ellipse_type <- match.arg(
    ellipse_type,
    c("none", "ellipse", "ellipse_CI", "ellipse_line", "ellipse_groups")
  )

  # variable
  if (!is.null(variable)) {
    colnames(plotdata)[which(colnames(plotdata) == variable)] <- "Compvar"
  } else {
    stop("There is no variable information, please check your inputdata")
  }

  # factorize Compvar
  if (!is.null(variable_name)) {
    plotdata$Compvar <- factor(plotdata$Compvar, levels = variable_name)
  } else {
    plotdata$Compvar <- factor(plotdata$Compvar, levels = unique(plotdata$Compvar))
    variable_name <- unique(plotdata$Compvar)
  }
  nlevels <- levels(plotdata$Compvar)

  # group color
  if (is.null(variable_color)) {
    gg_color_hue <- function(x) {
      hues <- seq(15, 375, length = x + 1)
      hcl(h = hues, l = 65, c = 100)[1:x]
    }
    group_cols <- gg_color_hue(length(unique(plotdata$Compvar)))
  } else {
    group_cols <- variable_color
  }

  # shape variable
  if (!is.null(var_shape)) {
    if (var_shape %in% colnames(plotdata)) {
      colnames(plotdata)[which(colnames(plotdata) == var_shape)] <- "ShapeVar"

      if (!is.null(var_shape_name)) {
        plotdata$ShapeVar <- factor(plotdata$ShapeVar, levels = var_shape_name)
      } else {
        var_shape_name <- unique(plotdata$ShapeVar)
        plotdata$ShapeVar <- factor(plotdata$ShapeVar, levels = var_shape_name)
      }
      # shape_value
      if (is.null(var_shape_value)) {
        var_shape_value <- c(6:(6+length(var_shape_name)))
      }
    } else {
      stop("There is no group information, please check your inputdata")
    }
  }

  # PERMANOVA
  if (display_test) {
    if ("PERMANOVA" %in% names(reslist)) {
      per_label <- reslist$PERMANOVA
      beta_label <- reslist$BETADISPER
    } else {
      stop("There is no PERMANOVA's result")
    }
  }

  # Axis label
  if ("explains" %in% names(reslist)) {
    xlabel <- reslist$explains[1]
    ylabel <- reslist$explains[2]
  } else {
    stop("There is no Axis' result")
  }

  # main plot
  if (is.null(var_shape)) {
    pmain <- ggplot(data = plotdata, aes(x = Axis1, y = Axis2)) +
      geom_point(aes(fill = Compvar),
                 size = point_size, shape = 21, stroke = .2, color = "black") +
      geom_hline(yintercept = 0, linetype = 1, linewidth = line_size, alpha = .5) +
      geom_vline(xintercept = 0, linetype = 1, linewidth = line_size, alpha = .5) +
      scale_fill_manual(name = variable,
                        values = group_cols,
                        labels = variable_name) +
      guides(fill = guide_legend(order = 1))
  } else {
    pmain <- ggplot(data = plotdata, aes(x = Axis1, y = Axis2, color = Compvar)) +
      geom_point(aes(shape = ShapeVar), size = point_size) +
      geom_hline(yintercept = 0, linetype = 1, linewidth = line_size, alpha = .5) +
      geom_vline(xintercept = 0, linetype = 1, linewidth = line_size, alpha = .5) +
      scale_color_manual(name = variable,
                        values = group_cols,
                        labels = variable_name) +
      guides(color = guide_legend(order = 1)) +
      scale_shape_manual(name = var_shape, values = var_shape_value) +
      guides(shape = guide_legend(order = 2))
  }

  # confidence interval
  if (ellipse_type == "ellipse") {
    p_circle <- pmain +
      stat_ellipse(aes(fill = Compvar), level = 0.95,
                   linetype = 2, size = ellipse_line_size, show.legend = FALSE) +
      scale_fill_manual(name = variable,
                        values = group_cols)
  } else if (ellipse_type == "ellipse_CI") {
    p_circle <- pmain +
      stat_ellipse(aes(fill = Compvar), geom = "polygon", level = 0.95,
                   size = ellipse_line_size, alpha = 0.2, show.legend = FALSE) +
      scale_fill_manual(name = variable,
                        values = group_cols)
  } else if (ellipse_type == "ellipse_line") {
    Compvar_position <- cbind(Axis1=tapply(plotdata$Axis1, plotdata$Compvar, mean),
                              Axis2=tapply(plotdata$Axis2, plotdata$Compvar, mean)) %>%
      data.frame() %>%
      rownames_to_column("Compvar")

    border_fun <- function(x) {
      dat <- plotdata[plotdata$Compvar %in% x, , F]
      res <- dat[chull(dat[[2]], dat[[3]]), ]
      return(res)
    }
    Compvar_border <- NULL
    for(i in 1:length(unique(plotdata$Compvar))) {
      Compvar_border <- rbind(Compvar_border,
                              border_fun(x = unique(plotdata$Compvar)[i]))
    }

    p_circle <- pmain +
      geom_text(data = Compvar_position,
                aes(x = Axis1, y = Axis2, label = Compvar, color = Compvar),
                size = geom_text_size) +
      geom_polygon(data = Compvar_border,
                   aes(fill = Compvar),
                   color = "black", alpha = 0.2, show.legend = FALSE) +
      scale_fill_manual(name = variable,
                        values = group_cols) +
      scale_color_manual(name = variable,
                         values = group_cols) +
      guides(color = "none", label = "none")
  } else if (ellipse_type == "ellipse_groups") {
  p_circle <- pmain +
    geom_enterotype(aes(color = Compvar, label = Compvar), show.legend = FALSE) +
    scale_fill_manual(name = variable,
                      values = group_cols) +
    scale_color_manual(name = variable,
                      values = group_cols) +
    guides(color = "none", label = "none")
  } else if (ellipse_type == "none") {
    p_circle <- pmain
  }

  p_theme <- p_circle +
    labs(x = xlabel, y = ylabel) +
    theme_bw() +
    theme(text = element_text(size = theme_text_size - 1, color = "black"),
          strip.text = element_text(size = theme_strip_size, color = "black"),
          panel.grid = element_blank(),
          # transparent legend (2023/12/1 update)
          legend.key = element_rect(color = NA, fill = NA),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent", color = NA),

          legend.position = legend_position,
          legend.justification = legend_justification,
          legend.title = element_text(size = legend_title_size, color = "black"),
          legend.text = element_text(size = legend_text_size, face = "bold", color = "black"))


  # PERMANOVA
  if (display_test) {
    p_test <- p_theme +
      annotate("text",
               x = max(plotdata$Axis1) - max(plotdata$Axis1) / 4,
               y = min(plotdata$Axis1) + min(plotdata$Axis1) / 4,
               label = per_label, size = test_annotate_size) +
      annotate("text",
               x = max(plotdata$Axis1) - max(plotdata$Axis1),
               y = min(plotdata$Axis1) + min(plotdata$Axis1) / 4,
               label = beta_label, size = test_annotate_size)
  } else {
    p_test <- p_theme
  }

  # Sample names
  if (sample_label) {
    p_sample <- p_test +
      ggrepel::geom_label_repel(aes(label = TempRowNames,
                                    color = Compvar),
                                size = geom_label_repel_size,
                                max.overlaps = Inf) +
      scale_color_manual(name = variable,
                         values = group_cols) +
      guides(color = "none")
  } else {
    p_sample <- p_test
  }

  # sidelinechart for taxa abundance and PCoA axis
  if (sideboxplot) {
    # comparisons
    cmp <- list()
    num <- utils::combn(length(unique(plotdata$Compvar)), 2)
    for (i in 1:ncol(num)) {
      cmp[[i]] <- num[, i]
    }
    xbp <- ggplot(data = plotdata, aes(x = Compvar, y = Axis2, fill = Compvar)) +
      stat_boxplot(geom = "errorbar", width = .12) +
      geom_boxplot(width = .2, outlier.shape = 3, outlier.size = 1) +
      ggpubr::stat_compare_means(comparisons = cmp,
                                 method = "wilcox.test",
                                 label = "p.signif") +
      labs(x = "", y = "") +
      scale_fill_manual(values = group_cols) +
      guides(fill = "none") +
      theme_void() +
      theme(axis.text = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())

    ybp <- ggplot(data = plotdata, aes(x = Compvar, y = Axis1, fill = Compvar)) +
      stat_boxplot(geom = "errorbar", width = .12) +
      geom_boxplot(width = .2, outlier.shape = 3, outlier.size = 1) +
      ggpubr::stat_compare_means(comparisons = cmp,
                                 method = "wilcox.test",
                                 label = "p.signif") +
      labs(x = "", y = "") +
      scale_fill_manual(values = group_cols) +
      guides(fill = "none") +
      coord_flip() +
      theme_void() +
      theme(axis.text = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())

    # merge three plot
    p1 <- cowplot::insert_xaxis_grob(p_sample, ybp,
                                     grid::unit(.2, "null"),
                                     position = "top")
    p2 <- cowplot::insert_yaxis_grob(p1, xbp,
                                     grid::unit(.2, "null"),
                                     position = "right")
    pl <- cowplot::ggdraw(p2)
  } else {
    pl <- p_sample
  }

  return(pl)
}


#' @references
#' https://stackoverflow.com/questions/42575769/how-to-modify-the-backgroup-color-of-label-in-the-multiple-ggproto-using-ggplot2
#'
#' multiple groups in enterotype plot
#'
#' @noRd
#'
#' @importFrom MASS cov.trob
#' @import ggplot2
#'
geom_enterotype <- function(
    mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    alpha = 0.2,
    prop = 0.5,
    lineend = "butt",
    linejoin = "round",
    linemitre = 1,
    arrow = NULL,
    na.rm = FALSE,
    parse = FALSE,
    nudge_x = 0,
    nudge_y = 0,
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.15, "lines"),
    label.size = 0.1,
    show.legend = TRUE,
    inherit.aes = TRUE,
    ...) {

  # create new stat and geom for scatterplot with ellipses
  StatEllipse <- ggproto(
    "StatEllipse",
    Stat,
    required_aes = c("x", "y"),
    compute_group = function(., data, scales, level = 0.75, segments = 51, ...) {
        #library(MASS)
        dfn <- 2
        dfd <- length(data$x) - 1
        if (dfd < 3) {
           ellipse <- rbind(c(NA, NA))
        } else {
          v <- MASS::cov.trob(cbind(data$x, data$y))
          shape <- v$cov
          center <- v$center
          radius <- sqrt(dfn * qf(level, dfn, dfd))
          angles <- (0:segments) * 2 * pi / segments
          unit.circle <- cbind(cos(angles), sin(angles))
          ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
        }
        ellipse <- as.data.frame(ellipse)
        colnames(ellipse) <- c("x", "y")

        return(ellipse)
      }
  )

  # write new ggproto
  GeomEllipse <- ggproto(
    "GeomEllipse",
    Geom,
    draw_group = function(data, panel_scales, coord) {
        n <- nrow(data)
        if (n == 1) {
          return(zeroGrob())
        }
        munched <- coord_munch(coord, data, panel_scales)
        munched <- munched[order(munched$group), ]
        first_idx <- !duplicated(munched$group)
        first_rows <- munched[first_idx, ]
        grid::pathGrob(munched$x, munched$y, default.units = "native",
                       id = munched$group,
                       gp = grid::gpar(col = first_rows$colour,
                                       fill = alpha(first_rows$fill, first_rows$alpha),
                                       lwd = first_rows$size * .pt,
                                       lty = first_rows$linetype))
        }, default_aes = aes(colour = "NA", fill = "grey20", size = 0.5,
                             linetype = 1, alpha = NA, prop = 0.5),
          handle_na = function(data, params) {data},
          required_aes = c("x", "y"),
          draw_key = draw_key_path
  )

  # create a new stat for PCA scatterplot with lines which totally directs to the center
  StatConline <- ggproto(
    "StatConline",
    Stat,
    compute_group = function(data, scales) {
        df <- data.frame(data$x, data$y)
        mat <- as.matrix(df)
        center <- MASS::cov.trob(df)$center
        names(center) <- NULL
        mat_insert <- insertRow(mat, 2, center)
        for(i in 1:nrow(mat)) {
            mat_insert <- insertRow(mat_insert, 2*i, center)
            next
        }
        mat_insert <- mat_insert[-c(2:3), ]
        rownames(mat_insert) <- NULL
        mat_insert <- as.data.frame(mat_insert, center)
        colnames(mat_insert) <- c("x", "y")
        return(mat_insert)
      }, required_aes = c("x", "y")
  )

  # create a new stat for PCA scatterplot with center labels
  StatLabel <- ggproto(
    "StatLabel",
    Stat,
    compute_group = function(data, scales) {
        df <- data.frame(data$x, data$y)
        center <- MASS::cov.trob(df)$center
        names(center) <- NULL
        center <- t(as.data.frame(center))
        center <- as.data.frame(cbind(center))
        colnames(center) <- c("x", "y")
        rownames(center) <- NULL
        return(center)
       }, required_aes = c("x", "y")
  )

  layer1 <- layer(data = data,
                  mapping = mapping,
                  stat = stat,
                  geom = GeomPoint,
                  position = position,
                  show.legend = show.legend,
                  inherit.aes = inherit.aes,
                  params = list(na.rm = na.rm, ...))

  layer2 <- layer(data = data,
                  mapping = mapping,
                  stat = StatEllipse,
                  geom = GeomEllipse,
                  position = position,
                  show.legend = FALSE,
                  inherit.aes = inherit.aes,
                  params = list(na.rm = na.rm, prop = prop, alpha = alpha, ...))

  layer3 <- layer(data = data,
                  mapping = mapping,
                  stat = StatConline,
                  geom = GeomPath,
                  position = position,
                  show.legend = show.legend,
                  inherit.aes = inherit.aes,
                  params = list(lineend = lineend, linejoin = linejoin,
                                linemitre = linemitre, arrow = arrow, na.rm = na.rm, ...))

  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      stop("Specify either `position` or `nudge_x`/`nudge_y`",
           call. = FALSE)
    }
    position <- position_nudge(nudge_x, nudge_y)
  }
  layer4 <- layer(data = data,
                  mapping = mapping,
                  stat = StatLabel,
                  geom = GeomLabel,
                  position = position,
                  show.legend = FALSE,
                  inherit.aes = inherit.aes,
                  params = list(parse = parse, label.padding = label.padding,
                                label.r = label.r, label.size = label.size, na.rm = na.rm, ...))

  return(list(layer1, layer2, layer3, layer4))
}
