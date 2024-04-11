

EnhancedVolcano_hai <-
function (toptable, lab, x, y, selectLab = NULL,
          xlim = c(min(toptable[[x]],na.rm = TRUE) - 1.5, max(toptable[[x]], na.rm = TRUE) + 1.5), 
          ylim = c(0, max(-log10(toptable[[y]]), na.rm = TRUE) + 5), 
          xlab = bquote(~Log[2] ~ "fold change"), 
          ylab = bquote(~-Log[10] ~ italic(P)), axisLabSize = 18, 
          title = "Volcano plot", 
          caption = paste0("total = ", nrow(toptable), " variables"), 
          titleLabSize = 18, 
          captionLabSize = 14, pCutoff = 1e-05, 
          FCcutoff = 1, cutoffLineType = "longdash", cutoffLineCol = "black", 
          cutoffLineWidth = 0.4, pointSize = 3, labSize = 5, labCol = "black", 
          labFace = "plain", labhjust = 0.5, labvjust = 1.5, 
          boxedLabels = FALSE, shape = 19, shapeCustom = NULL, col = c("grey30", 
                                                                       "forestgreen", "royalblue", "red2"), 
          colCustom = NULL, colAlpha = 1/2, colGradient = NULL, colGradientBreaks = c(pCutoff, 
                                                                                      1), colGradientLabels = c("0", "1.0"), colGradientLimits = c(0, 
                                                                                                                                                   1), legendLabels = c("NS", expression(Log[2] ~ 
                                                                                                                                                                                           FC), "p-value", expression(p - value ~ and ~ log[2] ~ 
                                                                                                                                                                                                                        FC)), legendPosition = "top", legendLabSize = 14, 
          legendIconSize = 5, legendDropLevels = TRUE, encircle = NULL, 
          encircleCol = "black", encircleFill = "pink", 
          encircleAlpha = 3/4, encircleSize = 2.5, shade = NULL, shadeFill = "grey", 
          shadeAlpha = 1/2, shadeSize = 0.01, shadeBins = 2, drawConnectors = FALSE, 
          widthConnectors = 0.5, typeConnectors = "closed", endsConnectors = "first", 
          lengthConnectors = unit(0.01, "npc"), colConnectors = "grey10", 
          arrowheads = TRUE, hline = NULL, hlineType = "longdash", 
          hlineCol = "black", hlineWidth = 0.4, vline = NULL, 
          vlineType = "longdash", vlineCol = "black", vlineWidth = 0.4, 
          gridlines.major = TRUE, gridlines.minor = TRUE, border = "partial", 
          borderWidth = 0.8, borderColour = "black", raster = FALSE) 
{
  if (!is.numeric(toptable[[x]])) {
    stop(paste(x, " is not numeric!", sep = ""))
  }
  if (!is.numeric(toptable[[y]])) {
    stop(paste(y, " is not numeric!", sep = ""))
  }
  if (raster) {
    geom_point <- geom_point_rast
  }
  i <- xvals <- yvals <- Sig <- NULL
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(abs(toptable[[x]]) > FCcutoff)] <- "FC"
  toptable$Sig[(toptable[[y]] < pCutoff)] <- "P"
  toptable$Sig[(toptable[[y]] < pCutoff) & (abs(toptable[[x]]) > 
                                              FCcutoff)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig, levels = c("NS", 
                                                  "FC", "P", "FC_P"))
  if (min(toptable[[y]], na.rm = TRUE) == 0) {
    warning(paste("One or more p-values is 0.", "Converting to 10^-1 * current", 
                  "lowest non-zero p-value..."), call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(toptable[which(toptable[[y]] != 
                                                                   0), y], na.rm = TRUE) * 10^-1
  }
  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(),
                                         plot.title = element_text(angle = 0, size = titleLabSize,
                                                                   face = "bold", vjust = 0.5,hjust=0.5),
                                         axis.text.x = element_text(angle = 0,
                                                                      size = axisLabSize, 
                                                                      vjust = 1), 
                                         axis.text.y = element_text(angle = 0,
                                                                    size = axisLabSize, 
                                                                    vjust = 0.5), 
                                         axis.title = element_text(size = axisLabSize), 
                                         legend.position = legendPosition,
                                         legend.key = element_blank(),
                                         legend.key.size = unit(0.5, "cm"),
                                         legend.text = element_text(size = legendLabSize),
                                         title = element_text(size = legendLabSize),
                                         legend.title = element_blank())
  if (!is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(names(shapeCustom))), alpha = colAlpha, 
                 size = pointSize, na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom))), 
                 alpha = colAlpha, shape = shape, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(guide = TRUE)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           4) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = Sig), alpha = colAlpha, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = c(NS = shape[1], FC = shape[2], 
                                    P = shape[3], FC_P = shape[4]), 
                         labels = c(NS = legendLabels[1],
                                    FC = legendLabels[2],
                                    P = legendLabels[3],
                                    FC_P = legendLabels[4]),
                         guide = TRUE, drop = legendDropLevels)
  }
  else if (is.null(colCustom) & !is.null(shapeCustom)) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + guides(colour = guide_legend(order = 1, 
                                          override.aes = list(size = legendIconSize)), 
                    shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
        geom_point(aes(color = Sig, shape = factor(names(shapeCustom))), 
                   alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
        scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                      P = col[3], FC_P = col[4]),
                           labels = c(NS = legendLabels[1],
                                      FC = legendLabels[2], P = legendLabels[3],
                                      FC_P = legendLabels[4]), drop = legendDropLevels) + 
        scale_shape_manual(values = shapeCustom)
    }
    else {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + guides(shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
        geom_point(aes(color = Sig, shape = factor(names(shapeCustom))), 
                   alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
        scale_colour_gradient(low = colGradient[1], high = colGradient[2], 
                              limits = colGradientLimits, breaks = colGradientBreaks, 
                              labels = colGradientLabels)
      scale_shape_manual(values = shapeCustom)
    }
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + guides(colour = guide_legend(order = 1, 
                                          override.aes = list(shape = shape, size = legendIconSize))) + 
        geom_point(aes(color = Sig), alpha = colAlpha, 
                   shape = shape, size = pointSize, na.rm = TRUE) + 
        scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                      P = col[3], FC_P = col[4]),
                           labels = c(NS = legendLabels[1],
                                      FC = legendLabels[2],
                                      P = legendLabels[3],
                                      FC_P = legendLabels[4]), drop = legendDropLevels)
    }
    else {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + geom_point(aes(color = yvals), alpha = colAlpha, 
                        shape = shape, size = pointSize, na.rm = TRUE) + 
        scale_colour_gradient(low = colGradient[1], high = colGradient[2], 
                              limits = colGradientLimits, breaks = colGradientBreaks, 
                              labels = colGradientLabels)
    }
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           4) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + 
        guides(colour = guide_legend(order = 1,
                                     override.aes = list(shape = c(NS = shape[1],
                                                                        FC = shape[2],
                                                                        P = shape[3],
                                                                        FC_P = shape[4]),
                                                         size = legendIconSize))) + 
        geom_point(aes(color = Sig,
                       shape = Sig), alpha = colAlpha, size = pointSize,
                   na.rm = TRUE) + scale_color_manual(values = c(NS = col[1],
                                                                 FC = col[2],
                                                                 P = col[3],
                                                                 FC_P = col[4]),
                                                      labels = c(NS = legendLabels[1],
                                                                 FC = legendLabels[2],
                                                                 P = legendLabels[3],
                                                                 FC_P = legendLabels[4]),
                                                    drop = legendDropLevels) +
        scale_shape_manual(values = c(NS = shape[1],
                                      FC = shape[2],
                                      P = shape[3],
                                      FC_P = shape[4]),
                           guide = FALSE, drop = legendDropLevels)
      }
    else {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + geom_point(aes(color = yvals, shape = Sig), 
                        alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
        scale_colour_gradient(low = colGradient[1], high = colGradient[2], 
                              limits = colGradientLimits, breaks = colGradientBreaks, 
                              labels = colGradientLabels) + 
        scale_shape_manual(values = c(NS = shape[1],
                                      FC = shape[2],
                                      P = shape[3],
                                      FC_P = shape[4]),
                           guide = FALSE, drop = legendDropLevels)
      }
  }
  plot <- plot + xlab(xlab) + ylab(ylab) + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) + 
    geom_vline(xintercept = c(-FCcutoff,FCcutoff), 
               linetype = cutoffLineType, 
               colour = cutoffLineCol,size = cutoffLineWidth) + 
    geom_hline(yintercept = -log10(pCutoff),linetype = cutoffLineType, 
               colour = cutoffLineCol, size = cutoffLineWidth)
  plot <- plot + labs(title = title, caption = caption)
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline, linetype = vlineType, 
                              colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline), 
                              linetype = hlineType, colour = hlineCol, size = hlineWidth)
  }
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = borderColour, 
                                                     fill = NA, size = borderWidth))
  }
  else if (border == "partial") {
    plot <- plot + 
      theme(axis.line = element_line(size = borderWidth,colour = borderColour),
            panel.border = element_blank(),
            panel.background = element_blank())
  }
  else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }
  if (gridlines.major) {
    plot <- plot + theme(panel.grid.major = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor) {
    plot <- plot + theme(panel.grid.minor = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  if (!boxedLabels) {
    if (drawConnectors && is.null(selectLab)) {
      if (arrowheads) {
        arr <- arrow(length = lengthConnectors, type = typeConnectors, 
                     ends = endsConnectors)
      }
      else {
        arr <- NULL
      }
      plot <- plot + geom_text_repel(data = 
                                       subset(toptable,toptable[[y]] < 
                                                pCutoff & abs(toptable[[x]]) > FCcutoff),
                                     aes(label = subset(toptable, toptable[[y]] <
                                                          pCutoff & abs(toptable[[x]]) > 
                                                          FCcutoff)[["lab"]]),
                                     size = labSize, segment.color = colConnectors, 
                                     segment.size = widthConnectors, arrow = arr, 
                                     hjust = labhjust, vjust = labvjust, colour = labCol, 
                                     fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors && !is.null(selectLab)) {
      if (arrowheads) {
        arr <- arrow(length = lengthConnectors, type = typeConnectors, 
                     ends = endsConnectors)
      }
      else {
        arr <- NULL
      }
      plot <- plot + geom_text_repel(data = 
                                       subset(toptable,!is.na(toptable[["lab"]])), 
                                     aes(label = subset(toptable,!is.na(toptable[["lab"]]))[["lab"]]), 
                                     size = labSize, segment.color = colConnectors, 
                                     segment.size = widthConnectors, arrow = arr, 
                                     hjust = labhjust, vjust = labvjust, colour = labCol, 
                                     fontface = labFace, na.rm = TRUE)
    }
    else if (!drawConnectors && !is.null(selectLab)) {
      plot <- plot + 
        geom_text(data = subset(toptable,!is.na(toptable[["lab"]])),
                  aes(label = subset(toptable,!is.na(toptable[["lab"]]))[["lab"]]), 
                               size = labSize, check_overlap = TRUE, hjust = labhjust, 
                               vjust = labvjust, colour = labCol, fontface = labFace, 
                               na.rm = TRUE)
    }
    else if (!drawConnectors && is.null(selectLab)) {
      plot <- plot + 
        geom_text(data = 
                    subset(toptable,toptable[[y]] < pCutoff & abs(toptable[[x]]) > FCcutoff),
                  aes(label = subset(toptable, toptable[[y]] < pCutoff & abs(toptable[[x]]) > FCcutoff)[["lab"]]), 
                               size = labSize, check_overlap = TRUE, hjust = labhjust, 
                               vjust = labvjust, colour = labCol, fontface = labFace, 
                               na.rm = TRUE)
    }
  }
  else {
    if (drawConnectors && is.null(selectLab)) {
      if (arrowheads) {
        arr <- arrow(length = lengthConnectors, type = typeConnectors, 
                     ends = endsConnectors)
      }
      else {
        arr <- NULL
      }
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    toptable[[y]] < pCutoff & abs(toptable[[x]]) > 
                                                      FCcutoff), aes(label = subset(toptable, toptable[[y]] < 
                                                                                      pCutoff & abs(toptable[[x]]) > FCcutoff)[["lab"]]), 
                                      size = labSize, segment.color = colConnectors, 
                                      segment.size = widthConnectors, arrow = arr, 
                                      hjust = labhjust, vjust = labvjust, colour = labCol, 
                                      fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors && !is.null(selectLab)) {
      if (arrowheads) {
        arr <- arrow(length = lengthConnectors, type = typeConnectors, 
                     ends = endsConnectors)
      }
      else {
        arr <- NULL
      }
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                                   !is.na(toptable[["lab"]]))[["lab"]]), 
                                      size = labSize, segment.color = colConnectors, 
                                      segment.size = widthConnectors, arrow = arr, 
                                      hjust = labhjust, vjust = labvjust, colour = labCol, 
                                      fontface = labFace, na.rm = TRUE)
    }
    else if (!drawConnectors && !is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                              !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                             !is.na(toptable[["lab"]]))[["lab"]]), 
                                size = labSize, hjust = labhjust, vjust = labvjust, 
                                colour = labCol, fontface = labFace, na.rm = TRUE)
    }
    else if (!drawConnectors && is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                              toptable[[y]] < pCutoff & abs(toptable[[x]]) > 
                                                FCcutoff), aes(label = subset(toptable, toptable[[y]] < 
                                                                                pCutoff & abs(toptable[[x]]) > FCcutoff)[["lab"]]), 
                                size = labSize, hjust = labhjust, vjust = labvjust, 
                                colour = labCol, fontface = labFace, na.rm = TRUE)
    }
  }
  if (!is.null(encircle)) {
    plot <- plot + geom_encircle(data = subset(toptable, 
                                               rownames(toptable) %in% encircle), colour = encircleCol, 
                                 fill = encircleFill, alpha = encircleAlpha, size = encircleSize, 
                                 show.legend = FALSE, na.rm = TRUE)
  }
  if (!is.null(shade)) {
    plot <- plot + stat_density2d(data = subset(toptable, 
                                                rownames(toptable) %in% shade), fill = shadeFill, 
                                  alpha = shadeAlpha, geom = "polygon", contour = TRUE, 
                                  size = shadeSize, bins = shadeBins, show.legend = FALSE, 
                                  na.rm = TRUE)
  }
  return(plot)
}

create_EnhancedVolcano_hai <- function(data,name,xlim,ylim,log2FoldChange=2 ){
  EnhancedVolcano_hai(data,
                      lab = rownames(data),
                      x = "log2FoldChange",
                      y = "padj",
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = c(ifelse(abs(data$log2FoldChange)>log2FoldChange, 4, 1)),
                      title = name,
                      ylab = bquote(~-Log[10] ~ italic(padj)),
                      titleLabSize=18,
                      caption = 'FC cutoff: 1.5; padj cutoff: 0.05',
                      legendPosition = "right",
                      legendLabSize = 18,
                      axisLabSize=18,
                      xlim=xlim,
                      ylim = ylim,
                      labSize = 6,
                      legendLabels = 
                        c("NS", expression(Log[2]~FC), "padj", "padj + FC"))
  
}


