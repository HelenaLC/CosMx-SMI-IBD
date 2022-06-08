suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrastr)
  library(ggrepel)
  library(patchwork)
  library(RColorBrewer)
  library(SingleCellExperiment)
})

# relative image resolution / scale factor
IMG_RES <- 0.05 

# parallelization
t <- ifelse(exists("THREADS"), THREADS, 1)
data.table::setDTthreads(t)
bp <- if (t == 1) {
  BiocParallel::SerialParam() 
} else BiocParallel::MulticoreParam(t)

# image, sample & clustering identifiers
names(iids) <- iids <- c("morphology", "segmentation")
names(sids) <- sids <- list.dirs("data", FALSE, FALSE)
names(kids) <- kids <- c("Louvain0", "Louvain", "SingleR1", "SingleR2")

# color palettes
gpal <- c("HC" = "royalblue", "UC" = "orange", "CD" = "tomato")
spal <- brewer.pal(9, "Paired")
pal <- \(obj, col) {
  n <- nlevels(factor(obj[[col]]))
  pal <- brewer.pal(12, "Paired")
  colorRampPalette(pal)(n)
}

# z-score normalization
.z <- \(x, mar = 2, th = 2.5) apply(x, mar, \(.) {
  sd <- stats::sd(., na.rm = TRUE)
  . <- . - mean(., na.rm = TRUE)
  if (sd != 0) 
    . <- ./sd
  .[. > th] <- th
  .[. < -th] <- -th
  return(.)
})

# 0-1 scaling using lower & upper 
# 1% quantiles as boundaries
.scale01 <- \(x, margin = 1, q = 0.01) {
    if (!is(x, "matrix")) 
        x <- as.matrix(x)
    qs <- c(rowQuantiles, colQuantiles)[[margin]]
    qs <- qs(x, probs = c(q, 1 - q))
    qs <- matrix(qs, ncol = 2)
    x <- switch(margin, `1` = (x - qs[, 1])/(qs[, 2] - qs[, 1]), 
        `2` = t((t(x) - qs[, 1])/(qs[, 2] - qs[, 1])))
    x[x < 0 | is.na(x)] <- 0
    x[x > 1] <- 1
    return(x)
}

# bar plot of relative cluster abundances
.plot_fq <- \(sce, xy, by = NULL, order = FALSE) {
  if (is.factor(sce[[xy[2]]]))
    sce[[xy[2]]] <- droplevels(sce[[xy[2]]])
  df <- table(x = sce[[xy[1]]], y = sce[[xy[[2]]]])
  xo <- if (order) {
    o <- hclust(dist(as.matrix(df)))$order
    scale_x_discrete(limits = \(.) .[o])
  }
  df <- as.data.frame(df)
  names(df) <- c(xy, "n")
  by <- if (!is.null(by)) {
    i <- match(df[[1]], sce[[xy[1]]])
    df[[by]] <- sce[[by]][i]
    list(facet_wrap(by, scales = "free_x"))
  }
  ggplot(df, 
    aes_string(xy[1], "n", fill = xy[2])) + 
    xo + by + geom_bar(
      col = "black", size = 0.2, width = 1,
      stat = "identity", position = "fill") +
    scale_fill_manual(NULL, values = pal(sce, xy[2])) +
    scale_y_continuous(labels = scales::percent_format()) +
    coord_cartesian(expand = FALSE) +
    theme_linedraw() + theme(
      aspect.ratio = 3/2,
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.key.size = unit(0.5, "lines"),
      strip.text = element_text(color = "black"),
      strip.background = element_rect(fill = NA),
      axis.text.x = element_text(angle = 45, hjust = 1))
}

# heatmap of mean expressions
# w/ hierarchical clustering 
# (x = genes, y = clusters)
.plot_hm <- \(sce, mgs, k, n) {
  mgs <- lapply(mgs, \(df) {
    df$gene <- rownames(df)
    data.frame(df, check.names = FALSE)
  })
  top <- mgs %>% 
    bind_rows(.id = "cluster") %>% 
    group_by(cluster) %>% 
    slice_min(Top, n = n)
  gs <- unique(top$gene)
  gs <- intersect(gs, rownames(sce))
  es <- as.matrix(t(logcounts(sce[gs, ])))
  ms <- aggregate(es, list(kid = sce[[k]]), mean)
  ms[, -1] <- scale(ms[, -1])
  df <- pivot_longer(ms, -kid)
  ggplot(df, 
    aes(name, kid, fill = value)) + 
    geom_tile() +
    scale_fill_gradient2(
      "scaled mean\nexpression",
      low = "red", high = "blue") +
    coord_cartesian(expand = FALSE) +
    scale_x_discrete(limits = gs) +
    scale_y_discrete(limits = \(.) rev(.)) +
    theme_linedraw() + theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(size = 8,
        angle = 90, hjust = 1, vjust = 0.5))
}

.plot_xy <- \(obj, xy = NULL, col = NULL, by = NULL, 
  size = NULL, trans = NULL, nas = FALSE, thm = c("black", "white"), ...) {
  thm <- match.arg(thm)
  if (is.null(xy))
    xy <- c("x", "y")
  if (!is.null(col)) {
    nai <- is.na(obj[[col]]) 
    .col <- ifelse(
      is.null(trans), col, 
      sprintf("%s(%s)", trans, col))
  } else {
    nai <- FALSE
    .col <- NULL
  }
  aes <- if (!is.null(size) && is.logical(obj[[size]]))
    list(
      guides(size = "none"),
      scale_size_manual(values = c("FALSE" = 0.25, "TRUE" = 1.5)))
  dr <- any(grepl("UMAP", xy))
  if (dr) {
    dots <- list(size = 0.1, alpha = 0.5, ...)
    thm <- "white"
  } else {
    dots <- list(...) 
  }
  if (length(dots) == 0) dots <- NULL
  if (is.data.frame(obj)) {
    fun <- ggplot
    if (!nas) obj <- obj[!nai, ]
  } else {
    fun <- scater::ggcells
    if (!nas) obj <- obj[, !nai]
  }
  p <- fun(obj,
    aes_string(xy[1], xy[2], col = .col, size = size)) + aes +
    do.call(geom_point_rast, c(list(shape = 16), dots)) + 
    ifelse(length(by) == 1, facet_wrap, facet_grid)(by) +
    theme_void() + theme(
      plot.title = element_text(hjust = 0.5),
      strip.text.y = element_text(angle = 90))
  if (!is.null(col)) {
    if (scale_type(obj[[col]]) == "continuous") {
      p <- p + theme(
        legend.position = "bottom",
        legend.key.height = unit(0.5, "lines"))
    } else {
      if (is.logical(obj[[col]])) {
        p <- p + scale_color_manual(NULL, values = c(
          "FALSE" = "lightgrey", "TRUE" = "orangered"))
      } else {
        p <- p + scale_color_manual(NULL, values = pal(obj, col))
      }
      p <- p + 
        theme(legend.key.size = unit(0.5, "lines")) +
        guides(col = guide_legend(override.aes = list(alpha = 1, size = 2)))
    }
  }
  switch(thm, white = p, p + theme(
    panel.background = element_rect(fill = "black")))
}
