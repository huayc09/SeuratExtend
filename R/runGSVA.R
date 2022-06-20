#' @title GSEA plot
#' @description Generate plot that mimic the Gene Set Enrichment computational analysis
#' published by the Broad Institute
#' @param seu Seurat object
#' @param group.by A variable name in meta.data to
#' group the violin plots by, or string with the same length of cells
#' @param geneset A list of genes
#' @param ident.1 Cell type name
#' @param ident.2 (Optional) Second cell type name to compare with 'ident.1'
#' @param cells Cell names to use, Default: all cells
#' @param slot Slot to pull feature data for, Default: 'data'
#' @param assay Name of assay to use, defaults to the active assay
#' @param logFC.cutoff Range of Log2 fold change of ranked genes, Default: c(-1, 1)
#' @param p.position Where to put ES and p value text, Default: c(0.8, 0.85)
#' @param sample Points used for drawing curves. Default: 500
#' @param title Title of geneset
#' @return Plot
#' @details P value is calculated using Kolmogorov-Smirnov Tests
#' @examples
#' GSEAplot(
#'   pbmc, ident.1 = "Naive CD4 T", title = "INTERFERON_GAMMA_RESPONSE",
#'   geneset = hall50$human$HALLMARK_INTERFERON_GAMMA_RESPONSE)
#' @rdname GSEAplot
#' @export

GSEAplot <-
  function(
    seu,
    group.by = NULL,
    geneset,
    ident.1 = NULL,
    ident.2 = NULL,
    cells = NULL,
    slot = "data",
    assay = NULL,
    logFC.cutoff = c(-1,1),
    p.position = c(0.8,0.85),
    sample = 500,
    title = deparse(substitute(geneset))
  ) {

    library(ggplot2)
    library(ggpubr)
    library(scales)

    std.matr <- Seu2Matr(
      seu = seu,
      features = rownames(seu),
      group.by = group.by,
      cells = cells,
      slot = "data",
      assay = assay,
      verbose = FALSE
    )
    std.matr$matr <-
      std.matr$matr[,apply(std.matr$matr, 2, function(x) any(x > 0))]

    ident.info <- CheckIdent(
      f = std.matr$f,
      ident.1 = ident.1,
      ident.2 = ident.2
    )

    geneset <- intersect(geneset, colnames(std.matr$matr))
    if(length(geneset) == 0) stop("No genes in 'geneset' found in expression matrix")

    title <- paste0("Enrichment Plot: ", title)
    title2 <- paste0("Rank in gene list (", ident.info$name.1,
                     " vs. ", ident.info$name.2, ")")

    cell.l <- (ident.info$cell.1 | ident.info$cell.2)
    scores <- expm1(std.matr$matr[cell.l, ])
    scores <- CalcTscore(scores, std.matr$f[cell.l], "logfc", idents = ident.info$name.1)
    scores <- scores[order(scores[[1]], decreasing = T),,drop = F]
    scores <- data.frame(
      genes = factor(rownames(scores), levels = rownames(scores)),
      LogFC = scores[[1]],
      geneset = (rownames(scores) %in% geneset),
      rank = 1:nrow(scores),
      centered_rank = min(which.min(abs(scores[[1]]))) - (1:nrow(scores))
    )
    step.T <- 1/sum(scores$geneset)
    step.F <- 1/sum(!scores$geneset)
    scores$`Enrichment Score (ES)` <-
      1:nrow(scores) %>%
      sapply(function(x) sum(scores$geneset[1:x]) * step.T - sum(!scores$geneset[1:x]) * step.F)

    es <- c(max(scores$`Enrichment Score (ES)`),
            min(scores$`Enrichment Score (ES)`))
    es <- es[which.max(abs(es))]
    p <- suppressWarnings(ks.test(which(scores$geneset), scores$rank))[["p.value"]]
    p.text <- paste0("ES = ",format(es, digits = 3),"\np ",
                     if(p == 0) "< 2.2e-16" else paste0("= ",format(p, digits = 3)))
    annotate_npc <- function(label, x, y, ...) {
      annotation_custom(grid::textGrob(
        x = unit(x, "npc"), y = unit(y, "npc"), label = label, ...))
    }

    samp <- c(1, nrow(scores),
              which.min(scores$`Enrichment Score (ES)`),
              which.max(scores$`Enrichment Score (ES)`),
              1:sample * ceiling(nrow(scores) / sample))
    score1 <- scores[scores$rank %in% samp,]
    p1 <- ggplot(score1, aes(x = rank, y = `Enrichment Score (ES)`)) +
      geom_hline(yintercept = 0, color = "#B5B5B5", size = 1) +
      geom_step(color = "#505050") +
      annotate_npc(p.text, p.position[1], p.position[2]) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(1,1,0,1)) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, nrow(scores))) +
      ggtitle(title)

    if(nrow(scores) > 256) {
      n <- 256
      width <- ceiling(nrow(scores) / n)
      width.half <- floor(width / 2)
      score2 <- scores[rep(1:width, length.out = nrow(scores)) == width.half,]
    } else score2 <- scores
    p2 <- ggplot() +
      geom_tile(
        mapping = aes(x = rank, y = 0.5, fill = centered_rank),
        data = score2, width = width + 1, color = NA) +
      geom_bar(
        mapping = aes(x = rank, y = 1),
        data = scores[scores$geneset == 1,], stat = "identity", width = 10) +
      theme_bw() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            plot.margin = margin(0,1,0,1)) +
      scale_fill_gradient2(low = muted("blue", 50), high = muted("red", 50)) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, nrow(scores)))

    p3 <- ggplot(score1, aes(x = rank, y = LogFC)) +
      geom_area() +
      theme_bw() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(0,1,1,1)) +
      coord_cartesian(ylim = logFC.cutoff) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, nrow(scores))) +
      xlab(title2)

    p <- ggarrange(p1, p2, p3, ncol = 1, align = "v", heights = c(6,1,3))

    return(p)
  }



