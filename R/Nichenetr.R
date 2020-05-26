load_Nichenetr_db <- function(db, local.path.nichenetr.db = "Nichenetr", saveRDS = T) {
  links = c(
    ligand_target_matrix = "https://zenodo.org/record/3260758/files/ligand_target_matrix.rds",
    lr_network = "https://zenodo.org/record/3260758/files/lr_network.rds",
    weighted_networks = "https://zenodo.org/record/3260758/files/weighted_networks.rds"
  )
  db.name <- paste0(db, ".rds")
  db.path <- paste0(local.path.nichenetr.db,"/",db.name)
  if(db.name %in% list.files(local.path.nichenetr.db)) {
    message(paste0("Loading ", db.path))
    d <- readRDS(db.path)
  }else{
    message(paste0("Download from: ", links[db]))
    d <- readRDS(url(links[db]))
    if(saveRDS) {
      if(!local.path.nichenetr.db %in% list.files()) dir.create(local.path.nichenetr.db)
      message(paste0("Save RDS to folder: ", local.path.nichenetr.db))
      saveRDS(d, db.path)
    }
  }
  return(d)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seu.sender PARAM_DESCRIPTION
#' @param seu.receiver PARAM_DESCRIPTION
#' @param ident.snd PARAM_DESCRIPTION, Default: NULL
#' @param ident.rcv PARAM_DESCRIPTION, Default: NULL
#' @param pct.snd PARAM_DESCRIPTION, Default: 0.25
#' @param pct.rcv PARAM_DESCRIPTION, Default: 0.1
#' @param fixed.ligand PARAM_DESCRIPTION, Default: NULL
#' @param feature.logfc.threshold.snd PARAM_DESCRIPTION, Default: 0
#' @param feature.logfc.threshold.rcv PARAM_DESCRIPTION, Default: 0.5
#' @param n_best_ligands PARAM_DESCRIPTION, Default: 30
#' @param n_ligand_target_links PARAM_DESCRIPTION, Default: 200
#' @param quantile_cutoff_target PARAM_DESCRIPTION, Default: 0
#' @param local.path.nichenetr.db PARAM_DESCRIPTION, Default: 'Nichenetr'
#' @param save.db PARAM_DESCRIPTION, Default: T
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[devtools]{remote-reexports}}
#'  \code{\link[dplyr]{desc}},\code{\link[dplyr]{select}}
#'  \code{\link[magrittr]{extract}}
#'  \code{\link[cowplot]{plot_grid}}
#'  \code{\link[ggpubr]{as_ggplot}},\code{\link[ggpubr]{get_legend}}
#' @rdname RunNichenetr
#' @export
#' @importFrom devtools install_github
#' @importFrom dplyr desc select
#' @importFrom magrittr set_rownames set_colnames
#' @importFrom cowplot plot_grid
#' @importFrom ggpubr as_ggplot get_legend

RunNichenetr <-
  function(seu.sender, seu.receiver, ident.snd = NULL, ident.rcv = NULL,
           pct.snd = 0.25, pct.rcv = 0.1, fixed.ligand = NULL,
           feature.logfc.threshold.snd = 0, feature.logfc.threshold.rcv = 0.5,
           n_best_ligands = 30, n_ligand_target_links = 200, quantile_cutoff_target = 0,
           local.path.nichenetr.db = "Nichenetr", save.db = T, spe = getOption("spe")) {
  check_spe(spe)
  if(!require(nichenetr)){
    message("Package \"nichenetr\" not installed")
    message("Please visit: https://github.com/saeyslab/nichenetr")
    devtools::install_github("saeyslab/nichenetr")
  }
  library(Seurat)
  library(nichenetr)
  library(tidyr)
  library(tibble)
  library(dplyr)
  library(rlang)
  library(openxlsx)
  library(ggplot2)

  ident.snd <- ident.snd %||% levels(Idents(seu.sender))
  ident.rcv <- ident.rcv %||% levels(Idents(seu.receiver))
  if(!"Nichenetr" %in% list.files()) dir.create("Nichenetr")

  # load and process necessary databases for nichenetr
  ligand_target_matrix <-
    load_Nichenetr_db("ligand_target_matrix",
                      local.path.nichenetr.db = local.path.nichenetr.db, saveRDS = save.db)
  lr_network <-
    load_Nichenetr_db("lr_network",
                      local.path.nichenetr.db = local.path.nichenetr.db, saveRDS = save.db)
  weighted_networks <-
    load_Nichenetr_db("weighted_networks",
                      local.path.nichenetr.db = local.path.nichenetr.db, saveRDS = save.db)
  weighted_networks_lr = weighted_networks$lr_sig %>%
    inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  if(spe == "mouse") {
    lr_network = lr_network %>%
      mutate(from = convert_human_to_mouse_symbols(from),
             to = convert_human_to_mouse_symbols(to)) %>%
      drop_na()
    colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>%
      convert_human_to_mouse_symbols()
    rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>%
      convert_human_to_mouse_symbols()
    ligand_target_matrix = ligand_target_matrix %>%
      .[!is.na(rownames(ligand_target_matrix)),
        !is.na(colnames(ligand_target_matrix))]
    weighted_networks_lr = weighted_networks_lr %>%
      mutate(from = convert_human_to_mouse_symbols(from),
             to = convert_human_to_mouse_symbols(to)) %>%
      drop_na()
  }

  # sender
  message("Check sender expressed genes")
  expressed_genes_sender_filename <-
    paste0("list_expressed_genes_sender_logFC_",feature.logfc.threshold.snd,".rds")
  if(expressed_genes_sender_filename %in% list.files("Nichenetr")) {
    message(paste0("Load list of sender expressed genes: ", expressed_genes_sender_filename))
    list_expressed_genes_sender <- readRDS(paste0("Nichenetr/",expressed_genes_sender_filename))
  }else{
    if(feature.logfc.threshold.snd > 0){
      list_expressed_genes_sender <-
        lapply(ident.snd, function(x){
          g <- FindMarkers(seu.sender, ident.1 = x,
                           logfc.threshold = feature.logfc.threshold.snd,
                           min.pct = pct.snd,
                           only.pos = T) %>% rownames()
          return(g)
        })
    }else{
      list_expressed_genes_sender <-
        lapply(ident.snd, gene_expressed, seu = seu.sender, pct = pct.snd)
    }
    message(paste0("Save list of sender expressed genes: Nichenetr/", expressed_genes_sender_filename))
    saveRDS(list_expressed_genes_sender, paste0("Nichenetr/",expressed_genes_sender_filename))
  }
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  ligands = lr_network %>% pull(from) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  nichenet_results <-
    list(list_expressed_genes_sender = list_expressed_genes_sender,
         expressed_ligands = expressed_ligands,
         par = list(
           ident.snd = ident.snd,
           ident.rcv = ident.rcv,
           fixed.ligand = fixed.ligand,
           pct.snd = pct.snd,
           pct.rcv = pct.rcv,
           feature.logfc.threshold.snd = feature.logfc.threshold.snd,
           feature.logfc.threshold.rcv = feature.logfc.threshold.rcv,
           n_best_ligands = n_best_ligands,
           n_ligand_target_links = n_ligand_target_links,
           quantile_cutoff_target = quantile_cutoff_target,
           local.path.nichenetr.db = local.path.nichenetr.db,
           save.db = save.db,
           spe = spe
         )
    )

  for (i in ident.rcv) {
    message(paste0("\nAnalyse start: ", i))
    feature_genes_folder <- paste0("Nichenetr/", i)
    if("feature_genes.rds" %in% list.files(feature_genes_folder)) {
      message(paste0("Load ", i, " feature genes: ", feature_genes_folder, "/feature_genes.rds"))
      feature_genes <- readRDS(paste0(feature_genes_folder, "/feature_genes.rds"))
    }else{
      dir.create(feature_genes_folder)
      feature_genes <- FindMarkers(seu.receiver, ident.1 = i,
                                   logfc.threshold = 0.25,
                                   only.pos = T)
      message(paste0("Save ", i, " feature genes: ", feature_genes_folder, "/feature_genes.rds"))
      saveRDS(feature_genes, paste0(feature_genes_folder, "/feature_genes.rds"))
    }
    expressed_genes_receiver = gene_expressed(i, seu.receiver, pct = pct.rcv)
    background_expressed_genes = intersect(expressed_genes_receiver, rownames(ligand_target_matrix))

    # DotPlot
    geneset_oi <- feature_genes %>%
      rownames_to_column("gene") %>%
      filter(p_val_adj <= 0.05 & abs(avg_logFC) >= feature.logfc.threshold.rcv) %>%
      pull(gene) %>%
      intersect(rownames(ligand_target_matrix))
    if(is_empty(geneset_oi)) {
      message(i, ": No feature genes after filtering. Decrease the logfc threshold?")
      next()
    }
    feature_genes_xlsx_path <-
      paste0(feature_genes_folder, "/filtered_feature_genes_logfc.thr_",
             feature.logfc.threshold.rcv, "_", i,".xlsx")
    message(paste0("Save filtered feature genes to Excel file: ", feature_genes_xlsx_path))
    write.xlsx(feature_genes[geneset_oi, ], rowNames = T, feature_genes_xlsx_path)

    receptors = lr_network %>% pull(to) %>% unique()
    expressed_receptors = intersect(receptors,expressed_genes_receiver)
    potential_ligands = lr_network %>%
      filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
      pull(from) %>% unique()
    ligand_activities = predict_ligand_activities(
      geneset = geneset_oi,
      background_expressed_genes = background_expressed_genes,
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = potential_ligands
    )
    ligand_activities = ligand_activities %>%
      arrange(-pearson) %>%
      mutate(rank = rank(dplyr::desc(pearson)))
    ligand_activities_path <- paste0(feature_genes_folder, "/ligand_activities.rds")
    message("Save ligand activities to: ", ligand_activities_path)
    saveRDS(ligand_activities, ligand_activities_path)

    best_upstream_ligands = ligand_activities %>%
      top_n(n_best_ligands, pearson) %>%
      arrange(-pearson) %>%
      pull(test_ligand) %>%
      unique() %>%
      union(intersect(ligand_activities$test_ligand, fixed.ligand))
    Ligands_dotplot <-
      DotPlot(seu.sender, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
    ggsave(paste0(feature_genes_folder, "/Ligands_dotplot.pdf"), Ligands_dotplot)

    # p_ligand_target_network
    active_ligand_target_links_df = best_upstream_ligands %>%
      lapply(get_weighted_ligand_target_links,
             geneset = geneset_oi,
             ligand_target_matrix = ligand_target_matrix,
             n = n_ligand_target_links) %>%
      bind_rows() %>% drop_na()
    active_ligand_target_links_df_path <- paste0(feature_genes_folder, "/active_ligand_target_links_df.rds")
    message("Save active ligand target links to: ", ligand_activities_path)
    saveRDS(active_ligand_target_links_df, active_ligand_target_links_df_path)

    active_ligand_target_links = prepare_ligand_target_visualization(
      ligand_target_df = active_ligand_target_links_df,
      ligand_target_matrix = ligand_target_matrix, cutoff = quantile_cutoff_target
    )
    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>%
      rev() %>% make.names()
    order_targets = active_ligand_target_links_df$target %>% unique() %>%
      intersect(rownames(active_ligand_target_links)) %>%
      make.names()
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

    p_ligand_target_network = vis_ligand_target %>%
      make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",
                          legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  +
      theme(axis.text.x = element_text(face = "italic")) +
      scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.006,0.012))
    ggsave(paste0(feature_genes_folder, "/p_ligand_target_network.pdf"), p_ligand_target_network)

    # p_ligand_receptor_network
    lr_network_top = lr_network %>%
      filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
      distinct(from,to)
    best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
    lr_network_top_df_large = weighted_networks_lr %>%
      filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
    lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
    lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>%
      magrittr::set_rownames(lr_network_top_df$to)

    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

    order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
    order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

    vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
    rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
    colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

    p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>%
      make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",
                          legend_title = "Prior interaction potential")
    ggsave(paste0(feature_genes_folder, "/p_ligand_receptor_network.pdf"), p_ligand_receptor_network)

    # p_ligand_receptor_network_strict
    lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
    ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
    receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

    lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>%
      inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
    lr_network_top_df_large_strict = lr_network_top_df_large_strict %>%
      inner_join(lr_network_top_df_large, by = c("from","to"))

    lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
    lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>%
      as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

    dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
    dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
    order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
    order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

    vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
    rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
    colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
    p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>%
      make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",
                          legend_title = "Prior interaction potential\n(bona fide)")
    ggsave(paste0(feature_genes_folder, "/p_ligand_receptor_network_strict.pdf"),
           p_ligand_receptor_network_strict)

    # combined heatmap: overlay ligand activities with target genes
    ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>%
      magrittr::set_rownames(ligand_activities$test_ligand)
    rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
    colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
    vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>%
      as.matrix(ncol = 1) %>%
      magrittr::set_colnames("Pearson")
    p_ligand_pearson = vis_ligand_pearson %>%
      make_heatmap_ggplot(
        "Prioritized ligands","Ligand activity",
        color = "darkorange",legend_position = "top",
        x_axis_position = "top",
        legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") +
      theme(legend.text = element_text(size = 9))

    figures_without_legend =
      cowplot::plot_grid(
        p_ligand_pearson +
          theme(legend.position = "none", axis.ticks = element_blank()) +
          theme(axis.title.x = element_text()), p_ligand_target_network +
          theme(legend.position = "none", axis.ticks = element_blank()) +
          ylab(""),
        align = "hv", nrow = 1,
        rel_widths = c(ncol(vis_ligand_pearson) + ncol(vis_ligand_target)/6, ncol(vis_ligand_target)))
    legends = cowplot::plot_grid(
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
      nrow = 1,
      align = "h")

    combined_plot = cowplot::plot_grid(
      figures_without_legend, legends,
      rel_heights = c(10,2), nrow = 2, align = "hv")
    ggsave(paste0(feature_genes_folder, "/combined_plot.pdf"), combined_plot)

    nichenet_results[[i]] <-
      list(feature_genes = feature_genes,
           ligand_activities = ligand_activities,
           active_ligand_target_links_df = active_ligand_target_links_df,
           best_upstream_ligands = best_upstream_ligands,
           best_upstream_ligands_strict = colnames(vis_ligand_receptor_network_strict),
           best_upstream_receptors = best_upstream_receptors,
           best_upstream_receptors_strict = rownames(vis_ligand_receptor_network_strict),
           ToPlot = list(
             vis_ligand_target = vis_ligand_target,
             vis_ligand_receptor_network = vis_ligand_receptor_network,
             vis_ligand_receptor_network_strict = vis_ligand_receptor_network_strict,
             vis_ligand_pearson = vis_ligand_pearson
           )
      )
  }
  saveRDS(nichenet_results, "Nichenetr/nichenet_results.rds")
  return(nichenet_results)
  }

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nichenet_results PARAM_DESCRIPTION
#' @param seu.sender PARAM_DESCRIPTION, Default: NULL
#' @param ident.rcv PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[cowplot]{plot_grid}}
#'  \code{\link[ggpubr]{as_ggplot}},\code{\link[ggpubr]{get_legend}}
#' @rdname Nichenetr_Plots
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom ggpubr as_ggplot get_legend

Nichenetr_Plots <- function(nichenet_results, seu.sender = NULL, ident.rcv = NULL) {
  library(nichenetr)
  library(Seurat)
  library(dplyr)
  library(rlang)
  library(ggplot2)

  Plots <- list()
  ident.rcv <- ident.rcv %||% nichenet_results$par$ident.rcv
  if(any(!ident.rcv %in% nichenet_results$par$ident.rcv)){
    warning(paste0("Ident(s) not found: ", setdiff(ident.rcv, nichenet_results$par$ident.rcv)))
    ident.rcv <- intersect(ident.rcv, nichenet_results$par$ident.rcv)
  }
  if(is.null(seu.sender)) {
    warning("Sender Seurat object not provided; Dot plot(s) not generated")
  }
  for (i in ident.rcv) {
    if(!is.null(seu.sender)) {
      best_upstream_ligands <- nichenet_results[[i]]$best_upstream_ligands
      Ligands_dotplot <-
        DotPlot(seu.sender, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") +
        RotatedAxis()
      Plots[[i]][["Ligands_dotplot"]] <- Ligands_dotplot
    }

    p_ligand_target_network =
      nichenet_results[[i]]$ToPlot$vis_ligand_target %>%
      make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",
                          legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  +
      theme(axis.text.x = element_text(face = "italic")) +
      scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.006,0.012))

    p_ligand_receptor_network =
      nichenet_results[[i]]$ToPlot$vis_ligand_receptor_network %>% t() %>%
      make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred",
                          x_axis_position = "top",
                          legend_title = "Prior interaction potential")

    p_ligand_receptor_network_strict =
      nichenet_results[[i]]$ToPlot$vis_ligand_receptor_network_strict %>% t() %>%
      make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",
                          legend_title = "Prior interaction potential\n(bona fide)")

    p_ligand_pearson = nichenet_results[[i]]$ToPlot$vis_ligand_pearson %>%
      make_heatmap_ggplot(
        "Prioritized ligands","Ligand activity",
        color = "darkorange",legend_position = "top",
        x_axis_position = "top",
        legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") +
      theme(legend.text = element_text(size = 9))

    figures_without_legend =
      cowplot::plot_grid(
        p_ligand_pearson +
          theme(legend.position = "none", axis.ticks = element_blank()) +
          theme(axis.title.x = element_text()), p_ligand_target_network +
          theme(legend.position = "none", axis.ticks = element_blank()) +
          ylab(""),
        align = "hv", nrow = 1,
        rel_widths = c(ncol(nichenet_results[[i]]$ToPlot$vis_ligand_pearson) +
                         ncol(nichenet_results[[i]]$ToPlot$vis_ligand_target) / 6,
                       ncol(nichenet_results[[i]]$ToPlot$vis_ligand_target)))
    legends = cowplot::plot_grid(
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
      nrow = 1,
      align = "h")

    combined_plot = cowplot::plot_grid(
      figures_without_legend, legends,
      rel_heights = c(10,2), nrow = 2, align = "hv")

    Plots[[i]] <-
      c(Plots[[i]],
        list(
          p_ligand_target_network = p_ligand_target_network,
          p_ligand_receptor_network = p_ligand_receptor_network,
          p_ligand_receptor_network_strict = p_ligand_receptor_network_strict,
          p_ligand_pearson = p_ligand_pearson,
          combined_plot = combined_plot
        ))
  }
  return(Plots)
}

# setwd("~/R documents/SeuratExtend_databases/2020-5-9 NicheNetr")
# options(max.print = 50, spe = "mouse")
# seu.receiver <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/PyMTEC_old.rds")
# seu.receiver <- NRAS13_newMalig
# seu.sender <- readRDS("~/R documents/2020-2-10 EC PyMT and E0771/rds/Tcell.rds")
# seu.sender <- NRAS13_all_2
# DefaultAssay(seu.sender) <- "RNA"
# Idents(seu.receiver) <- 'seurat_clusters'
# Idents(seu.sender) <- 'seurat_clusters'
# ident.snd <- levels(Idents(seu.sender))
# ident.snd <- "Endothelial cell"
# ident.rcv <- "LN_HEV"
# ident.rcv <- "5"
# local.path.nichenetr.db <- "Nichenetr"
# save.db = T
# spe = getOption("spe")
# pct.snd = 0.25
# pct.rcv = 0.1
# feature.logfc.threshold.snd = 0.25
# feature.logfc.threshold.rcv = 0.5
# quantile_cutoff_target = 0
# n_best_ligands = 30
# n_ligand_target_links = 200
# fixed.ligand = "Dll4"
# fixed.receptor = "Notch3"
# nichenet_results <- RunNichenetr(seu.sender = seu.sender, seu.receiver = seu.receiver)
# p <- Nichenetr_Plots(nichenet_results, seu.sender = seu.sender)
#
# # l-r interaction
# nichenet_results <- readRDS("~/R documents/2020-5-8 NicheNet for Florian/Nichenetr/nichenet_results.rds")
# best_upstream_ligands_strict <- nichenet_results$`5`$best_upstream_ligands_strict
# best_upstream_receptors_strict <- nichenet_results$`5`$best_upstream_receptors_strict
# lr <-
#   as.data.frame(lr_network_strict)[lr_network_strict$from %in% best_upstream_ligands_strict &
#                                      lr_network_strict$to %in% best_upstream_receptors_strict, ] %>%
#   distinct(from, to)
# lr$from <- lr$from %>% factor(levels = unique(.))
# lr$to <- lr$to %>% factor(levels = unique(.))
# l <- data.frame(gene = unique(lr$from)) %>% mutate(rank = rank(.$gene))
# r <- data.frame(gene = unique(lr$to)) %>% mutate(rank = rank(.$gene))
# lr$rank_l <- l$rank[lr$from]
# lr$rank_r <- r$rank[lr$to]
# p <-
#   ggplot() +
#   geom_tile(data = l, aes(x = 1, y = rank, fill = gene), color = "black") +
#   geom_text(data = l, aes(x = 0.4, y = rank, label = gene), hjust = 1) +
#   geom_tile(data = r, aes(x = 8, y = rank, fill = gene), color = "black") +
#   geom_text(data = r, aes(x = 8.6, y = rank, label = gene), hjust = 0) +
#   geom_segment(data = lr, aes(x = 1.5, xend = 7.5, y = rank_l, yend = rank_r),
#                arrow = arrow()) +
#   theme_classic() +
#   theme(axis.line = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         legend.position = "none") +
#   expand_limits(x = c(-1,10))
# ggsave("test.pdf", p)

# weighted_networks_lr = weighted_networks$lr_sig %>%
#   inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
# if(spe == "mouse") {
#   lr_network = lr_network %>%
#     mutate(from = convert_human_to_mouse_symbols(from),
#            to = convert_human_to_mouse_symbols(to)) %>%
#     drop_na()
#   colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>%
#     convert_human_to_mouse_symbols()
#   rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>%
#     convert_human_to_mouse_symbols()
#   ligand_target_matrix = ligand_target_matrix %>%
#     .[!is.na(rownames(ligand_target_matrix)),
#       !is.na(colnames(ligand_target_matrix))]
#   weighted_networks_lr = weighted_networks_lr %>%
#     mutate(from = convert_human_to_mouse_symbols(from),
#            to = convert_human_to_mouse_symbols(to)) %>%
#     drop_na()
# }

# download necessary databases for ligand-to-target signaling paths
# {
#   # weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
#   ligand_tf_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))
#
#   # lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#   sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
#   gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
# }
#
#
# ligands_all = "TNF" # this can be a list of multiple ligands if required
# ligands_all = "LTB" # this can be a list of multiple ligands if required
# targets_all = c("SELP","SELE","CXCL9","CXCL10","VCAM1")
# targets_all = c("NFKB2","RELB")
#
# active_signaling_network =
#   get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all,
#                             weighted_networks = weighted_networks)
#
# # For better visualization of edge weigths: normalize edge weights to make them comparable between
# # signaling and gene regulatory interactions
# active_signaling_network_min_max = active_signaling_network
# active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>%
#   mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
# active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>%
#   mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
#
# graph_min_max =
#   diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all,
#                                     targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")
#
# # To render the graph: uncomment following line of code
# DiagrammeR::render_graph(graph_min_max, layout = "kk")
