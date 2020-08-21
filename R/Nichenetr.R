load_Nichenetr_db <- function(db, local.path.nichenetr.db = "Nichenetr", saveRDS = T) {
  links = c(
    ligand_target_matrix = "https://zenodo.org/record/3260758/files/ligand_target_matrix.rds",
    lr_network = "https://zenodo.org/record/3260758/files/lr_network.rds",
    sig_network = "https://zenodo.org/record/3260758/files/signaling_network.rds",
    weighted_networks = "https://zenodo.org/record/3260758/files/weighted_networks.rds",
    gr_network = "https://zenodo.org/record/3260758/files/gr_network.rds"
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
#' @param object PARAM_DESCRIPTION
#' @param assay PARAM_DESCRIPTION, Default: NULL
#' @param features PARAM_DESCRIPTION
#' @param cols PARAM_DESCRIPTION, Default: 'RdYlBu'
#' @param col.min PARAM_DESCRIPTION, Default: -2.5
#' @param col.max PARAM_DESCRIPTION, Default: 2.5
#' @param dot.min PARAM_DESCRIPTION, Default: 0
#' @param dot.scale PARAM_DESCRIPTION, Default: 6
#' @param group.by PARAM_DESCRIPTION, Default: NULL
#' @param split.by PARAM_DESCRIPTION, Default: NULL
#' @param scale.by PARAM_DESCRIPTION, Default: 'radius'
#' @param scale.min PARAM_DESCRIPTION, Default: NA
#' @param scale.max PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetSeuratDotplotData
#' @export

GetSeuratDotplotData <-
  function(object, assay = NULL, features,
           cols = "RdYlBu", col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
           group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA, scale.max = NA) {
    library(rlang)
    PercentAbove <- function(x, threshold) {
      return(length(x = x[x > threshold]) / length(x = x))
    }
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                         radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    data.features$id <- if (is.null(x = group.by)) {
      Idents(object = object)
    } else {
      object[[group.by, drop = TRUE]]
    }
    if (!is.factor(x = data.features$id)) {
      data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
      splits <- object[[split.by, drop = TRUE]]
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enought colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
      data.features$id <- paste(data.features$id, splits, sep = "_")
      unique.splits <- unique(x = splits)
      id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
                          "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
      data.use <- data.features[data.features$id == ident,
                                1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
        return(mean(x = expm1(x = x)))
      })
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
                       threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
      data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                             FUN = function(x) {
                               data.use <- data.plot[data.plot$features.plot ==
                                                       x, "avg.exp"]
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min,
                                                  max = col.max)
                               return(data.use)
                             })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
      avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
                                           breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot,
                                      levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
      splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id),
                                        split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L),
                           2)
      data.plot$colors <- mapply(FUN = function(color, value) {
        return(colorRampPalette(colors = c("grey", color))(20)[value])
      }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
                       no = "colors")
    if (!is.na(x = scale.min)) {
      data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
      data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    return(data.plot)
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
#' @param feature.logfc.threshold.rcv PARAM_DESCRIPTION, Default: 0.25
#' @param feature.rcv PARAM_DESCRIPTION, Default: NULL
#' @param n_best_ligands PARAM_DESCRIPTION, Default: 30
#' @param n_ligand_target_links PARAM_DESCRIPTION, Default: 200
#' @param quantile_cutoff_target PARAM_DESCRIPTION, Default: 0
#' @param local.path.nichenetr.db PARAM_DESCRIPTION, Default: 'Nichenetr'
#' @param save.db PARAM_DESCRIPTION, Default: T
#' @param spe PARAM_DESCRIPTION, Default: getOption("spe")
#' @param custom_db PARAM_DESCRIPTION, Default: T
#' @param regulons PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{desc}},\code{\link[dplyr]{select}}
#'  \code{\link[magrittr]{extract}}
#'  \code{\link[cowplot]{plot_grid}}
#'  \code{\link[ggpubr]{as_ggplot}},\code{\link[ggpubr]{get_legend}}
#' @rdname RunNichenetr
#' @export
#' @importFrom dplyr desc select
#' @importFrom magrittr set_rownames set_colnames
#' @importFrom cowplot plot_grid
#' @importFrom ggpubr as_ggplot get_legend

RunNichenetr <-
  function(seu.sender, seu.receiver, ident.snd = NULL, ident.rcv = NULL,
           pct.snd = 0.25, pct.rcv = 0.1, fixed.ligand = NULL,
           feature.logfc.threshold.snd = 0, feature.logfc.threshold.rcv = 0.25, feature.rcv = NULL,
           n_best_ligands = 30, n_ligand_target_links = 200, quantile_cutoff_target = 0,
           local.path.nichenetr.db = "Nichenetr", save.db = T, spe = getOption("spe"),
           custom_db = T, regulons = NULL) {
    check_spe(spe)
    import("nichenetr")
    library(Seurat)
    library(tidyr)
    library(tibble)
    library(dplyr)
    library(rlang)
    library(openxlsx)
    library(ggplot2)

    ident.snd <- ident.snd %||% levels(Idents(seu.sender))
    ident.rcv <- ident.rcv %||% levels(Idents(seu.receiver))
    if(!"Nichenetr" %in% list.files()) dir.create("Nichenetr")

    # sender
    message("Check sender expressed genes")
    expressed_genes_sender_filename <-
      paste0("list_expressed_genes_sender_logFC_",feature.logfc.threshold.snd,".rds")
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
    expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
    names(list_expressed_genes_sender) <- ident.snd

    nichenet_results <-
      list(list_expressed_genes_sender = list_expressed_genes_sender,
           par = list(
             ident.snd = ident.snd,
             ident.rcv = ident.rcv,
             pct.snd = pct.snd,
             pct.rcv = pct.rcv,
             fixed.ligand = fixed.ligand,
             feature.logfc.threshold.snd = feature.logfc.threshold.snd,
             feature.logfc.threshold.rcv = feature.logfc.threshold.rcv,
             feature.rcv = feature.rcv,
             n_best_ligands = n_best_ligands,
             n_ligand_target_links = n_ligand_target_links,
             quantile_cutoff_target = quantile_cutoff_target,
             local.path.nichenetr.db = local.path.nichenetr.db,
             save.db = save.db,
             spe = spe,
             custom_db = custom_db,
             regulons = regulons)
      )

    for (i in ident.rcv) {
      # receiver
      message(paste0("\nAnalyse start: '", i,"'"))
      feature_genes_folder <- file.path("Nichenetr", sub("[/]","_",i))
      feature_genes_file <- paste0("feature_genes_", feature.logfc.threshold.rcv, ".rds")
      if(!is.null(feature.rcv)) {
        feature_genes <- feature.rcv[[i]]
      }else if(feature_genes_file %in% list.files(feature_genes_folder)) {
        message(paste0("Load '", i, "' feature genes: ", file.path(feature_genes_folder, feature_genes_file)))
        feature_genes <- readRDS(file.path(feature_genes_folder, feature_genes_file))
      }else{
        if(!dir.exists(feature_genes_folder)) dir.create(feature_genes_folder)
        feature_genes <- FindMarkers(seu.receiver, ident.1 = i, only.pos = T,
                                     logfc.threshold = feature.logfc.threshold.rcv)
        message(paste0("Save '", i, "' feature genes: ", file.path(feature_genes_folder, feature_genes_file)))
        saveRDS(feature_genes, file.path(feature_genes_folder, feature_genes_file))
      }
      expressed_genes_receiver = gene_expressed(i, seu.receiver, pct = pct.rcv)

      if(!custom_db){
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
      }else{
        message("Constructing customed database. This may take a while.\n",
                "If you want to use the default databases, please set 'custom_db' to FALSE")
        lr_network <- load_Nichenetr_db("lr_network")
        if(spe == "mouse"){
          lr_network = lr_network %>%
            mutate(from = convert_human_to_mouse_symbols(from),
                   to = convert_human_to_mouse_symbols(to)) %>%
            drop_na()
        }
        lr_network = lr_network %>%
          filter(database != "ppi_prediction_go" &
                   database != "ppi_prediction" &
                   from %in% expressed_genes_sender &
                   to %in% expressed_genes_receiver)

        sig_network = load_Nichenetr_db("sig_network")
        if(spe == "mouse"){
          sig_network = sig_network %>%
            mutate(from = convert_human_to_mouse_symbols(from),
                   to = convert_human_to_mouse_symbols(to)) %>%
            drop_na()
        }

        if(!is.null(regulons)) {
          library(AnnotationDbi)
          message("Using customed gene regulatory network")
          gr_network <-
            unlist2(regulons) %>%
            data.frame(from = names(.), to = ., source = "SCENIC",
                       database = "SCENIC", stringsAsFactors = F) %>%
            as_tibble() %>%
            filter(from %in% expressed_genes_receiver &
                     to %in% expressed_genes_receiver)
          source_weights_df_new <- rbind(source_weights_df, data.frame(source = "SCENIC", weight = 1))
        }else{
          gr_network = load_Nichenetr_db("gr_network")
          if(spe == "mouse"){
            gr_network = gr_network %>%
              mutate(from = convert_human_to_mouse_symbols(from),
                     to = convert_human_to_mouse_symbols(to)) %>%
              drop_na()
          }
          source_weights_df_new <- source_weights_df
        }

        message("Constructing weighted_networks")
        weighted_networks =
          construct_weighted_networks(
            lr_network = lr_network,
            sig_network = sig_network,
            gr_network = gr_network,
            source_weights_df = source_weights_df_new
          )

        weighted_networks_lr = weighted_networks$lr_sig %>%
          inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
        ligands = lr_network %>% pull(from) %>% unique()
        expressed_ligands = intersect(ligands,expressed_genes_sender)
        message("Constructing ligand_target_matrix")
        ligand_target_matrix =
          construct_ligand_target_matrix(
            weighted_networks = weighted_networks,
            ligands = as.list(expressed_ligands), algorithm = "PPR",
            damping_factor = hyperparameter_list$damping_factor,
            ltf_cutoff = hyperparameter_list$ltf_cutoff
          )
      }

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
        file.path(feature_genes_folder, paste0("filtered_feature_genes_logfc.thr_",feature.logfc.threshold.rcv,"_", i,".xlsx"))
      message(paste0("Save filtered feature genes to Excel file: ", feature_genes_xlsx_path))
      write.xlsx(feature_genes[geneset_oi, ], rowNames = T, feature_genes_xlsx_path)

      receptors = lr_network %>% pull(to) %>% unique()
      expressed_receptors = intersect(receptors,expressed_genes_receiver)
      potential_ligands = lr_network %>%
        filter(from %in% expressed_genes_sender & to %in% expressed_receptors) %>%
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
      ligand_activities_path <- file.path(feature_genes_folder, "ligand_activities.rds")
      message("Save ligand activities to: ", ligand_activities_path)
      saveRDS(ligand_activities, ligand_activities_path)

      best_upstream_ligands = ligand_activities %>%
        top_n(n_best_ligands, pearson) %>%
        arrange(-pearson) %>%
        pull(test_ligand) %>%
        unique() %>%
        union(intersect(ligand_activities$test_ligand, fixed.ligand))

      data.plot <- GetSeuratDotplotData(seu.sender, features = best_upstream_ligands %>% rev())
      makePlot <- function(plot.data){
        env <- new.env(parent = globalenv())
        env$subset <- plot.data

        my.plot <- with(env, {
          p <- ggplot(data = subset, mapping = aes(x = features.plot, y = id)) +
            geom_point(mapping = aes(size = pct.exp, color = avg.exp.scaled)) +
            guides(size = guide_legend(title = "Percent Expressed")) +
            labs(x = "Features", y = "Identity") +
            theme_classic() +
            scale_color_distiller(palette = "RdYlBu") +
            guides(color = guide_colorbar(title = "Average Expression")) +
            RotatedAxis()
          return(p)
        })
        return(my.plot)
      }
      Ligands_dotplot <- makePlot(data.plot)
      ggsave(file.path(feature_genes_folder, "Ligands_dotplot.pdf"), Ligands_dotplot)

      # p_ligand_target_network
      active_ligand_target_links_df = best_upstream_ligands %>%
        lapply(get_weighted_ligand_target_links,
               geneset = geneset_oi,
               ligand_target_matrix = ligand_target_matrix,
               n = n_ligand_target_links) %>%
        bind_rows() %>% drop_na()
      active_ligand_target_links_df_path <- file.path(feature_genes_folder, "active_ligand_target_links_df.rds")
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
      ggsave(file.path(feature_genes_folder, "p_ligand_target_network.pdf"), p_ligand_target_network)

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
      ggsave(file.path(feature_genes_folder, "p_ligand_receptor_network.pdf"), p_ligand_receptor_network)

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
      ggsave(file.path(feature_genes_folder, "combined_plot.pdf"), combined_plot)

      nichenet_results[[i]] <-
        list(feature_genes = feature_genes,
             potential_ligands = potential_ligands,
             ligand_activities = ligand_activities,
             active_ligand_target_links_df = active_ligand_target_links_df,
             best_upstream_ligands = best_upstream_ligands,
             best_upstream_receptors = best_upstream_receptors,
             Ligands_dotplot = Ligands_dotplot,
             combined_plot = combined_plot,
             p_ligand_receptor_network = p_ligand_receptor_network
        )
    }

    message("Save results to: Nichenetr/nichenet_results.rds")
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
