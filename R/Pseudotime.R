#' @title Visualize Gene Expression Dynamics Along Differentiation Trajectories
#' @description This function set provides tools to visualize smoothed gene expression trends along developmental trajectories derived from pseudotime analyses. GeneTrendCurve.Palantir() uses pseudotime results from Palantir, while GeneTrendCurve.Slingshot() utilizes outputs from Slingshot for similar purposes. These visualizations help elucidate the dynamics of gene expression as cells differentiate and develop.
#' @param seu A Seurat object that has undergone pseudotime analysis with either Palantir or Slingshot, containing necessary gene expression data and pseudotime calculations.
#' @param features A vector of gene names whose expression trends will be plotted along the pseudotime axis.
#' @param pseudotime.data The name of the slot within `seu@misc` where pseudotime data is stored or a data.frame containing pseudotime information; for Palantir, this is typically 'Pseudotime', and for Slingshot, 'PCA'.
#' @param magic Logical indicating whether to use denoised gene expression data for plotting, obtained via MAGIC. Requires prior execution of Palantir.Magic() if set to TRUE. Default: FALSE.
#' @param method The method used to smooth and visualize gene expression trends. Can be either "smooth" or "gam" to choose between simple moving averages or generalized additive models, respectively. Default: c("smooth", "gam").
#' @param point Logical indicating whether to plot raw data points on top of the smoothed expression curves. Default: FALSE.
#' @param line.size Size of the line used for plotting the gene expression trends. Default: NULL.
#' @param pt.alpha Opacity of the plotted points if `point = TRUE`. Default: 0.3.
#' @param pt.size Size of the plotted points if `point = TRUE`. Default: 1.
#' @param ncol Number of columns in the plot layout if multiple genes are plotted. Default: NULL.
#' @param se Logical indicating whether to include confidence intervals around the smoothed trends. Default: TRUE.
#' @param conda_env Name of the Conda environment where necessary Python dependencies are installed for operations that require Python. Default: 'seuratextend'.
#' @return A ggplot object displaying gene expression trends along computed pseudotime trajectories.
#' @details GeneTrendCurve.Palantir() and GeneTrendCurve.Slingshot() are designed for researchers looking to understand how gene expression changes over the course of cellular development, providing clear visualizations that highlight trends and variations in gene expression along calculated trajectories. Each function tailors its approach based on the source of pseudotime data, ensuring that the visualizations are closely aligned with the underlying data and pseudotime analysis results.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Load an example Seurat Object
#' mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))
#'
#' # Calculate diffusion map and pseudotime using Palantir
#' mye_small <- Palantir.RunDM(mye_small)
#' mye_small <- Palantir.Pseudotime(mye_small, start_cell = "sample1_GAGAGGTAGCAGTACG-1")
#'
#' # Basic gene trend curve visualization
#' GeneTrendCurve.Palantir(mye_small, features = c("CD14", "FCGR3A"))
#'
#' # Customizing gene trend curves with Palantir
#' ps <- mye_small@misc$Palantir$Pseudotime
#' colnames(ps)[3:4] <- c("fate1", "fate2")
#' GeneTrendCurve.Palantir(mye_small, pseudotime.data = ps, features = c("CD14", "FCGR3A"), point = TRUE, se = FALSE)
#'
#' # Run Slingshot
#' mye_small <- RunSlingshot(mye_small, group.by = "cluster", start.clus = "Mono CD14")
#'
#' # Using Slingshot for similar visualizations
#' GeneTrendCurve.Slingshot(mye_small, features = c("CD14", "FCGR3A"))
#'
#' @rdname GeneTrendCurve-Methods
#' @export

GeneTrendCurve.Palantir <- function(
    seu,
    features,
    pseudotime.data = "Pseudotime",
    magic = FALSE,
    method = c("smooth","gam"),
    point = FALSE,
    line.size = NULL,
    pt.alpha = 0.3,
    pt.size = 1,
    ncol = NULL,
    se = TRUE,
    conda_env = "seuratextend") {

  library(Seurat)
  library(reshape2)
  library(ggplot2)

  check_method(method)
  seu <- process_magic_and_set_assay(seu, magic, conda_env)

  # Check the pseudotime.data parameter
  pseudotime_df <- extract_pseudotime_data(seu, pseudotime.data, "Palantir", required_cols = c("Pseudotime", "Entropy"))
  fate_names <- setdiff(colnames(pseudotime_df), c("Pseudotime", "Entropy"))

  df <- cbind(pseudotime_df, FetchData(seu, vars = features))

  # Melt the dataframe to get the desired format
  df_melt_gene <- melt(
    df, id.vars = colnames(pseudotime_df),
    measure.vars = features,
    variable.name = "Gene",
    value.name = "Expression"
  )

  df_melt <- melt(
    df_melt_gene, id.vars = c("Pseudotime", "Entropy", "Gene", "Expression"),
    measure.vars = fate_names,
    variable.name = "cell_fate",
    value.name = "fate_prob"
  )
  df_melt$Gene <- factor(df_melt$Gene, levels = unique(df_melt$Gene))

  # Trim the dataframe
  df_melt <- trim_data_by_pseudotime(df, df_melt, fate_names)

  # Create the base plot
  p <- build_base_plot(point, df_melt, pt.alpha, pt.size, ncol)

  if (method[1] == "gam") {
    predictions <- perform_gam_analysis(df_melt, features, fate_names)
  } else {
    predictions <- df_melt
  }

  p <- add_plot_layer(p, method, line.size, se, predictions)

  return(p)
}

#' @title Visualize Gene Expression Dynamics Along Trajectories with Heatmaps
#' @description This function suite provides tools for visualizing gene expression dynamics along pseudotime trajectories using heatmaps, leveraging results from either Palantir or Slingshot pseudotime analyses. GeneTrendHeatmap.Palantir() and GeneTrendHeatmap.Slingshot() each cater to the specific formatting and requirements of the pseudotime data generated by their respective tools.
#' @param seu A Seurat object that has undergone pseudotime analysis with either Palantir or Slingshot, containing necessary gene expression data and pseudotime calculations.
#' @param features A vector of gene names for which expression dynamics will be visualized in the heatmap.
#' @param pseudotime.data The name of the slot within `seu@misc` where pseudotime data is stored or a data.frame containing pseudotime information; for Palantir, this is typically 'Pseudotime', and for Slingshot, 'PCA'.
#' @param lineage The specific lineage or trajectory to visualize in the heatmap, which corresponds to the clusters or trajectories defined within the pseudotime analysis.
#' @param magic Logical indicating whether to use denoised gene expression data for plotting, obtained via MAGIC. Requires prior execution of Palantir.Magic() if set to TRUE. Default: FALSE.
#' @param scale Boolean indicating whether to scale gene expression data between 0 and 1, which can help in normalizing expression levels for better visual contrast in the heatmap. Default: TRUE.
#' @param sort Boolean indicating whether to sort genes by expression peak along pseudotime from left to right. Default: TRUE.
#' @param color_scheme The color scheme to use for the heatmap. Default: 'D'.
#' @param conda_env Name of the Conda environment where necessary Python dependencies are installed for operations that require Python. Default: 'seuratextend'.
#' @return A ggplot object, displaying gene expression trends along computed pseudotime trajectories by heatmap.
#' @details These visualization tools are designed to highlight gene expression changes over pseudotime, providing insights into cellular behaviors and gene regulation during differentiation processes. Each function is tailored to the specific pseudotime analysis, ensuring that visualizations are closely aligned with the underlying data and pseudotime results.
#' @examples
#' library(Seurat)
#' library(SeuratExtend)
#'
#' # Load an example Seurat Object
#' mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))
#'
#' # Calculate diffusion map and pseudotime using Palantir
#' mye_small <- Palantir.RunDM(mye_small)
#' mye_small <- Palantir.Pseudotime(mye_small, start_cell = "sample1_GAGAGGTAGCAGTACG-1")
#'
#' # Customizing gene trend heatmap with Palantir
#' ps <- mye_small@misc$Palantir$Pseudotime
#' colnames(ps)[3:4] <- c("fate1", "fate2")
#' GeneTrendHeatmap.Palantir(
#'   mye_small,
#'   features = c("CD14", VariableFeatures(mye_small)[1:10]),
#'   pseudotime.data = ps,
#'   lineage = "fate1"
#' )
#'
#' # Run Slingshot
#' mye_small <- RunSlingshot(mye_small, group.by = "cluster", start.clus = "Mono CD14")
#'
#' # Using Slingshot for similar visualizations
#' GeneTrendHeatmap.Slingshot(
#'   mye_small,
#'   features = c("CD14", VariableFeatures(mye_small)[1:10]),
#'   lineage = "slingPseudotime_2"
#' )
#'
#' @rdname GeneTrendHeatmap-Methods
#' @export

GeneTrendHeatmap.Palantir <- function(
    seu,
    features,
    pseudotime.data = "Pseudotime",
    lineage,
    magic = FALSE,
    scale = TRUE,
    sort = TRUE,
    color_scheme = "D",
    conda_env = "seuratextend") {

  library(Seurat)
  library(reshape2)
  library(ggplot2)

  seu <- process_magic_and_set_assay(seu, magic, conda_env)

  # Check the pseudotime.data parameter
  pseudotime_df <- extract_pseudotime_data(seu, pseudotime.data, "Palantir", required_cols = c("Pseudotime", "Entropy"))
  fate_names <- setdiff(colnames(pseudotime_df), c("Pseudotime", "Entropy"))

  # If the lineage parameter is not provided, default to the first cell fate name
  if (is.null(lineage)) {
    lineage <- fate_names[1]
  }

  # Check if the provided lineage is one of the valid cell fate names
  if (!lineage %in% fate_names) {
    stop(paste("The provided lineage '", lineage, "' is not valid. Please choose one from the following cell fates:", paste(fate_names, collapse = ", "), "."))
  }

  pseudotime_df2 <- pseudotime_df[,c("Pseudotime",lineage),drop = F]
  df <- cbind(pseudotime_df2, FetchData(seu, vars = features))

  # Melt the dataframe to get the desired format
  df_melt_gene <- melt(
    df, id.vars = colnames(pseudotime_df2),
    measure.vars = features,
    variable.name = "Gene",
    value.name = "Expression"
  )

  df_melt <- melt(
    df_melt_gene, id.vars = c("Pseudotime", "Gene", "Expression"),
    measure.vars = lineage,
    variable.name = "cell_fate",
    value.name = "fate_prob"
  )
  df_melt$Gene <- factor(df_melt$Gene, levels = unique(df_melt$Gene))

  # Trim the dataframe
  df_melt <- trim_data_by_pseudotime(df, df_melt, fate_names)

  predictions <- perform_gam_analysis(df_melt, features, lineage)

  if (scale) {
    # Apply the rescale function to the Predicted values for each gene
    predictions$Predicted <- ave(predictions$Predicted, predictions$Gene, FUN = rescale)
  }

  # Reshape the predictions dataframe to create the matrix
  heatmap_matrix <- dcast(predictions, Gene ~ Pseudotime, value.var = "Predicted")

  # Convert the dataframe to a matrix
  heatmap_matrix <- as.matrix(heatmap_matrix[, -1])
  rownames(heatmap_matrix) <- levels(predictions$Gene)

  # Sort
  if(sort) {
    library(dplyr)

    # Calculate gradient scores for each row
    gradient_scores <- apply(heatmap_matrix, 1, compute_gradient)

    sort_df <- data.frame(
      row = rownames(heatmap_matrix),
      max = apply(heatmap_matrix, 1, which.max),
      gradient = gradient_scores
    )

    # Arrange by peak position first and then by gradient score
    sort_df <- sort_df %>%
      arrange(gradient) %>%
      arrange(max)

    heatmap_matrix <- heatmap_matrix[sort_df$row,,drop = F]
  }

  lab_fill <- ifelse(scale, "Relative\nExpression", "Expression")
  p <- Heatmap(heatmap_matrix, lab_fill = lab_fill, color_scheme = color_scheme, y_text_position = "left") +
    labs(title = lineage) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.2)))
  return(p)
}

#' @rdname GeneTrendCurve-Methods
#' @export

GeneTrendCurve.Slingshot <- function(
    seu,
    features,
    pseudotime.data = "PCA",
    magic = FALSE,
    method = c("smooth","gam"),
    point = FALSE,
    line.size = NULL,
    pt.alpha = 0.3,
    pt.size = 1,
    ncol = NULL,
    se = TRUE,
    conda_env = "seuratextend") {

  library(Seurat)
  library(reshape2)
  library(ggplot2)

  check_method(method)
  seu <- process_magic_and_set_assay(seu, magic, conda_env)

  # Check the pseudotime.data parameter
  pseudotime_df <- extract_pseudotime_data(seu, pseudotime.data, "Slingshot")
  pseudotime_df <- as.data.frame(pseudotime_df)
  fate_names <- colnames(pseudotime_df)

  df <- cbind(pseudotime_df, FetchData(seu, vars = features, clean = "none"))

  # Melt the dataframe to get the desired format
  df_melt_gene <- melt(
    df, id.vars = fate_names,
    measure.vars = features,
    variable.name = "Gene",
    value.name = "Expression"
  )

  df_melt <- melt(
    df_melt_gene, id.vars = c("Gene", "Expression"),
    measure.vars = fate_names,
    variable.name = "cell_fate",
    value.name = "Pseudotime"
  )
  df_melt$Gene <- factor(df_melt$Gene, levels = unique(df_melt$Gene))

  # Trim the dataframe
  df_melt <- na.omit(df_melt)

  # Create the base plot
  p <- build_base_plot(point, df_melt, pt.alpha, pt.size, ncol)

  if (method[1] == "gam") {
    predictions <- perform_gam_analysis(df_melt, features, fate_names)
  } else {
    predictions <- df_melt
  }

  p <- add_plot_layer(p, method, line.size, se, predictions)

  return(p)
}

#' @rdname GeneTrendHeatmap-Methods
#' @export

GeneTrendHeatmap.Slingshot <- function(
    seu,
    features,
    pseudotime.data = "PCA",
    lineage,
    magic = FALSE,
    scale = TRUE,
    sort = TRUE,
    color_scheme = "D",
    conda_env = "seuratextend") {

  library(Seurat)
  library(reshape2)
  library(ggplot2)

  seu <- process_magic_and_set_assay(seu, magic, conda_env)

  # Check the pseudotime.data parameter
  pseudotime_df <- extract_pseudotime_data(seu, pseudotime.data, "Slingshot")
  pseudotime_df <- as.data.frame(pseudotime_df)
  fate_names <- colnames(pseudotime_df)

  df <- cbind(pseudotime_df, FetchData(seu, vars = features))

  # If the lineage parameter is not provided, default to the first cell fate name
  if (is.null(lineage)) {
    lineage <- fate_names[1]
  }

  # Check if the provided lineage is one of the valid cell fate names
  if (!lineage %in% fate_names) {
    stop(paste("The provided lineage '", lineage, "' is not valid. Please choose one from the following cell fates:", paste(fate_names, collapse = ", "), "."))
  }

  pseudotime_df2 <- pseudotime_df[,lineage,drop = F]
  df <- cbind(pseudotime_df2, FetchData(seu, vars = features))

  # Melt the dataframe to get the desired format
  df_melt_gene <- melt(
    df, id.vars = colnames(pseudotime_df2),
    measure.vars = features,
    variable.name = "Gene",
    value.name = "Expression"
  )

  df_melt <- melt(
    df_melt_gene, id.vars = c("Gene", "Expression"),
    measure.vars = lineage,
    variable.name = "cell_fate",
    value.name = "Pseudotime"
  )
  df_melt$Gene <- factor(df_melt$Gene, levels = unique(df_melt$Gene))

  # Trim the dataframe
  df_melt <- na.omit(df_melt)

  predictions <- perform_gam_analysis(df_melt, features, lineage)

  if (scale) {
    # Apply the rescale function to the Predicted values for each gene
    predictions$Predicted <- ave(predictions$Predicted, predictions$Gene, FUN = rescale)
  }

  # Reshape the predictions dataframe to create the matrix
  heatmap_matrix <- dcast(predictions, Gene ~ Pseudotime, value.var = "Predicted")

  # Convert the dataframe to a matrix
  heatmap_matrix <- as.matrix(heatmap_matrix[, -1])
  rownames(heatmap_matrix) <- levels(predictions$Gene)

  # Sort
  if(sort) {
    library(dplyr)

    # Calculate gradient scores for each row
    gradient_scores <- apply(heatmap_matrix, 1, compute_gradient)

    sort_df <- data.frame(
      row = rownames(heatmap_matrix),
      max = apply(heatmap_matrix, 1, which.max),
      gradient = gradient_scores
    )

    # Arrange by peak position first and then by gradient score
    sort_df <- sort_df %>%
      arrange(gradient) %>%
      arrange(max)

    heatmap_matrix <- heatmap_matrix[sort_df$row,,drop = F]
  }

  lab_fill <- ifelse(scale, "Relative\nExpression", "Expression")
  p <- Heatmap(heatmap_matrix, lab_fill = lab_fill, color_scheme = color_scheme, y_text_position = "left") +
    labs(title = lineage) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.2)))
  return(p)
}


### Internal functions

check_method <- function(method) {
  if (!method[1] %in% c("gam", "smooth")) {
    stop("Invalid method. Please choose either 'gam' or 'smooth'.")
  }
}

process_magic_and_set_assay <- function(seu, magic, conda_env) {
  if (magic && !("magic" %in% names(seu@assays))) {
    seu <- Palantir.Magic(seu, conda_env = conda_env)  # Assuming Palantir.Magic modifies 'seu' by reference.
  }
  assay_to_use <- ifelse(magic, "magic", DefaultAssay(seu))
  DefaultAssay(seu) <- assay_to_use
  return(seu)
}

perform_gam_analysis <- function(df_melt, features, fate_names) {
  library(mgcv)

  # Create a dataframe to store predicted values
  predictions <- data.frame()

  # Loop through each feature and cell fate
  for (gene in features) {
    for (fate in fate_names) {
      subset_df <- df_melt[df_melt$Gene == gene & df_melt$cell_fate == fate,]

      # Fit the GAM model using the Pseudotime column
      model <- gam(Expression ~ s(Pseudotime), data = subset_df)

      # Create a grid for pseudotime based on the range of the current fate
      pseudo_seq <- seq(min(subset_df$Pseudotime), max(subset_df$Pseudotime), length.out = 200)
      grid_df <- data.frame(Pseudotime = pseudo_seq, Gene = gene, cell_fate = fate, Expression = NA)

      # Predict values and add to the predictions dataframe
      grid_df$Predicted <- predict(model, newdata = grid_df)
      predictions <- rbind(predictions, grid_df)
    }
  }
  predictions$Gene <- factor(predictions$Gene, levels = unique(predictions$Gene))

  return(predictions)
}

build_base_plot <- function(point, df_melt, pt.alpha, pt.size, ncol) {
  p <- ggplot() +
    labs(x = "Pseudotime", y = "Gene Expression", color = "Cell Fate") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(face = "bold", size = 10)) +
    facet_wrap(~ Gene, scales = "free", ncol = ncol) +
    scale_x_continuous(expand = c(0, 0))

  if (point) {
    p <- p + geom_point(
      data = df_melt,
      aes(Pseudotime, Expression, color = cell_fate),
      alpha = pt.alpha, size = pt.size,
      shape = 16)
  }

  return(p)
}

add_plot_layer <- function(p, method, line.size, se, data_for_plot) {
  if (method[1] == "gam") {
    p <- p + geom_line(
      data = data_for_plot,
      aes(Pseudotime, Predicted, color = cell_fate),
      linewidth = line.size)
  } else if (method[1] == "smooth") {
    p <- p + geom_smooth(
      data = data_for_plot,
      aes(x = Pseudotime, y = Expression, color = cell_fate),
      linewidth = line.size, se = se)
  }

  return(p)
}

extract_pseudotime_data <- function(seu, pseudotime.data, source, required_cols = NULL) {
  # Initialize pseudotime_df
  pseudotime_df <- NULL

  if (is.character(pseudotime.data)) {
    # Extract the data frame based on the source
    if (source == "Palantir") {
      if (!is.null(seu@misc$Palantir[[pseudotime.data]])) {
        pseudotime_df <- seu@misc$Palantir[[pseudotime.data]]
      } else {
        stop("The Palantir pseudotime output not found. Please run the Palantir.Pseudotime function to generate the Palantir pseudotime output in seu@misc$Palantir with the appropriate title.")
      }
    } else if (source == "Slingshot") {
      if (!is.null(seu@misc$slingshot[[pseudotime.data]]$SlingPseudotime)) {
        pseudotime_df <- seu@misc$slingshot[[pseudotime.data]]$SlingPseudotime
      } else {
        stop("The Slingshot pseudotime output not found. Please run the RunSlingshot function to generate the Slingshot pseudotime output in seu@misc$slingshot with the appropriate reducedDim.")
      }
    }
  } else if (is.data.frame(pseudotime.data) | is(pseudotime.data, "DFrame")) {
    # Validate the data frame if it's directly provided
    if (!is.null(required_cols) && !all(required_cols %in% colnames(pseudotime.data))) {
      stop(paste("The provided data frame does not contain the required columns:", paste(required_cols, collapse = ", "), "."))
    }
    if (!identical(rownames(pseudotime.data), colnames(seu))) {
      stop("The row names of the provided data frame do not match the cell names of the Seurat object.")
    }
    pseudotime_df <- pseudotime.data
  } else {
    stop(paste("The pseudotime.data parameter must either be a character specifying the title in seu@misc$", source, " or a data frame with the required columns and rows.", sep = ""))
  }

  return(pseudotime_df)
}

trim_data_by_pseudotime <- function(df, df_melt, fate_names) {
  # Determine the 5% cells with the lowest Pseudotime value
  threshold <- quantile(df$Pseudotime, 0.05)
  start_cells <- df[df$Pseudotime <= threshold, ]

  # Trim the dataframe based on the minimum fate values in the start cells
  for (fate in fate_names) {
    min_fate_value <- min(start_cells[[fate]])
    df_melt <- df_melt[!(df_melt$cell_fate == fate & df_melt$fate_prob < min_fate_value), ]
  }

  return(df_melt)
}

rescale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

compute_gradient <- function(row, n=10) {
  peak_position <- which.max(row)
  left_bound <- max(1, peak_position - n)
  right_bound <- min(length(row), peak_position + n)

  # Extract the sub-section of the row around the peak
  sub_section <- row[left_bound:right_bound]

  # Calculate the gradient score as the average difference between subsequent positions
  gradient_score <- mean(diff(sub_section))

  return(gradient_score)
}
