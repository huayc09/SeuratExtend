#' @title Default discrete color presets by 'I want hue'
#' @description Generate color presets from 'I want hue' online tool
#' @param n How many colors to generate
#' @param col.space Color space, Options: "default", "intense", "pastel",
#' "all" (k-Means) or "all_hard" (force vector)
#' @param set Several random presets, Default: 1
#' @return Vector of colors
#' @details
#' Random color presets generated from: https://medialab.github.io/iwanthue/
#' @examples
#' color_iwh(10)
#' @rdname color_iwh
#' @export

color_iwh <- function(
  n,
  col.space = c("default","intense","pastel","all","all_hard"),
  set = 1
) {
  message_only_once("color_iwh","The 'I want hue' color presets were generated from: https://medialab.github.io/iwanthue/")

  if(!col.space[1] %in% c("default","intense","pastel","all","all_hard", 1:5)) {
    message('"col.space" should be "default", "intense", "pastel", "all" or "all_hard"\n',
            'Using "default" preset.')
    col.space <- "default"
  }
  col.space <- col.space[1]
  if(!is.numeric(n)) stop("'n' should be integer (2-50)")
  if(n > 20 & n < 31 & col.space %in% c("default",1)) {
    message('Input color number is ', n, '. Using "intense" preset instead.')
    col.space <- "intense"
  }else if(n > 30 & n < 51 & col.space %in% c("default",1)) {
    message('Input color number is ', n, '. Using "all" preset instead.')
    col.space <- "all"
  }else if(n > 50 & col.space %in% c("default",1)) {
    stop("Input color number is ", n, ". Too many colors (n > 50).")
  }
  n.range <- names(presets_color_iwh[[col.space]])
  if(!n %in% n.range){
    stop('"n" should be integer and in range (',
         n.range[1],'~',tail(n.range, 1),') for preset "',col.space,'".')
  }
  set.range <- names(presets_color_iwh[[col.space]][[as.character(n)]])
  if(!set %in% set.range){
    stop('"set" should be in range (',set.range[1],'~',tail(set.range, 1),') for ',
         n,' colors in preset "', col.space,'".')
  }
  cols <- presets_color_iwh[[col.space]][[as.character(n)]][[set]]
  return(cols)
}

#' @title Professional discrete color presets by 'I want hue'
#' @description Generate professional color presets from 'I want hue' color tool
#' @param n How many colors to generate
#' @param col.space Color space, Options: "default", "light", "red", "yellow", "green", "blue", or "purple"
#' @param sort Sort colors by hue (default) or differentiation, Options: "hue", "diff"
#' @param set Several random presets, Default: 1
#' @return Vector of colors
#' @details
#' Random color presets generated from: https://medialab.github.io/iwanthue/
#' @examples
#' color_pro(10)
#' @rdname color_pro
#' @export

color_pro <- function(
    n,
    col.space = c("default","light","red","yellow","green","blue","purple"),
    sort = c("hue","diff"),
    set = 1
) {
  message_only_once("color_iwh","The 'I want hue' color presets were generated from: https://medialab.github.io/iwanthue/")

  if(!col.space[1] %in% c("default","light","red","yellow","green","blue","purple", 1:7)) {
    message('"col.space" should be "default", "light", "red", "yellow", "green", "blue", or "purple"\nUsing "default" preset.')
    col.space <- "default"
  }
  col.space <- col.space[1]

  if(!sort[1] %in% c("hue","diff", 1:2)) {
    message('"sort" should be "hue" or "diff". Sorting by "hue"')
    sort <- "hue"
  }
  sort <- sort[1]

  cols_set <- presets_color_pro[[col.space]][[sort]]

  n.range <- names(cols_set)
  if(!n %in% n.range | !is.numeric(n)){
    stop('"n" should be integer and in range (',
         n.range[1],'~',tail(n.range, 1),') for preset "',col.space,'".')
  }
  cols_set <- cols_set[[as.character(n)]]

  set.range <- length(cols_set)
  if(!set %in% 1:set.range){
    stop('"set" should be in range (1~',set.range,') for ',
         n,' colors in preset "', col.space,'".')
  }
  cols <- cols_set[[set]]
  return(cols)
}

#' @title Mix RYB color
#' @description Mix 3 primary colors (red, yellow and blue) and
#' return RGB code
#' @param ryb numeric vector range in [0,1],
#' e.g. c(r = 0.3, y = 0.5, b = 0.2)
#' @return RGB code
#' @details see above
#' @examples
#' library(scales)
#' library(dplyr)
#' data.frame(
#'   red = c(1,0,0),
#'   yellow = c(0,1,0),
#'   blue = c(0,0,1),
#'   orange = c(1,1,0),
#'   purple = c(1,0,1),
#'   green = c(0,1,1),
#'   black = c(1,1,1),
#'   grey = c(0.5,0.5,0.5),
#'   white = c(0,0,0)
#' ) %>%
#'   apply(2, ryb2rgb) %>%
#'   show_col()
#' @rdname ryb2rgb
#' @export

ryb2rgb <- function(ryb = c(0,0,0)){
  ryb2 <- sca(ryb)
  x <- ang(ryb2)
  tar <- seg.x(x)
  tar <- sca(c(tar,0,1), range = c(1-max(ryb), 1-min(ryb)))[1:3]
  col <- rgb(tar[1], tar[2], tar[3])
  return(col)
}

ryb_ref <- data.frame(
  # 0 = red, 2 = yellow, 4 = blue, 6 = red
  x = c(0,1,2,3,4,5,6),
  r = c(    1,  1,   1,   0,  0, 0.625,    1),
  y = c(0.125,0.5,0.75,0.75,0.5, 0.125,0.125),
  b = c(0.125,  0,   0,   0,  1,0.9375,0.125)
)
seg <- function(range, weight = c(1,1)){
  len <- range[,2] - range[,1]
  rat <- weight[2] / (weight[1] + weight[2])
  tar <- range[,1] + len * rat
  return(tar)
}
seg.x <- function(x){
  # x should be 0 to 6
  x2 <- c(floor(x),ceiling(x))
  if(x2[1] == x2[2]){
    res <- as.numeric(ryb_ref[ryb_ref$x == x, 2:4])
  }else{
    range <- t(ryb_ref[ryb_ref$x %in% x2, 2:4])
    res <- seg(range, weight = c(x2[2] - x, x - x2[1]))
  }
  return(res)
}
sca <- function(vec, range = c(0,1)){
  len.r <- range[2] - range[1]
  len.v <- max(vec) - min(vec)
  if(len.v == 0){
    res <- rep(mean(range), times = length(vec))
  }else{
    rat <- (vec - min(vec)) / len.v
    res <- range[1] + len.r * rat
  }
  return(res)
}
ang <- function(ryb = c(0,0.5,1)){
  # min should be 0 and max should be 1
  min.ryb <- which.min(ryb)
  range <- switch(
    min.ryb,
    "1" = c(2,4),
    "2" = c(6,4),
    "3" = c(0,2)
  )
  x <- seg(t(range), ryb[-min.ryb])
  return(x)
}

# internal functions for fill/color themes

scale_fill_cont_auto <- function(color_scheme) {
  if(is.null(color_scheme)) return(NULL)
  library(ggplot2)
  if(any(color_scheme %in% c("A","B","C","D","E"))) {
    import("viridis")
    cols <- scale_fill_viridis(option = color_scheme)
  }else if(!is.na(color_scheme["mid"])) {
    cols <- scale_fill_gradient2(low = color_scheme["low"],
                                 mid = color_scheme["mid"],
                                 high = color_scheme["high"])
  }else if(all(!is.na(color_scheme[c("low","high")]))) {
    cols <- scale_fill_gradient(low = color_scheme["low"],
                                high = color_scheme["high"])
  }else{
    cols <- scale_fill_gradientn(colors = color_scheme)
  }
  return(cols)
}

scale_color_cont_auto <- function(color_scheme) {
  if(is.null(color_scheme)) return(NULL)
  library(ggplot2)
  if(any(color_scheme %in% c("A","B","C","D","E"))) {
    import("viridis")
    cols <- scale_color_viridis(option = color_scheme)
  }else if(!is.na(color_scheme["mid"])) {
    cols <- scale_color_gradient2(low = color_scheme["low"],
                                  mid = color_scheme["mid"],
                                  high = color_scheme["high"])
  }else if(all(!is.na(color_scheme[c("low","high")]))) {
    cols <- scale_color_gradient(low = color_scheme["low"],
                                 high = color_scheme["high"])
  }else{
    cols <- scale_color_gradientn(colors = color_scheme)
  }
  return(cols)
}

scale_color_disc_auto <- function(color_scheme, n, labels = waiver()) {
  if(is.null(color_scheme)) return(NULL)
  library(ggplot2)
  pro_theme <- c("default","light","red","yellow","green","blue","purple")
  iwh_theme <- c("default","intense","pastel","all","all_hard")
  if(length(color_scheme) == 1 & n <= 50) {
    if(color_scheme %in% c("default","light",paste0("pro_",pro_theme))) {
      color_scheme <- sub("pro_","",color_scheme)
      if(n == 1 && color_scheme == "default") return(NULL)
      cols <- color_pro(n, col.space = color_scheme)
      cols <- scale_color_manual(values = cols, labels = labels)
    } else if(color_scheme %in% paste0("iwh_",iwh_theme)) {
      color_scheme <- sub("iwh_","",color_scheme)
      cols <- color_iwh(n, col.space = color_scheme)
      cols <- scale_color_manual(values = cols, labels = labels)
    } else {
      cols <- scale_color_brewer(palette = color_scheme, labels = labels)
    }
  } else {
    cols <- scale_color_manual(values = color_scheme, labels = labels)
  }
  return(cols)
}

scale_fill_disc_auto <- function(color_scheme, n) {
  if(is.null(color_scheme)) return(NULL)
  library(ggplot2)
  pro_theme <- c("default","light","red","yellow","green","blue","purple")
  iwh_theme <- c("default","intense","pastel","all","all_hard")
  if(length(color_scheme) == 1 & n <= 50) {
    if(color_scheme %in% c("default","light",paste0("pro_",pro_theme))) {
      color_scheme <- sub("pro_","",color_scheme)
      if(n == 1 && color_scheme == "default") return(NULL)
      cols <- color_pro(n, col.space = color_scheme)
      cols <- scale_fill_manual(values = cols)
    } else if(color_scheme %in% paste0("iwh_",iwh_theme)) {
      color_scheme <- sub("iwh_","",color_scheme)
      cols <- color_iwh(n, col.space = color_scheme)
      cols <- scale_fill_manual(values = cols)
    } else {
      cols <- scale_fill_brewer(palette = color_scheme)
    }
  } else {
    cols <- scale_fill_manual(values = color_scheme)
  }
  return(cols)
}

#' @title Save Custom Color Settings to a Seurat Object
#' @description This function stores custom color settings within a Seurat object, allowing these settings to be reused across various visualization functions within the Seurat environment.
#' @param seu A Seurat object to which the color settings will be saved.
#' @param col_list A list of color settings where names correspond to variable names and values are the colors to be associated with each variable.
#' @return The modified Seurat object with updated color settings in the `misc` slot.
#' @details `save_colors` enhances workflow efficiency by centralizing color management within the Seurat object. This enables consistent and coherent visualizations across multiple plots, reducing the need for repeated color specifications and ensuring that plots are both aesthetically pleasing and informative.
#' @examples
#' library(SeuratExtend)
#'
#' # Define custom colors for different clusters and genes
#' custom_colors <- list(
#'   "cluster" = "pro_blue",
#'   "CD14" = "D",
#'   "CD3D" = c("#EEEEEE", "black")
#' )
#' # Save these colors to the Seurat object
#' pbmc <- save_colors(pbmc, custom_colors)
#' # Use the colors in DimPlot2
#' DimPlot2(pbmc, features = c("cluster", "CD14", "CD3D"))
#'
#' @rdname save_colors
#' @export

save_colors <- function(seu, col_list) {
  # Ensure that col_list is a list before proceeding
  if (!is.list(col_list)) {
    stop("The 'col_list' parameter must be a list. Example format:\n",
         "list(cluster = 'pro_blue', CD14 = 'D', CD3D = c('#EEEEEE', 'black'))")
  }

  # Check if 'var_colors' slot exists and is not NULL; if not, initialize as an empty list
  if (is.null(seu@misc[["var_colors"]])) {
    seu@misc[["var_colors"]] <- list()
  }

  # Update 'var_colors' with new colors without overwriting existing entries
  # Only update the entries provided in col_list
  seu@misc[["var_colors"]][names(col_list)] <- col_list

  # Optionally return the Seurat object if you want to check the updated list
  return(seu)
}
