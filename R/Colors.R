#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname gg_color_hue
#' @export
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' @title color_iwh
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
  if(getOption("color_iwh.warning",TRUE)) {
    message("The 'I want hue' color presets were generated from: https://medialab.github.io/iwanthue/\n",
            "This message is shown once per session")
    options("color_iwh.warning"=FALSE)
  }

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
  n.range <- names(color_presets[[col.space]])
  if(!n %in% n.range){
    stop('"n" should be integer and in range (',
         n.range[1],'~',tail(n.range, 1),') for preset "',col.space,'".')
  }
  set.range <- names(color_presets[[col.space]][[as.character(n)]])
  if(!set %in% set.range){
    stop('"set" should be in range (',set.range[1],'~',tail(set.range, 1),') for ',
         n,' colors in preset "', col.space,'".')
  }
  cols <- color_presets[[col.space]][[as.character(n)]][[set]]
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
  tar <- sca(tar, range = c(1-max(ryb), 1-min(ryb)))
  col <- rgb(tar[1], tar[2], tar[3])
  return(col)
}

ryb_ref <- data.frame(
  # 0 = red, 2 = yellow, 4 = blue, 6 = red
  x = c(0,1,2,3,4,5,6),
  r = c(1,1,1,0.5,0,0.5,1),
  y = c(0,0.5,1,1,0.5,0,0),
  b = c(0,0,0,0,1,1,0)
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
