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
#' @param col.space Color space, Default: c("default", "intense", "pastel")
#' @param set Several random presets, Default: 1
#' @return Vector of colors
#' @details Random color presets generated from: https://medialab.github.io/iwanthue/
#' @examples
#' color_iwh(10)
#' @rdname color_iwh
#' @export

color_iwh <- function(n, col.space = c("default","intense","pastel"), set = 1) {
  if(getOption("color_iwh.warning",TRUE)) {
    message("The 'I want hue' color presets were generated from: https://medialab.github.io/iwanthue/\n",
            "This message is shown once per session")
    options("color_iwh.warning"=FALSE)
  }

  if(!col.space[1] %in% c("default","intense","pastel", 1:3)) {
    message('"col.space" should be "default", "intense" or "pastel"\n',
            'Using "default" preset.')
    col.space <- "default"
  }
  col.space <- col.space[1]
  if(n > 20 & n < 31 & col.space %in% c("default",1)) {
    message('Too many colors (n > 20). Using "intense" preset instead.')
    col.space <- "intense"
  }
  n.range <- names(color_presets[[col.space]])
  if(!n %in% n.range){
    stop('"n" should be in range (',n.range[1],'~',tail(n.range, 1),') for preset "',col.space,'".')
  }
  set.range <- names(color_presets[[col.space]][[as.character(n)]])
  if(!set %in% set.range){
    stop('"set" should be in range (',set.range[1],'~',tail(set.range, 1),') for ',
         n,' colors in preset "', col.space,'".')
  }
  cols <- color_presets[[col.space]][[as.character(n)]][[set]]
  return(cols)
}

