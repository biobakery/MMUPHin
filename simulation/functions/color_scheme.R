gg_color_hue <- function(values = NULL, n = NULL) {
  if(!is.null(values)) {
    if(anyDuplicated(values)) stop("Values should be unique if provided!")
    if(!is.character(values)) stop("Values should be of character class!")
    n <- length(values)
    hues <- seq(15, 375, length = n + 1)
    colors <- hcl(h = hues, l = 65, c = 100)[1:n]
    names(colors) <- values
    return(colors)
  }
  hues <- seq(15, 375, length = n + 1)
  colors <- hcl(h = hues, l = 65, c = 100)[1:n]
}
colors <- gg_color_hue(n=3)
colors_simulation_adjust.batch <- c("Original" = "black",
                                    "ComBat corrected" = colors[3],
                                    "MMUPHin corrected" = colors[1])
colors_simulation_unsupervised <- c("Original" = "black",
                                    "MMUPHin corrected" = colors[1])
