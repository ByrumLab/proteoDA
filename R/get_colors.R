#library(Polychrome)

#' Generate colors for Sample Groups
#' 
#' Assigns colors to each unique sample group in metadata. If the number of groups is 12 or fewer,
#' colors are assigned using a color-blind friendly palette. For more than 12 groups, the function
#' uses the Polychrome package to generate a distinct color palette.
#' 
#' @param group A vector representing the sample group classifications. 
#' If not a factor, it is converted into one. 
#' 
#' @return A named vector of colors, where names correspond to unique group labels. 
#' 
#' @details
#' - If there are 12 or fewer unique groups, the function assigns colors using a predefined
#'   color-blind friendly palette.
#' - If there are more than 12 groups, it leverages Polychrome to generate distinct colors.
#' - If grDevices is not installed, an error is thrown when handling more than 12 groups. 
#' - group must be a vector.
#' 
#' @importFrom Polychrome createPalette
#' @importFrom cli cli_abort 
#' 
#' @examples
#' 
#' group <- factor(c("A", "B", "C"))
#' get_colors(group)
#' 
#' @export

get_colors <- function(group) {
  
  if (!is.factor(group)) {
    group <- factor(group)
  }
  
  group <- droplevels(group)
  groups <- levels(group)
  num_groups <- length(groups)
  
  # Define a color-blind friendly palette (12 colors)
  cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                  "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", 
                  "#009E73", "#F0E442")
  
  if (num_groups <= 12L) {
    colors <- cb_palette[1:num_groups]
    names(colors) <- groups
  } else {
    
    if (!requireNamespace("grDevices", quietly = TRUE)) {
      cli::cli_abort("Package \"grDevices\" must be installed to make plots for more than 12 groups")
    }
    
    colors <- Polychrome::createPalette(N = num_groups, seedcolors = cb_palette)
    names(colors) <- groups
  }
  
  return(colors)
}
