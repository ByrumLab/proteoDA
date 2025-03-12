#' Generate colors for Sample Groups
#' 
#' Assigns colors to each unique sample group in metadata. If the number of groups is 12 or fewer,
#' colors are assigned using \code{\link{colorGroup}}. For more than 12 groups, the function
#' uses the \pkg{Polychrome} to generate a distinct color palette.
#' 
#' @param group A vector representing the sample group classifications. 
#' If not a factor, it is converted into one. 
#' 
#' @return A named vector of colors, where names correspond to unique group labels. 
#' 
#' @details
#' - If there are 12 or fewer unique groups, the function assigns colors using \code{\link{colorGroup}}.
#' - If there are more than 12 groups, it leverages \pkg{Polychrome} to generate distinct colors. The function 
#'  uses \code{\link{binfcolors}} as seed colors to create a new color palette 
#'  instead of creating rainbow color palette using colorGroup() function
#' - If \pkg{grDevices} is not installed, an error is thrown when handling more than 12 groups. 
#' - group must be a vector
#' 
#' @importFrom Polychrome createPalette
#' @importFrom cli cli_abort 
#'
#' @examples
#' \dontrun {
#' # example code
#' group <- factor(c("A","B","C"))
#' get_colors(group) 
#' 
#' # example2 
#' group_color <- get_colors(group = results$metadata$group)
#' 
#' }
#' 
#' @export

get_colors <- function(group){
 
 if(!is.factor(group)){
  group <- proteoDA:::make_factor(as.character(group))
 }
 group  <- droplevels(group)
 groups <- levels(group)
 num_groups <- length(groups)
 
 ## num_groups <= 12L
 ## if less than 12 groups use colorGroup
 if(num_groups <= 12L){
  
  colors <- proteoDA:::colorGroup(group = group)
  # Polychrome::swatch(colorset = colors,
  #                    main = paste("proteoDA:::colorGroup(N =", num_groups,
  #                                 ", seedcolors = binfcolors)"))
  names(colors) <- c(levels(group))
  
 } else {
  
  ## num_groups > 12L (by-pass rainbow colors)
  ## if > 12 groups use binfcolors as seed colors to create a new
  ## color palette using Polychrome package instead of creating rainbow
  ## color palette using colorGroup() function
  if (!requireNamespace("grDevices", quietly = TRUE)) {
   cli::cli_abort(c("Package \"grDevices\" must be installed to
                        make plots for more than 12 groups"))
  }
  
  colors <- Polychrome::createPalette(N=num_groups,
                                      seedcolors = proteoDA:::binfcolors)
  names(colors) <- c(levels(group))
  # Polychrome::swatch(colorset = colors,
  #                    main = paste("createPalette(N =", num_groups,
  #                                 ", seedcolors = binfcolors)"))
 }
 
 colors
 
}

