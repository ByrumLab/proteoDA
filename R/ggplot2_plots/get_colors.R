

## group must be a vector
## gets colors for each sample. if less than 12 groups 
## then standard colorGroup() used. PolyChrome used if 
## large number of groups
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

