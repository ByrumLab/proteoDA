
# ext2 <- proteoDA::normalize_data(ext, "log2")
# colors = mixOmics::color.mixo(length(unique(ext2$metadata$group)):1)
# norm2<- proteoDA::normalize_data(DAList = filt2,norm_method = "log2")


# qc_totInt_by_group_ggplot7(DAList = norm,
#                            label_column = "sample",
#                            grouping_column = "group",
#                            colors = mixOmics::color.mixo(length(unique(raw$metadata$group)):1)
# )
#
# qc_totInt_ggplot7(DAList = norm,
#                   label_column = "sample2",
#                   grouping_column = "group",
#                   colors =  mixOmics::color.mixo(length(unique(raw$metadata$group)):1),
#                   legend.position = "right")

# d <- qc_totInt_by_group_ggplot7(DAList = log2Raw,label_column = "sample3",grouping_column = "block3");
# pdf(file = "totInt.pdf", height=15, width=13); d; dev.off()




#######################################################################
## ADD THIS CODE TO CONTROL SIZE OF TEXT LABELS
##################################################################
          #   theme(axis.text.x  = element_text(angle = 45, vjust = 0.9, hjust = 1, size=text.sizes[3]),
          #         axis.title.x = element_blank(),
          #         axis.title.y = element_text(size=text.sizes[2]),
          #         axis.text.y  = element_text(size=text.sizes[3]),
          #         axis.text.y  = element_text(size=text.sizes[3]),
          #         plot.title   = element_text(size=text.sizes[1])) +
          #
          #   theme(legend.position = legend.position, legend.text = element_text(size=text.sizes[4]))
          #   geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
          #  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +




## individual barplots for each group
##
## percentile A double [0-1] used to identifiy the specified percentile of the data.
## takes the total intensity values for the samples, and identifies the intensity
## threshold for the say 10th percentile of the data. samples with values below
## the percentile are flagged as potential outliers. so it identifies samples with
## the lowest total intensities in the data. these samples may not have sequenced
## well and will likely have a lot of missing values. the horizontal line
## on the plots corresponds to the intensity threshold for the specified percentile
qc_totInt_by_group_ggplot7 <- function(DAList,
                                       label_column,
                                       grouping_column,
                                       percentile,
                                       colors = NULL, ## one per group
                                       nrow = NULL,
                                       ncol = NULL,
                                       legend.position = "none"){


  legend.position <- rlang::arg_match(arg = legend.position,
                                      values = c("none","right","top",
                                                 "left","bottom"),
                                      multiple = FALSE)

  label_column <- check_label_column(metadata = DAList$metadata,
                                     label_column = label_column)

  grouping_column <- check_grouping_column(metadata = DAList$metadata,
                                           grouping_column = grouping_column)

  ## check percentile parameter
  if(!check_num(x = percentile)){
    cli::cli_abort("percentile should be a numeric value between 0 and 1")
  }

  if(any(percentile < 0L,
         percentile > 1L)){
    cli::cli_abort("percentile should be a numeric value between 0 and 1")
  }


  ## create data.frame with sample labels, totalIntensity, and group column
  bar_data <- data.frame(ind = colnames(DAList$data),
                         label    = DAList$metadata[, label_column],
                         tot.int = as.numeric(colSums(DAList$data, na.rm = TRUE)),
                         tot.num = as.numeric(colSums(!is.na(DAList$data),
                                                      na.rm = FALSE)),
                         tot.na = as.numeric(colSums(is.na(DAList$data),
                                                     na.rm = FALSE)),
                         group    = DAList$metadata[, grouping_column],
                         check.names = FALSE,
                         fix.empty.names = FALSE,
                         stringsAsFactors = FALSE,
                         check.rows = FALSE)


  ## make group column a factor
  ## sort samples by group levels
  bar_data$group <- proteoDA:::make_factor(bar_data$group)
  bar_data <- bar_data[order(bar_data$group), ]

 ## number of groups
 num_groups <- length(unique(bar_data$group))

 ## make sample name a factor so they are not re-organized A-Z in plot
 bar_data$label <- proteoDA:::make_factor(bar_data$label)


 if(percentile == 0L | is.null(percentile)){
   perc_int_thresh <- 0
 } else {
  perc_int_thresh <- quantile(x = bar_data$tot.int, probs = percentile)
 }
 print(paste("{.val {percentile *100}}th percentile total intensity threshold = ",
                   perc_int_thresh))

 bar_data$tot.int.outlier <- bar_data$tot.int <= perc_int_thresh


 ## if colors is not supplied
 if(is.null(colors)){

   ## num_groups <= 12L
   ## if less than 12 groups use colorGroup
   if(num_groups <= 12L){

     colors <- proteoDA:::colorGroup(group = bar_data$group)
     # Polychrome::swatch(colorset = colors,
     #                    main = paste("proteoDA:::colorGroup(N =", num_groups,
     #                                 ", seedcolors = binfcolors)"))
     names(colors) <- c(levels(bar_data$group))

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
     names(colors) <- c(levels(bar_data$group))
     # Polychrome::swatch(colorset = colors,
     #                    main = paste("createPalette(N =", num_groups,
     #                                 ", seedcolors = binfcolors)"))
   }
 }



 ## check number colors equals number of unique groups
 if(length(colors) != length(levels(bar_data$group))){
  cli::cli_abort("number of colors does not match the number
                 of unique groups")
 }


 # colors <- as.character(colors)
 names(colors) <- c(levels(bar_data$group))


 ## total intensity
 p <- ggplot(bar_data, aes(x = label, y = tot.int, fill = group)) +
       geom_bar(stat = "identity") +
       scale_fill_manual(values = c(colors)) +
       scale_y_continuous() +
   geom_hline(yintercept = perc_int_thresh, linetype = "solid", color = "black", linewidth = 0.5) +
       labs(y = "Total Intensity", fill = "group") +
       theme_bw() +
       theme(axis.text.x = element_text(angle = 90)) +
       facet_wrap(~group, scales = "free_x", shrink = F, nrow = nrow, ncol = ncol) +
       theme(legend.position = legend.position) +
       guides(fill=guide_legend(title=grouping_column))
 p


 ## total number with intensity != NA
 p2 <- ggplot(bar_data, aes(x = label, y = tot.num, fill = group)) +
        geom_bar(stat = "identity", position="dodge") +
        ## adds numbers to top of each bar
        # geom_text(aes(label=tot.num), position=position_dodge(width=0.9), vjust=-0.5, size = 1.5) +
        scale_fill_manual(values = c(colors)) +
        scale_y_continuous() +
        labs(y = "Total Number (not NA)", fill = "group") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.9, hjust = 1, size = 10)) +
        facet_wrap(~group, scales = "free_x", shrink = F, nrow =nrow, ncol = ncol) +
        theme(legend.position = legend.position) +
        guides(fill=guide_legend(title=grouping_column))
 p2








   # p <- ggplot2::ggplot(data = input, aes(x = component, y = value, title=" ", fill = dummy)) +
   #   geom_bar(position = 'dodge', stat='identity') +
   #   geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
   #   scale_fill_manual(values = c("dodgerblue3")) +
   #   # labs(y = "Weighted average proportion variance", x="Effect", fill = "variable title") +
   #   labs(y = "Weighted average proportion variance", x="Effect", fill = "variable title") +
   #   theme_gray() +
   #   theme(axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 1,size=12),
   #         axis.text.y = element_text(size=12),
   #         plot.title = element_text(size=18),
   #         legend.position = "none") +
   #   ggtitle(label = title)






 ## total number proteins in each with NA
 p3 <- ggplot2::ggplot(bar_data, ggplot2::aes(x = label, y = tot.na, fill = group)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c(colors)) +
        scale_y_continuous() +
        labs(y = "Total Number (is NA)", fill = "group") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_wrap(~group, scales = "free_x", shrink = F, nrow =nrow, ncol = ncol) +
        theme(legend.position = legend.position) +
        guides(fill=guide_legend(title=grouping_column))
 p3


 out <- list(tot.int = p, tot.num=p2, tot.na=p3,
             bar_data = bar_data, colors=colors, percentile=percentile,
             perc_int_thresh=perc_int_thresh)
 out

}




## one big barplot samples organized/colored by group
##
## percentile A double [0-1] used to identifiy the specified percentile of the data.
## takes the total intensity values for the samples, and identifies the intensity
## threshold for the say 10th percentile of the data. samples with values below
## the percentile are flagged as potential outliers. so it identifies samples with
## the lowest total intensities in the data. these samples may not have sequenced
## well and will likely have a lot of missing values. the horizontal line
## on the plots corresponds to the intensity threshold for the specified percentile
qc_totInt_ggplot7 <- function(DAList,
                              label_column,
                              grouping_column,
                              percentile,
                              colors=NULL, ## one per group
                              legend.position = "right"){


 legend.position <- rlang::arg_match(arg = legend.position,
                                     values = c("none","right","top",
                                                "left","bottom"),
                                     multiple = FALSE)

 label_column <- check_label_column(metadata = DAList$metadata,
                                    label_column = label_column)

 grouping_column <- check_grouping_column(metadata = DAList$metadata,
                                          grouping_column = grouping_column)

  ## check percentile parameter
 if(!check_num(x = percentile)){
   cli::cli_abort("percentile should be a numeric value between 0 and 1")
 }

 if(any(percentile < 0L,
        percentile > 1L)){
   cli::cli_abort("percentile should be a numeric value between 0 and 1")
 }


 ## create data.frame with sample labels, totalIntensity, and group column
 bar_data <- data.frame(label    = DAList$metadata[, label_column],
                        tot.int = as.numeric(colSums(DAList$data, na.rm = TRUE)),
                        tot.num = as.numeric(colSums(!is.na(DAList$data),
                                                     na.rm = FALSE)),
                        group    = DAList$metadata[, grouping_column],
                        check.names = FALSE,
                        fix.empty.names = FALSE,
                        stringsAsFactors = FALSE,
                        check.rows = FALSE)


 ## make group column a factor
 ## sort samples by group levels
 bar_data$group <- proteoDA:::make_factor(bar_data$group)
 bar_data <- bar_data[order(bar_data$group), ]

 ## number of groups
 num_groups <- length(unique(bar_data$group))

 ## make sample name a factor so they are not re-organized A-Z in plot
 bar_data$label <- proteoDA:::make_factor(bar_data$label)

 ## percentile A double [0-1] used to identifiy the specified percentile of the data.
 ## takes the total intensity values for the samples, and identifies the intensity
 ## threshold for the say 10th percentile of the data. samples with values below
 ## the percentile are flagged as potential outliers. so it identifies samples with
 ## the lowest total intensities in the data. these samples may not have sequenced
 ## well and will likely have a lot of missing values. the horizontal line
 ## on the plots corresponds to the intensity threshold for the specified percentile

 if(percentile == 0L){
   perc_int_thresh <- 0
 } else {
   perc_int_thresh <- quantile(x = bar_data$tot.int, probs = percentile)
 }
  print(paste("intensity below 10th percentile = ", perc_int_thresh))

 ## samples with total intensity below input percentile threshold
 ## are flagged as outliers
 bar_data$tot.int.outlier <- bar_data$tot.int <= perc_int_thresh

 num_outliers <- sum(bar_data$tot.int.outlier)

 ## if colors is not supplied
 if(is.null(colors)){

  ## num_groups <= 12L
  ## if less than 12 groups use colorGroup
  if(num_groups <= 12L){

   colors <- proteoDA:::colorGroup(group = bar_data$group)
   # Polychrome::swatch(colorset = colors,
   #                    main = paste("proteoDA:::colorGroup(N =", num_groups,
   #                                 ", seedcolors = binfcolors)"))
   names(colors) <- c(levels(bar_data$group))

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
   names(colors) <- c(levels(bar_data$group))
   # Polychrome::swatch(colorset = colors,
   #                    main = paste("createPalette(N =", num_groups,
   #                                 ", seedcolors = binfcolors)"))
  }
 }

 ## check number colors equals number of unique groups
 if(length(colors) != length(levels(bar_data$group))){
  cli::cli_abort("number of colors does not match the number
                 of unique groups")
 }
 names(colors) <- c(levels(bar_data$group))


 p <- ggplot(bar_data, aes(x = label, y = tot.int, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(colors)) +
  geom_hline(yintercept = perc_int_thresh, linetype = "solid", color = "black", linewidth = 0.5) +
  labs(y = "Total Intensity", fill = "group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = legend.position) +
  guides(fill = guide_legend(title=grouping_column))
 p


  out <- list(p = p, bar_data = bar_data, colors=colors, percentile=percentile,
              perc_int_thresh=perc_int_thresh)
 out


}





