
## data = matrix or df of intensities
## groups = vector class membership
## sample_labels =  vector names ot display on plot
## title = character; plot title text
## text.sizes = title, x-axis, y-axis, x/y-axis labels, legend
qc_boxplot_plot5 <- function (data,
                            groups = NULL,
                            sample_labels = NULL,
                            title = NULL,
                            text.sizes = c(12,10,10,10),
                            # legend=TRUE,
                            legend.position="right"){


   ## check arg class types
  # testthat::expect_s3_class(object = data, class = c("matrix","data.frame"))
  # testthat::expect_is(object = groups, class = c("numeric","integer", "character","factor","NULL"))
  # testthat::expect_s3_class(object = sample_labels, class = c("character","NULL"))
  # testthat::expect_s3_class(object = title, class = c("numeric","integer","character","NULL"))
  # testthat::expect_s3_class(object = text.sizes, class = c("numeric","integer"))

   # Prep args
   if (is.null(groups)) { groups <- c(rep("1",ncol(data))) }
   groups <- proteoDA:::make_factor(as.character(groups), prefix = NULL)

   if (is.null(sample_labels)) { sample_labels <- colnames(data) }
   if (is.null(title)) { title <- " " }


   ## check number arg values match no cols data
   testthat::expect_equal(object = length(groups), expected = ncol(data))
   testthat::expect_equal(object = length(sample_labels), expected = ncol(data))
   testthat::expect_equal(object = length(text.sizes), expected = 4)

   ## groups is sorted according to group levels (index order returned)
   group_order <- sort.int(groups, index.return = T)$ix


   # Then, reorder the sample labels and groups, and data columns to match group levels
   sample_labels <- sample_labels[group_order]
   groups       <- groups[group_order]
   data         <- data[, group_order]

   ## targets
   plot.meta <- data.frame(ind = proteoDA:::make_factor(colnames(data)),
                           labels = proteoDA:::make_factor(sample_labels),
                           group = groups)

   #             ind labels       group
   # 1        DMSO_1     s1        DMSO
   # 2        DMSO_2     s5 inhibitor76
   # 3        DMSO_3     s9 inhibitor20
   # 4 inhibitor76_1     s2          KD


   ## long format (ind labels group   values)
   plot.data=merge(plot.meta, utils::stack(as.data.frame(data)), sort = F)

   #      ind labels group   values
   # 1 DMSO_1     s1  DMSO 17.93881
   # 2 DMSO_1     s1  DMSO 16.71879
   # 3 DMSO_1     s1  DMSO 19.53171
   # 4 DMSO_1     s1  DMSO 24.35849
   # 5 DMSO_1     s1  DMSO 14.98676
   # 6 DMSO_2     s5  DMSO 17.96687


   p1 <- ggplot(plot.data, aes(x = labels, y = values, fill = group)) +
          #  geom_violin(draw_quantiles = c(0.5), na.rm = T, col = "black") +
            geom_boxplot(width = 0.2, color="black", alpha=0.2) +
            scale_fill_manual(values = proteoDA:::colorGroup(plot.meta$group), name = NULL) +
            labs(y = "Density", x="", fill = "") +
            # ggtitle(title) +

            # theme_gray() +
            # theme(axis.title.x = element_blank(),
            #       axis.text.x  = element_text( size=text.sizes[3], angle = 45, vjust = 0.9, hjust = 1),
            #       axis.title.y = element_text(size=text.sizes[2], angle = 90, vjust=0.5, margin = unit(c(0,0,0,0), "line")),
            #       axis.text.y  = element_text(size=text.sizes[3]),
            #       plot.title   = element_text(size=text.sizes[1], face = "bold")) +
            # theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))

         # if(legend == TRUE){
         #    p1 <- p1 + theme(legend.position = legend.position, legend.text = element_text(size=text.sizes[4])) +
         #               theme(legend.key.size = unit(0.5, 'cm'))
         # }
         # if(legend == FALSE){
         #    p1 <- p1 + theme(legend.position = "none") +
         #               theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))
         # }
   theme_gray() +
     theme(axis.text.x  = element_text(angle = 45, vjust = 0.9, hjust = 1, size=text.sizes[3]),
           axis.title.x = element_blank(),
           axis.title.y = element_text(size=text.sizes[2]),
           axis.text.y  = element_text(size=text.sizes[3]),
           plot.title   = element_text(size=text.sizes[1])) +

     # geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
     theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
     theme(legend.position = legend.position, legend.text = element_text(size=text.sizes[4])) +
     theme(legend.key.size = unit(0.5, 'cm')) +
     ggtitle(title)

   # ggplotify::as.ggplot(p1)
  p1
}



