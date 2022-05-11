write_limma_results <- function(model_results,
                                annotation,
                                ilab,
                                out_dir = file.path(ilab, paste0(enrich, "_analysis")),
                                enrich = c("protein", "phospho"),
                                contrasts_subdir = "per_contrast_results",
                                summary_csv = "DE_summary.csv",
                                combined_file_csv = "combined_results.csv",
                                BQ_csv = paste0(ilab, "_results_BQ.csv"),
                                spreadsheet_xlsx = NULL) {

  ##########################
  ## Setup and arg checks ##
  ##########################

  enrich <- rlang::arg_match(enrich)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }

  # Add filename validation

  # Add checks of model result rownames in annotation rownames.

  statlist <- model_results$stats_by_contrast
  data <- model_results$data

  #######################
  ## Write summary csv ##
  #######################

  model_results <- results_reb

  summary <- lapply(X = names(statlist),
                    FUN = summarize_contrast_DE,
                    contrast_res_list = statlist) %>%
    do.call("rbind", .)
  summary$pval_thresh <- model_results$min.pval
  summary$lfc_thresh <- model_results$min.lfc
  summary$p_adj_method <- model_results$adj.method

  summary_ouput_file <- file.path(out_dir, summary_csv)

  utils::write.csv(x = summary,
                   file = summary_ouput_file,
                   row.names = F)

  if (!file.exists(summary_ouput_file)) {
    cli::cli_abort(c("Failed to write summary {.path .csv} to {.path {summary_ouput_file}}"))
  }


  #############################
  ## Write per-contrast csvs ##
  #############################

  contrast_csv_success <- write_per_contrast_csvs(annotation_df = annotation,
                                                  data = data,
                                                  results_statlist = statlist,
                                                  output_dir = file.path(out_dir, contrasts_subdir))
  if (!all(contrast_csv_success)) {
    failed <- names(contrast_csv_success)[!contrast_csv_success]
    cli::cli_abort(c("Failed to write {.path .csv} results for {cli::qty(sum(!contrast_csv_success))} contrast{?s}:",
                     "!" = "{.val {failed}"))
  }

  ################################
  ## Write combined results csv ##
  ################################

  # Get common row order from the first element of the stat list
  row_order <- rownames(statlist[[1]])
  # Reorder rows and rename cols for each element of the results statlist
  results_for_combine <- lapply(X = names(statlist),
                                function(x) {
                                  tmp <- statlist[[x]][row_order, ,drop = F]
                                  colnames(tmp) <- paste(colnames(tmp), x, sep = "_")
                                  tmp
                                })
  combined_output_file <- file.path(out_dir, combined_file_csv)
  combined_results <- cbind(annotation[row_order, ],
                            data[row_order, ],
                            results_for_combine)
  utils::write.csv(x = combined_results,
                   file = combined_output_file,
                   row.names = F)

  if (!file.exists(combined_output_file)) {
    cli::cli_abort(c("Failed to write combined results {.path .csv} to {.path {combined_output_file}}"))
  }

  #################################
  ## Write combined BiqQuery csv ##
  #################################

  BQ_output_file <- file.path(out_dir, BQ_csv)
  BQ_output <- combined_results
  BQ_output[is.na(BQ_output)] <- 0
  BQ_output[BQ_output == ""] <- 0

  # colnames can't begin with numbers?
  col_beings_number <- stringr::str_detect(colnames(BQ_output), "^[:digit:]")
  if (any(col_beings_number)) {
    colnames(BQ_output)[col_beings_number] <- paste0("X", colnames(BQ_output)[col_beings_number])
  }
  utils::write.csv(BQ_output, file = BQ_output_file, row.names = F)

  if (!file.exists(BQ_output_file)) {
    cli::cli_abort(c("Failed to write BigQuery results {.path .csv} to {.path {BQ_output_file}}"))
  }

  #############################
  ## Write excel spreadsheet ##
  #############################

  # TODO: deal with this.

  ############
  ## Finish ##
  ############

  # If everything works, return combined results
  invisible(combined_results)
}


summarize_contrast_DE <- function(contrast_name, contrast_res_list) {
  tmp <- contrast_res_list[[contrast_name]][,c("sig.PVal", "sig.FDR")]

  data.frame(cbind(contrast = contrast_name,
                   type = c("down", "nonsig", "up"),
                   rbind(colSums(tmp == -1, na.rm = T),
                         colSums(tmp == 0, na.rm = T),
                         colSums(tmp == 1, na.rm = T))))

}


write_per_contrast_csvs <- function(annotation_df,
                                    data,
                                    results_statlist,
                                    output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }

  # make list combining
  # annotation, intensity data, and statistical results
  # reordering everything by rownames to match and subset
  # Note: technically, these could be in different orders across contrasts,
  # If they differed across contrasts.
  per_contrast_results <- lapply(X = results_statlist,
                                 FUN = function(x) cbind(annotation_df[rownames(x), ],
                                                         data[rownames(x), ],
                                                         x))
  # Make corresponding filenames
  filenames <- file.path(output_dir, paste0(names(per_contrast_results), ".csv"))
  # write
  tmp <- lapply(seq_along(per_contrast_results),
                function(x)
                  {
                  utils::write.csv(per_contrast_results[[x]],
                                   file = filenames[x],
                                   row.names = F)
                  }
                )

  # check
  write_success <- file.exists(filenames)
  names(write_success) <- names(per_contrast_results)

  write_success
}

