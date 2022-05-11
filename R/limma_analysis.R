## ---------------------------
##  LIMMA EBAYES (NORMAL)
## ---------------------------
## paired==FALSE is used for group comparisons,
## group comparisons correcting for e.g. batch or gender effect
## group comparisons for paired samples etc.



## ---------------------------------
##  LIMMA EBAYES (MIXED EFFECTS)
## ---------------------------------
## paired==TRUE is used for group comparisons,
## when comparing within (paired samples) and across (groups)
## subjects

fit_limma_model <- function(data,
                            targets,
                            design,
                            contrasts,
                            paired = FALSE) {
  # Original fxn had a check that rownames in the data equaled rownames in the
  # annotation. Took that out for now, and removed annotation as an argument,
  # since we don't actually need the annotation to fit the model. But may want
  # to put that back in. We definitely need to check that at some point (maybe in
  # later functions?), but it may be worth doing it now to just not waste time running the
  # model on something that is incorrect.

  # On the other hand, depending on the various outputs, may be able to
  # reorder and match things back up if they ever are actually out of sync.

  # TODO:
  # NEED TO ADD SOME SORT OF CHECK FOR THE CONTRASTS


  # Do some checks of the input data
  # All samples in data have targets
  if (!all(colnames(data) %in% rownames(targets))) {
    problemCols <- colnames(data)[colnames(data) %notin% rownames(targets)]
    cli::cli_abort(c("Not all column names in {.arg data} have a matching rowname in {.arg targets}",
                     "!" = "{cli::qty(length(problemCols))} Column{?s} without a match: {.val {problemCols}}"))
  }
  # All targets have corresponding data
  if (!all(rownames(targets) %in% colnames(data))) {
    problemTargets <- rownames(targets)[rownames(targets) %notin% colnames(data)]
    cli::cli_abort(c("Not all target rownames in {.arg targets} have a matching colname in {.arg data}",
                     "!" = "{cli::qty(length(problemTargets))} Row{?s} without a match: {.val {problemTargets}}"))
  }

  # Ensure data cols are in same order as rows in targets
  data <- data[, rownames(targets)]
  groups <- targets$group

  # Double-check that the order is all the same
  # TODO: Shouldn't really ever trigger this, right? maybe remove?
  stopifnot(identical(colnames(data), rownames(targets)))


  cli::cli_rule()

  # On to model fitting
  # First step of fitting differs between paired and unpaired
  if (paired) {

    cli::cli_inform("Performing paired analysis with mixed effects model")
    # check for paired col
    if ("paired" %notin% colnames(targets)) {
      cli::cli_abort(c("{.arg targets} does not contain required column names {.val paired}",
                       "i" = "Add {.val paired} column, or using {.fun make_design}",
                       "i" = "Ensure that design formula does not include {.val paired}",
                       "i" = "e.g. ~0+group OR ~0+group+batch NOT ~0+group+paired"))
    }

    # Coerce paired col to factor, if it isn't
    if (!is.factor(targets$paired)) {
      targets$paired <- make_factor(targets$paired, prefix = "")
    }

    corfit <- limma::duplicateCorrelation(object = data, design = design, block = targets$paired)
    corfit_display <- round(corfit$consensus.correlation, 3)
    cli::cli_inform("Estimated inter-duplicate correlation = {.val {corfit_display}}")

    if (corfit$consensus.correlation < 0.1) {
      cli::cli_inform(cli::col_yellow("Estimated inter-duplicate correlation is low,
                                      which may indicate little or no paired influence"))
    }

    fit <- limma::lmFit(object = data, design = design,
                        block = targets$paired,
                        correlation = corfit$consensus.correlation)
  } else {
    cli::cli_inform("Performing standard un-paired model")
    fit <- limma::lmFit(object = data, design = design)
    }

  # contrasts fit and eBayes are the same across paired and not paired
  con.fit <- limma::contrasts.fit(fit = fit, contrasts = contrasts)
  efit <- limma::eBayes(fit = con.fit, robust = TRUE)
  cli::cli_inform("limma DE analysis with {.arg paired} == {paired} complete")


  # Set up return
  model <- list(eBayes_fit = efit,
                contrasts_fit = con.fit,
                lm_fit = fit,
                data = data,
                design = design,
                contrasts = contrasts)
  if (paired) model[["corr_fit"]] <- corfit


  cli::cli_rule()
  cli::cli_inform(c("v" = "Success!"))

  model
}





extract_limma_DE_results <- function(limma_fit, min.pval = 0.055, min.lfc = 1, adj.method = "BH") {

  # TODO:
  # ADD a check on adjustment method?

  # Just get it going initially, without too much testing.
  # Really, shouldn't need to do much testing if we're just passing in the whole limma fit object

  efit <- limma_fit$eBayes_fit # This should fail if we're not working with our data
  contrast_names <- colnames(efit$coefficients)

  # Get dfs where 0 = insig, -1 is sig downregulated, and 1 is sig upregulated
  outcomes_table_rawp <- limma::decideTests(efit, adjust.method = "none", p.value = min.pval, lfc = min.lfc) %>%
    as.data.frame(.)
  outcomes_table_adjp <- limma::decideTests(efit, adjust.method = adj.method, p.value = min.pval, lfc = min.lfc) %>%
    as.data.frame(.)

  perc_sig_rawp <- check_DE_perc(outcomes_table_rawp, min.pval = min.pval, min.lfc = min.lfc, adj.method = "none")
  perc_sig_adjp <- check_DE_perc(outcomes_table_adjp, min.pval = min.pval, min.lfc = min.lfc, adj.method = adj.method)

  # Make statlist, for use in later functions
  # A list where length is number of contrasts
  # And each element is a dataframe of results for that contrast
  results_per_contrast <- list()
  limmaStatColums <- c("logFC", "CI.L", "CI.R", "AveExpr", "t", "B", "P.Value", "adj.P.Val")
  results_per_contrast <- base::lapply(contrast_names, function(x) {
    one_contrast_table <- limma::topTable(efit,
                             coef = x, number = Inf, adjust.method = adj.method,
                             sort.by = "none", p.value = 1, lfc = 0, confint = TRUE
    )
    outcomes <- cbind(outcomes_table_rawp[, x, drop = F], outcomes_table_adjp[, x, drop = F])
    colnames(outcomes) <- c("sig.PVal", "sig.FDR")
    one_contrast_results <- cbind(one_contrast_table[, limmaStatColums], outcomes[rownames(outcomes), ])
  })
  names(results_per_contrast) <- contrast_names

  list(stats_by_contrast = results_per_contrast,
       data = limma_fit$data,
       min.pval = min.pval,
       min.lfc = min.lfc,
       adj.method = adj.method) # just passing data through
}



# utility function used in extract_limma_DE_results
check_DE_perc <- function(DE_outcomes_table, threshold = 0.1, min.pval, min.lfc, adj.method) {
  perc_sig <- colSums(DE_outcomes_table != 0, na.rm = T)/colSums(!is.na(DE_outcomes_table))

  if (any(perc_sig > threshold)) {
    above_thresh <- names(perc_sig)[perc_sig > threshold]
    thresh_perc <- threshold*100
    cli::cli_inform(c("!" = "Warning: more than {.perc {thresh_perc}}% of the data is DE in {cli::qty(length(above_thresh))} {?A/some} contrast{?s}",
                      "!" = "Criteria for DE: min.lfc > {.val {min.lfc}}, min.pval < {.val {min.pval}}, p.value adjustment = {.val {adj.method}}",
                      "!" = "{cli::qty(length(above_thresh))} Problematic contrast{?s}: {.val {above_thresh}}",
                      "!" = "Assumption that most genes/proteins/phospho are not DE may be violated"))
  }

  perc_sig
}


## TODO:
## NEED TO REVISIT LOGGING FOR ALL THIS.

# dummy_function <- function(data,
#                            annot,
#                            targets,
#                            design,
#                            contrasts,
#                            min.pval = 0.055,
#                            min.lfc = 1,
#                            adj.method = "BH",
#                            paired = FALSE,
#                            pipe = "DIA",
#                            enrich = "protein",
#                            dir = NULL,
#                            save = TRUE,
#                            ilab = "PI_DATE") {
#   enrich <- match.arg(arg = enrich, choices = c("protein", "phospho"), several.ok = FALSE)
#
#
#
#
#   pipe <- match.arg(arg = pipe, choices = c("DIA", "TMT", "phosphoTMT", "LF"), several.ok = FALSE)
#
#   adj.method <- match.arg(
#     arg = adj.method,
#     choices = c("none", "BH", "BY", "holm"), several.ok = FALSE
#   )
#
#   ## MATCH TARGETS, DATA, ANNOTATION
#   # from original limma analysis
#   if (all(rownames(data) %in% rownames(annot))) {
#     annot <- annot[rownames(data), ]
#   } else {
#     stop("Error! Row names of norm. data and row names of annotation
#                     do not match.")
#   }
#
#   stopifnot(identical(rownames(data), rownames(annot)))
#
#   ## DE stat results are extracted in 3 formats. statList is list object, where each
#   ## item in the list is a data.frame of the stat results for a particular contrast.
#   ## This list object is used to create/save DE plots and individual stat result files.
#   ## comboStats = combined stat results in wide format. This data.frame is used to
#   ## create results file for Big Query upload. This data.frame is used to create results file
#   ## returned to the investigator.
#   res <- extract_limma_results(
#     efit = efit, annot = annot, data = data, min.pval = min.pval,
#     min.lfc = min.lfc, adj.method = adj.method, dir = dir, save = save,
#     enrich = enrich, ilab = ilab
#   )
#
#   ## SAVE DE PLOTS
#   if (save) {
#     ## SUMMARY OF DIFF EXPRESSION
#     sumfile <- paste0("./", dir, "/summary.txt")
#     sink(file = sumfile)
#     cat(paste0("\n##", paste(rep("-", 40), collapse = "")))
#     cat("\n##  Summary of Differential Expression")
#     cat(paste0("\n##", paste(rep("-", 40), collapse = "")))
#     cat("\n\n")
#     cat(paste0("Significance Criteria: (|logFC| >= ", min.lfc, " & p-value <= ", min.pval, ")"))
#     cat("\n\n")
#     print(res$sum.dtp)
#     cat(paste0("Significance Criteria: (|logFC| >= ", min.lfc, " & adj. p-value <= ", min.pval, ")"))
#     cat("\n\n")
#     print(res$sum.dt)
#     cat("\n\n")
#     sink()
#   } ## SAVE
#
#
#   ## save summary DE results to log file
#   if (!dir.exists("logs")) {
#     dir.create("logs", recursive = TRUE)
#   }
#   sink(file = "./logs/processing.log", append = TRUE)
#   title <- "LIMMA DE SUMMARY (P-VALUE)"
#   cat(paste0("\n##", paste(rep("-", 40), collapse = "")))
#   cat(paste0("\n##  ", title, "\n"))
#   cat(paste0("##", paste(rep("-", 40), collapse = "")))
#   cat("\n\n")
#   print(res$sum.dtp)
#   cat("\n\n")
#   sink()
#
#   sink(file = "./logs/processing.log", append = TRUE)
#   title <- "LIMMA DE SUMMARY (ADJ. P-VALUE)"
#   cat(paste0("\n##", paste(rep("-", 40), collapse = "")))
#   cat(paste0("\n##  ", title, "\n"))
#   cat(paste0("##", paste(rep("-", 40), collapse = "")))
#   cat("\n\n")
#   print(res$sum.dt)
#   cat("\n\n")
#   sink()
#
#   param <- stats <- list()
#   param[["min.lfc"]] <- min.lfc
#   param[["min.pval"]] <- min.pval
#   param[["adj.method"]] <- adj.method
#   param[["paired"]] <- paired
#   param[["robust"]] <- TRUE
#   param[["pipe"]] <- pipe
#   param[["enrich"]] <- enrich
#   param[["ilab"]] <- ilab
#
#   ## save DE parameters to log file
#   logs <- make_log(param = param, stats = stats, title = "LIMMA DE ANALYSIS", save = TRUE)
#
#
#   data2 <- list(
#     statList = res$statList, comboStats = res$comboStats, comboStats.BQ = res$comboStats.BQ,
#     de = res$de, dep = res$dep, dt = res$dt, dtp = res$dtp, sum.dt = res$sum.dt, sum.dtp = res$sum.dtp,
#     contrastNames = res$contrastNames, model = model, targets = res$targets, groups = groups,
#     data = res$data, annot = res$annot, param = logs$param, stats = logs$stats
#   )
#
#   return(data2)
# }



