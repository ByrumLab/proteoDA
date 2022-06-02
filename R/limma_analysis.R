#' Fit the limma differential expression model
#'
#' Fits the limma differential expresison model to the expression data, following
#' the specified design and contrasts matrices. When paired = T, first estimates
#' the inter-duplicate correlation across paired samples with
#' \code{\link[limma:duplicateCorrelation]{limma::duplicateCorrelation}}. Then,
#' uses \code{\link[limma:lmFit]{limma::lmFit}} to fit the initial model,
#' \code{\link[limma:contrasts.fit]{limma::contrasts.fit}} to re-parameterize
#' the results in terms of your desired contrasts, and then recomputes moderated
#' statistic following limma's empirical Bayes  model with
#' \code{\link[limma:eBayes]{limma::eBayes}}.
#'
#' @param data The data on which to fit the limma model. In our pipeline,
#'   this should be a single normalized dataset from the normalized data.
#' @param design_obj A list object output from \code{\link{make_design}} which
#'   specifies the statistical design.
#' @param contrasts_obj A list object output from \code{\link{make_contrasts}}
#'   which specifies
#' @param paired Is the statistical design paired (TRUE) or not (FALSE)? Default
#'   is FALSE.
#'
#' @return A list with either six or seven slots: \enumerate{
#'   \item "eBayes_fit"- The results returned from \code{\link[limma:eBayes]{limma::eBayes}}.
#'   \item "contrasts_fit"- The results returned from
#'     \code{\link[limma:contrasts.fit]{limma::contrasts.fit}}.
#'   \item "lm_fit"- The results returned from
#'     \code{\link[limma:lmFit]{limma::lmFit}}.
#'   \item "data"- A data frame of the data that were input into the limma DE model.
#'   \item "design" - The design matrix used for model fitting.
#'   \item "contrasts"- The contrasts matrix used for model fitting.
#'   \item "corr_fit"- When paired == TRUE, output includes the results returned
#'     by \code{\link[limma:duplicateCorrelation]{limma::duplicateCorrelation}}.
#' }
#' @export
#'
#' @examples
#' # No examples yet
#'

fit_limma_model <- function(data,
                            design_obj,
                            contrasts_obj,
                            paired = FALSE) {
  # Original fxn had a check that rownames in the data equalled rownames in the
  # annotation. Took that out for now, and removed annotation as an argument,
  # since we don't actually need the annotation to fit the model. But may want
  # to put that back in. We definitely need to check that at some point (maybe in
  # later functions?), but it may be worth doing it now to just not waste time running the
  # model on something that is incorrect.

  # On the other hand, depending on the various outputs, may be able to
  # reorder and match things back up if they ever are actually out of sync.

  # Extract elements out of their objects
  design <- design_obj$design
  targets <- design_obj$targets
  contrasts <- contrasts_obj$contrasts

  #TODO: could possibly add some checks for the above, to make sure the objects are what we're expecting.

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
  # Levels in contrasts and levels in design match
  # TODO: need to figure out the best way to deal with this check in conjunction with
  # issue #14: these may not always match perfectly if, e.g., there's a
  # factor we're controlling for in the design matrix but not comparing in the
  # contrasts. So, I think everything in contrasts needs to be present in design,
  # but not necessarily vice versa (at least, with the way things are now)
  if (!all(rownames(contrasts) %in% colnames(design))) {
    problemContrasts <- rownames(contrasts)[rownames(contrasts) %notin% colnames(design)]
    cli::cli_abort(c("{cli::qty(length(problemContrasts))} {?A/Some} level{?s} in the contrast matrix {?is/are} not present in the design matrix",
                     "!" = "{cli::qty(length(problemContrasts))} Level{?s} without a match: {.val {problemContrasts}}"))
  }


  # Ensure data cols are in same order as rows in targets
  data <- data[, rownames(targets)]


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




#' Extract differential expression results from a model fit
#'
#' Extracts statistical results describing differential expression from a results
#' object created by \code{\link{fit_limma_model}}.
#'
#' @param limma_fit A results object from the \code{\link{fit_limma_model}}
#'   function.
#' @param pval.thresh The p-value threshold used to determine significance
#'   (significant when p < pval.thresh). Default is 0.055.
#' @param lfc.thresh The logFC threshold used to determine significance
#'   (significant when |logFC| > lfc.tresh). Default is 1. LogFC are base 2.
#' @param adj.method The method used for adjusting P-values. Default is "BH",
#'   for the Benjamini-Hochberg correction
#'
#' @return A list with five slots \enumerate{
#'   \item "stats_by_contrast"- A list with length equal to the number of contrasts.
#'     Each element in the list is a data frame giving the DE results for that contrast.
#'   \item "data"- A data frame of the data that were input into the limma DE model.
#'   \item "pval.thresh"- The p-value threshold used to determine significance.
#'   \item "lfc.thresh"- The logFC threshold used to determine significance.
#'   \item "adj.method"- The p-value adjustment method used.
#' }
#' @export
#'
#' @examples
#' # No examples yet
extract_limma_DE_results <- function(limma_fit, pval.thresh = 0.055, lfc.thresh = 1, adj.method = "BH") {

  # check args
  adj.method <- rlang::arg_match(
    arg = adj.method,
    values = c("none", "BH", "BY", "holm"),
    multiple = FALSE
  )


  # Just get it going initially, without too much testing.
  # Really, shouldn't need to do much testing if we're just passing in the whole limma fit object

  efit <- limma_fit$eBayes_fit # This should fail if we're not working with our data
  contrast_names <- colnames(efit$coefficients)

  # Get dfs where 0 = insig, -1 is sig downregulated, and 1 is sig upregulated
  outcomes_table_rawp <- as.data.frame(limma::decideTests(efit, adjust.method = "none", p.value = pval.thresh, lfc = lfc.thresh))
  outcomes_table_adjp <- as.data.frame(limma::decideTests(efit, adjust.method = adj.method, p.value = pval.thresh, lfc = lfc.thresh))

  perc_sig_rawp <- check_DE_perc(outcomes_table_rawp, pval.thresh = pval.thresh, lfc.thresh = lfc.thresh, adj.method = "none")
  perc_sig_adjp <- check_DE_perc(outcomes_table_adjp, pval.thresh = pval.thresh, lfc.thresh = lfc.thresh, adj.method = adj.method)

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
       pval.thresh = pval.thresh,
       lfc.thresh = lfc.thresh,
       adj.method = adj.method) # just passing data through
}


#' Check percentage of DE genes
#'
#' Internal utility function, used in \code{\link{extract_limma_DE_results}} to
#' check if assumptions are met.
#'
#' @param DE_outcomes_table DE results dataframe. Should be the output of
#' \code{\link[limma:decideTests]{limma::decideTests}}, coerced to a dataframe.
#' @param DE_warn_threshold Proporion of DE genes at which we warn user.
#' @param pval.thresh P-value threshold used.
#' @param lfc.thresh logFC threshold used.
#' @param adj.method P-value adjustment method used.
#'
#' @return A vector of numeric values giving the % of significant DE proteins within
#'   each contrast.
#'
#' @examples
#' # No examples yet
check_DE_perc <- function(DE_outcomes_table, DE_warn_threshold = 0.2, pval.thresh, lfc.thresh, adj.method) {
  perc_sig <- colSums(DE_outcomes_table != 0, na.rm = T)/colSums(!is.na(DE_outcomes_table))

  if (any(perc_sig > DE_warn_threshold)) {
    above_thresh <- names(perc_sig)[perc_sig > DE_warn_threshold]
    thresh_perc <- DE_warn_threshold*100
    cli::cli_inform(c("!" = "Warning: more than {.perc {thresh_perc}}% of the data is DE in {cli::qty(length(above_thresh))} {?a/some} contrast{?s}",
                      "!" = "Criteria for DE: |logFC| > {.val {lfc.thresh}}, p-value < {.val {pval.thresh}}, p.value adjustment = {.val {adj.method}}",
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



