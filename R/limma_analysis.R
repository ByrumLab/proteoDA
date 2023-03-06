#' Fit the limma differential abundance model
#'
#' Fits the limma differential abundance model to the intensity data, following
#' the specified design and (optional) contrast matrices. When a random factor is included,
#' uses \code{\link[limma:duplicateCorrelation]{limma::duplicateCorrelation}} to estimate
#' the intra-block correlation within groups. Uses \code{\link[limma:lmFit]{limma::lmFit}}
#' to fit the initial model, optionally re-parameterizes the results in terms of
#' contrasts with \code{\link[limma:contrasts.fit]{limma::contrasts.fit}},
#' and then recomputes moderated statistics following limma's empirical Bayes
#' model with \code{\link[limma:eBayes]{limma::eBayes}}.
#'
#' @param DAList A DAList, which must contain a statistical design.
#'
#' @return A DAList object, with the model fit added in the eBayes_fit slot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' model_fit <- fit_limma_model(DAList)
#' }
#'

fit_limma_model <- function(DAList) {

  # Make sure there's a design matrix present already,
  # tell user to set it first if not
  # put these checks above DAList validation, to try to get better errors for endusers
  if (is.null(DAList$design)) {
    cli::cli_abort(c("Input DAList does not have a statistical design",
                     "i" = "Run {.code DAList <- add_design(DAList, ~ formula)} before fitting model"))
  }

  # Warn user if data aren't normalized, unlikely
  # to want to do stats on non-normalized data
  if (is.null(DAList$tags$normalized)) {
    cli::cli_inform("Data in DAList are not normalized. You may wish to normalize before fitting model.")
  }
  if (!is.null(DAList$tags$normalized)) {
    if (!DAList$tags$normalized) {
      cli::cli_inform("Data in DAList are not normalized. You may wish to normalize before fitting model.")
    }
  }


  # Check input arguments generally
  validate_DAList(DAList)

  # With a random factor
  if (!is.null(DAList$design$random_factor)) {

    if (!requireNamespace("statmod", quietly = TRUE)) {
      cli::cli_abort(c("Package \"statmod\" must be installed to model a random effect"))
    }

    block <- DAList$metadata[, DAList$design$random_factor, drop = T]
    corfit <- limma::duplicateCorrelation(object = DAList$data,
                                          design = DAList$design$design_matrix,
                                          block = block)

    corfit_display <- round(corfit$consensus.correlation, 3)
    cli::cli_inform("Estimated intra-block correlation = {.val {corfit_display}}")

    if (corfit$consensus.correlation < 0) {
      cli::cli_abort(c("Estimated intra-block correlation is negative.",
                       "Rerun model without the random effect."))
    } else if (corfit$consensus.correlation < 0.05) {
      cli::cli_inform(cli::col_yellow(c("Estimated intra-block correlation is low.",
                                        "Consider using a model with no random effect.")))
    } else {
      intra_block_cor <- corfit$consensus.correlation
    }
  } else { # No random factor
    intra_block_cor <- NULL
    block <- NULL
  }

  # Fit the initial model
  fit <- limma::lmFit(object = DAList$data,
                      design = DAList$design$design_matrix,
                      block = block,
                      correlation = intra_block_cor)

  # If contrasts are specified, re-fit
  if (!is.null(DAList$design$contrast_matrix)) {
    fit <- limma::contrasts.fit(fit = fit, contrasts = DAList$design$contrast_matrix)
  }

  # Annoying statmod issue still...
  if (!requireNamespace("statmod", quietly = TRUE)) {
    cli::cli_abort(c("Package \"statmod\" must be installed to perform empirical Bayes moderation of test statistics"))
  }

  # Fit empirical bayes model
  efit <- limma::eBayes(fit = fit, robust = TRUE)

  # If there are already results, warn about overwriting
  if (!is.null(DAList$eBayes_fit) | !is.null(DAList$results)) {
    cli::cli_inform("DAList already contains statistical results. Overwriting.")
    # Get rid of any old stuff
    DAList["eBayes_fit"] <- list(NULL)
    DAList["results"] <- list(NULL)
  }

  # Add results
  DAList$eBayes_fit <- efit

  # Validate and return
  validate_DAList(DAList)
}

#' Extract differential abundance results from a model fit
#'
#' Extracts tables of statistical results from the model fit created by
#' \code{\link{fit_limma_model}} and add them to the results slot of the DAList.
#' The results are a list of tables, one for each contrast/term specific in the
#' statistical design. In models with intercepts these are ignored by default,
#' but can be output by setting extract_intercept to TRUE. See
#' \code{\link[limma:decideTests]{limma::decideTests}} and
#' \code{\link[limma:topTable]{limma::topTable}} for information on the statistical results.
#'
#' @param DAList A DAList, which must contain a model fit.
#' @param pval_thresh The p-value threshold used to determine significance
#'   (significant when p < pval_thresh). Default is 0.05.
#' @param lfc_thresh The logFC threshold used to determine significance
#'   (significant when |logFC| > lfc.tresh). Default is 1. LogFC are base 2.
#' @param adj_method The method used for adjusting P-values. Possible values are
#'   "none", "BH", "BY", and "holm". Default is "BH", for the Benjamini-Hochberg
#'   correction. See \code{\link[stats:p.adjust]{stats::p.adjust}} for details.
#' @param extract_intercept For models with an intercept term, should results for
#'   the intercept be extracted? Default if FALSE.
#'
#' @return A DAList object, with differential abundance results added
#'   to the results slot.
#' @export
#'
#' @examples
#' \dontrun{
#' # Using default thresholds and p-value adjustment
#' results <- extract_DA_results(DAList)
#'
#' # Relax significance and log fold-change thresholds
#' results <- extract_DA_results(DAList,
#'                               pval_thresh = 0.1,
#'                               lfc_thresh = 0)
#'
#' # Include intercept term in results
#' results <- extract_DA_results(DAList,
#'                               extract_intercept = T)
#' }
extract_DA_results <- function(DAList, pval_thresh = 0.05, lfc_thresh = 1, adj_method = "BH", extract_intercept = F) {

  # check args
  adj_method <- rlang::arg_match(
    arg = adj_method,
    values = c("none", "BH", "BY", "holm"),
    multiple = FALSE
  )

  if (any(!is.numeric(pval_thresh),
          !(pval_thresh > 0),
          !(pval_thresh <= 1))) {
    cli::cli_abort(c("{.arg pval_thresh} must be a numeric value greater > 0 and <= 1."))
  }

  if (any(!is.numeric(lfc_thresh),
          lfc_thresh <= 0)) {
    cli::cli_abort(c("{.arg lfc_thresh} must be a numeric value greater >= 0."))
  }

  # validate the input dataset
  validate_DAList(DAList)

  # Make sure there's an eBayes fit already
  # tell user to set it first if not
  if (is.null(DAList$eBayes_fit)) {
    cli::cli_abort(c("Input DAList does not have a model fit",
                     "i" = "Run {.code DAList <- fit_limma_model(DAList)} to fit a model."))
  }


  # check if there are already results, warn about overwriting if so
  if (!is.null(DAList$results)) {
    cli::cli_inform("DAList already contains DA results. Overwriting.")
    # Get rid of any old stuff
    DAList["results"] <- list(NULL)
  }

  # Grab fit object and extract results from it
  efit <- DAList$eBayes_fit
  contrast_names <- colnames(efit$coefficients)

  # Get dfs where 0 = insig, -1 is sig downregulated, and 1 is sig upregulated
  outcomes_table_rawp <- as.data.frame(limma::decideTests(efit, adjust.method = "none", p.value = pval_thresh, lfc = lfc_thresh))
  outcomes_table_adjp <- as.data.frame(limma::decideTests(efit, adjust.method = adj_method, p.value = pval_thresh, lfc = lfc_thresh))

  # Remove intercept term
  # if we're not extracting it
  if (!extract_intercept) {
    # Get the non-intercept indices, based on contrast names
    # use string searching instead of index, so it shouldn't for no-intercept models
    non_intercept_terms <- which(!stringr::str_detect(stringr::str_to_lower(contrast_names), "intercept"))

    outcomes_table_rawp <- outcomes_table_rawp[,non_intercept_terms, drop = F]
    outcomes_table_adjp <- outcomes_table_adjp[,non_intercept_terms, drop = F]
    contrast_names <- contrast_names[non_intercept_terms]
  }


  perc_sig_rawp <- check_DA_perc(outcomes_table_rawp, pval_thresh = pval_thresh, lfc_thresh = lfc_thresh, adj_method = "none")
  perc_sig_adjp <- check_DA_perc(outcomes_table_adjp, pval_thresh = pval_thresh, lfc_thresh = lfc_thresh, adj_method = adj_method)

  # make a list of statistical results for each contrast/term/comparison
  results_per_contrast <- list()
  limmaStatColums <- c("logFC", "CI.L", "CI.R", "AveExpr", "t", "B", "P.Value", "adj.P.Val")
  results_per_contrast <- base::lapply(contrast_names, function(x) {
    one_contrast_table <- limma::topTable(efit,
                             coef = x, number = Inf, adjust.method = adj_method,
                             sort.by = "none", p.value = 1, lfc = 0, confint = TRUE
    )
    outcomes <- cbind(outcomes_table_rawp[, x, drop = F], outcomes_table_adjp[, x, drop = F])
    colnames(outcomes) <- c("sig.PVal", "sig.FDR")
    one_contrast_results <- cbind(one_contrast_table[, limmaStatColums], outcomes[rownames(outcomes), ])
    # limma names for expression/array studies, rename to average_intensity for protein data
    colnames(one_contrast_results)[colnames(one_contrast_results) == "AveExpr"] <- "average_intensity"
    one_contrast_results
  })
  names(results_per_contrast) <- contrast_names

  # Add results
  DAList$results <- results_per_contrast
  DAList$tags$DA_criteria$pval_thresh <- pval_thresh
  DAList$tags$DA_criteria$lfc_thresh <- lfc_thresh
  DAList$tags$DA_criteria$adj_method <- adj_method
  DAList$tags$extract_intercept <- extract_intercept

  # Validate and return
  validate_DAList(DAList)
}


#' Check percentage of DA genes
#'
#' Internal utility function, used in \code{\link{extract_DA_results}} to
#' check if assumptions are met.
#'
#' @param DA_outcomes_table DA results data frame. Should be the output of
#' \code{\link[limma:decideTests]{limma::decideTests}}, coerced to a data frame.
#' @param DA_warn_threshold Proportion of DA genes at which we warn user.
#' @param pval_thresh P-value threshold used.
#' @param lfc_thresh logFC threshold used.
#' @param adj_method P-value adjustment method used.
#'
#' @return A vector of numeric values giving the % of significant DA proteins within
#'   each contrast.
#'
check_DA_perc <- function(DA_outcomes_table, DA_warn_threshold = 0.2, pval_thresh, lfc_thresh, adj_method) {
  perc_sig <- colSums(DA_outcomes_table != 0, na.rm = T)/colSums(!is.na(DA_outcomes_table))

  # Don't check the intercept column, if it exists
  perc_sig <- perc_sig[names(perc_sig) %notin% c("Intercept", "intercept")]

  if (any(perc_sig > DA_warn_threshold)) {
    above_thresh <- names(perc_sig)[perc_sig > DA_warn_threshold]
    thresh_perc <- DA_warn_threshold*100
    cli::cli_inform(c("!" = "Warning: more than {.perc {thresh_perc}}% of the data is differentially adundant in {cli::qty(length(above_thresh))} {?a/some} term{?s}",
                      "!" = "Criteria for DA: |logFC| > {.val {lfc_thresh}}, p-value < {.val {pval_thresh}}, p.value adjustment = {.val {adj_method}}",
                      "!" = "{cli::qty(length(above_thresh))} Problematic term{?s}: {.val {above_thresh}}",
                      "!" = "Assumption that most proteins are not DA may be violated"))
  }

  perc_sig
}
