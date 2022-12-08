#' Fit the limma differential expression model
#'
#' Fits the limma differential expresison model to the expression data, following
#' the specified design and (optional) contrast matrices. When a random factor is included,
#' uses \code{\link[limma:duplicateCorrelation]{limma::duplicateCorrelation}} to estimate
#' the intra-block correlation within groups. Uses \code{\link[limma:lmFit]{limma::lmFit}}
#' to fit the initial model, optionally re-parameterizes the results in terms of
#' contrasts with \code{\link[limma:contrasts.fit]{limma::contrasts.fit}},
#' and then recomputes moderated statistics following limma's empirical Bayes
#' model with \code{\link[limma:eBayes]{limma::eBayes}}.
#'
#' @param DIAlist A DIAlist, which must contain a statistical design.
#'
#' @return A DIAlist object, with the model fit.
#'
#' @export
#'
#' @examples
#' # No examples yet
#'

fit_limma_model <- function(DIAlist) {

  # Check input arguments generally
  validate_DIAlist(DIAlist)

  # Make sure there's a design matrix present already,
  # tell user to set it first if not
  if (is.null(DIAlist$design)) {
    cli::cli_abort(c("Input DIAlist does not have a statistical design",
                     "i" = "Run {.code DIAlist <- add_design(DIAlist, ~ formula)} before fitting model"))
  }

  # Warn user if data aren't normalized, unlikely
  # to want to do stats on non-normalized data
  if (is.null(DIAlist$tags$normalized)) {
    cli::cli_inform("Data in DIAlist are not normalized. You may wish to normalize before fitting model.")
  }
  if (!is.null(DIAlist$tags$normalized)) {
    if (!DIAlist$tags$normalized) {
      cli::cli_inform("Data in DIAlist are not normalized. You may wish to normalize before fitting model.")
    }
  }

  # With a random factor
  if (!is.null(DIAlist$design$random_factor)) {

    if (!requireNamespace("statmod", quietly = TRUE)) {
      cli::cli_abort(c("Package \"statmod\" must be installed to model a random effect"))
    }

    block <- DIAlist$metadata[, DIAlist$design$random_factor, drop = T]
    corfit <- limma::duplicateCorrelation(object = DIAlist$data,
                                          design = DIAlist$design$design_matrix,
                                          block = block)

    corfit_display <- round(corfit$consensus.correlation, 3)
    cli::cli_inform("Estimated intra-block correlation = {.val {corfit_display}}")

    if (corfit$consensus.correlation < 0.1) {
      cli::cli_inform(cli::col_yellow(c("Estimated intra-block correlation is low.",
                                        "Consider using a model with no random effect.")))
    }
    intra_block_cor <- corfit$consensus.correlation
  } else { # No random factor
    intra_block_cor <- NULL
    block <- NULL
  }

  # Fit the initial model
  fit <- limma::lmFit(object = DIAlist$data,
                      design = DIAlist$design$design_matrix,
                      block = block,
                      correlation = intra_block_cor)

  # If contrasts are specified, re-fit
  if (!is.null(DIAlist$design$contrast_matrix)) {
    fit <- limma::contrasts.fit(fit = fit, contrasts = DIAlist$design$contrast_matrix)
  }

  # Annoying statmod issue still...
  if (!requireNamespace("statmod", quietly = TRUE)) {
    cli::cli_abort(c("Package \"statmod\" must be installed to perform empirical Bayes moderation of test statistics"))
  }

  # Fit empirical bayes model
  efit <- limma::eBayes(fit = fit, robust = TRUE)

  # If there are already results, warn about overwriting
  if (!is.null(DIAlist$eBayes_fit) | !is.null(DIAlist$results)) {
    cli::cli_inform("DIAlist already contains statistical results. Overwriting.")
    # Get rid of any old stuff
    DIAlist$eBayes_fit <- NULL
    DIAlist$results <- NULL
  }

  # Add results
  DIAlist$eBayes_fit <- efit

  # Validate and return
  validate_DIAlist(DIAlist)
}

#' Extract differential expression results from a model fit
#'
#' Extracts statistical results describing differential expression from the.
#'
#'
#' @param DIAlist A DIAlist, which must contain a statistical design.
#' @param pval_thresh The p-value threshold used to determine significance
#'   (significant when p < pval_thresh). Default is 0.055.
#' @param lfc_thresh The logFC threshold used to determine significance
#'   (significant when |logFC| > lfc.tresh). Default is 1. LogFC are base 2.
#' @param adj_method The method used for adjusting P-values. Default is "BH",
#'   for the Benjamini-Hochberg correction
#'
#' @param DIAlist A DIAlist, which must contain a statistical design.
#'
#' @return A DIAlist object, with differential expression results.
#' @export
#'
#' @examples
#' # No examples yet
extract_DE_results <- function(DIAlist, pval_thresh = 0.055, lfc_thresh = 1, adj_method = "BH") {

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
          !(pval_thresh >= 0))) {
    cli::cli_abort(c("{.arg lfc_thresh} must be a numeric value greater >= 0."))
  }

  # validate the input dataset
  validate_DIAlist(DIAlist)

  # Make sure there's an eBayes fit already
  # tell user to set it first if not
  if (is.null(DIAlist$eBayes_fit)) {
    cli::cli_abort(c("Input DIAlist does not have a model fit",
                     "i" = "Run {.code DIAlist <- fit_limma_model(DIAlist)} to fit a model."))
  }


  # check if there are already results, warn about overwriting if so
  if (!is.null(DIAlist$results)) {
    cli::cli_inform("DIAlist already contains DE results. Overwriting.")
    # Get rid of any old stuff
    DIAlist$results <- NULL
  }

  # Grab fit object and extract results from it
  efit <- DIAlist$eBayes_fit
  contrast_names <- colnames(efit$coefficients)

  # Get dfs where 0 = insig, -1 is sig downregulated, and 1 is sig upregulated
  outcomes_table_rawp <- as.data.frame(limma::decideTests(efit, adjust.method = "none", p.value = pval_thresh, lfc = lfc_thresh))
  outcomes_table_adjp <- as.data.frame(limma::decideTests(efit, adjust.method = adj_method, p.value = pval_thresh, lfc = lfc_thresh))

  perc_sig_rawp <- check_DE_perc(outcomes_table_rawp, pval_thresh = pval_thresh, lfc_thresh = lfc_thresh, adj_method = "none")
  perc_sig_adjp <- check_DE_perc(outcomes_table_adjp, pval_thresh = pval_thresh, lfc_thresh = lfc_thresh, adj_method = adj_method)

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
  })
  names(results_per_contrast) <- contrast_names

  # Add results
  DIAlist$results <- results_per_contrast
  DIAlist$tags$DE_criteria$pval_thresh <- pval_thresh
  DIAlist$tags$DE_criteria$lfc_thresh <- lfc_thresh
  DIAlist$tags$DE_criteria$adj_method <- adj_method

  # Validate and return
  validate_DIAlist(DIAlist)
}


#' Check percentage of DE genes
#'
#' Internal utility function, used in \code{\link{extract_DE_results}} to
#' check if assumptions are met.
#'
#' @param DE_outcomes_table DE results dataframe. Should be the output of
#' \code{\link[limma:decideTests]{limma::decideTests}}, coerced to a dataframe.
#' @param DE_warn_threshold Proporion of DE genes at which we warn user.
#' @param pval_thresh P-value threshold used.
#' @param lfc_thresh logFC threshold used.
#' @param adj_method P-value adjustment method used.
#'
#' @return A vector of numeric values giving the % of significant DE proteins within
#'   each contrast.
#'
#' @examples
#' # No examples yet
check_DE_perc <- function(DE_outcomes_table, DE_warn_threshold = 0.2, pval_thresh, lfc_thresh, adj_method) {
  perc_sig <- colSums(DE_outcomes_table != 0, na.rm = T)/colSums(!is.na(DE_outcomes_table))

  # Don't check the intercept column, if it exists
  perc_sig <- perc_sig[names(perc_sig) %notin% c("Intercept")]

  if (any(perc_sig > DE_warn_threshold)) {
    above_thresh <- names(perc_sig)[perc_sig > DE_warn_threshold]
    thresh_perc <- DE_warn_threshold*100
    cli::cli_inform(c("!" = "Warning: more than {.perc {thresh_perc}}% of the data is DE in {cli::qty(length(above_thresh))} {?a/some} term{?s}",
                      "!" = "Criteria for DE: |logFC| > {.val {lfc_thresh}}, p-value < {.val {pval_thresh}}, p.value adjustment = {.val {adj_method}}",
                      "!" = "{cli::qty(length(above_thresh))} Problematic term{?s}: {.val {above_thresh}}",
                      "!" = "Assumption that most genes/proteins/phospho are not DE may be violated"))
  }

  perc_sig
}
