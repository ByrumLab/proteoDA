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

  # With a random factor
  if (!is.null(DIAlist$design$random_factor)) {
    block <- DIAlist$metadata[, DIAlist$design$random_factor, drop = T]
    corfit <- limma::duplicateCorrelation(object = DIAlist$data,
                                          design = DIAlist$design$design_matrix,
                                          block = block)

    corfit_display <- round(corfit$consensus.correlation, 3)
    cli::cli_inform("Estimated intra-block correlation = {.val {corfit_display}}")

    if (corfit$consensus.correlation < 0.1) {
      cli::cli_inform(cli::col_yellow(c("Estimated intra-block correlation is low.",
                                        "Consider using a model withno random effect.")))
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
