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
    cli::cli_abort(c(
      "Input DAList does not have a statistical design",
      "i" = "Run {.code DAList <- add_design(DAList, ~ formula)} before fitting model"
    ))
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
  
  # ------------------------------------------------------------------
  # Random factor handling
  # ------------------------------------------------------------------
  intra_block_cor <- NULL
  block           <- NULL
  rf              <- DAList$design$random_factor
  
  # With a random factor
  if (!is.null(rf)) {
    
    if (!requireNamespace("statmod", quietly = TRUE)) {
      cli::cli_abort(c(
        "Package \"statmod\" must be installed to model a random effect"
      )) # nocov
    }
    
    block <- DAList$metadata[, rf, drop = TRUE]
    
    corfit <- limma::duplicateCorrelation(
      object = DAList$data,
      design = DAList$design$design_matrix,
      block  = block
    )
    
    corfit_display <- round(corfit$consensus.correlation, 3)
    cli::cli_inform("Estimated intra-block correlation = {.val {corfit_display}}")
    
    if (corfit$consensus.correlation < 0) {
      cli::cli_abort(c(
        "Estimated intra-block correlation is negative.",
        "Rerun model without the random effect."
      ))
    } else if (corfit$consensus.correlation < 0.05) {
      cli::cli_inform(cli::col_yellow(c(
        "Estimated intra-block correlation is low.",         # nocov
        "Consider using a model with no random effect."      # nocov
      )))
    } else {
      intra_block_cor <- corfit$consensus.correlation
    }
  }
  
  # ------------------------------------------------------------------
  # Fit the initial model
  # ------------------------------------------------------------------
  fit <- limma::lmFit(
    object      = DAList$data,
    design      = DAList$design$design_matrix,
    block       = block,
    correlation = intra_block_cor
  )
  
  # If contrasts are specified, re-fit
  if (!is.null(DAList$design$contrast_matrix)) {
    fit <- limma::contrasts.fit(
      fit       = fit,
      contrasts = DAList$design$contrast_matrix
    )
  }
  
  if (!requireNamespace("statmod", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package \"statmod\" must be installed to perform empirical Bayes moderation of test statistics"
    )) # nocov
  }
  
  # Fit empirical Bayes model
  efit <- limma::eBayes(fit = fit, robust = TRUE)
  
  # ------------------------------------------------------------------
  # Tag the fit with random-factor information (for validate_DAList)
  # ------------------------------------------------------------------
  if (!is.null(rf)) {
    attr(efit, "random_factor_name")  <- rf
    attr(efit, "random_factor_block") <- block
  }
  
  # If there are already results, warn about overwriting
  if (!is.null(DAList$eBayes_fit) | !is.null(DAList$results)) {
    cli::cli_inform("DAList already contains statistical results. Overwriting.")
    # Get rid of any old stuff
    DAList["eBayes_fit"] <- list(NULL)
    DAList["results"]    <- list(NULL)
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
          lfc_thresh < 0)) {
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


# Run filtered limma analysis per contrast
# What changed (and why)
# No _vs_ splitting or sample subsetting. We no longer assume the label encodes two groups. This lets you use labels like Diff_in_effect_CHLA_SKNF1 freely.
# Contrast expressions drive everything. Each loop pulls the full "Label = ..." expression and passes it to add_contrasts(). As long as the RHS uses valid design column names (e.g., group-level columns from your ~ 0 + group design, or whatever your design_formula yields), the label on the LHS can be any descriptive string.
# Per-contrast rows (proteins) still supported. If you've precomputed DAList$data_per_contrast[[label]] or DAList$filtered_proteins_per_contrast[[label]], the function uses them. It does not drop samples based on label guesses.
# This should resolve the error you hit and let you define interaction contrasts like:
#  contrast_info$involved_levels is conservative: it only includes values from metadata[[group_col]] whose exact strings also appear as design columns and in the RHS expression. This avoids false positives and works for both simple contrasts and interaction-type expressions whose terms are design columns.
#  If your add_design() stores the design matrix under a different field than $design$X or tracks group_col elsewhere, adjust the two lookups accordingly.
#  With this in place, write_per_contrast_csvs() (called by write_limma_tables()) can read results[[c]]$tags$contrast_info and stop guessing from labels.  

#' Run limma analysis per contrast on filtered proteins
#'
#' This function subsets the DAList for each contrast according to
#' `filtered_proteins_per_contrast` and fits the limma model separately.
#' It is useful when proteins are filtered differently per contrast.
#' After model fitting, the function computes rolling standard deviations (moving SDs)
#' of log fold changes based on intensity-sorted proteins, and computes z-scores.
#'
#' @param DAList A DAList object containing full data and per-contrast filtered proteins.
#' @param design_formula A formula for the model design (e.g., `~0 + group`, or `~0 + cell:treatment`).
#' @param pval_thresh P-value threshold used for defining significance. Default = 0.05.
#' @param lfc_thresh Log2 fold change threshold for significance. Default = 1.
#' @param adj_method Adjustment method for multiple testing ("BH", "BY", "none", etc). Default = "BH".
#' @param contrasts_file Optional CSV file with contrast definitions if not present in the DAList.
#'   Each row should contain a full contrast statement of the form `Label = expression`.
#'   The label (left-hand side) may be any string; it does not need to encode group names.
#' @param binsize Either an integer for the moving window size, or
#'   \code{"auto"} (default) to select an appropriate bin size automatically.
#' @param plot_movingSD Logical. If TRUE (default), plot moving SD curves for each contrast.
#'
#' @return The input DAList, updated with per-contrast model fits and results,
#'   moving SDs, and logFC z-scores.
#'
#' @section Contrast metadata (for downstream writers):
#' For each contrast, a sidecar list is saved at
#' \code{DAList$tags$per_contrast[[label]]$contrast_info} so downstream
#' functions (e.g., \code{write_limma_tables()}) do not need to parse labels.
#' \code{contrast_info} contains:
#' \itemize{
#'   \item \code{label}: contrast label (left-hand side)
#'   \item \code{contrast_expression_raw}: the RHS as written in the file
#'   \item \code{contrast_expression}: the RHS actually passed to limma (after translation)
#'   \item \code{design_formula}: the design formula as a character string
#'   \item \code{group_col}: grouping column if applicable (or \code{NULL} for pure interaction designs)
#'   \item \code{factors}: factor names participating in a simple two-factor interaction (e.g., \code{c("cell","treatment")})
#'   \item \code{design_columns_involved}: model-matrix columns referenced by the translated expression
#'   \item \code{involved_levels}: levels from \code{group_col} appearing in the expression (if applicable)
#' }
#'
#' @examples
#' \dontrun{
#' # Using default binsize auto-selection
#' filtered_DAList <- run_filtered_limma_analysis(filtered_DAList)
#'
#' # With specified bin size and no plots
#' filtered_DAList <- run_filtered_limma_analysis(
#'     filtered_DAList, 
#'     binsize = 200, 
#'     plot_movingSD = FALSE)
#'}
#' @export
run_filtered_limma_analysis <- function(
    DAList,
    design_formula = ~ 0 + group,
    pval_thresh = 0.05,
    lfc_thresh = 1,
    adj_method = "BH",
    contrasts_file = NULL,
    binsize = "auto",
  #  binsize_range = c(50, 100, 200, 400, 500, 1000),
    plot_movingSD = FALSE
) {
  validate_DAList(DAList)
  if (!exists("compute_movingSD_zscores")) {
    stop("Function compute_movingSD_zscores() must be sourced before running this.")
  }
  
  # helpers -------------------------------------------------------------
  .re_escape <- function(x) gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", x)
  
  .replace_tokens <- function(expr, map) {
    if (is.null(map) || length(map) == 0) return(expr)
    keys <- names(map)[order(nchar(names(map)), decreasing = TRUE)]
    out <- expr
    for (k in keys) {
      patt <- paste0("(?<![A-Za-z0-9_.])", .re_escape(k), "(?![A-Za-z0-9_.])")
      out <- gsub(patt, map[[k]], out, perl = TRUE)
    }
    out
  }
  
  # 1) ensure a design exists ------------------------------------------
  DAList <- add_design(DAList = DAList, design_formula = design_formula)
  
  # 2) build contrast map (names = labels; values = "Label = expr") ----
  contrast_map <- NULL
  if (!is.null(DAList$design$contrast_vector)) {
    cv  <- DAList$design$contrast_vector
    lhs <- stringr::str_trim(stringr::str_split_fixed(cv, "=", 2)[, 1])
    contrast_map <- stats::setNames(cv, lhs)
  } else if (!is.null(contrasts_file)) {
    contrast_table <- utils::read.csv(contrasts_file, header = FALSE, stringsAsFactors = FALSE)
    cv  <- as.character(contrast_table[[1]])
    lhs <- stringr::str_trim(stringr::str_split_fixed(cv, "=", 2)[, 1])
    contrast_map <- stats::setNames(cv, lhs)
    cli::cli_alert_info("Loaded contrasts from file: {.file {contrasts_file}}")
  } else if (!is.null(DAList$filtered_proteins_per_contrast)) {
    cli::cli_abort(c(
      "Could not determine contrast expressions for the provided contrast names",
      "i" = "Provide {.arg contrasts_file} or set {.code DAList$design$contrast_vector} first."
    ))
  } else {
    cli::cli_abort(c(
      "Could not determine contrasts",
      "i" = "Provide {.arg contrasts_file} or set {.code DAList$design$contrast_vector} in DAList."
    ))
  }
  
  # containers ----------------------------------------------------------
  labels <- names(contrast_map)
  filtered_results <- list()
  filtered_ebayes  <- list()
  if (is.null(DAList$tags)) DAList$tags <- list()
  if (is.null(DAList$tags$per_contrast)) DAList$tags$per_contrast <- list()
  
  # main loop -----------------------------------------------------------
  for (label in labels) {
    cli::cli_inform(paste("Processing contrast:", label))
    
    # extract RHS expr from "Label = expr"
    cont_string <- contrast_map[[label]]
    split2 <- stringr::str_split_fixed(cont_string, "=", 2)
    rhs_expr <- if (ncol(split2) == 2) stringr::str_trim(split2[, 2]) else cont_string
    orig_rhs <- rhs_expr
    
    # subset rows (proteins) for this contrast if provided
    if (!is.null(DAList$data_per_contrast) && !is.null(DAList$data_per_contrast[[label]])) {
      sub_data <- DAList$data_per_contrast[[label]]
    } else if (!is.null(DAList$filtered_proteins_per_contrast) && !is.null(DAList$filtered_proteins_per_contrast[[label]])) {
      keep_proteins <- DAList$filtered_proteins_per_contrast[[label]]
      sub_data <- DAList$data[keep_proteins, , drop = FALSE]
    } else {
      sub_data <- DAList$data
    }
    
    # per-contrast view: same samples/metadata; only rows change
    sub_DAList <- DAList
    sub_DAList$filtered_proteins_per_contrast <- NULL
    sub_DAList$data <- sub_data
    sub_DAList$annotation <- DAList$annotation[rownames(sub_DAList$data), , drop = FALSE]
    sub_DAList$metadata <- DAList$metadata[colnames(sub_DAList$data), , drop = FALSE]
    rownames(sub_DAList$metadata) <- colnames(sub_DAList$data)
    
    # rebuild design for the current sample set
    sub_DAList <- add_design(DAList = sub_DAList, design_formula = design_formula)
    
    # robust token -- design-column translation for 2-factor interactions
    X <- tryCatch(sub_DAList$design$X, error = function(e) NULL)
    if (!is.null(X)) sub_DAList$design$design_matrix <- X  # ensure slot for add_contrasts()
    dcols <- if (!is.null(sub_DAList$design$design_matrix)) colnames(sub_DAList$design$design_matrix) else character(0)
    
    rhs_tokens <- unique(stringr::str_extract_all(rhs_expr, "[A-Za-z0-9_.]+")[[1]])
    rhs_tokens <- rhs_tokens[!is.na(rhs_tokens) & nzchar(rhs_tokens)]
    term_labels <- tryCatch(attr(terms(design_formula), "term.labels"), error = function(e) character(0))
    has_two_factor_interaction_only <- length(term_labels) == 1 && grepl(":", term_labels)
    
    token_map <- list()
    factors <- NULL
    
    if (has_two_factor_interaction_only && length(dcols) > 0) {
      fvars <- strsplit(term_labels, ":", fixed = TRUE)[[1]]
      fvars <- trimws(fvars)
      f1 <- fvars[1]; f2 <- fvars[2]
      factors <- c(f1, f2)
      
      f1_lvls <- if (f1 %in% colnames(sub_DAList$metadata)) levels(droplevels(factor(sub_DAList$metadata[[f1]]))) else character(0)
      f2_lvls <- if (f2 %in% colnames(sub_DAList$metadata)) levels(droplevels(factor(sub_DAList$metadata[[f2]]))) else character(0)
      
      find_interaction_col_variants <- function(lv1, lv2) {
        variants <- c(
          paste0(f1, lv1, ":", f2, lv2),
          paste0(       lv1, ":", f2, lv2),
          paste0(f1, lv1, ".", f2, lv2),
          paste0(       lv1, ".", f2, lv2)
        )
        hit <- dcols[dcols %in% variants]
        if (length(hit)) return(hit[1])
        # last resort fuzzy
        hit <- dcols[
          grepl(paste0("(^|[:.])", .re_escape(lv1), "([:.]|$)"), dcols) &
            grepl(paste0("(^|[:.])", .re_escape(lv2), "([:.]|$)"), dcols) &
            grepl(.re_escape(f2), dcols)
        ]
        if (length(hit)) return(hit[1])
        character(0)
      }
      
      for (tok in rhs_tokens) {
        if (tok %in% dcols) next
        if (!grepl("_", tok, fixed = TRUE)) next
        parts <- strsplit(tok, "_", fixed = TRUE)[[1]]
        if (length(parts) != 2) next
        
        cand <- character(0)
        if (parts[1] %in% f1_lvls && parts[2] %in% f2_lvls) {
          cand <- find_interaction_col_variants(parts[1], parts[2])
        } else if (parts[1] %in% f2_lvls && parts[2] %in% f1_lvls) {
          cand <- find_interaction_col_variants(parts[2], parts[1])
        }
        if (length(cand)) token_map[[tok]] <- cand
      }
      
      if (length(token_map)) {
        rhs_expr <- .replace_tokens(rhs_expr, token_map)
      }
    }
    
    # optional: print translated RHS & design columns
    cli::cli_inform(c(
      "i" = "Contrast '{label}' RHS (post-translation): {rhs_expr}",
      "i" = "Design columns: {paste(dcols, collapse = ', ')}"
    ))
    
    # add this (possibly translated) contrast
    cont_string_translated <- paste(label, "=", rhs_expr)
    sub_DAList <- add_contrasts(sub_DAList, contrasts_vector = cont_string_translated)
    
    # fit & extract
    sub_DAList <- fit_limma_model(sub_DAList)
    sub_DAList <- extract_DA_results(
      DAList = sub_DAList,
      pval_thresh = pval_thresh,
      lfc_thresh  = lfc_thresh,
      adj_method  = adj_method
    )
    
    if (!is.null(sub_DAList$results[[label]])) {
      filtered_results[[label]] <- sub_DAList$results[[label]]
      filtered_ebayes[[label]]  <- sub_DAList$eBayes_fit
      
      # contrast metadata sidecar --------------------------------------
      design_cols <- colnames(sub_DAList$design$design_matrix)
      tokens_for_cols <- unique(stringr::str_extract_all(rhs_expr, "[A-Za-z0-9_.:]+")[[1]])
      cols_in_expr <- intersect(tokens_for_cols, design_cols)
      
      # group_col if available (mainly for ~0 + group models)
      group_col <- NULL
      if (!is.null(sub_DAList$design$group_col) && is.character(sub_DAList$design$group_col)) {
        group_col <- sub_DAList$design$group_col
      } else if ("group" %in% colnames(sub_DAList$metadata)) {
        group_col <- "group"
      }
      
      involved_levels <- character(0)
      if (!is.null(group_col) && group_col %in% colnames(sub_DAList$metadata)) {
        grp_lvls <- levels(droplevels(factor(sub_DAList$metadata[[group_col]])))
        involved_levels <- intersect(rhs_tokens, grp_lvls)
      }
      
      DAList$tags$per_contrast[[label]] <- list(
        contrast_info = list(
          label = label,
          contrast_expression_raw = orig_rhs,
          contrast_expression = rhs_expr,
          design_formula = paste(deparse(design_formula), collapse = ""),
          group_col = group_col,
          factors = factors,
          design_columns_involved = cols_in_expr,
          involved_levels = involved_levels
        )
      )
    } else {
      cli::cli_alert_warning("No results found for contrast '{label}' -- skipping.")
    }
  }
  
  DAList$results    <- filtered_results
  DAList$eBayes_fit <- filtered_ebayes
  
  # moving SD z-scores --------------------------------------------------
  DAList <- compute_movingSD_zscores(
    DAList = DAList,
    binsize = binsize,
  #  binsize_range = binsize_range,
    plot = plot_movingSD,
    contrasts_file = contrasts_file
  )
  
  DAList$tags$DA_criteria <- list(
    pval_thresh = pval_thresh,
    lfc_thresh  = lfc_thresh,
    adj_method  = adj_method
  )
  
  new_DAList(DAList)
}

