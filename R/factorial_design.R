# ---- helpers ---------------------------------------------------------------


# Rescale a contrast's logFC by `factor` (default 0.5) in your DAList-like object.
# Keeps t / p / q / B unchanged. Also rescales CI.L/CI.R if present.
# Optionally renames the contrast and updates design + tags consistently.
#' Rescale log-fold changes for a single contrast
#'
#' Rescales the \code{logFC} (and optional \code{CI.L}/\code{CI.R}) for one
#' contrast inside a DAList-like object, leaving inferential statistics
#' (\code{t}, \code{P.Value}, \code{adj.P.Val}, \code{B}, etc.) unchanged.
#' Optionally renames the contrast and updates contrast metadata/tags.
#'
#' The function is robust to several result table shapes:
#' \itemize{
#'   \item a plain \code{data.frame} with a \code{logFC} column;
#'   \item a list with a \code{$table} \code{data.frame};
#'   \item a list with a numeric \code{$logFC} vector;
#'   \item a plain numeric vector (assumed to be \code{logFC}).
#' }
#'
#' When \code{new_label} is provided, the function also:
#' \itemize{
#'   \item renames the entry in \code{x$results};
#'   \item moves any per-contrast tags from the old label to the new one;
#'   \item renames the corresponding column in \code{x$design$contrast_matrix}, if present.
#' }
#'
#' A cumulative \code{rescale_factor} is tracked in
#' \code{x$tags$per_contrast[[label]]$rescale_factor}.
#'
#' @param x A DAList-like object containing \code{$results}, and optionally
#'   \code{$tags} and \code{$design$contrast_matrix}.
#' @param label Character. Name of the contrast to rescale (must exist in \code{x$results}).
#' @param factor Numeric scalar. Multiplicative factor applied to \code{logFC}
#'   (and \code{CI.L}/\code{CI.R} if present). Default \code{0.5}.
#' @param new_label Optional character. If supplied and different from \code{label},
#'   the contrast is renamed consistently across \code{$results}, \code{$tags},
#'   and \code{$design$contrast_matrix}.
#'
#' @return The modified object \code{x}.
#'
#' @examples
#' \dontrun{
#' # Halve logFCs for a contrast and rename it
#' results <- rescale_contrast_logFC(results, label = "Avg_Treat_Bio_vs_DMSO",
#'                                   factor = 0.5, new_label = "Avg_Treat_x0.5")
#' }
#'
#' @export

rescale_contrast_logFC <- function(x, label, factor = 0.5, new_label = NULL) {
  if (!("results" %in% names(x)) || is.null(x$results[[label]])) {
    warning(sprintf("Contrast '%s' not found; nothing rescaled.", label))
    return(x)
  }
  
  # Helper to rescale columns in a data.frame without touching stats
  .rescale_df <- function(df) {
    if (!is.data.frame(df)) return(df)
    if ("logFC" %in% names(df)) df$logFC <- df$logFC * factor
    # Rescale CIs if present
    if ("CI.L" %in% names(df)) df$CI.L <- df$CI.L * factor
    if ("CI.R" %in% names(df)) df$CI.R <- df$CI.R * factor
    # DO NOT change: t, P.Value, adj.P.Val, B, AveExpr, etc.
    df
  }
  
  # Warn if already rescaled before (based on tags)
  already <- tryCatch({
    z <- x$tags$per_contrast[[label]]$rescale_factor
    is.numeric(z) && is.finite(z) && z != 1
  }, error = function(e) FALSE)
  if (already) {
    prev <- x$tags$per_contrast[[label]]$rescale_factor
    warning(sprintf("Contrast '%s' appears already rescaled (factor=%.3g). Applying an additional factor=%.3g.",
                    label, prev, factor))
  }
  
  res <- x$results[[label]]
  
  # Case A: plain data.frame with logFC
  if (is.data.frame(res)) {
    res <- .rescale_df(res)
    
    # Case B: list with $table data.frame
  } else if (is.list(res) && !is.null(res$table) && is.data.frame(res$table)) {
    res$table <- .rescale_df(res$table)
    
    # Case C: list with a numeric $logFC vector
  } else if (is.list(res) && !is.null(res$logFC) && is.numeric(res$logFC)) {
    res$logFC <- res$logFC * factor
    
    # Case D: plain numeric vector (assume it's logFC)
  } else if (is.numeric(res)) {
    res <- res * factor
    
  } else {
    warning(sprintf("Contrast '%s' structure not recognized; no changes made.", label))
    x$results[[label]] <- res
    return(x)
  }
  
  # Write back
  x$results[[label]] <- res
  
  # Update tags: mark the rescale
  if (is.null(x$tags)) x$tags <- list()
  if (is.null(x$tags$per_contrast)) x$tags$per_contrast <- list()
  if (is.null(x$tags$per_contrast[[label]])) x$tags$per_contrast[[label]] <- list()
  # accumulate factor if repeated
  old_f <- x$tags$per_contrast[[label]]$rescale_factor
  x$tags$per_contrast[[label]]$rescale_factor <- if (is.null(old_f)) factor else as.numeric(old_f) * factor
  x$tags$per_contrast[[label]]$rescaled <- TRUE
  
  # Optional: rename the contrast key and keep design/metadata in sync
  if (!is.null(new_label) && nzchar(new_label) && new_label != label) {
    # Rename in results
    pos <- match(label, names(x$results))
    names(x$results)[pos] <- new_label
    x$results[[new_label]] <- x$results[[pos]]
    x$results[[label]] <- NULL
    
    # Move per-contrast tags
    if (!is.null(x$tags$per_contrast[[label]])) {
      x$tags$per_contrast[[new_label]] <- x$tags$per_contrast[[label]]
      x$tags$per_contrast[[label]] <- NULL
      # store the label & factor in the moved tag
      x$tags$per_contrast[[new_label]]$contrast_info$label <- new_label
    }
    
    # If a contrast_matrix exists, rename its column too
    if (!is.null(x$design) && !is.null(x$design$contrast_matrix)) {
      cm <- x$design$contrast_matrix
      if (!is.null(colnames(cm))) {
        colnames(cm)[colnames(cm) == label] <- new_label
        x$design$contrast_matrix <- cm
      }
    }
  }
  
  x
}

#####################  
# Extract a (logFC, adj.P.Val, P.Value, t, B) row for one protein & one contrast
#' @keywords internal
#' @noRd
.get_row_for_protein <- function(res_obj, protein, protein_col, DA_results) {
  tab <- res_obj
  if (is.null(tab) || NROW(tab) == 0) return(NULL)
  
  # If the table already has an explicit ID column
  if (!is.null(protein_col) && protein_col %in% colnames(tab)) {
    hit <- tab[tab[[protein_col]] == protein, , drop = FALSE]
    if (NROW(hit)) return(hit[1, , drop = FALSE]) else return(NULL)
  }
  # Otherwise resolve via rownames/annotation
  rix <- .resolve_row_index(
    DA_results = DA_results,
    res_tab    = tab,
    protein    = protein,
    ann_id_col = "uniprot_id",
    ann_row_col = NULL
  )
  if (is.na(rix)) return(NULL)
  tab[rix, , drop = FALSE]
}

# Nice sign + significance label
#' @keywords internal
#' @noRd
.effect_label <- function(logFC, adjp, alpha = 0.05) {
  if (is.na(logFC)) return("NA")
  dir <- if (logFC > 0) "\u2191" else if (logFC < 0) "\u2193" else "0"
  sig <- if (!is.na(adjp) && adjp < alpha) "*" else ""
  paste0(dir, sig)
}

# Try hard to map `protein` (e.g., "Q16666") to a row in `res_tab`
# using DA_results$annotation as needed.
# Normalize any single-row data.frame to these columns in this order
#' @keywords internal
#' @noRd
.norm_cols <- function(df) {
  needed <- c("logFC","adj.P.Val","P.Value","t","B")
  for (nm in needed) if (!(nm %in% names(df))) df[[nm]] <- NA_real_
  # Keep only needed (but don't drop if already single row)
  df[, needed, drop = FALSE]
}

#' @keywords internal
#' @noRd
.resolve_row_index <- function(DA_results, res_tab, protein,
                               ann_id_col = "uniprot_id",
                               ann_row_col = NULL) {
  ann <- DA_results$annotation
  rn  <- rownames(res_tab)
  
  # Case A: rownames equal UniProt IDs
  if (!is.null(rn)) {
    if (protein %in% rn) return(protein)
    if (!is.null(ann) && ann_id_col %in% names(ann) &&
        sum(rn %in% ann[[ann_id_col]]) > 0L && protein %in% ann[[ann_id_col]]) {
      if (protein %in% rn) return(protein)
    }
  }
  
  # Case B: rownames are some other feature key; auto-detect
  if (is.null(ann_row_col) && !is.null(ann)) {
    overlaps <- vapply(names(ann), function(col) sum(ann[[col]] %in% rn, na.rm = TRUE), numeric(1))
    ann_row_col <- if (any(overlaps > 0)) names(ann)[which.max(overlaps)] else NULL
  }
  
  if (!is.null(ann) && !is.null(ann_row_col) &&
      ann_id_col %in% names(ann) && ann_row_col %in% names(ann)) {
    idx_ann <- match(protein, ann[[ann_id_col]])
    if (!is.na(idx_ann)) {
      candidate_rn <- ann[[ann_row_col]][idx_ann]
      if (!is.na(candidate_rn) && candidate_rn %in% rn) return(candidate_rn)
    }
  }
  
  # Case C: last resort. look for an ID column inside res_tab
  possible_cols <- intersect(c("uniprot_id","Protein","Gene","Accession"), colnames(res_tab))
  for (col in possible_cols) {
    hit <- which(res_tab[[col]] == protein)
    if (length(hit) == 1L) return(rownames(res_tab)[hit])
  }
  
  NA_character_
}
# ---- main function ---------------------------------------------------------
#' Interpret a protein across factorial contrasts
#'
#' Extracts per-contrast statistics for a single protein from a DAList-like
#' object fitted with a factorial (e.g., BioID) design, and returns a compact
#' table plus a human-readable summary.
#'
#' Row matching is flexible:
#' \itemize{
#'   \item If \code{protein_col} names a column in the per-contrast table,
#'         it is used directly.
#'   \item Otherwise, the function attempts to resolve the row using
#'         \code{rownames()} or \code{DA_results$annotation} (preferring
#'         \code{annotation$uniprot_id} if available).
#' }
#'
#'
#' The returned table includes an \code{effect} label indicating direction
#' (up/down) and significance (\code{"*"} if \code{adj.P.Val < alpha}).
#' 
#' @param DA_results A DAList-like object with \code{$results} (per-contrast
#'   result tables) and optionally \code{$annotation}.
#' @param protein Character. Protein identifier (e.g., UniProt ID) to look up.
#' @param protein_col Optional character. Name of an explicit ID column in the
#'   per-contrast result tables. If \code{NULL}, matching falls back to row names
#'   and/or \code{DA_results$annotation}.
#' @param alpha Numeric. FDR threshold used to annotate the \code{effect} column.
#'   Default \code{0.05}.
#' @param contrast_map Named character vector mapping display labels (names)
#'   to keys inside \code{DA_results$results}. Defaults assume contrasts like
#'   average treatment effect, interaction, and per-cell effects.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{table}}{A \code{data.frame} with columns
#'         \code{contrast}, \code{logFC}, \code{adj.P.Val}, \code{effect},
#'         \code{t}, \code{B}, \code{P.Value}.}
#'   \item{\code{summary}}{A single character string summarizing the key contrasts.}
#' }
#'
#' @examples
#' \dontrun{
#' out <- interpret_protein_factorial(DA_results, protein = "Q16666")
#' out$summary
#' out$table
#' }
#'
#' @export

interpret_protein_factorial <- function(
    DA_results,
    protein,
    protein_col = NULL,
    alpha = 0.05,
    contrast_map = c(
      Avg_Treat_Bio_vs_DMSO      = "Avg_Treat_Bio_vs_DMSO",
      Interaction_CHLA_vs_SKNF1  = "Interaction_CHLA_vs_SKNF1",
      TreatEffect_CHLA90         = "TreatEffect_CHLA90",
      TreatEffect_SKNF1          = "TreatEffect_SKNF1",
      CHLA90_Bio_vs_SKNF1_Bio    = "CHLA90_Bio_vs_SKNF1_Bio",
      CHLA90_DMSO_vs_SKNF1_DMSO  = "CHLA90_DMSO_vs_SKNF1_DMSO"
    )
) {
  rows <- lapply(names(contrast_map), function(lbl) {
    key <- contrast_map[[lbl]]
    res_obj <- DA_results$results[[key]]
    
    if (is.null(res_obj)) {
      out <- data.frame(contrast = lbl, logFC = NA_real_, adj.P.Val = NA_real_,
                        P.Value = NA_real_, t = NA_real_, B = NA_real_, check.names = FALSE)
      return(out[, c("contrast","logFC","adj.P.Val","P.Value","t","B"), drop = FALSE])
    }
    
    row <- .get_row_for_protein(res_obj, protein, protein_col, DA_results)
    if (is.null(row)) {
      out <- data.frame(contrast = lbl, logFC = NA_real_, adj.P.Val = NA_real_,
                        P.Value = NA_real_, t = NA_real_, B = NA_real_, check.names = FALSE)
      return(out[, c("contrast","logFC","adj.P.Val","P.Value","t","B"), drop = FALSE])
    } else {
      row_std <- .norm_cols(row)
      cbind(contrast = lbl, row_std)
    }
  })
  
  tab <- do.call(rbind, rows)
  
  .effect_label <- function(lfc, q, alpha = 0.05) {
    if (is.na(lfc) || is.na(q)) return("NA")
    if (q >= alpha) return("NS")
    if (lfc > 0) "\u2191*" else if (lfc < 0) "\u2193*" else "NS"
  }
  
  tab$effect <- mapply(.effect_label, tab$logFC, tab$adj.P.Val, MoreArgs = list(alpha = alpha))
  
  # ----- replace the block where avg/int/ch/sk and texts are built -----
  
  # length-1 safe row getter from tab
  
  .safe_row <- function(lbl) {
    r <- tab[tab$contrast == lbl, , drop = FALSE]
    if (nrow(r) == 0) {
      data.frame(logFC = NA_real_, adj.P.Val = NA_real_, P.Value = NA_real_,
                 t = NA_real_, B = NA_real_)
    } else {
      r[1, c("logFC","adj.P.Val","P.Value","t","B"), drop = FALSE]
    }
  }
  fmtp <- function(q) ifelse(is.na(q), "NA", sprintf("%.2g", q))
  
  avg <- .safe_row("Avg_Treat_Bio_vs_DMSO")
  int <- .safe_row("Interaction_CHLA_vs_SKNF1")
  ch  <- .safe_row("TreatEffect_CHLA90")
  sk  <- .safe_row("TreatEffect_SKNF1")
  
  avg_txt <- if (!is.na(avg$logFC[1])) {
    sprintf("Average Bio vs DMSO: log2FC=%.2f (adjP=%s)", avg$logFC[1], fmtp(avg$adj.P.Val[1]))
  } else "Average: NA"
  
  int_txt <- {
    if (!is.na(int$logFC[1])) {
      int_dir <- if (int$logFC[1] > 0) "stronger in CHLA90"
      else if (int$logFC[1] < 0) "stronger in SKNF1"
      else "similar across cells"
      sprintf("Interaction (CHLA90-SKNF1): log2FC=%.2f (%s; adjP=%s)",
              int$logFC[1], int_dir, fmtp(int$adj.P.Val[1]))
    } else "Interaction: NA"
  }
  
  ch_txt <- if (!is.na(ch$logFC[1])) {
    sprintf("CHLA90 Bio-DMSO: log2FC=%.2f (adjP=%s)", ch$logFC[1], fmtp(ch$adj.P.Val[1]))
  } else "CHLA90: NA"
  
  sk_txt <- if (!is.na(sk$logFC[1])) {
    sprintf("SKNF1  Bio-DMSO: log2FC=%.2f (adjP=%s)", sk$logFC[1], fmtp(sk$adj.P.Val[1]))
  } else "SKNF1: NA"
  
  summary <- paste(sprintf("[%s] %s.", protein, avg_txt), ch_txt, sk_txt, int_txt, sep = " ")
  
  list(table = tab[, c("contrast","logFC","adj.P.Val","effect","t","B","P.Value")],
       summary = summary)
}

### Batch ---------------
#' Classify proteins under a factorial design
#'
#' Applies a simple rule-based classifier to each protein using per-contrast
#' statistics extracted by \code{\link{interpret_protein_factorial}}.
#' Intended for BioID-style designs with average treatment, interaction,
#' and baseline (between-cell) contrasts.
#'
#' Classes returned include (non-exhaustive):
#' \itemize{
#'   \item \code{"Shared_Tx"} - significant average treatment effect; non-significant interaction.
#'   \item \code{"WTPreferential_Tx"} - significant average treatment; interaction significantly
#'         negative; baseline near zero.
#'   \item \code{"WTPreferential_Tx_withBaselineShift"} - as above but with significant baseline shift.
#'   \item \code{"TruncationGain"} - significant average treatment; interaction significantly
#'         positive; baseline near zero.
#'   \item \code{"TruncationGain_withBaselineShift"} - as above but with significant baseline shift.
#'   \item \code{"BaselineDifferenceOnly"} - significant baseline difference; average treatment not significant.
#'   \item \code{"Unclassified"} - none of the above rules matched.
#' }
#'
#' @param DA_results A DAList-like object with \code{$results} and optionally
#'   \code{$annotation}.
#' @param id_vector Optional character vector of protein IDs to classify. If
#'   \code{NULL}, the function tries \code{unique(DA_results$annotation$uniprot_id)}.
#' @param contrast_map Named character vector mapping display labels to
#'   keys in \code{DA_results$results}. Defaults match the typical factorial setup.
#' @param p_thresh Numeric. FDR threshold for significance tests. Default \code{0.05}.
#' @param lfc_thresh Numeric. Absolute log2 fold-change threshold used to decide
#'   whether the baseline contrast is \emph{near zero}. Default \code{1}.
#'
#' @return A \code{data.frame} with columns:
#' \itemize{
#'   \item \code{id} - protein ID,
#'   \item \code{class} - assigned class label,
#'   \item \code{base_logFC}, \code{base_q} - baseline contrast statistics,
#'   \item \code{treat_logFC}, \code{treat_q} - average treatment statistics,
#'   \item \code{int_logFC}, \code{int_q} - interaction statistics.
#' }
#'
#' @examples
#' \dontrun{
#' cls <- classify_factorial_batch(DA_results)
#' table(cls$class)
#'
#' # Use stricter thresholds
#' cls2 <- classify_factorial_batch(DA_results, p_thresh = 0.01, lfc_thresh = 1.5)
#' }
#'
#' @seealso \code{\link{interpret_protein_factorial}}
#' @export

classify_factorial_batch <- function(
    DA_results,
    id_vector = NULL,
    contrast_map = c(
      Avg_Treat_Bio_vs_DMSO      = "Avg_Treat_Bio_vs_DMSO",
      Interaction_CHLA_vs_SKNF1  = "Interaction_CHLA_vs_SKNF1",
      CHLA90_Bio_vs_SKNF1_Bio    = "CHLA90_Bio_vs_SKNF1_Bio",
      CHLA90_DMSO_vs_SKNF1_DMSO  = "CHLA90_DMSO_vs_SKNF1_DMSO"
    ),
    p_thresh = 0.05,
    lfc_thresh = 1
) {
  if (is.null(id_vector)) {
    if (!is.null(DA_results$annotation) && "uniprot_id" %in% names(DA_results$annotation)) {
      id_vector <- unique(DA_results$annotation$uniprot_id)
    } else stop("Provide id_vector or ensure annotation$uniprot_id exists.")
  }
  
  classify_one <- function(prot) {
    x <- interpret_protein_factorial(
      DA_results  = DA_results,
      protein     = prot,
      protein_col = NULL,
      contrast_map = contrast_map
    )$table
    
    getrow <- function(lbl) x[x$contrast == lbl, , drop = FALSE]
    avg  <- getrow("Avg_Treat_Bio_vs_DMSO")
    int  <- getrow("Interaction_CHLA_vs_SKNF1")
    base <- getrow("CHLA90_DMSO_vs_SKNF1_DMSO")
    
    sig   <- function(row) isTRUE(!is.na(row$adj.P.Val) && row$adj.P.Val < p_thresh)
    near0 <- function(row) isTRUE(!is.na(row$logFC) && abs(row$logFC) < lfc_thresh)
    
    class <- "Unclassified"
    if (sig(avg)) {
      if (sig(int) && isTRUE(int$logFC > 0)) {
        class <- if (near0(base)) "TruncationGain" else "TruncationGain_withBaselineShift"
      } else if (sig(int) && isTRUE(int$logFC < 0)) {
        class <- if (near0(base)) "WTPreferential_Tx" else "WTPreferential_Tx_withBaselineShift"
      } else {
        class <- "Shared_Tx"
      }
    } else if (sig(base)) {
      class <- "BaselineDifferenceOnly"
    }
    
    data.frame(
      id = prot,
      class = class,
      base_logFC = base$logFC, base_q = base$adj.P.Val,
      treat_logFC = avg$logFC, treat_q = avg$adj.P.Val,
      int_logFC = int$logFC, int_q = int$adj.P.Val,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, lapply(id_vector, classify_one))
}
