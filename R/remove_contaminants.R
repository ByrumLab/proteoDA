#' Remove protein contaminant rows
#'
#' The second subfunction called by \code{\link{read_DIA_data}}. Removes
#' contaminant rows from the data.
#'
#' @param data Raw Maxquant intensity data, from which contaminants need to be
#'   removed. Within the \code{\link{read_DIA_data}} function, this is "diaData",
#'   which is the data object from the \code{\link{read_maxquant_delim}} function that has
#'   gone through additional processing within the body of the
#'   \code{\link{read_DIA_data}} function
#'
#' @return A list of data with contaminants removed, with three elements:
#'   \enumerate{
#'     \item "data", a data frame with the contaminants removed. Rows are
#'       proteins, columns include both protein info and individual intensities.
#'     \item "num_contam", an integer giving the number of contaminant rows removed
#'     \item "num_qfilter", an integer giving the number of retained rows
#'       (ie, nrow() of the "data" slot).
#'   }
#'
#' @export
#'
#' @examples
#' # No examples yet
#'

remove_contaminants <- function(data) {
  ## vector of required contaminant column names
  contamColums <- diaContamColums

  num_rows <- nrow(data)
  data <- data[!grepl("DECOY", data[, contamColums]),]
  data <- data[!grepl("Group of", data[, contamColums]),]

  num_contam  <- num_rows - nrow(data)
  num_qfilter <- nrow(data)

  cli::cli_inform("{num_contam} contaminantes removed")
  cli::cli_inform("{num_qfilter} DIA protein entries retained")

  list(data = data,
       num_contam = num_contam,
       num_qfilter = num_qfilter)
}
