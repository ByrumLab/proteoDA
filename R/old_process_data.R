
#' Process raw protein data
#'
#' A wrapper function, which calls some fubfunctions to filter and normalize
#' raw data. First, calls \code{\link{filter_data}}, which does some filtering
#' that I do not understand yet. Then, calls \code{\link{normalize_data}}, which
#' runs multiple normalization methods on the the filtered data, so that later
#' you can pick the best normalization method. Also calls \code{\link{make_logs}}
#' to make some log files.
#'
#' @param data Raw data to be processed. In the pipeline, this is generally the
#'   "data" slot of the list output by \code{\link{extract_data}}.
#' @param targets An inputs target dataframe. In the pipeline, this is generally
#'   the subsetted targets data frame (that is, with pool or other non-desired
#'   samples removed).
#' @param group_col Optional. The name of the column in the targets dataframe which
#'   lists the experimental group names. Default is "group"
#' @param min.reps Optional. NEED TO FIGURE OUT
#' @param min.grps Optional. NEED TO FIGURE OUT
#'
#' @return A list with five elements:
#'   \enumerate{
#'     \item "normList"- a list of length 8, where each item in the list is a
#'       named dataframe. Names give the normalization method, and the dataframe
#'       gives the normalized intensity data
#'     \item "targets" - A dataframe of targets. Should be nearly the same as
#'       the targets dataframe that was input, possibly with new rownames.
#'     \item "filt"- The output of the \code{\link{filter_data}} function: a list
#'       of length 5 containing filtered data and other stats.
#'       See \code{\link{filter_data}}.
#'     \item "param"- A dataframe giving the parameters and arguments used for
#'       data processing
#'     \item "stats"- A dataframe giving statistics on data processing and
#'       filtering.
#'   }
#' @export
#'
#' @examples
#' # No examples yet
#'
#'
process_data <- function(data,
                         targets,
                         group_col = "group",
                         min.reps = 3,
                         min.grps = 1) {

  cli::cli_rule()

  ## conditional filter applied to the data
  filt <- filter_data(data = data,
                      targets = targets,
                      group_col = group_col,
                      min.reps = min.reps,
                      min.grps = min.grps)

  ## apply 8 norm. methods to the filtered data set.
  ## output returned includes list object containing df for each norm. tech.
  ## filtered unnormalized dataset (data) is also returned

  # TODO: As far as I can tell, the above statement is not true:
  # normalize data only returns the 8 normalized data frames,
  # it does not return a raw data frame. Not sure if it needs to?
  norm <- normalize_data(data = filt$data, targets = filt$targets)

  param <- filt$param
  colnames(param) <- NULL

  stats <- filt$stats
  colnames(stats) <- NULL

  cli::cli_inform("Writing logs")
  logs <- make_log(param = as.list(unlist(param)),
                   stats = as.list(unlist(stats)),
                   title = "DATA PROCESSING",
                   save = TRUE)

  cli::cli_rule()
  cli::cli_inform(c("v" = "Data processing complete"))

  # Return output
  list(normList = norm$normList,
       targets = norm$targets,
       filt = filt,
       param = logs$param,
       stats = logs$stats)

}


#' Filter data
#'
#' NEED TO WRITE UP
#' subfunctions: make_logs
#'
#' @inheritParams process_data
#'
#' @return NEED TO WRITE UP
#' @export
#'
#' @examples
#' # No examples yet
#'
filter_data <- function(data,
                        targets,
                        group_col = "group",
                        min.reps = 3,
                        min.grps = 1) {

  # Check args
  if (any(rownames(targets) %notin% colnames(data))) {
    missing_targets <- rowname(targets)[rowname(targets) %notin% colnames(data)]

    cli::cli_abort(c("Some targets are not present in input data.",
                     "x" = "Missing target{?s}: {missing_targets}",
                     "i" = "Check {.code colnames({data})"))

  }

  if (group_col %notin% colnames(targets)) {
    cli::cli_abort(c("Column {.arg {group_col}} not found in {.arg targets} data frame"))
    }


  param <- stats <- list()

  stats[["total_input_samples"]] <- ncol(data)



  data <- data[, rownames(targets)]
  stopifnot(rownames(targets) == colnames(data))

  stats[["no_removed_samples"]] <- stats[["total_input_samples"]] - ncol(data)
  stats[["no_filtered_samples"]] <- ncol(data)
  stats[["total_input_rows"]] = nrow(data)

  print(paste(stats[["no_removed_samples"]], "samples were removed from the data matrix..."))

  ## change data column names and targets row names
  ## to sample name i.e. sample column in targets
  if ("sample" %in% colnames(targets)) {
    stopifnot(rownames(targets) == colnames(data))
    print("FILTER_DATA: column names of data and row names of targets converted to targets sample names...");cat("\n")
    rownames(targets) <- targets$sample
    colnames(data)    <- rownames(targets)
    stopifnot(colnames(data) == rownames(targets))
  }


  ## extract group column as character vector and make a factor.
  groups <- make_factor(as.character(targets[, group_col]))

  param[["group"]] <- group_col
  param[["groups"]] <- paste(paste0(names(table(groups)), "=" , table(groups)) ,collapse = ", ")

  nreps <- table(groups)              ##  number of samples in each group
  ngrps <- length(unique(groups))      ##  number of groups
  # print(group.names);print(nreps);print(ngrps)



  ## if min.reps is NULL set min.reps to the smallest group size
  ## if min.grps is NULL set the min.grps to 1
  if (is.null(min.reps)) {min.reps <- min(nreps)}
  if (is.null(min.grps)) {min.grps <- 1}

  ## if input no min.reps exceeds the max no. reps / group
  ## min.reps value is lowered to the smallest samle group.
  repsCutoff <- (min.reps <= nreps)
  if (all(repsCutoff) == FALSE) {
    warning(paste0("The min.reps threshold min.reps = ", min.reps, " exceeds the max. ",
                   "number of replicates per group: ", paste(names(nreps),nreps,sep="=",collapse=", "),
                   ". \nThe min.reps ",
                   "threshold was lowered to equal the smallest sample group.\n",
                   "min.reps = ", min(nreps)))
    min.reps <- min(nreps)
    repsCutoff <- (min.reps <= nreps)
  }

  ## repsCutoff is used to determine how many groups meet the min.reps criteria.
  ## if min.grps parameter is greater than grpsCutoff, then the threshold is lowered
  ## to one (i.e. )
  grpsCutoff <- sum(as.numeric(repsCutoff))

  if (min.grps > grpsCutoff) {
    message(paste0("Based on the min.reps parameter ",min.reps,". The min.grps threshold : ",
                   min.grps," exceeds what is allowed by the dataset: ", grpsCutoff,".\n",
                   "The min.grps threshold has been lowered to one ",
                   "group. \nmin.grps = ", 1))
    min.grps <- 1
  }

  param[["min.reps"]] <- min.reps
  param[["min.grps"]] <- min.grps
  print(paste0("extracting entries with intensity > 0 in at least ", min.reps,
               " of the samples in ", min.grps," or more groups..."))

  ## FILTERING
  ## zeros are replaced with NA
  tmpData <- data
  tmpData[,][tmpData[,] == 0] <- NA

  ## calc. no samples in each group with intensities > 0
  noSamplesPerGroup <- NULL
  for (x in unique(groups)) {
    keep<-groups %in% x
    grpData<- data.frame(tmpData[,keep])
    ## no samples in each group with intensities > 0
    noSamplesPerGroup <- cbind(noSamplesPerGroup, apply(grpData, 1, FUN = function(x) {sum(!is.na(x))}))
  }
  noSamplesPerGroup <- data.frame(noSamplesPerGroup)
  colnames(noSamplesPerGroup) <- unique(groups)
  rownames(noSamplesPerGroup) <- rownames(grpData)
  head(noSamplesPerGroup)

  ## min X samples in at least Y groups
  aboveCutoff <- apply(noSamplesPerGroup, 1 , FUN = function(x) {sum(x >= min.reps) >= min.grps })
  noSamplesPerGroup$aboveCutoff <-aboveCutoff
  filterData <- tmpData[aboveCutoff == TRUE, ]
  removeData <- tmpData[aboveCutoff == FALSE, ]

  stats[["no_removed_rows"]] = nrow(removeData)
  stats[["no_filtered_rows"]] = nrow(filterData)

  ## filtering processing stats
  print(paste("A total of", nrow(removeData), "entries were removed from the data set leaving",
              nrow(filterData), "entries for further analysis. Success!!"))

  logs <- make_log(param, stats, title="FILTERING X in Y", save = TRUE)

  # Output results
  list(data = filterData,
      targets = targets,
      noSamplesPerGroup = noSamplesPerGroup,
      rm.data = removeData,
      param = logs$param,
      stats = logs$stats)
}


