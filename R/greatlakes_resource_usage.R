# ------------------------------------------------------------------------------
# GetResourceUsage
# 2020/01/29
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Take a job id and return the Job Wall-clock time in hours and the Memory
# Utilized in GB as a numeric vector
# ------------------------------------------------------------------------------

# 2020/03/24 Katie Saund added two catches for edge cases:
# (1) When the memory usage is basically 0MB
# (2) When the run time is more than 1 day

#' Get memory usage and run time for a Great Lakes sbatch job given the job id
#'
#' @description This function needs to be run from the command line, not within
#'   an Rstudio session due to a call of system(). T
#'
#' @param jobId Character or numeric
#'
#' @return Numeric Vector. Memory in GB and time in hours.
#' @export
GetResourceUsage = function( jobId )
{
  # Run the system command to get the seff result for this job ID
  seff = system( paste("seff", jobId ), intern = TRUE)

  # Find the position of the time and the memory utilized in the
  # seff result
  wallClockStr = "Job Wall-clock time: "
  wallTime = gsub(wallClockStr, '', seff[ grep(wallClockStr, seff) ])
  memoryUseStr   = "Memory Utilized: "
  memoryUsage = gsub(memoryUseStr, '', seff[ grep(memoryUseStr, seff) ])

  # Catch the case where, if memory usage from seff is essentially 0 it says
  #   "0.00 MB (estimated maximum)"
  if (grepl("(estimated maximum)", memoryUsage)) {
    memoryUsage <- "0.00 MB"
  }

  # If memory usage from seff if returned in MB, convert to GB
  if ( grepl("MB", memoryUsage) )
  {
    memoryUsage = gsub("MB", '', memoryUsage)
    memoryUsage = as.numeric( trimws(memoryUsage) ) / 1000
  } else if ( grepl("GB", memoryUsage) ) {
    memoryUsage = gsub("GB", '', memoryUsage)
    memoryUsage = as.numeric( trimws(memoryUsage) )
  } else {
    memoryUsage = gsub("KB", '', memoryUsage)
    memoryUsage = as.numeric( trimws(memoryUsage) ) * 1000
  }

  return( c(ConvertTimeToHours(wallTime), memoryUsage) )
}




#' Convert sbatch job run time into hours
#'
#' @description the `$ seff <jobid>` command on Great Lakes returns a string
#'   that includes the run time of an sbatch job. These strings look like:
#'   * Job Wall-clock time: 10:30:00
#'   or
#'   * Job Wall-clock time: 1-10:30:00
#'   and indicate 10.5 hours or 34.5 hours, respectively. Given such a string
#'   this function converts the number portion into just hours. The time is
#'   recorded as hours:minutes:seconds.
#'
#' @param timeStr String of the form: 1-10:30:00 or 10:30:00. The number before
#'   the dash indicates days.
#'
#' @return Number of hours rounded to 3 digits. Numeric.
#' @export
#'
#' @examples
#' ConvertTimeToHours("1-10:30:00")
#' ConvertTimeToHours("10:30:00")
ConvertTimeToHours = function(timeStr)
{
  # Check function input -------------------------------------------------------
  if (!is.character(timeStr)) {
    stop("Input must be a character string.")
  }
  if (!grepl("[0-9][0-9][:][0-9][0-9][:][0-9][0-9]", timeStr)) {
    stop("Input must be of the format 10:00:00 or 1-10:00:00")
  }

  # Function -------------------------------------------------------------------
  days <- 0
  # When there are 1 or more days of run time
  if (grepl("-", timeStr)) {
    days <- as.numeric(unlist(strsplit(timeStr, "-"))[1])
  }
  days_in_hours <- days * 24

  # Time ignoring days
  timeStr_without_days <- gsub(".*[-]", "", timeStr)
  times <- as.numeric( unlist(strsplit(timeStr_without_days, ':') ) )
  expVals = c(0, 1, 2)
  hours   = sum( sapply( seq(3), function(i) times[i] / 60**expVals[i] ) )

  # Add hours from days plus hours without days together
  total <- sum(hours, days_in_hours)
  return(round(total, 3))
}

# ------------------------------------------------------------------------------
