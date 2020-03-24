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

# Take the seff time result in hours:minutes:seconds and convert the 
# results to hours
# KS modified to take data of the form: 1-15:04:51, where 1 == 1 day, so it's actually 24+16=39 hours. 
ConvertTimeToHours = function(timeStr)
{
  days <- 0
  if (grepl("-", timeStr)) {
    days <- as.numeric(unlist(strsplit(timeStr, "-"))[1])
    
  }
  timeStr_without_days <- gsub(".*[-]", "", timeStr)
  times   = as.numeric( unlist(strsplit(timeStr_without_days, ':') ) )
  expVals = c(0, 1, 2)
  hours   = sum( sapply( seq(3), function(i) times[i] / 60**expVals[i] ) )
  days_in_hours <- days * 24
  total <- sum(hours, days_in_hours)
  return(  round(total, 3) )
}

# ------------------------------------------------------------------------------
