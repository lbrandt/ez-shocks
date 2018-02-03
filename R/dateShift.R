dateShift <- function(x, unit, rule) {
  # Shifts elements of date vector according to rule like MATLAB's dateshift
  # NOT YET ALL FEATURES
  
  x = as.POSIXlt(x)
  
  if(unit == 'month'){
    if(rule == 'start'){
      x$mday = 1
    }else if(rule == 'mid'){
      x$mday = 15
    }else if(rule == 'end'){
      x$mday = 1 # Go to start of month
      x$mon  = x$mon + 1 # Go to next month
      x$mday = x$mday - 1 # Subtract a day
    }
  }else if(unit == 'year'){
    if(rule == 'start'){
      x$mday = 1
      x$mon  = 0
    }else if(rule == 'mid'){
      x$mday = 2
      x$mon  = 6
    }else if(rule == 'end'){
      x$mday = 31
      x$mon  = 12
    }
  }else{warning("Inputs not valid.")}
  
  # Convert to R date
  x = as.Date(x)
  return(x)
}
