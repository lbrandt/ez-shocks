##################################################################################
#                                                                                #
#   fun.chain                                                                    #
#   Chain linking two overlapping series                                         #
#                                                                                #
#   Lennart Brandt                                                               #
#   Nov 2017                                                                     #
#                                                                                #
##################################################################################


##
# !!! ONLY BACKWARDS LINKING FROM FIRST POSSIBLE OVERLAP IMPLEMENTED AT THE MOMENT ! ! !
##


# This function chain links two overlapping time series of equal length.
# x is the "target" series, which is being extended via growth rates of y.


fun.chain = function(x, y, overlap, forwards){
  
  if(missing(forwards)){
    
    warning('If direction of link is not supplied, function defaults to backwards linking, i.e. from most recent to oldest value.')
    forwards = 0
  }
  
  if(missing(overlap)){
    
    warning('If overlap is not supplied, function defaults to overlap = 0 and will chose first possible overlap in chosen direction.')
    overlap = 0
  }
  
  if(length(x) == length(y)){
    
    series.length = length(x)
  }else{
    
    stop('Series must be of equal length. Must be padded with NAs such that they align.')
  }
  
  
  # Linking backwards
  if(forwards == 0){
    
    for(i in series.length:1){
      
      if(overlap == 0){
        
        if(is.na(x[i]) & is.numeric(y[i]) & is.numeric(y[i+1]) == TRUE){
          
          x[i] = x[i+1]* y[i]/y[i+1]
          
        }
      }
    }
  }
  
  out = x
  return(out)
}
