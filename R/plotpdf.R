# Function which plots a dataframe to pdf
plotpdf = function(data, filepath, mfrow, type){
  
  path = normalizePath(filepath, winslash = "/")
  pdf(file = path)
  
  columns = dim(data)[2]
  
  par(mfrow = mfrow)
  for(i in columns){
    
    plot(data[[i]], type = type, col = "blue", xlab = "Time", ylab = "Variable",
         main = colnames(data[i]), lwd = 2)
  }
  dev.off()
  par(mfrow = c(1,1))
  
}
