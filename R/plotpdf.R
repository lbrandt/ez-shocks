# Function which plots a dataframe to pdf
plotpdf = function(data, rows, columns, filepath, mfrow, type){
  
  path = normalizePath(filepath, winslash = "/")
  pdf(file = path)
  
  par(mfrow = mfrow)
  for(i in columns){
    
    plot(data[[i]], type = type, col = "blue", xlab = "Time", ylab = "Variable",
         main = colnames(data[i]), lwd = 2)
  }
  dev.off()
  par(mfrow = c(1,1))
  
}

ez_data[2]
ez_data[[3:10, 2]]

plot(ez_data[[2]])

plot(ez_data[[3:10, 2])
