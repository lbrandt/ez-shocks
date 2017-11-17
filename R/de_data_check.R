# Check de_data2

#rm(list=ls())

require(tidyverse)

load("de_data2.RData")


de.T = dim(data)[1]
de.N = dim(data)[2] 

qplot(1:de.T, data$BDPOPTOTP) + geom_line()


# From stats pc tut
path = file.path(getwd(), paste("de_data_plot.pdf"))
pdf(file = path)

par(mfrow = c(3,2))

for(i in 1:de.N){

  plot(data[[i]], type = "l", col = "blue", xlab = "Time", ylab = "Variable",
       main = colnames(data[i]), lwd = 2)
  #lines(n, Y2, col = "green", lwd = 2)
  #legend("topright", legend = c("T_n", "T'_n"), fill = c("blue", "green"))
  #box()
  
}
dev.off()
par(mfrow = c(1,1))








plot(0, type="n", xlab="Time", ylab="Series", xlim = c(1, de.T), ylim = c(0, 100))

plot(data$BDBQ9059A)

lines(data$BDDB_GDP , col = "blue")
lines(data$BDBQ9059A, col = "red")

