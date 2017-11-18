# Check de_data2

#rm(list=ls())

require(tidyverse)

# Load cleaned German data
load("de_data2.RData")
load("pw_data_final.RData")

# Read Pirschel-Wolters data for comparison and save as RData

#pw_data_final <- read_excel("C:/Dateien/My Dropbox/Lennart/Thesis/Data/pw_data_final.xls", 
#                            sheet = "datafinal", range = "C90:EL189", col_names = FALSE, na = "NA")
#
#save(pw_data_final, file = "pw_data_final.RData")


de.T = dim(data)[1]
de.N = dim(data)[2] 


# Check if columns contain observations that fulfil some condition
test1 = is.na(x)
test2 = colSums(test1) # Count cases
test3 = test2[test2>0] # Extract cases

# Restrict sample to 1991Q1-2016Q4
x = slice(data, 1:(de.T-4))





# Plot raw series
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


# Plot logdiffs
path = file.path(getwd(), paste("de_data2_plot.pdf"))
pdf(file = path)

par(mfrow = c(3,2))

for(i in 1:de.N){
  
  plot(dlndata[[i]], type = "l", col = "blue", xlab = "Time", ylab = "Variable",
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

