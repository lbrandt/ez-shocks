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


de.T = dim(dlndata)[1]
de.N = dim(dlndata)[2] 


# Check if rows/columns contain observations that fulfil some condition
test1 = is.na(dlndata)
rowSums(test1) # Count cases


#### Compare PW and GS
# Which series are in both datesets?
var1 = colnames(data)
var2 = colnames(data2)

test2 = cbind(c(var1, var2))

test1 = matrix(NA,100,105)

for(i in 1:100){
  
  for(j in 1:105){
    
    test1[i,j] = var1[i] == var2[j]
  }
}

sum(test1)

tindex = which(test1, arr.ind = TRUE)

# There are 8 series which are in both data sets. Merging should result in a 323x197 df
data3 = merge(data, data2)


# Restrict sample to 1991Q1-2016Q4
x = slice(data2, 1:(de.T-5))





# Plot raw series
path = file.path(getwd(), paste("de_gsdata_plot.pdf"))
pdf(file = path)

par(mfrow = c(3,2))

for(i in 1:de.N){

  plot(data2[[i]], type = "l", col = "blue", xlab = "Time", ylab = "Variable",
       main = colnames(data2[i]), lwd = 2)
  #lines(n, Y2, col = "green", lwd = 2)
  #legend("topright", legend = c("T_n", "T'_n"), fill = c("blue", "green"))
  #box()
  
}
dev.off()
par(mfrow = c(1,1))


# Plot logdiffs
path = file.path(getwd(), paste("de_dlndata_plot.pdf"))
pdf(file = path)

par(mfrow = c(3,2))

for(i in 1:de.N){
  
  plot(dlndata[[i]], type = "l", col = "blue", xlab = "Time", ylab = "Variable",
       main = colnames(dlndata[i]), lwd = 2)
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

