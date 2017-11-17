# Check de_data2

#rm(list=ls())


load("de_data2.RData")


de.T = dim(data)[1]
de.N = dim(data)[2] 




plot(0, type="n", xlab="Time", ylab="Series", xlim = c(1, de.T), ylim = c(0, 100))

plot(data$BDBQ9059A)

lines(data$BDDB_GDP , col = "blue")
lines(data$BDBQ9059A, col = "red")

lines(diff(log(data$BDUSMB36B)))
lines(diff(log(data$BDUSMC36B)))

lines(pw_data_final$X__45)


test3 = diff(log(data$BDUSMB36B))/diff(log(data$BDUSMC36B))