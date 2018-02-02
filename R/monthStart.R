# Sets day of date vector to 1, like Matlab dateshift(start of month)
monthStart <- function(x) {
  x = as.POSIXlt(x)
  x$mday = 1
  as.Date(x)
}
