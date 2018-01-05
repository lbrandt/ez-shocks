# Tut4 Histogram


x = rnorm(15,2,1)

hist(x)
hist(x,probability = TRUE) # histogram with scaled y-axis

H = hist(x,probability = TRUE)



################
# 1. Histogram function

set.seed(1)
x = rnorm(10,3,1)

# 1.a)
tut.hist.1 = function(x,a,b,jnum){
  
  # Calculate bin width h
  h = (b-a)/jnum
  n = length(x)
  f = rep(0,times=n)
  
  for(j in 1:jnum){
    
    # Interval borders
    int.l = a + (j-1)*h
    int.u = a + j*h
    
    # Count obs in bin
    count = sum((x>=int.l)&(x<int.u))
    
    # Assign values to density estimator f
    for(i in 1:n){
      if((x[i]>=int.l)&(x[i]<int.u)){f[i] = count/(n*h)}
    }
    
    
  }
  return(f)
}


# Check
hist1 = tut.hist.1(x,min(x),max(x),5)

plot(x,hist1)


# 1.b)
tut.hist.2 = function(x,jnum,h){
  
  n = length(x)
  f = rep(0,times=n)
  
  for(i in 1:n){
    f[i] = 1/(n*h)*sum(abs((x-x[i])/h) <= 0.5)
  }
  return(f)
}

# Check
hist2 = tut.hist.2(x,3,0.6)


################
# 2.a) Exponential function

y = rexp(n=10,rate=2)

hist(y,5,0.5)

plot(y,tut.hist.2(y,5,0.5))
lines(density(y))

# Calculate bias and variance
f.y = tut.hist.2(y,5,0.5)
f.true = dexp(y,rate=2)

# bias
mean(f.y - f.true)

# variance
mean((f.y - f.true)^2)


# Monte-Carlo for a grid of n
MC = 100
n = 5:100
bins = 3
h = 0.5

var.est = c()
bias.est = c()

for(i in 1:length(n)){
  
  var.mc = c()
  bias.mc = c()
  
  for(imc in 1:MC){
    
    y = rexp(n[i],rate=2)
    f.y = tut.hist.2(y,bins,h)
    f.true = dexp(y,rate=2)
    
    var.mc[imc] = mean((f.y - f.true)^2)
    bias.mc[imc] = mean(f.y - f.true)
    
  }
  var.est[i] = mean(var.mc)
  bias.est[i] = mean(bias.mc)
}

plot(n,var.est,type="l")
plot(n,bias.est,type="l")


# Monte Carlo for a grid of h
MC = 100
n = 10
bins = 4
h = seq(from=0.01,to=1,length.out = 100)

var.est = c()
bias.est = c()

for(i in 1:length(h)){
  
  var.mc = c()
  bias.mc = c()
  
  for(imc in 1:MC){
    
    y = rexp(h[i],rate=2)
    f.y = tut.hist.2(y,bins,h[i])
    f.true = dexp(y,rate=2)
    
    var.mc[imc] = mean((f.y - f.true)^2)
    bias.mc[imc] = mean(f.y - f.true)
    
  }
  var.est[i] = mean(var.mc) # Does not compute!!
  bias.est[i] = mean(bias.mc)
}

plot(h,var.est,type="l")
plot(h,bias.est,type="l")


