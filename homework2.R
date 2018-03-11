# https://www.r-project.org/conferences/useR-2009/slides/Thorin+Mallem+Noireaud+Desfontis.pdf
#https://heuristicandrew.blogspot.com/2018/01/type-i-error-rates-in-two-sample-t-test.html
# http://rstudio-pubs-static.s3.amazonaws.com/3683_9693836f3ec34563a59bb63f17756086.html
# http://support.sas.com/resources/papers/proceedings13/425-2013.pdf
vitaminD <- read.csv2("http://staff.pubhealth.ku.dk/~linearpredictors/datafiles/VitaminD.csv",
                      sep = ";",
                      dec = ".",
                      header = TRUE,
                      colClasses =
                        ## "country","category","vitd", "age","bmi","sunexp","vitdintake"
                        c("factor", "factor","numeric","numeric","numeric","factor","numeric"),
                      na.strings="."
)
library("lubridate")
library("tidyverse")
library("dplyr")
data<-as_tibble(vitaminD)
  View(data)
head(data)
## omitting individuals considered underweight
vitaminD <- vitaminD[vitaminD$bmi >= 18.5,]

## Irish women
irlwomen <- subset(vitaminD, vitaminD$country == 4 & vitaminD$category == 2)



###





source('http://192.38.117.59/~linearpredictors/datafiles/readVitaminD.R') 

## Summary statistics of Vitamin D within different weight groups

irlwomen$overweight <- (irlwomen$bmi >= 25) 

label1 <- c('normal', 'overweight')
n1 <- tapply(irlwomen$vitd, irlwomen$overweight, length)
xbar1 <- tapply(irlwomen$vitd, irlwomen$overweight, mean)
tableContent1<-cbind('Overweight group'=label1,
                     n=n1,
                     mean=round(xbar1,3))

print(tableContent1,quote=FALSE)

irlwomen$weightgroup <- (irlwomen$bmi >= 25) + (irlwomen$bmi >= 30) 

label2 <- c('normal', 'slight overweight', 'obese')
n2 <- tapply(irlwomen$vitd, irlwomen$weightgroup, length)
xbar2 <- tapply(irlwomen$vitd, irlwomen$weightgroup, mean)
tableContent2<-cbind('Overweight group'=label2,
                     n=n2,
                     mean=round(xbar2,3))

print(tableContent2,quote=FALSE)



####  


source("http://192.38.117.59/~linearpredictors/datafiles/readVitaminD.R")

## 3 groups of bmi...
irlwomen$bmigroup3 <- (irlwomen$bmi > 0) + (irlwomen$bmi > 25) + (irlwomen$bmi > 30)

## ...which we turn into a factor with midpoints as levels
irlwomen$bmigroup3 <- factor(irlwomen$bmigroup3, labels = c(21.75,27.5,32.5))

## Level of Vitamin D plotted against midpoint scores.

stripchart(irlwomen$vitd ~ irlwomen$bmigroup3,
           pch = 1,
           method = "jitter", jitter = 0.03,
           vertical = TRUE,
           ylab = "Vitamin D",
           xlab = "BMI score"
)



################


########################
# R code to conduct: Monte Carlo chi-square on 
# correlated data, 12x12 grid of locations,
# exponential covariance function.
#
# Discussion in Waller, L.A., Smith, D., Childs, J.E., and Real, L.A. (2003, to appear)
# "Monte Carlo Assessments of Goodness of Fit for Ecological Simulation Models". 
# Ecological Modelling.
########################
##  http://web1.sph.emory.edu/users/lwaller/chisqMC.r

library(MASS)
x <- rep(1:12,12)
y <- c(rep(1,12),
       rep(2,12),
       rep(3,12),
       rep(4,12),
       rep(5,12),
       rep(6,12),
       rep(7,12),
       rep(8,12),
       rep(9,12),
       rep(10,12),
       rep(11,12),
       rep(12,12))

dist <- matrix(0,length(x),length(x))

for (i in 1:length(x)){
  dist[i,] <- sqrt( (x[i]-x)^2 + (y[i]-y)^2 )
}

numsim <- 100
sig2 <- 1.0
gamma <- 0.9
Sigma <- sig2*exp(-gamma*dist)
mu <- rep(1,length(x))

simdata <- mvrnorm(numsim,mu,Sigma)
chisqMC <- rep(0,length(numsim))

for (i in 1:numsim) {
  O <- simdata[i,]
  E <- apply(simdata[-i,],2,mean)
  V <- apply(simdata[-i,],2,var)
  chisqMC[i] <- sum( (O-E)^2/V )
}

hist(chisqMC,probability=T)
xvals <- 0:200
lines( xvals, dchisq(xvals,df=144) )


############## T-test  for type one error

# https://stats.stackexchange.com/questions/148526/how-to-simulate-type-i-error-and-type-ii-error

#To observe the type I error of a test we need to generate/simulate data from the same
#distribution that follows null hypothesis.
n=10000 # testing 10,000 times
t1err=0
for (i in 1:n){
  #x=rnorm(100, 0, 1)
  #http://rstudio-pubs-static.s3.amazonaws.com/633_048ca606ec374db2b7edc3a6fd397cef.html
  x <-  arima.sim(n = 100, list(ar = 0.5), innov=rnorm(100,0,25))
  if (((t.test(x, mu=0))$p.value)<=0.05) (t1err=t1err+1) 
}

cat("Type I error rate in percentage is", (t1err/n)*100,"%")

# To test Type II error we have to generate/simulate data from another
# distribution than that followed by null hypothesis.

n=10000 # testing 10,000 times
t2err=0
for (i in 1:n){
  x=rnorm(100, 2, 0)
  if (((t.test(x, mu=0))$p.value)>0.05) (t2err=t2err+1) 
}
cat("Type II error rate in percentage is", (t2err/n)*100,"%")

# You will see 0.0%. As the variance is really low. 
# If you increase variance to 5, you will see about 2% error as Type II error.


## Ques No: 8
sigma[i,i]<-25
cov[i,j]<-rho*Sigma



#
rho <- 0.5 
sigma <- 25 
x <- diag(5) 
x <- sigma * rho^abs(row(x)-col(x))


# create a variance covariance matrix



x<-function(i,var,rho){
  sigma<-var
  corr<-rho
  x<-diag(i)
  x<-sigma * corr^abs(row(x)-col(x))
  return(x)
  }

yy<-x(50,25,-.5)

library(MASS)

#If n = 1 a vector of the same length as mu, 
# otherwise an n by length(mu) matrix with one sample in each row
#The matrix decomposition is done via eigen; although a Choleski 
# decomposition might be faster, the eigendecomposition is stabler
# https://heuristicandrew.blogspot.com/2018/01/type-i-error-rates-in-two-sample-t-test.html
y<- mvrnorm(n=1,mu=rep(0,50),yy)

# type I error rate analysis

n=10000 # testing 10,000 times
t1err=0
for (i in 1:n){
  #x=rnorm(100, 0, 1)
  x<- mvrnorm(n=1,mu=rep(0,50),yy)
  if (((t.test(x, mu=0))$p.value)<=0.05) (t1err=t1err+1) 
}

cat("Type I error rate in percentage is", (t1err/n)*100,"%")


library(dplyr)
setwd("E:/KUMC SPRING/SPRING 2018/Design of Experiment/Lecture/HOMEWORK")
hyperten<-read_csv("Hypertension.csv")

head(hyperten)

x<-ifelse(hyperten$Treatment,"A","B")






