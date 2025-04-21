# Load Rcpp package
library(Rcpp)
library(lamW)      ## for lambert W
library(VGAMdata)  ## for dbell()
library(multicool) ## for Bell()
library(glmmTMB)

# Load the functions from the C++ file
sourceCpp("bell_functions.cpp")

## not really clear to me what these are.
## do we really not need the Bell numbers in order to compute deviances?
## (maybe Bell numbers are in the normalization constant/data-dependent component?)
## dev_bell() doesn't depend on the shape parameter -- how can that be?
source("bell_funs.R")

## copied from the web (https://www.statisticshowto.com/bells-numbers-bell-triangle/)
bell_nums <- c(1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975)
## Test the bell function against these and existing implementation
as.integer(sapply(1:8, bell))  ## uh-oh
multicool::Bell(1:8)

## Test the dbell() function against the VGAMdata version
d1 <- VGAMdata::dbell(0:15, shape = 1.5)
d2 <- dbell(0:15, theta = 1.5)
all.equal(d1[1:3], d2[1:3])
## things start to go wrong at x = 4
d1[4]
d2[4]

## Test the lambertW function in R (against an independent implementation)

z <- 1.0
all.equal(lambertW(z), lamW::lambertW0(1.0)) ## TRUE


