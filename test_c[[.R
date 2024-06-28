library(RcppArmadillo)
library(Rcpp)



############
#rpois <- Vectorize(rpois,"lambda")
rt <- readRDS('rt.rds') ## Rt input
## Set up the serial interval
si <- function(ndays,alpha,beta){
  prob <- numeric(ndays)
  for (i in 1:ndays){
    prob[i] <- pgamma(i,shape = alpha,rate = beta) - pgamma(i-1,shape = alpha,rate = beta)
  }
  result <- prob/sum(prob)
  return(result)
}

sip <- si(14,4.29,1.18)
#
# simdata <- function(n0,r0,rt,r,m,maxdelay,ld){
n0 = 100
r0 = 1.8
rt = rt
r = 3
m = 9
maxdelay = 20
ld = 60

  n1=rpois(1,n0*r0)
  inf1=rgamma(n1,shape = 4.29,rate = 1.18) 
  i1=ceiling(inf1) ## Date of infection
  inc1=exp(rnorm(n1,1.621,0.418)) 
  o1=i1+ceiling(inc1) ## Date of onset sympotom
  p=dnbinom(0:maxdelay,size = r, mu = m)
  d1=o1+sample(0:maxdelay,size = n1,prob = p,replace = TRUE) ## Date of report
  dat=cbind(i1,o1,d1)
  while(any(i1<=ld)){
    date=as.numeric(names(table(i1))) ## Infection dates of previous generation
    no=as.numeric(table(i1)) ## Number of cases on each date
    mu=no*rt[date]
    n1=rpois(1,mu)
    inf1=rgamma(sum(n1),shape = 4.29,rate = 1.18)
    i1=rep(date,times=n1)+ceiling(inf1)
    inc1=exp(rnorm(sum(n1),1.621,0.418)) 
    o1=i1+ceiling(inc1) ## Date of onset sympotom
    d1=o1+sample(0:maxdelay,size = sum(n1),prob = p,replace = TRUE) ## Date of report
    dat1=cbind(i1,o1,d1)
    dat=rbind(dat,dat1)
  }
  ## 3 outputs
  ## 1-Report Curve
  ## 2-Epidemic Curve
  ## 3-Data used for analysis
  dat=data.frame(dat)
  colnames(dat)=c("infect","onset","report")
  epic=table(dat[dat$onset<=ld,]$onset)
  dat <- dat[dat$report<=ld,c(2,3)]
  dat[sample(nrow(dat),round(nrow(dat)*0.6)),1] <- NA
  repc=table(dat$report) ## Report curve
  out = list(repc=repc,epic=epic,dat=dat)
#}

library(tidyverse)
#out <- simdata(100,1.8,rt,3,9,maxdelay,60)
d <- out$dat
d <- d%>%arrange(report)
d <- d%>%mutate(weekend=ifelse(report%%7==0|report%%7==6,1,0))
d <- d%>%mutate(delay=report-onset,minday=min(report))
d <- d%>%mutate(report=report-minday+1,week=ceiling(report/7))
epic <- out$epic
repc <- out$repc
minday <- unique(d$minday)

#source("backnow.cpp")
sourceCpp("backnow.cpp")

out1 <- backnow(outcome=d$delay,
                days=d$report,
                week=d$week,
                weekend=d$weekend,
                iter=2,
                sigma=0.2,
                maxdelay=20,si=sip,size=6)

# NumericVector prop(NumericVector x, NumericVector onset, int maxdelay, int cd) {
#   
#   LogicalVector v = (x <= maxdelay) & (onset >= cd);
#   NumericVector x1 = x[v];
#   int dem = x1.size();
#   NumericVector p1 (maxdelay);
#   
# 
#     for (int i = 0; i < maxdelay; ++i){
#       p1[i] = sum(x1 == (maxdelay - i));
#     }
#   
#   NumericVector p = p1 / dem;
#   NumericVector p2 = cumsum(p);
#   NumericVector result = 1 - p2;
#   
#   return result;
# }

prop_r <- function(x, onset, maxdelay, cday) {
  v = (x <= maxdelay) & (onset >= cd)
  x1 = x[v]
  dem = length(x1)
  p1 <- vector("numeric", maxdelay)
  for(i in 0:(maxdelay-1)) {
    p1[i] = sum((x1 == (maxdelay - i))*1)
  }
  p = p1 / dem
  p2 = cumsum(p)
  1 - p2
}

prop_r2 <- function(x, onset, maxdelay, cd) {
  # Logical vector to select valid entries
  v <- (x <= maxdelay) & (onset >= cd)
  
  # Filter x based on logical vector
  x1 <- x[v]
  
  # Number of valid entries
  dem <- length(x1)
  
  # Initialize p1
  p1 <- numeric(maxdelay)
  
  # Populate p1
  for (i in 1:maxdelay) {
    p1[i] <- sum(x1 == (maxdelay - i + 1))
  }
  
  # Calculate proportion
  p <- p1 / dem
  
  # Calculate cumulative sum and final result
  p2 <- cumsum(p)
  result <- 1 - p2
  
  return(result)
}


# so weight happens for the last 20

# ok so outcome1 = reporting delay
# backd = d$report - out1$outcome1 so this is onset time

weights = prop(out1$outcome1, out1$backd, out1$maxdelay, out1$cday)
weights

prop_r2(out1$outcome1, out1$backd, out1$maxdelay, out1$cday)

# mapply
sapply(1:20, function(i) rnb(out1$back2[i], weights[i]))

#Eh, prop can be done outside of stan right?

##
# (1) get weights

# (2) get the tail of backd (> back1)

# (3) if any tail values are 0, set to 1

# (4) do mapply for rnb between weights and existing values

# (5) add to tail, subtract 1, and set floor to 0

# (6) overwride existing tail

# other things that tenglong did
# - to fit the log-likelihood in the first place, he replaces any missing value
#   with a random draw from sample from 1:maxdelay. 
#   > this would bias the distribution towards the tail 
# - he also calculates a standard error for everyone, not just the missing value
#  > this overestimates the SE, because you have true data for some people
#  so the true confidence intervals will be smaller



##
# NumericVector weights = prop(outcome1,backd,maxdelay,cday);

# back1 = backc[seq(nd,nd+maxdelay-1)];
# NumericVector check00 (back1.size());
# check0 = check00;
# 
# LogicalVector ll = (back1==0);
# l = ll;
# check0[l] = 1;
# 
# back2 = back1 + check0;  
# trunc = mapply(back2,weights,rnb);
# now = back1+trunc;
# now1 = now - check0;
# check = (now1<0);
# now1[check] = 0;
# backc[seq(nd,nd+maxdelay-1)] = now1;
