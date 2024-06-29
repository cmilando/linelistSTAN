library(rstan)
library(parallel)
library(lubridate)
library(tidyverse)
library(linelistBayes)

options(mc.cores = 4)

#############
caseCounts <- create_caseCounts(sample_dates, sample_location, sample_cases)[1:80,]
set.seed(123)
ll <- convert_to_linelist(caseCounts, reportF_missP = 0.60)

which(is.na(ll$onset_date))
sip <- si(14, 4.29, 1.18)
out_list_demo <- run_backnow(ll, 
                             MAX_ITER = as.integer(2000), 
                             norm_sigma = 0.2,
                             sip = sip,
                             NB_maxdelay = as.integer(20),
                             NB_size = as.integer(6),
                             printProgress = 1)
plot(out_list_demo, 'est')
#############

# ll <- data.frame(ll) %>% filter(!is.na(delay_int))

head(ll)

ll$value <- 1
ll$id <- 1:nrow(ll)
ll <- data.frame(ll)
Y <- ll$delay_int

## THIS SHOULD BE DEFINED BY REPORTING DATE
reference_date = as.Date("2020-01-01")
ll$report_date_int <- as.vector(ll$report_date - reference_date)
ll$onset_date_int <- as.vector(ll$onset_date - reference_date)

###
rr <- ll %>% group_by(report_date_int) %>% tally()

dt_wide <- ll %>% arrange(onset_date) %>%
  pivot_wider(id_cols = c(id, onset_date_int, report_date_int, is_weekend),
              names_from = week_int, names_prefix = 'week',
              values_fill = 0) %>% arrange(id)

n_weeks = max(ll$week_int)

##
miss_rows <- is.na(dt_wide$onset_date_int)
miss_rows2 <- is.na(Y)
identical(miss_rows, miss_rows2)

## CHECK THAT YOU HAVE AT LEAST ONE PERSON PER WEEK
## This is the reason he does that is that it ensures
## that you have estimates for every week
ll %>%
  group_by(week_int) %>%
  summarize(
    .groups = 'keep',
    n = n(),
    n_not_na = sum(!is.na(onset_date) * 1)
  )



########

stan_data <- list(
  ##
  J = as.integer(n_weeks + 1),
  sipN = as.integer(length(sip)),
  sip = sip,
  maxdelay = as.integer(20),
  missvector = as.integer(1*miss_rows),
  ndays = max(dt_wide$report_date_int),
  windowsize = as.integer(7),
  ##
  N_obs = as.integer(nrow(dt_wide[!miss_rows, ])),
  dum_obs = as.matrix(dt_wide[!miss_rows, -c(1:3)]),
  Y_obs = as.integer(Y[!miss_rows]), ## DELAY
  ReportOnset = as.integer(unlist(dt_wide[!miss_rows, 2])), ## ONSET
  ##
  N_miss = as.integer(nrow(dt_wide[miss_rows, ])),
  dum_miss = as.matrix(dt_wide[miss_rows, -c(1:3)]),
  ReportDays = as.integer(unlist(dt_wide[miss_rows, 3]))
)

########

mod1 <- stan(file = "linelistBayes.stan", data = stan_data, chains = 1)

########

out <- rstan::extract(mod1)

# any(is.na(out$mu_miss))

# any(out$day_onset_tally_tail < 1)
# tt <- as.matrix(dt_wide[miss_rows, -c(1:3)]) %*% out$betas[6,]
# dim(tt)
# head(tt)

########

out_df <- data.frame(
  x = reference_date + out$day_onset_tally_x[1, ],
  med = apply(out$day_onset_tally, 2, quantile, probs = 0.5), 
  lb = apply(out$day_onset_tally, 2, quantile, probs = 0.025), 
  ub = apply(out$day_onset_tally, 2, quantile, probs = 0.975)
)

plot(out_list_demo, 'est')
lines(x = out_df$x, y = out_df$med, col='blue')
lines(x = out_df$x, y = out_df$lb, col='green')
lines(x = out_df$x, y = out_df$ub, col='green')

legend("topright", 
       legend = c("Reported cases", "Predicted Onset_new", "Empircal CI", 
                  "Predicted Onset_old"), 
       col = c("black", "blue", "green", "red"), 
       lty = c(NA, 1, 1, 1), # Line types
       pch = c(1, NA, NA, NA), # Point types (1 is a default point type)
       cex = 0.8) # Text size


########
length(out$day_onset_tally_x[1, ]) # ndays + maxdelay
length(out$rt[1,])
# ndays + maxdelay - windowsize - 1

rt_df <- data.frame(
  x = reference_date + out$day_onset_tally_x[1, ],
  med = apply(out$rt, 2, quantile, probs = 0.5), 
  lb = apply(out$rt, 2, quantile, probs = 0.025), 
  ub = apply(out$rt, 2, quantile, probs = 0.975)
) 
head(rt_df)

plot(out_list_demo, 'rt')
lines(x = rt_df$x, y = rt_df$med, col='blue')
lines(x = rt_df$x, y = rt_df$lb, col='green')
lines(x = rt_df$x, y = rt_df$ub, col='green')

# legend("topright", 
#        legend = c("Reported cases", "Predicted Onset_new", "Empircal CI", 
#                   "Predicted Onset_old"), 
#        col = c("black", "blue", "green", "red"), 
#        lty = c(NA, 1, 1, 1), # Line types
#        pch = c(1, NA, NA, NA), # Point types (1 is a default point type)
#        cex = 0.8) # Text size

########
# install.packages("EpiEstim")
library(EpiEstim)
library(pbapply)
dim(out_df)

windowsize = 6
t_start <- seq(2, length(out$day_onset_tally_x[1,])-(windowsize - 1))
t_end <- t_start + (windowsize  - 1)

allRts <- pblapply(1:nrow(out$day_onset_tally), function(i) {
  epiR <- estimate_R(data.frame(dates = out$day_onset_tally_x[1, ], 
                                I = out$day_onset_tally[i, ]), 
             method = 'non',
             config = make_config(list(si_distr = c(0, sip),
                                       t_start = t_start,
                                       t_end = t_end)))
  unlist(epiR$R$`Median(R)`)
})
allRts_mat <- as.matrix(do.call(rbind, allRts))

Epirt_df <- data.frame(
  x = reference_date - 20 + 1 + t_end,
  med = apply(allRts_mat, 2, quantile, probs = 0.5), 
  lb = apply(allRts_mat, 2, quantile, probs = 0.025), 
  ub = apply(allRts_mat, 2, quantile, probs = 0.975)
) 

head(Epirt_df)
plot(out_list_demo, 'rt')
lines(Epirt_df$x, Epirt_df$med)
lines(Epirt_df$x, Epirt_df$lb)
lines(Epirt_df$x, Epirt_df$ub)

legend("topright",
       legend = c("EpiEstim", "Chad's STAN", "backnow.cpp"),
       col = c("black", "blue", "red"),
       lty = c(1, 1, 1), # Line types
       pch = c(NA, NA, NA), # Point types (1 is a default point type)
       cex = 0.8) # Text size

########
## try and calculate it with matrices

ndays = max(dt_wide$report_date_int)
maxdelay = as.integer(20)
tau = 6


boot_i = 1
sipN = 14
library(pbapply)
# maybe I'm doing lambda incorrectly
boot_rt <- pblapply(1:nrow(out$day_onset_tally), function(boot_i) {
  
  rt = vector("numeric", ndays + maxdelay)
  
  for(t in 1:(ndays + maxdelay)) {
    
    k_init = t - tau
    
    if(k_init < 1) {
      
      rt[t] = NA
      
    } else {
      
      numerator = 0
      denominator = 0
      
      for(k in k_init:t) {
      
        ## numerator
        numerator = numerator + out$day_onset_tally[boot_i, k]
        
        ## denominator
        lambda = 0
        
        # maybe this is supposed to be k + 1?
        for(j in 1:min(sipN, k)) {
          # if you don't do k - j + 1, then you get a 0 index when j = k
          lambda = lambda + out$day_onset_tally[boot_i, k - j + 1] * sip[j] 
          #print(paste("lambda:", lambda))
        }
        #print(paste("lambda:", lambda))
        denominator = denominator + lambda
      
      }
      
      rt[t] = (numerator + 1) / (denominator + 0.2)
  
    }
    
  }
  
  rt
})
  
rt_mat <- do.call(rbind, boot_rt)

Calc2epirt_df <- data.frame(
  #x = reference_date - 20 + 1 + t_end,
  x = reference_date + out$day_onset_tally_x[1,],
  med = apply(rt_mat, 2, quantile, probs = 0.5, na.rm = T), 
  lb = apply(rt_mat, 2, quantile, probs = 0.025, na.rm = T), 
  ub = apply(rt_mat, 2, quantile, probs = 0.975, na.rm = T)
) 

dim(Calc2epirt_df)

head(Epirt_df)
plot(out_list_demo, 'rt')
lines(Calc2epirt_df$x, Calc2epirt_df$med)

lines(Epirt_df$x, Epirt_df$lb)
lines(Epirt_df$x, Epirt_df$ub)

