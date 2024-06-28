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

# which(is.na(ll$onset_date))
# sip <- si(14, 4.29, 1.18)
# out_list_demo <- run_backnow(ll, 
#                              MAX_ITER = as.integer(2000), 
#                              norm_sigma = 0.2,
#                              sip = sip,
#                              NB_maxdelay = as.integer(20),
#                              NB_size = as.integer(6),
#                              printProgress = 1)
# plot(out_list_demo, 'est')
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

########

stan_data <- list(
  ##
  J = as.integer(n_weeks + 1),
  sipN = as.integer(length(sip)),
  sip = sip,
  maxdelay = as.integer(30),
  missvector = as.integer(1*miss_rows),
  ndays = max(dt_wide$report_date_int),
  windowsize = as.integer(6),
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

## RIGHT so the reason he does that is that it ensures
## that you have estimates for every week

########

mod1 <- stan(file = "linelistBayes.stan", data = stan_data, chains = 1)

########

out <- rstan::extract(mod1)
any(is.na(out$mu_miss))

med <- apply(out$day_onset_tally, 2, quantile, probs = 0.5)
lb <- apply(out$day_onset_tally, 2, quantile, probs = 0.025)
ub <- apply(out$day_onset_tally, 2, quantile, probs = 0.975)
out_df <- data.frame(
  x = out$day_onset_tally_x[1, ],
  med, lb, ub
)

########

plot(out_list_demo, 'est')
lines(x = reference_date+out_df$x, y = out_df$med, col='blue')
lines(x = reference_date+out_df$x, y = out_df$lb, col='green')
lines(x = reference_date+out_df$x, y = out_df$ub, col='green')

legend("topright", 
       legend = c("Reported cases", "Predicted Onset_new", "Empircal CI", 
                  "Predicted Onset_old"), 
       col = c("black", "blue", "green", "red"), 
       lty = c(NA, 1, 1, 1), # Line types
       pch = c(1, NA, NA, NA), # Point types (1 is a default point type)
       cex = 0.8) # Text size



