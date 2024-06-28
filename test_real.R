library(rstan)
library(parallel)
library(lubridate)
library(tidyverse)
library(linelistBayes)

options(mc.cores = 4)

# get onset dates

#############
# caseCounts <- create_caseCounts(sample_dates, sample_location, sample_cases)[1:80,]
# ll <- convert_to_linelist(caseCounts, reportF_missP = 0.75)
# which(is.na(ll$onset_date))
# sip <- si(14, 4.29, 1.18)
# out_list_demo <- run_backnow(caseCounts, 
#                              MAX_ITER = as.integer(2000), 
#                              norm_sigma = 0.2,
#                              sip = sip,
#                              NB_maxdelay = as.integer(20),
#                              NB_size = as.integer(6),
#                              printProgress = 1,
#                              reportF_missP = 0.6)
# plot(out_list_demo, 'est')
#############
# how are the confidence intervals filled in here
# shouldn't it get wider if there are fewer data points?
# 

#############
caseCounts <- create_caseCounts(sample_dates, sample_location, sample_cases)[1:80,]
set.seed(123)
ll <- convert_to_linelist(caseCounts, reportF_missP = 0.5)

which(is.na(ll$onset_date))
sip <- si(14, 4.29, 1.18)
out_list_demo <- run_backnow(caseCounts, 
                             MAX_ITER = as.integer(2000), 
                             norm_sigma = 0.2,
                             sip = sip,
                             NB_maxdelay = as.integer(20),
                             NB_size = as.integer(6),
                             printProgress = 1,
                             reportF_missP = 0.6)
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
head(ll)
dim(ll)

rr <- ll %>% group_by(report_date_int) %>% tally()
# oo <- ll %>% group_by(onset_date_int) %>% tally()

dt_wide <- ll %>% arrange(onset_date) %>%
  pivot_wider(id_cols = c(id, onset_date_int, report_date_int, is_weekend),
              names_from = week_int, names_prefix = 'week',
              values_fill = 0) %>% arrange(id)
head(dt_wide)
head(ll)

n_weeks = max(ll$week_int)

x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
onset <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
maxdelay <- 10
cd <- 1
plot(prop(x, onset, maxdelay, cd))

## R:: gives a scalar, so you can assume that n = 1
## R::rnbinom(    size, prob )

rnbinom()

##
miss_rows <- is.na(dt_wide$onset_date_int)
miss_rows2 <- is.na(Y)
identical(miss_rows, miss_rows2)

# //
#   int<lower=1>     J;         // Number of Betas: aka n weeks + 1 is_weekend
# int<lower=1>     sipN;      // SIP length
# vector[sipN]     sip;       // SIP
# 
# // OBSERVED
# int<lower=1>     N_obs;           // Number of individuals
# matrix[N_obs, J] dum_obs;         // matrix of indicator values
# int              Y_obs[N_obs];    // observed reporting delays
# int              ReportOnset[N_obs]; // which days did these occur on ...
# 
# // MISSING
# int<lower=1>      N_miss;        // Number of individuals with missing data
# matrix[N_miss, J] dum_miss;      // matrix of indicator values
# int               ReportDays[N_miss]; // which days were things reported on

stan_data <- list(
  ##
  J = as.integer(n_weeks + 1),
  sipN = as.integer(length(sip)),
  sip = sip,
  maxdelay = as.integer(20),
  missvector = as.integer(1*miss_rows),
  ndays = max(dt_wide$report_date_int),
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

mod1 <- stan(file = "test.stan", data = stan_data, chains = 1)

out <- rstan::extract(mod1)

########

summary(out$phi)

dim(out$guessOnset) # 4000 iterations, 1025 missing
out$guessOnset[1:5, 1:5]

dim(out$day_onset_tally)

plot(out_list_demo, 'est')

med <- apply(out$day_onset_tally, 2, quantile, probs = 0.5)
lb <- apply(out$day_onset_tally, 2, quantile, probs = 0.025)
ub <- apply(out$day_onset_tally, 2, quantile, probs = 0.975)
out_df <- data.frame(
  x = out$day_onset_tally_x[1, ],
  med, lb, ub
)
head(out_df)

out$day_onset_tally_x[1, ]

# ok so get the sums that result from this
# library(pbapply)
# out1 <- pblapply(1:nrow(out$guessOnset), function(i) {
#   xx1 <- dt_wide[, c('id', 'onset_date_int')]
#   miss3 <- is.na(xx1$onset_date_int)
#   identical(miss_rows, miss3)
#   xx1[miss3, 'onset_date_int'] <- out$guessOnset[i, ]
#   xx1 %>% group_by(onset_date_int) %>% tally() %>% mutate(iter = i)
# })

# out1 <- pblapply(1:nrow(out$guessOnset), function(i) {
#   xx1 <- dt_wide[, c('id', 'onset_date_int')]
#   miss3 <- is.na(xx1$onset_date_int)
#   identical(miss_rows, miss3)
#   xx1[miss3, 'onset_date_int'] <- out$guessOnset[i, ]
#   xx1 %>% group_by(onset_date_int) %>% tally() %>% mutate(iter = i)
# })
# 
# out_df <- do.call(rbind, out1)
# max(out_df$iter)
# head(out_df)
# 
# summaryx <- out$day_onset_tally %>%
#   group_by(onset_date_int) %>%
#   summarize(
#     nx = n(),
#     lb = quantile(n, probs = c(0.025)),
#     ex = quantile(n, probs = c(0.5)),
#     ub = quantile(n, probs = c(0.975))
#   )
# head(summaryx)
# which(summaryx$lb != summaryx$ex)

plot(out_list_demo, 'est')
# abline(v = as.Date('2020-01-08'))
lines(x = as.Date('2020-01-01')+out_df$x, y = out_df$med, col='blue')
lines(x = as.Date('2020-01-01')+out_df$x, y = out_df$lb, col='green')
lines(x = as.Date('2020-01-01')+out_df$x, y = out_df$ub, col='green')

legend("topright", 
       legend = c("Reported cases", "Predicted Onset_new", "Empircal CI", 
                  "Predicted Onset_old"), 
       col = c("black", "blue", "green", "red"), 
       lty = c(NA, 1, 1, 1), # Line types
       pch = c(1, NA, NA, NA), # Point types (1 is a default point type)
       cex = 0.8) # Text size

# head(out_list_demo$report_cases)
# head(rr, 10)
# 
# plot(rr)
# head(rr)
# #plot(out_list_demo$est_back_date, out_list_demo$est_back)
# #lines(oo, col='blue')
# dt_seq <- seq.Date(from = as.Date('2020-01-01'))
# lines(x = summaryx$onset_date_int, y = summaryx$lb, col='green')
# lines(x = summaryx$onset_date_int, y = summaryx$ub, col='green')
# lines(x = summaryx$onset_date_int, y = summaryx$ex, col='red')
# legend("topright", 
#        legend = c("Reported cases", "Observed Onset", "Empircal CI", 
#                   "Predicted Onset"), 
#        col = c("black", "blue", "green", "red"), 
#        lty = c(NA, 1, 1, 1), # Line types
#        pch = c(1, NA, NA, NA), # Point types (1 is a default point type)
#        cex = 0.8) # Text size

