library(rstan)
library(parallel)
library(lubridate)
library(tidyverse)
library(linelistBayes)

options(mc.cores = 4)

# get onset dates
caseCounts <- create_caseCounts(sample_dates, sample_location, sample_cases)[1:80,]
ll <- convert_to_linelist(caseCounts, reportF_missP = 0.01)
ll <- data.frame(ll) %>% filter(!is.na(delay_int))

head(ll)

ll$value <- 1
ll$id <- 1:nrow(ll)
ll <- data.frame(ll)
Y <- ll$delay_int

## THIS SHOULD BE DEFINED BY REPORTING DATE
ll$start_dt <- as.vector(ll$report_date - min(ll$report_date))

xx <- ll %>%
  mutate(actual_onset = start_dt,
         actual_report = start_dt + delay_int)

rr <- xx %>% group_by(actual_report) %>% tally()
oo <- xx %>% group_by(actual_onset) %>% tally()

dt_wide <- ll %>% arrange(onset_date) %>%
  pivot_wider(id_cols = c(id, onset_date, is_weekend),
              names_from = week_int, names_prefix = 'week',
              values_fill = 0)
n_weeks = 12
head(dt_wide)

dim(dt_wide)
length(Y)

stan_data <- list(
  dum = as.matrix(dt_wide[, c(-1, -2)]),
  N = as.integer(nrow(dt_wide)),
  J = as.integer(n_weeks + 1),
  Y = as.integer(Y)
)

mod1 <- stan(file = "test.stan", data = stan_data,
             cores = 4)

out <- rstan::extract(mod1)

summary(out$phi)

dim(out$y_rep)
out$y_rep[1:5, 1:5]

out1 <- lapply(1:nrow(out$y_rep), function(i) {
  xx1 <- xx
  xx1$pred_delay = out$y_rep[i, ]
  xx1 <- xx1 %>%
    mutate(pred_onset = actual_report - pred_delay)
  xx1 %>% group_by(pred_onset) %>% tally() %>% mutate(iter = i)
})

out_df <- do.call(rbind, out1)

summaryx <- out_df %>%
  group_by(pred_onset) %>%
  summarize(
    nx = n(),
    lb = quantile(n, probs = c(0.025)),
    ex = quantile(n, probs = c(0.5)),
    ub = quantile(n, probs = c(0.975))
  )
head(summaryx)

plot(rr)
lines(oo, col='blue')
lines(x = summaryx$pred_onset, y = summaryx$lb, col='green')
lines(x = summaryx$pred_onset, y = summaryx$ub, col='green')
lines(x = summaryx$pred_onset, y = summaryx$ex, col='red')
legend("topright", 
       legend = c("Reported cases", "Observed Onset", "Empircal CI", 
                  "Predicted Onset"), 
       col = c("black", "blue", "green", "red"), 
       lty = c(NA, 1, 1, 1), # Line types
       pch = c(1, NA, NA, NA), # Point types (1 is a default point type)
       cex = 0.8) # Text size
