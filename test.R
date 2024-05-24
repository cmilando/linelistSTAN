library(rstan)
library(parallel)
library(lubridate)
library(tidyverse)

options(mc.cores = 4)

# get onset dates
dt1 <- make_date(2020, 1, 1)
dt2 <- make_date(2020, 3, 20)
seqdt <- seq.Date(dt1, dt2, by = 'day')

set.seed(123)
pp <- sample(c(1, 2, 3), size = length(seqdt), replace = T)

onset <- sapply(1:length(pp), function(i) rep(seqdt[i], times = pp[i]) )
onset <- do.call(c, onset)

dt <- data.frame(onset_dates = onset)

dt <- dt %>%
  mutate(week = week(onset_dates),
         id = 1:nrow(dt),
         is_weekend = 1*(wday(onset_dates) %in% c(6,7)))
n_weeks <- length(unique(dt$week))

dt_wide <- dt %>%
  mutate(value = 1) %>%
  pivot_wider(id_cols = c(id, onset_dates, is_weekend),
              names_from = week, names_prefix = 'week',
              values_fill = 0)


# simulate reporting delay
## this is slightly different
## https://stats.oarc.ucla.edu/r/dae/negative-binomial-regression/
## the r parameter theta is the inverse of the stan parameter
## 
r_true <- 1
set.seed(123)
betas <- c(rnorm(1), rnorm(n_weeks, mean = 1, sd = 1))
betas

mu <- exp(as.matrix(dt_wide[, c(-1, -2)]) %*% betas)

## right now to sample this, you need to get 
## 1000 at each timestep and take the median
N_REPS <- 10000
all_samples <- sapply(1:nrow(dt), function(i) {
  round(median(rnbinom(N_REPS, size = r_true, mu = mu[i])))
})
all_samples

stan_data <- list(
  dum = as.matrix(dt_wide[, c(-1, -2)]),
  N = as.integer(nrow(dt)),
  J = as.integer(n_weeks + 1),
  Y = as.integer(all_samples)
)

mod1 <- stan(file = "test.stan", data = stan_data,
             cores = 4)

out <- rstan::extract(mod1)

summary(out$phi)

betas
summary(out$betas)
dim(out$y_rep)
y_pred <- t(apply(out$y_rep, 2, function(j) 
  quantile(j, probs = c(0.025, 0.5, 0.975))))
y_pred <- data.frame(y_pred)
y_pred$real_y <- all_samples
y_pred$id <- 1:nrow(y_pred)
head(y_pred)

ggplot(y_pred) +
  geom_pointrange(aes(x = id, ymin = X2.5., y = X50.,
                 ymax = X97.5.), size = .1) +
  geom_point(aes(x = id, y = real_y), color = 'red')
