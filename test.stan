#include prop.stan
  
data {
  //
  int<lower=1>     J;         // Number of Betas: aka n weeks + 1 is_weekend
  int<lower=1>     sipN;      // SIP length
  vector[sipN]     sip;       // SIP
  int<lower=1>  maxdelay;
  int<lower=1>  ndays;
  
  // OBSERVED
  int<lower=1>     N_obs;           // Number of individuals
  matrix[N_obs, J] dum_obs;         // matrix of indicator values
  int              Y_obs[N_obs];    // observed reporting delays
  int              ReportOnset[N_obs]; // which days did these occur on ...
  
  // MISSING
  int<lower=1>      N_miss;        // Number of individuals with missing data
  matrix[N_miss, J] dum_miss;      // matrix of indicator values
  int               ReportDays[N_miss]; // which days were things reported on
  
  // ADD A FUNCTION FOR CALCULATING R(T) BY SIP AND THE OTHERS
  int  missvector[N_obs + N_miss];
}


parameters {
  vector[J] betas;         // one param for each week + 1 is_weekend
  real<lower=0.01> phi;    // a single dispersion param
}


transformed parameters {
  vector[N_obs] mu_obs;          // each person has their own mu
  mu_obs = exp(dum_obs * betas); // dot-product, gives mu_vector
}

model {
  
  // prior for beta and size
  betas ~ normal(0, 1);
  phi ~ normal(0, 1);
  
  // likelihood
  // -- UPDATE THIS TO BE JUST THE DATA YOU HAVE
  Y_obs ~ neg_binomial_2(mu_obs, phi);
  
}

generated quantities {

 // -- THE DATA WITH EXISTING ONSET INFORMATION
 // -- I guess in theory you don't need to do this
 // -- because you know when the onset is happening
 // vector[N_obs] y_rep_obs;
 // for(n in 1:N_obs) {
 //    y_rep_obs[n] = neg_binomial_2_rng(mu_obs[n], phi); 
 //}

 // -- THE DATA THAT WERE MISSING ONSET INFORMATION
  vector[N_miss] mu_miss;          // each person has their own mu
  vector[N_miss] y_rep_miss;
  vector[N_miss] guessOnset;
  mu_miss = exp(dum_miss * betas); // dot-product, gives mu_vector
  for(n in 1:N_miss) {
    y_rep_miss[n] = neg_binomial_2_rng(mu_miss[n], phi); 
    guessOnset[n] = ReportDays[n] - y_rep_miss[n];
  }
  
 // -- CALCULATE RT --> could happen outside of this as well right?
 // or ... you could pass in the days vector as well
 // get summed counts by day
 // then it should be fine ...
 vector[N_obs + N_miss] allOnset;
 vector[N_obs + N_miss] allY;
 int i_miss = 1;
 int i_true = 1;
 for(n in 1:(N_obs + N_miss)) {
   if(missvector[n] == 0) {
     allOnset[n] = ReportOnset[i_true];
     allY[n] = Y_obs[i_true];
     i_true += 1;
   } else {
    allOnset[n] = guessOnset[i_miss];
    allY[n] = ReportDays[i_miss] - guessOnset[i_miss];
    i_miss += 1;
   }
 }

 // so now you need to summarize the values in allOnset
 vector[maxdelay] weights;
 weights = prop(allY, allOnset, maxdelay, -maxdelay);
 
 // ok now create the full summed out-matrix
 vector[ndays + maxdelay] day_onset_tally;
 vector[ndays + maxdelay] day_onset_tally_x;
 for (j in 1:(ndays + maxdelay)){
    int local_sum = 0;
    for(n in 1:(N_obs + N_miss)) {
      if(allOnset[n] == (j - maxdelay + 1)) {
        local_sum += 1;
      }
    }
    day_onset_tally[j] = local_sum;
    day_onset_tally_x[j] = j - maxdelay + 1;
 }
 
 // Now do now-casting on the daily tallys of the tail
 vector[maxdelay] day_onset_tally_tail;
 day_onset_tally_tail = day_onset_tally[(ndays + 1):(ndays + maxdelay)];
 vector[maxdelay] check;
 for(i in 1:maxdelay) {
   check[i] = 0;
   if(day_onset_tally_tail[i] == 0) {
     check[i] = 1;
     day_onset_tally_tail[i] = 1;
   }
 }
 
 // here's where you calculate what to add
 // rnb(out1$back2[i], weights[i]))
 //     size, prob
 // so the size is back2, which is the onset
 // and the prob is from weights
 vector[maxdelay] trunc;
 
 for(i in 1:maxdelay) {
   real mu_local;
   real phi_local;
   mu_local = day_onset_tally_tail[i] * (1 - weights[i]) / weights[i];
   phi_local = day_onset_tally_tail[i];
   trunc[i] =  neg_binomial_2_rng(mu_local, phi_local);
   day_onset_tally_tail[i] = day_onset_tally_tail[i] + trunc[i];
   if(check[i] == 1) {
     day_onset_tally_tail[i] -= 1;
   }
   if(day_onset_tally_tail[i] < 0) {
     day_onset_tally_tail[i] = 0;
   }
   day_onset_tally[ndays + i] = day_onset_tally_tail[i];
 }
 
}

