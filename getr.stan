// -- CALCULATE R(t) ---
 // getr(backc,si,size)

 int n = ndays + maxdelay;

 // calculate lambda
 vector[ndays + maxdelay - 1] lambda;
 
 for(i in 1:(n-1)) {
   if(i <= sipN) {
     vector[i] c = day_onset_tally[1:i];
     vector[i] sic = sip[1:i];
     vector[i] rev_sic;
     int i_backwards = i;
     for(ii in 1:i) {
       rev_sic[ii] = sic[i_backwards];
       i_backwards -= 1;
     }
     lambda[i] = dot_product(to_vector(c), to_vector(rev_sic));
     
   } else {
     vector[sipN] c = day_onset_tally[(i - sipN + 1):i];
     vector[sipN] sic = sip;
     lambda[i] = dot_product(to_vector(c), to_vector(sic));
   }
   
 }
 
 // finish RT calc
 int nr = ndays + maxdelay - windowsize - 1;
 vector[ndays + maxdelay - windowsize - 1] rt;
 real shape;
 real scale;
 vector[windowsize + 1] incid;
 vector[windowsize + 1] dem1;

  for (i in 1:nr) {
      incid = day_onset_tally[(i+1):(i + windowsize + 1)];
      shape = sum(incid) + 1;
      dem1 = lambda[(i):(i + windowsize)];
      scale = 1 / (sum(dem1) + 0.2);
      rt[i] = shape * scale;
  }