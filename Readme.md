# linelistSTAN

The obvious To-Dos are to implement the right-truncated NB distribution function
in the likelihood and rng forms

See here: https://discourse.mc-stan.org/t/rng-for-truncated-distributions/3122/14 
for starters

# next step 
See here: https://discourse.mc-stan.org/t/rng-for-truncated-distributions/3122/14 
for starters

# notes from TengLong

[x] 1-Impute the missing report delays (I did not see you do that in the Stan code; but that’s not hard). 

[x] 2-Back calculation (based on the missing report delays; harder than we thought in Stan).

3-Nowcasting (getting harder; see the technical details in the original PLOS CB paper).

4-Calculate the reproductive number (getting even harder, because it depends on 2, 3 and the serial interval). 

5-Assemble everything together in every single iteration (one should depends on another, I would say this is probably the hardest part, how you design it). 