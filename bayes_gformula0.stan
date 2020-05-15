//////////////////////////////////////////////////////////////////////////////////////////
// 2) Bayesian g-formula (using data generated in SAS simulation)
//////////////////////////////////////////////////////////////////////////////////////////
data {
   int<lower=0> N;
   int<lower=0,upper=1> y[N];
   row_vector[N] x1;
   row_vector[N] x2;
   int<lower=0,upper=1> l2[N];
}
transformed data{
// intevention variables
   int<lower=0,upper=1> g1[N];
   int<lower=0,upper=1> g0[N];
   g1 <- rep_array(1, N);   
   g0 <- rep_array(0, N);   
}
parameters {
   real inta;
   real intb;
   vector[2] a;
   vector[4] b;
}
model {
// priors
   inta ~ normal(0, 1000);
   inta ~ normal(0, 1000);
   b ~ normal(0, 15);
   a ~ normal(0, 15);
// joint model for observed data
for(i in 1:N){
   l2[i] ~ bernoulli_logit(inta + x1[i]*a[1]);
   y[i]  ~ bernoulli_logit(intb + b[1]*x1[i] + b[2]*x2[i] + b[3]*l2[i]);
   }
}

generated quantities {
   real rd;
   real bias;
   real rdi[N];
   // posterior predictive distribution (potential outcomes)
   // note: there is no baseline covariate distribution in this example, so no baseline population to draw from - the
   //  sample size of the simulated data need not equal N and should often be larger in complex problems;
  for(i in 1:N){
    // individual potential risk difference 
    // (more efficient than generating potential outcomes)
    rdi[i] <- bernoulli_rng(inv_logit(intb + b[1]*g1[i] + b[2]*g1[i] + 
                            b[3]*bernoulli_rng(inv_logit(inta + g1[i]*a[1])))) - 
              bernoulli_rng(inv_logit(intb + b[1]*g0[i] + b[2]*g0[i] + 
                            b[3]*bernoulli_rng(inv_logit(inta + g0[i]*a[1]))));
   }
   rd <- mean(rdi);
   bias <- rd-0.2;
} 