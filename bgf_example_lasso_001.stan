/*---------------------------------------------------------------------------------------------------------------------
 Author: Alex Keil
 Program: bgf_example_lasso_001.stan
 Language: STAN 2.6.0
 Date: Monday, June 1, 2015 at 2:43:12 PM
 Project: Bayes g-formula
 Tasks: Bayesian g-formula analysis using a real world example
 Data in: 
 Data out:
 Description: Examining the effects of time-varying environmental tobacco smoke on childhood BMI z-scores
   This program fits a model with all null centered priors under no hierarchy.
 Difference from previous versions: simplified the model
 Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
---------------------------------------------------------------------------------------------------------------------*/
//

functions{
//expexpon function from sas (exponenential-exponential distribution)
// note that this is distinct from the double exponenential distribution
// call using expexpon(scale);
 real expexpon_log(real y, real scale){
  return log(1 / scale * exp(y) * exp( -exp(y) / scale));
  }
}
data{
 int<lower=0> N;                 // number of observations (num IDs * 3 time points each)
 int<lower=0> Jb;                 // number of parameters, model for physical activity
 int<lower=0> Ja;                 // number of parameters, outcome model
 row_vector[N] id;
 row_vector[N] mage_bl;
 row_vector[N] mbmi_bl;
 row_vector[N] mheight_bl;
 row_vector[N] ets_bl;
 row_vector[N] raceth_bl;
 row_vector[N] educat_bl;
 row_vector[N] male_bl;
 row_vector[N] breastfed_bl;
 row_vector[N] mwhite_bl;
 row_vector[N] mblack_bl;
 row_vector[N] hsgrad_bl;
 row_vector[N] somecoll_bl;
 row_vector[N] collgrad_bl;
 row_vector[N] magesq_bl;
 row_vector[N] magecu_bl;
 row_vector[N] bmisq_bl;
 row_vector[N] bmicu_bl;
 row_vector[N] agemos;
 row_vector[N] agemossq;
 row_vector[N] agemoscu;
 row_vector[N] bmiz;              // primary outcome
 row_vector[N] bmizsq;
 row_vector[N] bmizcu;
 row_vector[N] visit;
 row_vector[N] ets;               // primary exposure
 int pa[N];                       // time varying confounder
}
transformed data{
 int visit2[N];
 //cumulative variables, lagged variables
 row_vector[N] cumpa;
 row_vector[N] cumets;
 row_vector[N] bmizl1;
 row_vector[N] bmizsql1;
 row_vector[N] bmizcul1;
 row_vector[N] pa_vec; // need vector form of PA for model vectorization
 real mu0[Jb];

 //prior mean for beta 
  mu0 <- rep_array(0.0, Jb);

 
 for(i in 1:N){
  pa_vec[i] <- pa[i];
  if(visit[i]==4.0){
  visit2[i] <- 1;
  }
  if(visit[i]==6.0){
  visit2[i] <- 2;
  }
  if(visit[i]==7.0){
  visit2[i] <- 3;
  }
 }

 cumpa[1] <- 0;
 cumets[1] <- 0;
 bmizl1[1] <- 0;
 bmizsql1[1] <- bmizl1[1]^2;
 bmizcul1[1] <- bmizl1[1]^3; 
 for(i in 2:N){
  if(visit2[i]==1) {
   cumpa[i] <- 0;
   cumets[i] <- 0;
   bmizl1[i] <- 0;
  }
  else{
   cumpa[i] <- cumpa[i-1] + pa[i-1];
   cumets[i] <- cumets[i-1] + ets[i-1];
   bmizl1[i] <- bmiz[i-1];
  }
  bmizsql1[i] <- bmizl1[i]^2;
  bmizcul1[i] <- bmizl1[i]^3;
 }
 

 
}
parameters{
 real a[Ja];
 real b[Jb];
 real a0;
 real b0;
 vector<lower=0>[Jb] omega;
 real<lower=0> sigma2;
 real<lower=0> lambda;

}
transformed parameters{
}
model{
 vector[Jb] tau;
 row_vector[N] py_;
 row_vector[N] mua_; 
 
 //hyper-hyperpriors for shrinkage
//  lambda ~ gamma(1.0, 10.0); 
  lambda ~ gamma(1.0, 0.1); // more shrinkage



 //hyperpriors for beta
  sigma2 ~ inv_gamma(0.1, 0.1); //scale parameter is switched in stan vs. sas (stan scale parameter = iscale parameter in sas)
  for(j in 1:Jb){
   //unable to vectorize this due to limitations with creating custom sampling distributions
   omega[j] ~ expexpon(2.0/lambda); //custom sampler function created above
   tau[j] <- sqrt(sigma2*exp(omega[j]));
  }
 
 
 b ~ normal(mu0,tau);
 a ~ normal(0,.3);
 b0 ~ normal(0, 2000);
 a0 ~ normal(0, 2000);

 mua_ <- a0 + a[1]*cumets + a[2]*cumpa + a[3]*mage_bl + a[4]*mbmi_bl + a[5]*mheight_bl + a[6]*ets_bl + a[7]*male_bl + a[8]*breastfed_bl + a[9]*mwhite_bl + a[10]*mblack_bl + a[11]*hsgrad_bl + a[12]*somecoll_bl + a[13]*collgrad_bl + a[14]*magesq_bl + a[15]*magecu_bl + a[16]*bmisq_bl + a[17]*bmizl1 + a[18]*agemos + a[19]*cumets .* mage_bl + a[20]*cumets .* mbmi_bl + a[21]*cumets .* mheight_bl + a[22]*cumets .* ets_bl + a[23]*cumets .* male_bl + a[24]*breastfed_bl .* cumets + a[25]*cumets .* mwhite_bl + a[26]*cumets .* mblack_bl + a[27]*cumets .* hsgrad_bl + a[28]*cumets .* somecoll_bl + a[29]*collgrad_bl .* cumets + a[30]*cumets .* magesq_bl + a[31]*cumets .* magecu_bl + a[32]*bmisq_bl .* cumets + a[33]*bmizl1 .* cumets + a[34]*agemos .* cumets;
 py_ <- b0 + b[1]*ets + b[2]*cumets + b[3]*cumpa + b[4]*pa_vec + b[5]*mage_bl + b[6]*mbmi_bl + b[7]*mheight_bl + b[8]*ets_bl + b[9]*male_bl + b[10]*breastfed_bl + b[11]*mwhite_bl + b[12]*mblack_bl + b[13]*hsgrad_bl + b[14]*somecoll_bl + b[15]*collgrad_bl + b[16]*magesq_bl + b[17]*magecu_bl + b[18]*bmisq_bl + b[19]*bmizl1 + b[20]*agemos + b[21]*ets .* mage_bl + b[22]*ets .* mbmi_bl + b[23]*ets .* mheight_bl + b[24]*ets .* ets_bl + b[25]*ets .* male_bl + b[26]*breastfed_bl .* ets + b[27]*ets .* mwhite_bl + b[28]*ets .* mblack_bl + b[29]*ets .* hsgrad_bl + b[30]*ets .* somecoll_bl + b[31]*collgrad_bl .* ets + b[32]*ets .* magesq_bl + b[33]*ets .* magecu_bl + b[34]*bmisq_bl .* ets + b[35]*bmizl1 .* ets + b[36]*agemos .* ets;

 pa ~ bernoulli_logit(mua_); 
 bmiz ~ normal(py_, sqrt(sigma2));
 
}
generated quantities{
 real meanpa[2,3];
 real meanbmi[2, 3];
 real cnt[2,3];
   // note: sample size of the simulated data need not equal N and should often be larger in complex problems;
   // this will generally require sampling baseline covariates either through a non-parametric sampling scheme (e.g. Bayesian bootstrap) or
   // via a parametric model for baseline covariates. Here, the baseline covariates are fixed across samples, making this a 'conditional' or finite-sample causal effect
 row_vector[N] pa_int;
 row_vector[N] ets_int;
 row_vector[N] cumpa_int;
 row_vector[N] cumets_int;
 row_vector[N] bmiz_int;
 row_vector[N] bmizl1_int;
 
 meanpa <- rep_array(0, 2, 3);
 meanbmi <- rep_array(0, 2, 3);
 cnt <- rep_array(0, 2, 3);

 //interventions
 for(set in 1:2){
  ets_int <- rep_row_vector(set-1.0, N);
  for(i in 1:N){
   if(i==1)  cumpa_int[1] <- 0;
   if(i==1)  cumets_int[1] <- 0;
   if(i==1)  bmizl1_int[1] <- 0;
   if(visit2[i]==1) {
    cumpa_int[i] <- 0;
    cumets_int[i] <- 0;
    bmizl1_int[i] <- 0;
   }
   else{
    cumpa_int[i] <- cumpa_int[i-1] + pa_int[i-1];
    cumets_int[i] <- cumets_int[i-1] + ets_int[i-1];
    bmizl1_int[i] <- bmiz_int[i-1];
   }

  pa_int[i] <- bernoulli_rng(inv_logit(
   a0 + a[1]*cumets_int[i] + a[2]*cumpa_int[i] + a[3]*mage_bl[i] + a[4]*mbmi_bl[i] + a[5]*mheight_bl[i] + a[6]*ets_bl[i] + a[7]*male_bl[i] + a[8]*breastfed_bl[i] + a[9]*mwhite_bl[i] + a[10]*mblack_bl[i] + a[11]*hsgrad_bl[i] + a[12]*somecoll_bl[i] + a[13]*collgrad_bl[i] + a[14]*magesq_bl[i] + a[15]*magecu_bl[i] + a[16]*bmisq_bl[i] + a[17]*bmizl1_int[i] + a[18]*agemos[i] + a[19]*cumets_int[i]*mage_bl[i] + a[20]*cumets_int[i]*mbmi_bl[i] + a[21]*cumets_int[i]*mheight_bl[i] + a[22]*cumets_int[i]*ets_bl[i] + a[23]*cumets_int[i]*male_bl[i] + a[24]*breastfed_bl[i]*cumets_int[i] + a[25]*cumets_int[i]*mwhite_bl[i] + a[26]*cumets_int[i]*mblack_bl[i] + a[27]*cumets_int[i]*hsgrad_bl[i] + a[28]*cumets_int[i]*somecoll_bl[i] + a[29]*collgrad_bl[i]*cumets_int[i] + a[30]*cumets_int[i]*magesq_bl[i] + a[31]*cumets_int[i]*magecu_bl[i] + a[32]*bmisq_bl[i]*cumets_int[i] + a[33]*bmizl1_int[i]*cumets_int[i] + a[34]*agemos[i]*cumets_int[i]
   ));
  
  bmiz_int[i] <- 
   b0 + b[1]*ets_int[i] + b[2]*cumets_int[i] + b[3]*cumpa_int[i] + b[4]*pa_int[i] + b[5]*mage_bl[i] + b[6]*mbmi_bl[i] + b[7]*mheight_bl[i] + b[8]*ets_bl[i] + b[9]*male_bl[i] + b[10]*breastfed_bl[i] + b[11]*mwhite_bl[i] + b[12]*mblack_bl[i] + b[13]*hsgrad_bl[i] + b[14]*somecoll_bl[i] + b[15]*collgrad_bl[i] + b[16]*magesq_bl[i] + b[17]*magecu_bl[i] + b[18]*bmisq_bl[i] + b[19]*bmizl1_int[i] + b[20]*agemos[i] + b[21]*ets_int[i]*mage_bl[i] + b[22]*ets_int[i]*mbmi_bl[i] + b[23]*ets_int[i]*mheight_bl[i] + b[24]*ets_int[i]*ets_bl[i] + b[25]*ets_int[i]*male_bl[i] + b[26]*breastfed_bl[i]*ets_int[i] + b[27]*ets_int[i]*mwhite_bl[i] + b[28]*ets_int[i]*mblack_bl[i] + b[29]*ets_int[i]*hsgrad_bl[i] + b[30]*ets_int[i]*somecoll_bl[i] + b[31]*collgrad_bl[i]*ets_int[i] + b[32]*ets_int[i]*magesq_bl[i] + b[33]*ets_int[i]*magecu_bl[i] + b[34]*bmisq_bl[i]*ets_int[i] + b[35]*bmizl1_int[i]*ets_int[i] + b[36]*agemos[i]*ets_int[i]; 

  for(v in 1:3){   
    if(visit2[i]==v){
     meanbmi[set, v] <- meanbmi[set, v]+bmiz_int[i];
     meanpa[set, v] <- meanpa[set, v]+pa_int[i];
     cnt[set,v] <- cnt[set,v] + 1;
    }
   }
  }//end i loop
  
  
 }//end set loop
 
 //kludgy way to get average bmi under each intervention at each time point
 for(set in 1:2){
  for(v in 1:3){
   meanbmi[set,v] <- meanbmi[set,v] / cnt[set,v];
   meanpa[set,v] <- meanpa[set,v] / cnt[set,v];
   }
  }
 /*
 */
}