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
 Difference from previous versions: simplified the model, added covariate resampling, more stan friendly parameterization
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
  mu0 = rep_array(0.0, Jb);

 
 for(i in 1:N){
  pa_vec[i] = pa[i];
  if(visit[i]==4.0){
  visit2[i] = 1;
  }
  if(visit[i]==6.0){
  visit2[i] = 2;
  }
  if(visit[i]==7.0){
  visit2[i] = 3;
  }
 }

 cumpa[1] = 0;
 cumets[1] = 0;
 bmizl1[1] = 0;
 bmizsql1[1] = bmizl1[1]^2;
 bmizcul1[1] = bmizl1[1]^3; 
 for(i in 2:N){
  if(visit2[i]==1) {
   cumpa[i] = 0;
   cumets[i] = 0;
   bmizl1[i] = 0;
  }
  else{
   cumpa[i] = cumpa[i-1] + pa[i-1];
   cumets[i] = cumets[i-1] + ets[i-1];
   bmizl1[i] = bmiz[i-1];
  }
  bmizsql1[i] = bmizl1[i]^2;
  bmizcul1[i] = bmizl1[i]^3;
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
 real sigma = sqrt(sigma2);
 vector[Jb] tau = sigma*sqrt(exp(omega));
 real il2 = 2.0/lambda;
}
model{
 row_vector[N] py_;
 row_vector[N] mua_; 
 
 //hyper-hyperpriors for shrinkage
//  lambda ~ gamma(1.0, 10.0); 
  lambda ~ gamma(1.0, 0.1); // more shrinkage

 //hyperpriors for beta
  sigma2 ~ inv_gamma(0.1, 0.1); //scale parameter is switched in stan vs. sas (stan scale parameter = iscale parameter in sas)
  for(j in 1:Jb){
   //unable to vectorize this due to limitations with creating custom sampling distributions
   //tau[j] = sqrt(sigma2*exp(omega[j]));
   omega[j] ~ expexpon(il2); //custom sampler function created above
  }
 
 
 b ~ normal(mu0,tau);
 a ~ normal(0,.3);
 b0 ~ normal(0, 2000);
 a0 ~ normal(0, 2000);

 mua_ = a0 + a[1]*cumets + a[2]*cumpa + a[3]*mage_bl + a[4]*mbmi_bl + a[5]*mheight_bl + a[6]*ets_bl + a[7]*male_bl + a[8]*breastfed_bl + a[9]*mwhite_bl + a[10]*mblack_bl + a[11]*hsgrad_bl + a[12]*somecoll_bl + a[13]*collgrad_bl + a[14]*magesq_bl + a[15]*magecu_bl + a[16]*bmisq_bl + a[17]*bmizl1 + a[18]*agemos + a[19]*cumets .* mage_bl + a[20]*cumets .* mbmi_bl + a[21]*cumets .* mheight_bl + a[22]*cumets .* ets_bl + a[23]*cumets .* male_bl + a[24]*breastfed_bl .* cumets + a[25]*cumets .* mwhite_bl + a[26]*cumets .* mblack_bl + a[27]*cumets .* hsgrad_bl + a[28]*cumets .* somecoll_bl + a[29]*collgrad_bl .* cumets + a[30]*cumets .* magesq_bl + a[31]*cumets .* magecu_bl + a[32]*bmisq_bl .* cumets + a[33]*bmizl1 .* cumets + a[34]*agemos .* cumets;
 py_ = b0 + b[1]*ets + b[2]*cumets + b[3]*cumpa + b[4]*pa_vec + b[5]*mage_bl + b[6]*mbmi_bl + b[7]*mheight_bl + b[8]*ets_bl + b[9]*male_bl + b[10]*breastfed_bl + b[11]*mwhite_bl + b[12]*mblack_bl + b[13]*hsgrad_bl + b[14]*somecoll_bl + b[15]*collgrad_bl + b[16]*magesq_bl + b[17]*magecu_bl + b[18]*bmisq_bl + b[19]*bmizl1 + b[20]*agemos + b[21]*ets .* mage_bl + b[22]*ets .* mbmi_bl + b[23]*ets .* mheight_bl + b[24]*ets .* ets_bl + b[25]*ets .* male_bl + b[26]*breastfed_bl .* ets + b[27]*ets .* mwhite_bl + b[28]*ets .* mblack_bl + b[29]*ets .* hsgrad_bl + b[30]*ets .* somecoll_bl + b[31]*collgrad_bl .* ets + b[32]*ets .* magesq_bl + b[33]*ets .* magecu_bl + b[34]*bmisq_bl .* ets + b[35]*bmizl1 .* ets + b[36]*agemos .* ets;

 pa ~ bernoulli_logit(mua_); 
 bmiz ~ normal(py_, sqrt(sigma2));
 
}
generated quantities{
 real meanpa[2,3] = rep_array(0., 2, 3);
 real meanbmi[2,3] = rep_array(0., 2, 3);
 real meandiff[3];
   {
      // note: sample size of the simulated data need not equal N and should often be larger in complex problems;
      // this will generally require sampling baseline covariates either through a non-parametric sampling scheme (e.g. Bayesian bootstrap) or
      // via a parametric model for baseline covariates. Here, the baseline covariates are fixed across samples, making this a 'conditional' or finite-sample causal effect
      int M = 1000; // resample size
      real Mr = M; // resample size as a real number
      int sampsize =  0;
      real amos;
      row_vector[3] ages = [-1.1956522, 0.1388889, 1.3429952];
      real pa_int;
      real ets_int;
      real cumpa_int;
      real cumets_int;
      real bmizl1_int;
      real bmiz_int;
      vector[N] psample = rep_vector(inv(N), N); // note that this only works when all participants are observed at every time, otherwise need to sample on individual basis!
      
      for(i in 1:N){
        if(i==1)
          sampsize += 1;
        else if(id[i] != id[i-1])
          sampsize += 1;
        
      }
     
      //interventions
      for(set in 1:2){
       ets_int = set-1.0;
       for(i in 1:M){
         int iidx = categorical_rng(psample);// draw random individual from population (just draws a random number from 1:N - works in this case due to balanced data!)
         cumpa_int = 0.;
         cumets_int = 0.;
         bmizl1_int = 0;
         amos = 0.; // -1, 0, 1
         for(time in 1:3){
         amos = ages[time];
         if(time>1){
           cumpa_int += pa_int; // auto lagged from previous time
           cumets_int += (set-1.0);
           bmizl1_int = bmiz_int;
         }
     
         pa_int = bernoulli_rng(inv_logit(
          a0 + a[1]*cumets_int + a[2]*cumpa_int + a[3]*mage_bl[iidx] + a[4]*mbmi_bl[iidx] + a[5]*mheight_bl[iidx] + a[6]*ets_bl[iidx] + a[7]*male_bl[iidx] + a[8]*breastfed_bl[iidx] + a[9]*mwhite_bl[iidx] + a[10]*mblack_bl[iidx] + a[11]*hsgrad_bl[iidx] + a[12]*somecoll_bl[iidx] + a[13]*collgrad_bl[iidx] + a[14]*magesq_bl[iidx] + a[15]*magecu_bl[iidx] + a[16]*bmisq_bl[iidx] + a[17]*bmizl1_int + a[18]*amos + a[19]*cumets_int*mage_bl[iidx] + a[20]*cumets_int*mbmi_bl[iidx] + a[21]*cumets_int*mheight_bl[iidx] + a[22]*cumets_int*ets_bl[iidx] + a[23]*cumets_int*male_bl[iidx] + a[24]*breastfed_bl[iidx]*cumets_int + a[25]*cumets_int*mwhite_bl[iidx] + a[26]*cumets_int*mblack_bl[iidx] + a[27]*cumets_int*hsgrad_bl[iidx] + a[28]*cumets_int*somecoll_bl[iidx] + a[29]*collgrad_bl[iidx]*cumets_int + a[30]*cumets_int*magesq_bl[iidx] + a[31]*cumets_int*magecu_bl[iidx] + a[32]*bmisq_bl[iidx]*cumets_int + a[33]*bmizl1_int*cumets_int + a[34]*amos*cumets_int
          ));
       
         bmiz_int = normal_rng(
          b0 + b[1]*ets_int + b[2]*cumets_int + b[3]*cumpa_int + b[4]*pa_int + b[5]*mage_bl[iidx] + b[6]*mbmi_bl[iidx] + b[7]*mheight_bl[iidx] + b[8]*ets_bl[iidx] + b[9]*male_bl[iidx] + b[10]*breastfed_bl[iidx] + b[11]*mwhite_bl[iidx] + b[12]*mblack_bl[iidx] + b[13]*hsgrad_bl[iidx] + b[14]*somecoll_bl[iidx] + b[15]*collgrad_bl[iidx] + b[16]*magesq_bl[iidx] + b[17]*magecu_bl[iidx] + b[18]*bmisq_bl[iidx] + b[19]*bmizl1_int + b[20]*amos + b[21]*ets_int*mage_bl[iidx] + b[22]*ets_int*mbmi_bl[iidx] + b[23]*ets_int*mheight_bl[iidx] + b[24]*ets_int*ets_bl[iidx] + b[25]*ets_int*male_bl[iidx] + b[26]*breastfed_bl[iidx]*ets_int + b[27]*ets_int*mwhite_bl[iidx] + b[28]*ets_int*mblack_bl[iidx] + b[29]*ets_int*hsgrad_bl[iidx] + b[30]*ets_int*somecoll_bl[iidx] + b[31]*collgrad_bl[iidx]*ets_int + b[32]*ets_int*magesq_bl[iidx] + b[33]*ets_int*magecu_bl[iidx] + b[34]*bmisq_bl[iidx]*ets_int + b[35]*bmizl1_int*ets_int + b[36]*amos*ets_int
          , sigma);
     
          meanbmi[set, time] += bmiz_int;
          meanpa[set, time] +=  pa_int;
        } // end time loop
       }//end i loop
       
      }//end set loop
      
      //kludgy way to get average bmi under each intervention at each time point
      for(set in 1:2){
       for(time in 1:3){
        meanbmi[set,time] = meanbmi[set,time] / Mr;
        meanpa[set,time] = meanpa[set,time] / Mr;
        }
       }
       for(time in 1:3){
         meandiff[time] = meanbmi[2,time] - meanbmi[1,time];
       }
      /*
      */
   }
}

