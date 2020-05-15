DM LOG 'clear;' CONTINUE; DM OUT 'clear;' CONTINUE; 
/**********************************************************************************************************************
* Author: Alex Keil
* Program: bayes_gformula
* System notes: SAS 9.4 (SAS/STAT 13.1), Windows 8.1 [does not work in lower versions of SAS due to MCMC changes]
* Date: Wednesday, July 2, 2014 at 11:41:42 AM
* Project: miscellaneous
* Tasks: Implement bayesian g-formula algorithm in PROC MCMC
* Data in: NA
* Data out: NA
* Description: program simulates data for a small longitudinal analysis
      1) frequentist g-formula
      2) bayesian g-formula under vague null centered priors for every coefficient
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

* Next steps:
 1) explore a way to place priors directly on the causal risk difference 
 2) compare bootstrap SE in frequentist analysis to MCSD from single run of Bayes g-formula
 3) generalize for counting process style data
**********************************************************************************************************************/
*clear the log window and the output window;
OPTIONS MERGENOBY = warn NODATE NONUMBER LINESIZE = 120  PAGESIZE=80 SKIP = 2 FORMDLIM = '-' MPRINT NOCENTER;
OPTIONS FORMCHAR = '|----|+|---+=|-/\<>*';



%LET n=60; *sample size - does Bayes g-formula outperform frequentist g-formula in small samples?;
%LET trueRD = .2;


************************************************************************
1) simulate data
************************************************************************;
DATA d (keep = id X: L: Y: py00);
CALL STREAMINIT(1192887);
*potential variables under interventions - set all to missing;
/* y_0=.; y_1=.; l2_0=.; l2_1=.;
LABEL  y0="Potential outcome - never exposed"  y1 = "Potential outcome - always exposed" 
 l0="Potential covariate - never exposed"  l1 = "Potential covariate - always exposed" ;
*/
*observed data;
 DO id = 1 TO &N;
  x1 = RAND("bernoulli", 0.5);
  py00 = RAND("uniform")*0.1 + 0.4;  
  l2 = RAND("bernoulli", 1/(1+exp(-1 + x1 + py00)));
  x2 = RAND("bernoulli", 1/(1+exp(-1 + x1 + l2)));
  py = py00 + &trueRD*((x1 + x2)/2); *true risk difference per unit exposure;
  y = RAND("bernoulli", py);
  OUTPUT;
 END;
RUN;

*OUTPUT FOR STAN;
PROC EXPORT DATA = d OUTFILE="Z:/Documents/Papers/2015_bayes_gformula/Simulations/simdata.csv" DBMS=CSV;RUN;

*0.1) Example data;
PROC PRINT DATA = d (OBS=10); TITLE "example data";
PROC MEANS DATA = d MAXDEC=3 FW=5; TITLE "observed covariate distributions"; VAR x1 x2 l2 y py00; RUN;
************************************************************************
2) Frequentist g-formula
************************************************************************;
PROC LOGISTIC NOPRINT DATA = d DESCENDING OUTEST=lm (KEEP=intercept x1 RENAME=(intercept=a0 x1=a1));
 MODEL l2 = x1;
PROC LOGISTIC NOPRINT DATA = d DESCENDING OUTEST=ym (KEEP=intercept x1 x2 l2 RENAME=(intercept=b0 x1=b1 x2=b2 l2=b3));
 MODEL y = x1 x2 l2;
DATA parms;
IF _N_=1 THEN SET lm;
IF _N_=1 THEN SET ym;
IF _n_=1 THEN OUTPUT;
DATA e;
 SET d;
 IF _n_=1 THEN SET parms;
PROC SURVEYSELECT DATA = e OUT=dmc N=10000 SEED=1223 METHOD=URS OUTHITS NOPRINT;  
DATA interv;
 CALL STREAMINIT(121232);
  SET dmc(KEEP=a0 a1 b0 b1 b2 b3);
  DO setx = 0,1;
   lpl = a0 + a1*setx;
   l2 = RAND("bernoulli", 1/(1+exp(-lpl)));
   lpy = b0 + (b1+b2)*setx + b3*l2;
   py=1/(1+exp(-lpy));
   OUTPUT;
  END;
PROC SORT DATA = interv; BY setx; 
PROC MEANS DATA = interv MEAN NOPRINT;
  BY setx;
  VAR py;
  OUTPUT OUT=my MEAN=my;
 DATA my;
  SET my END=eof;
  RETAIN rd py0 py1 0;
  LABEL rd = "risk difference (frequentist g-formula)"
        py0 = "Potential risk given no exposure"
        py1 = "Potential risk given always exposed";
  IF setx=1 THEN py1=my; 
  ELSE IF setx=0 THEN py0=my;
  rd = py1-py0;
  bias = rd-&trueRD;
 IF eof THEN OUTPUT;
 PROC PRINT DATA = my;
  TITLE "frequentist g-formula risk difference (truth=&truerd)";
  VAR rd py0 py1 bias;
 RUN;
 
 
************************************************************************
3) Bayes g-formula
************************************************************************;
DATA d;
 SET d;
 ind = _n_;
RUN;

%LET m=1000;
PROC MCMC DATA = d SEED = 12123 NMC = 100000 NBI=50000 OUTPOST=dpost 
MONITOR=(rd py_1 py_0 bias);
 TITLE "Bayesian g-formula risk difference (truth=&truerd)";
 ODS SELECT PostSumInt ESS TADPanel;
*priors - weak shrinkage;
 PARMS a0 a1 b0 b1 b2 b3;
 PRIOR a1 b1 b2 b3 ~ NORMAL(0, var=15); PRIOR a0 b0~ NORMAL(LOG(0.5), var=1000);

*joint model for observed data;
  pl2 = LOGISTIC(a0+a1*x1);
  MODEL L2 ~ BINARY(pl2);
  py = LOGISTIC(b0 + b1*x1 + b2*x2 + b3*l2);
  MODEL y ~ BINARY(py);

*posterior predictive risk difference();
* note: there is no baseline covariate distribution in this example, so no baseline population to draw from - the
*  sample size of the simulated data need not equal N and should often be larger in complex problems;
  ARRAY rdi[&M];
  *rdi[IND] = y_1-y_0;
  rdi[IND] = RAND("bernoulli",  LOGISTIC(b0 + b1 + b2 + 
				  b3*RAND("bernoulli", LOGISTIC(a0 + a1)))) - 
             RAND("bernoulli", LOGISTIC(b0  + 
                  b3*RAND("bernoulli", LOGISTIC(a0))));
  IF ind=&n THEN DO;
   rd =  MEAN(of rdi:);
   bias = rd-&trueRD;
  END;
  LABEL rd="Average difference of potential outcome ";
RUN;






