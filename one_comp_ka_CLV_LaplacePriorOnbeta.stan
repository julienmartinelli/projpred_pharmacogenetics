//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions{

  real OneCptModel1(real time, int dose, real CL, real V1, real ka){
    real k;
    k=CL/V1;

    return (dose / V1) * ka/(ka-k)*(exp(-k*time)/(1-exp(-k*12))+exp(-ka*time)/(1-exp(-ka*12)));
  }

  vector OneCptModel(real[] time, int dose, real CL, real V1, real ka){
    vector[size(time)] conc;

    for(i in 1:size(time))
      conc[i] = OneCptModel1(time[i], dose, CL, V1, ka);

    return conc;
  }
  
}

data{
  int<lower = 1> ntot; // total number of observations nbr tot
  int<lower=1> N; // number of patients  
  int off_data[N]; // index at which the ith individual data starts
  real time[ntot]; // vector of observed times
  real cObs[ntot];    // vector of observed concentrations
  int dose;          // dose (here same dose for all)
  int nSNP; // number of SNPs explored
  vector[nSNP] tabSNP[N]; // array of SNPs
}



parameters{
  real<lower = 0> kaHat;
  real<lower = 0> CLHat;
  real<lower = 0> V1Hat;
  vector[nSNP] beta;  // coefficient d'effet associé aux autres snp
  real <lower =0> lambda;   // paramètre de la prior de Laplace
  
  real<lower = 0> omegaka;
  real<lower = 0> omegaCL;
  real<lower = 0> omegaV1;
  
  real sigma1;
  
  vector[N] etaka;
  vector[N] etaCL;
  vector[N] etaV1; 
}

transformed parameters{
  real<lower = 0> ka[N];
  real<lower = 0> CL[N];
  real<lower = 0> V1[N];
  vector<lower = 0>[ntot] cHat;


  for(i in 1:N){
    ka[i] = kaHat*exp(omegaka*etaka[i]);
    CL[i] = CLHat*exp(omegaCL*etaCL[i])*exp(to_row_vector(tabSNP[i])*beta);
    V1[i] = V1Hat*exp(omegaV1*etaV1[i]);
  
    cHat[off_data[i]:(off_data[i]+2)] = OneCptModel(time[off_data[i]:(off_data[i]+2)],
                                        dose, CL[i], V1[i], ka[i]);
  }
}

model{
    kaHat ~ lognormal(.01, 0.5);
    CLHat ~ lognormal(.9, 0.5);
    V1Hat ~ lognormal(5.3, 0.5);
    
    omegaka ~ cauchy(0, 2.5);
    omegaCL ~ cauchy(0, 2.5);
    omegaV1 ~ cauchy(0, 2.5);
    
    sigma1 ~ cauchy(0, 2.5);
    
    lambda ~ uniform(0, 40);// prior on parameter of the Laplace prior
    
    for(iSNP in 1:nSNP)
      beta[iSNP] ~ double_exponential(0, 1.0 / lambda); //Laplace prior for all SNP
    
    // Inter-individual variability
    etaka ~ normal(0, 1);
    etaCL ~ normal(0, 1);
    etaV1 ~ normal(0, 1);
    
    cObs ~ normal(cHat, cHat*sigma1); 
}

 generated quantities{
  real CLPred[N];
  real kaPred[N];
  real V1Pred[N];
  vector[ntot] cHatPred;
   vector[ntot] cObsCond;
  vector[ntot] cObsPred;

  vector[N] etakaPred;
  vector[N] etaCLPred;
  vector[N] etaV1Pred;

  for(i in 1:N){
    etakaPred[i] = normal_rng(0, 1);
    etaCLPred[i] = normal_rng(0, 1);
    etaV1Pred[i] = normal_rng(0, 1);

    kaPred[i] = kaHat*exp(omegaka*etakaPred[i]);
    CLPred[i] = CLHat*exp(omegaCL*etaCLPred[i])*exp(to_row_vector(tabSNP[i])*beta);//
    V1Pred[i] = V1Hat*exp(omegaV1*etaV1Pred[i]);

    cHatPred[off_data[i]:(off_data[i]+2)] = OneCptModel(time[off_data[i]:(off_data[i]+2)],
                                        dose, CLPred[i], V1Pred[i], kaPred[i]);
  }

   for(j in 1:ntot){
     cObsCond[j] = normal_rng(cHat[j], sigma1*cHat[j]);
     cObsPred[j] = normal_rng(cHatPred[j], sigma1*cHatPred[j]);
   }
 }


