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
  real < lower =0> scale_global ; // scale for the half -t prior for tau
                                  // ( tau0 = scale_global * sigma )
  real < lower =1> nu_global ; // degrees of freedom for the half -t prior for tau
  real < lower =1> nu_local ; // degrees of freedom for the half -t priors for lambdas
                              // ( nu_local = 1 corresponds to the horseshoe )
  real <lower =0> slab_scale ; // slab scale for the regularized horseshoe
  real <lower =0> slab_df ; // slab degrees of freedom for the regularized horseshoe
}



parameters{
  real<lower = 0> kaHat;
  real<lower = 0> CLHat;
  real<lower = 0> V1Hat;
  
  real<lower = 0> omegaka;
  real<lower = 0> omegaCL;
  real<lower = 0> omegaV1;
  
  real sigma1;
  
  vector[N] etaka;
  vector[N] etaCL;
  vector[N] etaV1; 

  // auxiliary variables that define the global and local parameters
  vector [nSNP] z;
  real < lower =0> r1_global ;
  real < lower =0> r2_global ;
  vector < lower =0 >[nSNP] r1_local ;
  vector < lower =0 >[nSNP] r2_local ;
  real < lower = 0 > caux ;
}

transformed parameters{
  real<lower = 0> ka[N];
  real<lower = 0> CL[N];
  real<lower = 0> V1[N];
  vector<lower = 0>[ntot] cHat;
  real < lower =0> tau; // global shrinkage parameter
  real < lower = 0 > c ; // slab scale
  vector < lower =0 >[nSNP] lambda ; // local shrinkage parameters
  vector < lower =0 >[nSNP] lambda_tilde ; // truncated local shrinkage parameters
  vector[nSNP] beta;  // coefficient d'effet associé aux autres snp

  lambda = r1_local .* sqrt ( r2_local );
  tau = r1_global * sqrt ( r2_global ) * scale_global * omegaCL;
  c = slab_scale * sqrt ( caux );
  lambda_tilde = sqrt ( c ^2 * square ( lambda ) ./ (c ^2 + tau ^2* square ( lambda )) );
  beta = z .* lambda_tilde * tau;

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
    
    // Shrinkage

    z ~ normal(0, 1);
    r1_local ~ normal (0.0 , 1.0);
    r2_local ~ inv_gamma (0.5* nu_local , 0.5* nu_local );
    // half -t prior for tau
    r1_global ~ normal (0.0 , 1.0);
    r2_global ~ inv_gamma (0.5* nu_global , 0.5* nu_global );
    caux ~ inv_gamma (0.5* slab_df, 0.5* slab_df );
    
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


