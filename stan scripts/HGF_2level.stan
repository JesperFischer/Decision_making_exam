data {
  int<lower=0> ntrials;
  int<lower=0> nsubs;
  int<lower=0,upper=1> y[ntrials, nsubs];
  int<lower=0,upper=1> u[ntrials, nsubs];
}

parameters {
  real mu_omega;
  real<lower=0> sigma_omega;
  real<lower=0> mu_beta;
  real<lower=0> sigma_beta;
  real omega[nsubs];
  real<lower=0> beta[nsubs];
}



transformed parameters {
  real mu2[ntrials,nsubs];
  real<lower=0> sa2[ntrials,nsubs];
  real sa2hat[ntrials,nsubs];
  real mu1hat[ntrials,nsubs];
  real sa1hat[ntrials,nsubs];
  real da[ntrials,nsubs];
  real exp_p[ntrials,nsubs];
  

  for (s in 1:nsubs) {
    mu2[1,s] = 0;
    sa2[1,s] = 2;
    
    for (t in 2:ntrials) {
  
      sa2hat[t,s] = sa2[t-1,s]+exp(omega[s]);
      mu1hat[t,s] = inv_logit(mu2[t-1,s]);
      
      sa1hat[t,s] = mu1hat[t,s]*(1-mu1hat[t,s]);
      
      da[t,s] = u[t,s]-mu1hat[t,s];
      
      sa2[t,s] = 1/((1/sa2hat[t,s])+sa1hat[t,s]);
      
      mu2[t,s] = mu2[t-1,s]+da[t,s]*sa2[t,s];
      
      exp_p[t,s] = (mu1hat[t,s]^beta[s])/((mu1hat[t,s]^beta[s])+(1-mu1hat[t,s])^(beta[s]));
  }
}
}



model {
  mu_omega ~ normal(-1,0.3);
  sigma_omega ~ normal(1,0.2);
  sigma_beta~ normal(0.3,0.2);
  mu_beta ~ gamma(4,2);
  
  for (s in 1:nsubs) {
    omega[s] ~ normal(mu_omega,sigma_omega);
    beta[s] ~ normal(mu_beta,sigma_beta);
    
    for (t in 2:ntrials) {
      
      y[t,s] ~ bernoulli_logit(exp_p[t,s]);
  }
}
}

