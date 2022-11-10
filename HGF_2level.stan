functions{
  
  real unitsquare(x, beta){
    return x^beta/((x^beta)+(1-x)^beta)
  }
}

data {
  int<lower=0> N;   //trials (306)
  int u[N];         //inputs
  int y[N];         //responses
}
parameters {
  real omega;
  real<lower=0> beta;
}



transformed parameters {
  
  real<lower=0,upper=1> mu1hat[N];
  real<lower=0,upper=0.25> sa1hat[N];
  real<lower=0> sa2hat[N];
  real<lower=-1,upper=1> da[N];
  real mu2[N];
  real<lower=0> sa2[N];
  
  
  for (k in 2:N){
    mu1hat[k] = inv_logit(mu2[k]);       //first level belief (prediction belief);
    sa1hat[k] = mu1hat[k]*(1-mu1hat[k]); //first level prediction uncertainty;
    da[k] = u[k]-mu1hat[k];              //prediction error;
    //updates
    sa2hat[k] = sa2[k-1]+exp(omega);
    sa2[k] = 1/((1/sa2hat[k])+sa1hat[k]);
    mu2[k] = mu2[k-1]+da[k]*sa2[k];
    
    
  }
  
  
  //response model:
  
  
}


model {
  
  for (j in 1:N){
    y[j] ~ unitsquare(mu1hat[j],beta);
  }
  
  mu2 ~ normal(0,sa2);
  sa2 ~ normal(0,0);
  omega ~ normal(-3.6,3);
  beta ~ normal(2,5);
  
  

    
}
  