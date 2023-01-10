
sigmoid = function(x) {
  1 / (1 + exp(-x))
}


hgf_agent = function(u,omega,beta){



ntrials = 307

mu2 = array(NA,ntrials)
sa2 = array(NA,ntrials)
mu1hat = array(NA,ntrials)
da = array(NA,ntrials)
sa1hat = array(NA,ntrials)
sa2hat = array(NA,ntrials)
yhat= array(NA,ntrials)
exp_p= array(NA,ntrials)


mu2[1] = dnorm(0,1)
sa2[1] = dnorm(0.4,5)
mu1hat[1] = 0
sa1hat[1] = 0
sa2hat[1] = 0
da[1] = 0
yhat[1] = 0
exp_p[1] = 0


for  (t in 2:ntrials){
  
  sa2hat[t] = sa2[t-1]+exp(omega)
  mu1hat[t] = sigmoid(mu2[t-1])
  
  
  sa1hat[t] = mu1hat[t]*(1-mu1hat[t])
  
  da[t] = u[t]-mu1hat[t]
  
  
  sa2[t] = 1/((1/sa2hat[t])+sa1hat[t])
  
  mu2[t] = mu2[t-1]+da[t]*sa2[t]
  
  
  
  exp_p[t] = mu1hat[t]^beta/((mu1hat[t]^beta)+(1-mu1hat[t])^(beta))
  
  
  
  yhat[t] = rbinom(n = 1,size = 1, prob = exp_p[t])
  
}

data = data.frame(u,mu1hat,da,sa1hat,mu2,sa2,sa2hat,exp_p,yhat)

}

