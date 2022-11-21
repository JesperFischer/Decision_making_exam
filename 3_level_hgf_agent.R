
sigmoid = function(x) {
  1 / (1 + exp(-x))
}

hgf_agent_3level = function(u,omega,theta,kappa,beta){



ntrials = 307

mu2 = array(NA,ntrials)
sa2 = array(NA,ntrials)
mu1hat = array(NA,ntrials)
da = array(NA,ntrials)
sa1hat = array(NA,ntrials)
sa2hat = array(NA,ntrials)
yhat= array(NA,ntrials)
exp_p= array(NA,ntrials)

mu3 = array(NA,ntrials)
sa3 = array(NA,ntrials)
pi3hat = array(NA,ntrials)
pi3 = array(NA,ntrials)
da2 = array(NA,ntrials)
r2 = array(NA,ntrials)
w2 = array(NA,ntrials)



mu2[1] = 0
sa2[1] = 2
mu1hat[1] = 0
sa1hat[1] = 0
sa2hat[1] = 0
da[1] = 0
yhat[1] = 0
exp_p[1] = 0

mu3[1] <- 0
sa3[1] <- 2
pi3hat[1] <- 0
pi3[1] <- 1/sa3[1]



da2[1] <- 0
r2[1] <- 0
w2[1] <- 0


for  (t in 2:ntrials){
  
  sa2hat[t] <- sa2[t-1]+exp(kappa*mu3[t-1]+omega)
  mu1hat[t] = sigmoid(mu2[t-1])
  
  
  sa1hat[t] = mu1hat[t]*(1-mu1hat[t])
  
  da[t] = u[t]-mu1hat[t]
  
  
  sa2[t] = 1/((1/sa2hat[t])+sa1hat[t])
  
  mu2[t] = mu2[t-1]+da[t]*sa2[t]
  
  
  da2[t] <- ((sa2[t]+(mu2[t]-mu2[t-1])^2)/(sa2[t-1]+exp(kappa*mu3[t-1]+omega)))-1
  
  r2[t] <-(exp(kappa*mu3[t-1]+omega)-sa2[t-1])/(sa2[t-1]+exp(kappa*mu3[t-1]+omega))
  
  w2[t] <-exp(kappa*mu3[t-1]+omega)/(sa2[t-1]+exp(kappa*mu3[t-1]+omega))
  
  pi3hat[t] <- 1/(sa3[t-1]+theta)
  
  pi3[t] <- pi3hat[t]+(kappa^2/2)*w2[t]*(w2[t]+r2[t]*da2[t])
  
  
  sa3[t] = 1/pi3[t]
  
  mu3[t] <- mu3[t-1]+sa3[t]*(kappa/2)*w2[t]*da2[t]
  
  
  exp_p[t] = mu1hat[t]^beta/((mu1hat[t]^beta)+(1-mu1hat[t])^(beta))
  
  
  
  yhat[t] = rbinom(n = 1,size = 1, prob = exp_p[t])
  
}

data = data.frame(u,mu1hat,da,sa1hat,mu2,sa2,sa2hat,mu3,sa3,pi3hat,w2,r2,da2,exp_p,yhat)

}

