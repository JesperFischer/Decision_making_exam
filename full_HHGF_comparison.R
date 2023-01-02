#Full comparison HGF
rtrunc <- function(n, distr, lower = -Inf, upper = Inf, ...){
  makefun <- function(prefix, FUN, ...){
    txt <- paste(prefix, FUN, "(x, ...)", sep = "")
    function(x, ...) eval(parse(text = txt))
  }
  if(length(n) > 1) n <- length(n)
  pfun <- makefun("p", distr, ...)
  qfun <- makefun("q", distr, ...)
  lo <- pfun(lower, ...)
  up <- pfun(upper, ...)
  u <- runif(n, lo, up)
  qfun(u, ...)
}


sigmoid = function(x) {
  1 / (1 + exp(-x))
}




hgf_agent_c_hier_compar = function(u,u4,mu_omega,sa_omega,mu_beta,sa_beta,mu_c,sa_c,nsubs,mu_alpha, sa_alpha){
  
  u1 = u
  
  for (i in 1:((nsubs/2)-1)){
    
    u1 = cbind(u1,u)
  }
  
  u = u1
  
  
  uu = u4
  
  for (i in 1:((nsubs/2)-1)){
    
    uu = cbind(uu,u4)
  }
  
  u4 = uu
  
  
  ntrials = 307
  
  mu2 = array(NA,c(ntrials[1],nsubs/2))
  sa2 = array(NA,c(ntrials[1],nsubs/2))
  mu1hat = array(NA,c(ntrials[1],nsubs/2))
  da = array(NA,c(ntrials[1],nsubs/2))
  sa1hat = array(NA,c(ntrials[1],nsubs/2))
  sa2hat = array(NA,c(ntrials[1],nsubs/2))
  yhat= array(NA,c(ntrials[1],nsubs/2))
  exp_p= array(NA,c(ntrials[1],nsubs/2))
  
  data1 = as.data.frame(NULL)
  
  for (s in 1:(nsubs/2)){
    alpha1 = rnorm(1, mean = mu_alpha, sd = sa_alpha)
    c1 = rnorm(1,mu_c+(alpha1/2),sa_c)
    omega1 = rnorm(1,mu_omega,sa_omega)
    beta1 = rtrunc(1, "norm", lower = 0, mean = mu_beta,sd = sa_beta)
    
    
    mu2[1,s] = 0
    sa2[1,s] = 4
    mu1hat[1,s] = 0
    sa1hat[1,s] = 0
    sa2hat[1,s] = 0
    da[1,s] = 0
    yhat[1,s] = 0
    exp_p[1,s] = 0
    
    
    for(t in 2:ntrials){
      
      
      sa2hat[t,s] = sa2[(t-1),s]+exp(omega1)
      mu1hat[t,s] = sigmoid(mu2[(t-1),s])
      
      
      sa1hat[t,s] = mu1hat[t,s]*(1-mu1hat[t,s])
      
      da[t,s] = u[t,s]-mu1hat[t,s]
      
      
      sa2[t,s] = 1/((1/sa2hat[t,s])+sa1hat[t,s])
      
      mu2[t,s] = mu2[t-1,s]+da[t,s]*sa2[t,s]
      
      
      
      exp_p[t,s] = mu1hat[t,s]^beta1/((mu1hat[t,s]^beta1)+(1-mu1hat[t,s])^(beta1+c1))
      
      
      yhat[t,s] = rbinom(n = 1,size = 1, prob = exp_p[t,s])
      
    }
    data = data.frame(u = u[,s], mu1hat = mu1hat[,s],yhat = yhat[,s],sub = rep(s,307),beta = rep(beta1,307),omega = rep(omega1, 307),c = rep(c1,307),alpha = rep(alpha1,307), seq = rep(1,307))
    colnames(data) = c("u","muhat","yhat","sub","beta","omega","c","alpha","seq")
    data1 = rbind(data1,data)
  }
  
  
  
  mu2 = array(NA,c(ntrials[1],nsubs/2))
  sa2 = array(NA,c(ntrials[1],nsubs/2))
  mu1hat = array(NA,c(ntrials[1],nsubs/2))
  da = array(NA,c(ntrials[1],nsubs/2))
  sa1hat = array(NA,c(ntrials[1],nsubs/2))
  sa2hat = array(NA,c(ntrials[1],nsubs/2))
  yhat= array(NA,c(ntrials[1],nsubs/2))
  exp_p= array(NA,c(ntrials[1],nsubs/2))
  
  
  data4 = as.data.frame(NULL)
  
  
  for (s in 1:(nsubs/2)){
    
    alpha = rtrunc(1,"norm",lower = 0, mean = mu_alpha, sd = sa_alpha)
    c4 = rnorm(1,mu_c-(alpha/2),sa_c)
    omega4 = rnorm(1,mu_omega,sa_omega)
    beta4 = rtrunc(1, "norm", lower = 0, mean = mu_beta,sd = sa_beta)
    
    
    mu2[1,s] = 0
    sa2[1,s] = 4
    mu1hat[1,s] = 0
    sa1hat[1,s] = 0
    sa2hat[1,s] = 0
    da[1,s] = 0
    yhat[1,s] = 0
    exp_p[1,s] = 0
    
    
    for(t in 2:ntrials){
      
      
      sa2hat[t,s] = sa2[(t-1),s]+exp(omega4)
      mu1hat[t,s] = sigmoid(mu2[(t-1),s])
      
      
      sa1hat[t,s] = mu1hat[t,s]*(1-mu1hat[t,s])
      
      da[t,s] = u4[t,s]-mu1hat[t,s]
      
      
      sa2[t,s] = 1/((1/sa2hat[t,s])+sa1hat[t,s])
      
      mu2[t,s] = mu2[t-1,s]+da[t,s]*sa2[t,s]
      
      
      
      exp_p[t,s] = mu1hat[t,s]^beta4/((mu1hat[t,s]^beta4)+(1-mu1hat[t,s])^(beta4+c4))
      
      
      yhat[t,s] = rbinom(n = 1,size = 1, prob = exp_p[t,s])
      
    }
    data = data.frame(u4 = u4[,s], mu1hat = mu1hat[,s],yhat = yhat[,s],sub = rep(s,307),beta = rep(beta4,307),omega = rep(omega4, 307),c = rep(c4,307),alpha = rep(alpha4,307), seq = rep(4,307))
    
    colnames(data) = c("u","muhat","yhat","sub","beta","omega","c","alpha","seq")
    data4 = rbind(data4,data)
  }
  
  data = rbind(data1,data4)
  return(data)
  
}