set.seed(12345)
#require(mvtnorm)

load(file="disease.rda")

f1 <- function(I, par){
  ((par[4]/(par[1]*I[1] + par[2]*I[2] + par[3])) + 1)*I[1]
}

f2 <- function(I,par,md){
  ((par[1]-par[2])*I + par[3]*md)/par[1]
}

makeAV_mcmc <- function(nsamp, pJ1, pA1, lamJ1, lamA1, V=1){
  temp <- data.frame(matrix(rep(0,5*nsamp), ncol = 5))
  names(temp) <- c("J","A","DJ","DA","V")
  for(i in 1:nsamp){
    # Sample for environment
    # Sample for age bins
    nJ <- rpois(1,lamJ1)
    nA <- rpois(1,lamA1)
    # Sample for diseased in each bin
    DJ <- rbinom(1,nJ, pJ1)
    DA <- rbinom(1,nA, pA1)
    temp[i,] <- c(nJ - DJ, nA - DA, DJ, DA,V)
  }
  return(temp)
}

makeMeans <- function(nsim, pJ1, pA1, lamJ1, lamA1){
  means <- matrix(rep(0,nsim*5),ncol=5)
  for(i in 1:nsim){
    temp <- makeAV_mcmc(1000, pJ1, pA1, lamJ1, lamA1)
    means[i,] <- colMeans(temp)
  }
  colnames(means) <- c("J","A","DJ","DA","V")
  return(means)
}

runMCMC <- function(nrepeat,nhouse){
  set.seed(12345)
  V=1
  hh1 = subset(disease, V == 1)
  
  pJ1 = mean(hh1$DJ) / mean(hh1$J + hh1$DJ)
  pA1 = mean(hh1$DA) / mean(hh1$A + hh1$DA)
  
  # Poisson params
  lamJ1 <- mean(hh1$J + hh1$DJ)
  lamA1 <- mean(hh1$A + hh1$DA)
  
  dat <- data.frame(makeMeans(nhouse,pJ1, pA1, lamJ1, lamA1))
  mean.DJ=lamJ1*pJ1
  mean.DA=lamA1*pA1
  nJ = (dat$J + dat$DJ)[1]
  nA = (dat$A + dat$DA)[1]
  
  d1 = dat$J 
  d2 = dat$A 
 
  #assigning hyper parameters and tunning constants
  tunning=2.4
  a_para=1.5; 
  b_para=3; 
  mean.d1=mean(d1); mean.d2=mean(d2)
  max.d1=max(d1); max.d2=max(d2)
  n.d1=length(d1); n.d2=length(d2)
  
  # propal distribution of multivariate normal for parameters set beta1, beta2, phi, and gamma 
  proposal.1 <-function(current,cov){
    res=rmvnorm(1, current, sigma=tunning*cov/4, method=c("chol"))
    min.res=min(res)
    while(min.res<0){
      res=rmvnorm(1, current, sigma=tunning*cov/4, method=c("chol"))
      min.res=min(res)  
    }
    return(res)
  }
  
  # propal distribution of multivariate normal for parameters set alpha, nu, and delta  
  proposal.2 <-function(current,cov){
    res=rmvnorm(1, current, sigma=tunning*cov/(3), method=c("chol"))
    min.res=min(res)
    while(min.res<0){
      res=rmvnorm(1, current, sigma=tunning*cov/(3), method=c("chol"))
      min.res=min(res)  
    }
    return(res)
  }
  
  # propal distribution of univariate normal for I1 or I2  
  proposal.3 <-function(current,tun){
    res=rnorm(1, current, tun*1)
    while(res<0){
      res=rnorm(1, current, tun*1)
    }
    return(res)
  }
  
  # log likelihood function for f1 and f2 in one likelihood function 
  l.likelihood<-function(m,I1,I2,par1,par2,m.D){
    mean=(f1(c(I1,I2),par1)+f2(I1,par2,m.D))/2
    sd=sigma/sqrt(n.d1)
    res=log(dnorm(m,mean,sd)+1e-320)
    return(res)
  }
  
  # condtional posterior of I1 and I2
  gibbs.I <-function(par,md1,md2,sig,mean.p,sig.p,upper){
    mean.x=par[1]/(par[1]-par[2])*(md1-par[3]/par[1]*md2)
    sig.x=((par[1]/(par[1]-par[2]))*sig)**2
    sig=1/(1/sig.p+1/sig.x)
    mean=sig*(mean.p/sig.p+mean.x/sig.x)
    res=rnorm(1,mean, sqrt(sig))
    repeat{
      res=rnorm(1,mean, sqrt(sig))
      if((res>0)&&(res<upper))break
    }
    return(res)
  }
  
  # choose every inc sample for summary 
  pick.data<-function(data,inc){
    nrepeat=nrow(data)
    new.data=NULL
    for(i in 1:nrepeat){
      if(i%%inc==0){
        i1=i/inc
        new.data=rbind(new.data,data[i,])
      }
    }
    return(new.data)
  }
  
  # summary statistics of MCMC: only use left half samples 
  stat<-function(data){
    nobs=nrow(data)
    id=1:nobs
    data=cbind(data,id)
    data=data[which(id>nobs/2),]
    MC_mean=apply(data[,1:16],2,mean)
    MC_sd=apply(data[,1:16],2,sd)
    CI=NULL
    for(l in 1:16) CI=rbind(CI,quantile(data[,l],c(0.025,0.975)))
    res=cbind(MC_mean, MC_sd, CI)
    rownames(res)=c("beta11","beta12","phi1","gamma1","beta11","beta21","phi2","gamma2","alpha1","nu1","delta1","alpha2","nu2","delta2","I1","I2")
    
    return(res)
  }
  
  
  tunning=2.4
  tun_I1=0.6
  tun_I2=0.6
  a_para=1.5; 
  b_para=3; 
  a_I1=2; b_I1=2.3;
  a_I2=2; b_I2=1;
  sigma=2.7234;
  cov.1=diag(1,4)*0.6
  cov.2=diag(1,4)*0.4
  cov.3=diag(1,3)*0.3
  cov.4=diag(1,3)*0.1

  
  
  ###################################################################################
  # MCMC procedure                                                                  # 
  ###################################################################################
  init.param=rep(0.5,14)
  init_a = 1; init_b = 2.5;
  #MCMC_gen=cbind(init.param,init_a,init_b)
  MCMC_gen=NULL
  count.1=0
  count.2=0
  count.3=0
  count.4=0
  count.I1=0
  count.I2=0
  
  current.param=init.param;
  current.param.1=init.param[1:4]
  current.param.2=init.param[5:8]
  current.param.3=init.param[9:11];current.param.3[1]=0.6;
  current.param.4=init.param[12:14];current.param.4[1]=0.6;
  current.I1=init_a
  current.I2=init_b
  
  for(k in 1:nrepeat){
    if (k%%100==0)   if (k%%1000==0) cat("|") else cat("*");
    cand.param.1=proposal.1(current.param.1,cov.1)
    l.current.param.prior.1=sum(log(dgamma(current.param.1,shape=a_para,rate=b_para)))
    l.cand.param.prior.1=sum(log(dgamma(cand.param.1,shape=a_para,rate=b_para)))
    l.current.likelihood.1=l.likelihood(mean.d1,current.I1,current.I2,current.param.1,current.param.3,mean.DJ)
    l.cand.likelihood.1   =l.likelihood(mean.d1,current.I1,current.I2,cand.param.1   ,current.param.3,mean.DJ)
    l.MH.1=l.cand.param.prior.1+l.cand.likelihood.1-l.current.param.prior.1-l.current.likelihood.1
    if(l.MH.1<0){
      MH=exp(l.MH.1)
      uni=runif(1)
      if(uni<MH){
        current.param.1=cand.param.1
        count.1=count.1+1
      }
    }else{
      current.param.1=cand.param.1
      count.1=count.1+1
    }
    
    #Gibbs for theta2
    cand.param.2=proposal.1(current.param.2,cov.2)  
    l.current.param.prior.2=sum(log(dgamma(current.param.2,shape=a_para,rate=b_para)))
    l.cand.param.prior.2=sum(log(dgamma(cand.param.2,shape=a_para,rate=b_para)))
    l.current.likelihood.2=l.likelihood(mean.d2,current.I2,current.I1,current.param.2,current.param.4,mean.DA)
    l.cand.likelihood.2   =l.likelihood(mean.d2,current.I2,current.I1,cand.param.2   ,current.param.4,mean.DA)
    l.MH.2=l.cand.param.prior.2+l.cand.likelihood.2-l.current.param.prior.2-l.current.likelihood.2
    if(l.MH.2<0){
      MH=exp(l.MH.2)
      uni=runif(1)
      if(uni<MH){
        current.param.2=cand.param.2
        count.2=count.2+1
      }
    }else{
      current.param.2=cand.param.2
      count.2=count.2+1
    }
    
    
    #Gibbs for theta3
    cand.param.3=proposal.2(current.param.3,cov.3)
    l.current.param.prior.3=sum(log(dgamma(current.param.3,shape=a_para,rate=b_para)))
    l.cand.param.prior.3=sum(log(dgamma(cand.param.3,shape=a_para,rate=b_para)))
    l.current.likelihood.3=l.likelihood(mean.d1,current.I1,current.I2,current.param.1,current.param.3,mean.DJ)
    l.cand.likelihood.3   =l.likelihood(mean.d1,current.I1,current.I2,current.param.1,   cand.param.3,mean.DJ)
    l.MH.3=l.cand.param.prior.3+l.cand.likelihood.3-l.current.param.prior.3-l.current.likelihood.3
    if(cand.param.3[1]-cand.param.3[2]>0){
      if(l.MH.3<0){
        MH=exp(l.MH.3)
        uni=runif(1)
        if(uni<MH){
          current.param.3=cand.param.3
          count.3=count.3+1
        }
      }else{
        current.param.3=cand.param.3
        count.3=count.3+1
      }      
    }
    
    #Gibbs for theta4
    cand.param.4=proposal.2(current.param.4,cov.4)
    l.current.param.prior.4=sum(log(dgamma(current.param.4,shape=a_para,rate=b_para)))
    l.cand.param.prior.4=sum(log(dgamma(cand.param.4,shape=a_para,rate=b_para)))
    l.current.likelihood.4=l.likelihood(mean.d2,current.I2,current.I1,current.param.2,current.param.4,mean.DA)
    l.cand.likelihood.4   =l.likelihood(mean.d2,current.I2,current.I1,current.param.2,   cand.param.4,mean.DA)
    l.MH.4=l.cand.param.prior.4+l.cand.likelihood.4-l.current.param.prior.4-l.current.likelihood.4
    if(cand.param.4[1]-cand.param.4[2]>0){
      if(l.MH.4<0){
        MH=exp(l.MH.4)
        uni=runif(1)
        if(uni<MH){
          current.param.4=cand.param.4
          count.4=count.4+1
        }
      }else{
        current.param.4=cand.param.4
        count.4=count.4+1
      }
    }
    
    
    #Gibbs for I1
    cand.I1=proposal.3(current.I1,tun_I1)
    l.current.prior.I1=log(dgamma(current.I1,shape=a_I1,rate=b_I1))
    l.cand.prior.I1   =log(dgamma(cand.I1   ,shape=a_I1,rate=b_I1))
    l.current.likelihood.I1=l.likelihood(mean.d1,current.I1,current.I2,current.param.1,current.param.3,mean.DJ)
    l.cand.likelihood.I1   =l.likelihood(mean.d1,   cand.I1,current.I2,current.param.1,current.param.3,mean.DJ)
    l.MH.I1=l.cand.prior.I1+l.cand.likelihood.I1-l.current.prior.I1-l.current.likelihood.I1
    if(cand.I1<max.d1){
      if(l.MH.I1<0){
        MH=exp(l.MH.I1)
        uni=runif(1)
        if(uni<MH){
          current.I1=cand.I1
          count.I1=count.I1+1
        }
      }else{
        current.I1=cand.I1
        count.I1=count.I1+1
      }
    }
    
    #Gibbs for I1
    cand.I2=proposal.3(current.I2,tun_I2)
    l.current.prior.I2=log(dgamma(current.I2,shape=a_I2,rate=b_I2))
    l.cand.prior.I2   =log(dgamma(cand.I2   ,shape=a_I2,rate=b_I2))
    l.current.likelihood.I2=l.likelihood(mean.d2,current.I2,current.I1,current.param.2,current.param.4,mean.DA)
    l.cand.likelihood.I2   =l.likelihood(mean.d2,   cand.I2,current.I1,current.param.2,current.param.4,mean.DA)
    l.MH.I2=l.cand.prior.I2+l.cand.likelihood.I2-l.current.prior.I2-l.current.likelihood.I2
    if(cand.I2<max.d2){
      if(l.MH.I2<0){
        MH=exp(l.MH.I2)
        uni=runif(1)
        if(uni<MH){
          current.I2=cand.I2
          count.I2=count.I2+1
        }
      }else{
        current.I2=cand.I2
        count.I2=count.I2+1
      }
    }
    
    MCMC_gen=rbind(MCMC_gen,c(current.param.1,current.param.2,current.param.3,current.param.4,current.I1,current.I2))  
  }
  colnames(MCMC_gen)=c("beta11","beta12","phi1","gamma1","beta11","beta21","phi2","gamma2","alpha1","nu1","delta1","alpha2","nu2","delta2","I1","I2")
  count=c(count.1, count.2, count.3, count.4, count.I1, count.I2)/nrepeat
  cat("\n")
  print(count)
  return(MCMC_gen)
}

result <- runMCMC(nrepeat = 40000, nhouse = 100)

