

sim.step=function(threshold.type=c("sigmoid","step","quadratic"),X.ditr=c('norm',"unif"),thres='NA',shape="NA",mu.x,mu.z,sd.x,sd.z,coef.X,eps.sd,seed,n){
  mu.X<-as.matrix(c(mu.z,mu.x)) # mean X =4.7, mean z=0
  cov.X <- diag(c(sd.z^2,sd.x^2)) # sd.z=1, sd.x=1.6
  npar=length(mu.X) # X including Z and x, it is a vector. If npar=8, x is the 8th, z are 1-7
  beta=coef.X
  set.seed(seed)
  if(X.ditr=='norm'){
    allcovariates=mvrnorm(n, mu=mu.X, Sigma=cov.X, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    z <- allcovariates[,-npar] # returning all but he npar th column==all zs
    x <- allcovariates[,npar] # returning all x
    if(threshold.type=='step'){
      e <- thres
      xe<-ifelse(x>e,1,0)
      Xe <- cbind(1,z,xe) # Z including 1 and z # Xe including (1,z,xe) 
      Y=Xe %*% beta+rnorm(n,0,eps.sd)
    } else if (threshold.type=="sigmoid"){
      x.new <- exp(shape * (x-thres))/(1+exp(shape * (x-thres)))
      Xe <- cbind(1,z,x.new)
      Y <- Xe %*% beta+rnorm(n,0,eps.sd)
    }
    dat=data.frame(Y,z,x)
    return(dat)
  } else if (X.ditr=="unif"){
    x=runif(n)*4*sd.x + mu.x-2*sd.x # (0,6.4)+1.5=uniform (1.5,7.9), mid: 4.7 sd.x=1.6
    z=rnorm(n, mean=mu.z, sd=sd.z)
    if(threshold.type=='step'){
      e <- thres
      xe<-ifelse(x>e,1,0)
      Xe <- cbind(1,z,xe) # Z including 1 and z # Xe including (1,z,xe)
      Y=Xe %*% beta+rnorm(n,0,eps.sd)
    } else if (threshold.type=="sigmoid"){
      x.new <- exp(shape * (x-thres))/(1+exp(shape * (x-thres)))
      Xe <- cbind(1,z,x.new)
      Y <- Xe %*% beta+rnorm(n,0,eps.sd)
    } else if (threshold.type=="quadratic"){
      x.quad <- x*x
      Xe <- cbind(1,z,x,x.quad)
      beta.quad <- as.matrix(c(-1,log(1.4),-1,0.3))
      Y <- Xe %*% beta.quad+rnorm(n,0,eps.sd)
    }
    dat=data.frame(Y,z,x)
    return(dat)
  }
}


#coef.X="null"
#dat <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=n)
