getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(cubature)
library(lava) 
library(ggplot2)
library(mcGlobaloptim) 
library(DiceKriging)
library(nloptr)
library(MASS)
library(mcmc)
library(geoR)
library(RobustCalibration)
library(lhs)
library(RobustGaSP)
library(numDeriv)
library(magrittr)
library(FRegSigCom)
library(proxy)

load("s2-sams.RData")
load("s2-ini-data.RData")

################################################################################
##Setting#########################################################################
x.sp <- function(y_row,s0){
  fit <- smooth.spline(s.sp,y_row)
  predict(fit,s0)$y
}

lower.s=0; upper.s=1; lower.t=0; upper.t=1; lower.x = c(0.5,0.5); upper.x = c(1.5,1.5)

n.sam=100; n0=10 
t.mc = sams.re$t.mc
t.l2 = sams.re$t.l2

s.sp = sams.re$s.sp
s.test = t.test = sams.re$s.test

s.dis = sams.re$s.dis
t.dis = sams.re$t.dis

################################################################################
##Model#########################################################################
x.para <- function(s0,a,b) a*cos(pi/(b+sin(exp(s0)+pi*s0)^2))
x.star <- function(s0) cos(pi/(1+sin(exp(s0)+pi*s0)^2))

f.obj <- function(x,t0){
  it1 = adaptIntegrate(function(s0) (x(s0)-x.star(s0))^2,lower.s,upper.s)$integral
  return(-3*exp(t0*it1))
} 

err <- function(t0) {
  set.seed(as.integer(t0))
  gp.sam <- rnorm(1,0,1e-4)
  return(gp.sam)
}

weight <- function(t0) 1
g.obj <- function(y) adaptIntegrate(y,lower.t,upper.t)$integral

opar<-par(no.readonly = TRUE)
par(mfrow = c(1, 2))
par(mar = c(4.5,4.5,2,2),xpd=TRUE)
plot(s.test, apply(s.test,1,x.star),type="l", lwd=3, ylab="x")
plot(t.test, apply(t.test,1,function(t0) f.obj(function(s0) x.para(s0,1/3,1/4),t0)),
     type="l", lwd=3, ylab="f")
par(opar)

x.star.dis = t(as.matrix(apply(s.dis,1,x.star)))
g.obj.dis <- function(y) mean(y)


################################################################################
## training mtgp, optimizing parameters
lower.th=c(1e-3,1e-2,1e-2,1e-10); upper.th=c(30,10,2,1e-4)

dis.mtgp <- function(x1,x2) as.matrix((dist(x1,x2,method = "Euclidean"))^2)

k.x.mtgp <- function(dis.x,tau.x) matern(sqrt(dis.x),tau.x,5/2)
k.y.mtgp <- function(tau.y) outer(as.vector(t.dis),as.vector(t.dis),function(t1,t2) exp(-abs(t1-t2)/tau.y))

like.mtgp <- function(dis.obs,y.obs,n,th){
  sigma2=th[1]; tau.x=th[2]; tau.y=th[3]; lambda=th[4]
  
  R1 = matrix(sapply(dis.obs,function(dis.x) k.x.mtgp(dis.x,tau.x)),n,n)
  eigen.x <- eigen(R1); ev.x <- eigen.x$values; evec.x <- eigen.x$vectors
  
  K.y = k.y.mtgp(tau.y)
  eigen.y <- eigen(K.y); ev.y <- eigen.y$values; evec.y <- eigen.y$vectors
  
  ev.sig = kronecker(evec.x,evec.y)
  
  Sigma = ev.sig%*%diag((kronecker(ev.x,ev.y)+lambda/sigma2))%*%t(ev.sig)
  Sol.Sigma = ev.sig%*%diag(1/(kronecker(ev.x,ev.y)+lambda/sigma2))%*%t(ev.sig)
  
  likeli = determinant(Sigma,logarithm=TRUE)$modulus+t(y.obs)%*%Sol.Sigma%*%y.obs
  return(likeli)
}


## predict mtgp
pre.f.mtgp <- function(dis.obs,dis.pre,y.obs,n,n.test,th){
  mu.prior = mean(y.obs)
  
  sigma2=th[1]; tau.x=th[2]; tau.y=th[3]; lambda=th[4]
  
  R1 = matrix(sapply(dis.obs,function(dis.x) k.x.mtgp(dis.x,tau.x)),n,n)
  eigen.x <- eigen(R1); ev.x <- eigen.x$values; evec.x <- eigen.x$vectors
  
  K.y = k.y.mtgp(tau.y)
  eigen.y <- eigen(K.y); ev.y <- eigen.y$values; evec.y <- eigen.y$vectors
  
  ev.sig = kronecker(evec.x,evec.y)
  
  Sigma = ev.sig%*%diag((kronecker(ev.x,ev.y)+lambda/sigma2))%*%t(ev.sig)
  Sol.Sigma = ev.sig%*%diag(1/(kronecker(ev.x,ev.y)+lambda/sigma2))%*%t(ev.sig)
  
  R.pre = matrix(sapply(dis.pre,function(dis.x) k.x.mtgp(dis.x,tau.x)),n,n.test)
  K.pre = kronecker(R.pre,K.y)
  
  mu.f = mu.prior+t(K.pre)%*%Sol.Sigma%*%(y.obs-mu.prior)
  cov.f = sigma2*(K.y-t(K.pre)%*%Sol.Sigma%*%K.pre)
  
  mu.g = mean(mu.f)
  cov.g = abs(tr(cov.f))
  
  result = list(mu.f=mu.f, cov.f=cov.f, mu.g=mu.g, cov.g=cov.g)
  return(result)
}



############## FFBO ##################
J.mtgp = 10
y.for.mtgp = x.for.mtgp = g.for.mtgp = list()
f.for.mtgp = lapply(1:J.mtgp, function(x) list())
hyper.for.mtgp = lapply(1:J.mtgp, function(x) list())
beta.mtgp = ucb.for.mtgp = lapply(1:J.mtgp, function(x) list())
time.for.mtgp = list()

n0=10
x.ini.dis = ini.data$x.ini.dis
y.ini.for.dis = ini.data$y.ini.for.dis
g.ini.for.dis = ini.data$g.ini.for.dis


# j.mtgp = 1;
delta=0.05
for(j.mtgp in 1:J.mtgp){
  ptm = proc.time()
  n=10; m=20
  
  x.for.mtgp[[j.mtgp]] = x.ini.dis
  g.for.mtgp[[j.mtgp]] = g.ini.for.dis[[j.mtgp]]
  y.for.mtgp[[j.mtgp]] = y.ini.for.dis[[j.mtgp]]
  
  opar<-par(no.readonly = TRUE)
  par(mfrow = c(2, 2))
  par(mar = c(4.5,4.5,2,2),xpd=TRUE)
  for(i in 1:n){
    plot(s.dis, x.for.mtgp[[j.mtgp]][i,],type="l", lwd=3, ylab="x")
    plot(t.dis, y.for.mtgp[[j.mtgp]][i,],type="l", lwd=3, ylab="f")
  }
  
  dis.obs = dis.mtgp(x.for.mtgp[[j.mtgp]],x.for.mtgp[[j.mtgp]])
  
  hyper.for.mtgp.old = directL(function(th) like.mtgp(dis.obs,
                                                      as.matrix(vec(t(y.for.mtgp[[j.mtgp]]))),n,th),lower.th,upper.th,control=list(xtol_rel=1e-8, maxeval=1000))$par
  hyper.for.mtgp[[j.mtgp]][[1]] = bobyqa(hyper.for.mtgp.old, function(th) like.mtgp(dis.obs,
                                                                                    as.matrix(vec(t(y.for.mtgp[[j.mtgp]]))),n,th),lower=lower.th,upper=upper.th)$par
  hyper.mtgp.new = hyper.for.mtgp[[j.mtgp]][[1]]
  
  x.old.mtgp = x.for.mtgp[[j.mtgp]][which.max(g.for.mtgp[[j.mtgp]]),]
  
  
  for(i.mtgp in 1:m){
    beta.mtgp[[j.mtgp]][[i.mtgp]] = 2*log(i.mtgp^2*pi^2/3)*sqrt(log(2/delta))
    
    ## construct ucb criteria
    ucb.mtgp <- function(x){
      x = t(as.matrix(x))
      dis.pre = dis.mtgp(x.for.mtgp[[j.mtgp]],x)
      result = pre.f.mtgp(dis.obs,dis.pre,as.matrix(vec(t(y.for.mtgp[[j.mtgp]]))),n,1,hyper.mtgp.new)
      mu.g=result$mu.g; cov.g=result$cov.g
      ucb = mu.g+sqrt(beta.mtgp[[j.mtgp]][[i.mtgp]]*cov.g/n)
      return(ucb)
    }
    
    #########L-BFGS Algorithm
    x.new.mtgp = bobyqa(x.for.mtgp[[j.mtgp]][which.max(g.for.mtgp[[j.mtgp]]),],
                        ucb.mtgp,lower=rep(-2,20),upper=rep(2,20))$par
    
    plot(s.sp,apply(s.sp,1,x.star),type="l",lwd=2,col=1)
    lines(s.dis,x.new.mtgp,type="l",lwd=2,col=2)
    
    x.old.mtgp = x.new.mtgp
    fit.mtgp <- smooth.spline(s.dis,x.new.mtgp)
    x.sp.mtgp <- function(s0) predict(fit.mtgp,s0)$y
    
    ## Update design
    y.new.sp.mtgp = apply(t.dis,1,function(t0) f.obj(x.sp.mtgp,t0)+err(t0*j.mtgp*i.mtgp))
    g.new.sp.mtgp = g.obj.dis(y.new.sp.mtgp)
    
    x.for.mtgp[[j.mtgp]] = rbind(x.for.mtgp[[j.mtgp]],x.new.mtgp)
    y.for.mtgp[[j.mtgp]] = rbind(y.for.mtgp[[j.mtgp]],y.new.sp.mtgp)
    g.for.mtgp[[j.mtgp]] = c(g.for.mtgp[[j.mtgp]],g.new.sp.mtgp)
    
    n=n+1
    dis.obs = dis.mtgp(x.for.mtgp[[j.mtgp]],x.for.mtgp[[j.mtgp]])
    
    if(i.mtgp %% 10 == 0){
      hyper.mtgp.new = bobyqa(hyper.mtgp.new, function(th) like.mtgp(dis.obs,
                                                                     as.matrix(vec(t(y.for.mtgp[[j.mtgp]]))),n,th),lower=lower.th,upper=upper.th)$par
    }else{
      hyper.mtgp.new = hyper.mtgp.new
    }
    hyper.for.mtgp[[j.mtgp]][[i.mtgp+1]] = hyper.mtgp.new
    
    print(i.mtgp)
  }
  
  time.mtgp = (proc.time()-ptm)[1]
  time.for.mtgp[[j.mtgp]] = time.mtgp
}


re.mtgp = list(x.for.mtgp=x.for.mtgp, y.for.mtgp=y.for.mtgp, g.for.mtgp=g.for.mtgp,
               hyper.for.mtgp=hyper.for.mtgp, beta.mtgp=beta.mtgp, ucb.for.mtgp=ucb.for.mtgp,
               time.for.mtgp=time.for.mtgp) 
save(re.mtgp, file="s2-re-mtgp.RData")  







