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

load("s3-sams.RData")
load("s3-ini-data.RData")

################################################################################
##Setting#########################################################################
x.sp <- function(y_row,s0){
  fit <- smooth.spline(s.sp,y_row)
  predict(fit,s0)$y
}

lower.s=0; upper.s=1; lower.t=0; upper.t=1; lower.x = c(0.1,0.1,0.1); upper.x = c(0.9,0.9,0.9)

n.sam=100; n0=10 
t.mc = sams.re$t.mc
t.l2 = sams.re$t.l2

s.sp = sams.re$s.sp
s.test = t.test = sams.re$s.test

s.dis = sams.re$s.dis
t.dis = sams.re$t.dis


################################################################################
##Model#########################################################################
x.para <- function(s0, theta) {
  theta[1] * sin(2 * pi * s0) +
    theta[2] * cos(2 * pi * s0) +
    theta[3] * exp(-5*(s0-0.5)^2)
}
theta.star <- c(1/2,1/3,1/4)
x.star <- function(s0) x.para(s0, theta.star)

f.obj <- function(x, t0) {
  dist <- adaptIntegrate(function(s0) (x(s0) - x.star(s0))^2, 0, 1)$integral
  asym <- adaptIntegrate(function(s0) x(s0) * sin(3 * pi * s0), 0, 1)$integral
  return(20 * exp(-5 * dist) + 10 * sin(3 * pi * t0) * asym)
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
plot(t.test, apply(t.test,1,function(t0) f.obj(function(s0) x.para(s0,c(0.1,0.2,0.3)),t0)),
     type="l", lwd=3, ylab="f")
par(opar)

x.star.dis = t(as.matrix(apply(s.dis,1,x.star)))
g.obj.dis <- function(y) mean(y)
g.obj.dis.figp <- function(x) mean(apply(t.dis,1,function(t0) f.obj(x,t0)))



################################################################################
## training figp, optimizing parameters
lower.th.figp=c(1e-3,1e-2,1e-10); upper.th.figp=c(30,10,1)

# G.figp <- function(x1,x2,n1,n2,s) (matrix(x1(s),n1,n2,byrow=F)-matrix(x2(s),n1,n2,byrow=T))^2
# dis.figp <- function(x1,x2,n1,n2) matrix(adaptIntegrate(function(s0) G.figp(x1,x2,n1,n2,s0),lower.s,upper.s,fDim = n1*n2)$integral,n1,n2)
dis.figp <- function(x1,x2) as.matrix((dist(x1,x2,method = "Euclidean"))^2)

k.x.figp <- function(dis.x,tau.x) matern(sqrt(dis.x),tau.x,5/2)

like.figp <- function(dis.obs,g.obs,n,th){
  mu.prior = mean(g.obs)
  sigma2=th[1]; tau.x=th[2]; lambda=th[3]
  
  R1 = matrix(sapply(dis.obs,function(dis.x) k.x.figp(dis.x,tau.x)),n,n)
  Sigma.g = sigma2*R1+lambda*diag(n)
  
  i11 = determinant(Sigma.g,logarithm=TRUE)$modulus
  i22 = t(g.obs-mu.prior)%*%solve(Sigma.g)%*%(g.obs-mu.prior)
  
  likeli = i11+i22
  return(likeli)
}

## predict figp
pre.g.figp <- function(dis.obs,dis.pre,g.obs,n,n.test,th){
  mu.prior = mean(g.obs)
  sigma2=th[1]; tau.x=th[2]; lambda=th[3]
  
  R1 = matrix(sapply(dis.obs,function(dis.x) k.x.figp(dis.x,tau.x)),n,n)
  R.pre = matrix(sapply(dis.pre,function(dis.x) k.x.figp(dis.x,tau.x)),n,n.test)
  
  Sigma.g = solve(R1+lambda/sigma2*diag(n))
  
  mu.g = mu.prior+t(R.pre)%*%Sigma.g%*%(g.obs-mu.prior)
  cov.g = sigma2*(1-t(R.pre)%*%Sigma.g%*%R.pre)
  
  result = list(mu.g=mu.g, cov.g=cov.g)
  return(result)
}


############## FFBO ##################
J.figp = 10
x.for.figp = g.for.figp = list()
hyper.for.figp = lapply(1:J.figp, function(x) list())
beta.figp = ucb.for.figp = lapply(1:J.figp, function(x) list())
time.for.figp = list()

###
n0=10
b.ini <- ini.data$b.ini 
x.ini <- ini.data$x.ini 
g.ini.for.dis <- ini.data$g.ini.for.dis

delta=0.05
for(j.figp in 1:J.figp){
  ptm = proc.time()
  n=10; m=20
  
  x.for.figp[[j.figp]] = apply(s.sp,1,function(s0) unlist(x.ini(s0)))
  g.for.figp[[j.figp]] = g.ini.for.dis[[j.figp]]

  opar<-par(no.readonly = TRUE)
  par(mfrow = c(2, 2))
  par(mar = c(4.5,4.5,2,2),xpd=TRUE)
  for(i in 1:n){
    plot(s.sp, x.for.figp[[j.figp]][i,],type="l", lwd=3, ylab="x")
  }
  plot(g.for.figp[[j.figp]],type="l", lwd=3, ylab="f")
  
  dis.obs = dis.figp(x.for.figp[[j.figp]],x.for.figp[[j.figp]])
  hyper.for.figp.old = directL(function(th) like.figp(dis.obs,
                                                      g.for.figp[[j.figp]],n,th),lower.th.figp,upper.th.figp,control=list(xtol_rel=1e-8, maxeval=1000))$par
  hyper.for.figp[[j.figp]][[1]] = bobyqa(hyper.for.figp.old, function(th) like.figp(dis.obs,
                                                                                    g.for.figp[[j.figp]],n,th),lower=lower.th.figp,upper=upper.th.figp)$par
  hyper.figp.new = hyper.for.figp[[j.figp]][[1]]
  
  for(i.figp in 1:m){
    beta.figp[[j.figp]][[i.figp]] = 2*log(i.figp^2*pi^2/3)*sqrt(log(2/delta))
    
    ## construct ucb criteria
    ucb.figp <- function(dis.pre){
      sigma2=hyper.figp.new[1]; tau.x=hyper.figp.new[2]; lambda=hyper.figp.new[3]
      
      R1 = matrix(sapply(dis.obs,function(dis.x) k.x.figp(dis.x,tau.x)),n,n)
      R.pre = matrix(sapply(dis.pre,function(dis.x) k.x.figp(dis.x,tau.x)),n,1)
      
      sol.Sigma.g = solve(R1+lambda/sigma2*diag(n))
      al.1 = sol.Sigma.g%*%g.for.figp[[j.figp]]
      
      mu.g = t(R.pre)%*%al.1
      cov.g = sigma2*(1-t(R.pre)%*%sol.Sigma.g%*%R.pre)+1e-8
      
      ucb = mu.g+sqrt(beta.figp[[j.figp]][[i.figp]]*cov.g/n^2)
      d.ucb.k = al.1-as.numeric(sqrt(beta.figp[[j.figp]][[i.figp]]/(n^2*cov.g)))*sol.Sigma.g%*%R.pre
      
      result = list(ucb=ucb, d.ucb.k=d.ucb.k)
      return(result)
    }
    
    Ta.figp = 100; eps.figp=1e-8; ucb.val.figp = vector()
    x.old.figp = t(as.matrix(x.for.figp[[j.figp]][which.max(g.for.figp[[j.figp]]),]-1e-2))
    
    dis.pre.0 = dis.figp(x.for.figp[[j.figp]],x.old.figp)
    re.ucb.val.figp = ucb.figp(dis.pre.0)
    ucb.val.figp[1] = re.ucb.val.figp$ucb
    
    d.k.dis.figp = jacobian(function(dis) k.x.figp(dis,hyper.figp.new[2]), dis.pre.0)
    d.ucb.dis.figp = d.k.dis.figp%*%re.ucb.val.figp$d.ucb.k
    
    d.dis.x.figp <- function(x) 2*(matrix(x,n,50)-x.for.figp[[j.figp]])
    d.ucb.x.figp <- function(x) t(d.ucb.dis.figp)%*%as.matrix(d.dis.x.figp(x))
    
    #########FGA Algorithm
    for(l.figp in 1:Ta.figp){
      # grad.old = apply(s.sp,1,function(s0) d.ucb.x.figp(x.old.figp,s0)^2)
      sr.figp = 0.00005/sqrt(l.figp) #/sqrt(mean(grad.old))
      
      x.new.figp = x.old.figp+sr.figp*d.ucb.x.figp(x.old.figp)
      
      ##sparsify
      x.new.hat.figp <- function(s0) x.sp(x.new.figp,s0)
      plot(s.sp,apply(s.sp,1,x.new.hat.figp),type="l",lwd=2,col=1)
      lines(s.sp,apply(s.sp,1,x.star),type="l",lwd=2,col=3)
      
      dis.pre.0 = dis.figp(x.for.figp[[j.figp]],x.new.figp)
      
      re.ucb.val.figp = ucb.figp(dis.pre.0)
      ucb.val.figp[l.figp+1] = re.ucb.val.figp$ucb
      d.k.dis.figp = jacobian(function(dis) k.x.figp(dis,hyper.figp.new[2]), dis.pre.0)
      d.ucb.dis.figp = d.k.dis.figp%*%re.ucb.val.figp$d.ucb.k
      print(paste(l.figp,ucb.val.figp[l.figp+1]))
      
      x.old.figp <- x.new.figp
      if((ucb.val.figp[l.figp+1]-ucb.val.figp[l.figp])^2<eps.figp) break
    }
    plot(ucb.val.figp,lwd=2,type="l",col=4,ylab="UCB",pch=3)
    
    ## Update design
    g.new.sp.figp = g.obj.dis.figp(x.new.hat.figp)
    
    x.for.figp[[j.figp]] = rbind(x.for.figp[[j.figp]],x.new.figp)
    g.for.figp[[j.figp]] = c(g.for.figp[[j.figp]],g.new.sp.figp)
    
    n=n+1
    dis.obs = dis.figp(x.for.figp[[j.figp]],x.for.figp[[j.figp]])
    hyper.figp.new = bobyqa(hyper.figp.new, function(th) like.figp(dis.obs,
                                                                   g.for.figp[[j.figp]],n,th),lower=lower.th.figp,upper=upper.th.figp)$par
    hyper.for.figp[[j.figp]][[i.figp+1]] = hyper.figp.new
    
    print(paste("rep:", j.figp, ", sequential:", i.figp))
  }
  time.figp = (proc.time()-ptm)[1]
  time.for.figp[[j.figp]] = time.figp
}


re.figp = list(x.for.figp=x.for.figp, g.for.figp=g.for.figp,
               hyper.for.figp=hyper.for.figp, beta.figp=beta.figp, ucb.for.figp=ucb.for.figp,
               time.for.figp=time.for.figp) 
save(re.figp, file="s3-re-figp.RData")  





