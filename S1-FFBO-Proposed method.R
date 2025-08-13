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

load("s1-sams.RData")
load("s1-ini-data.RData")

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
  return(-(3+sin(2*pi*t0))*(it1))
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


################################################################################
## training ffgp, optimizing parameters
k = 17
lower.th=c(1e-3,1e-2,1e-2,1e-10); upper.th=c(30,10,2,1)

G.ffgp <- function(x1,x2,n1,n2,s) (matrix(x1(s),n1,n2,byrow=F)-matrix(x2(s),n1,n2,byrow=T))^2
dis.ffgp <- function(x1,x2,n1,n2) matrix(adaptIntegrate(function(s0) G.ffgp(x1,x2,n1,n2,s0),lower.s,upper.s,fDim = n1*n2)$integral,n1,n2)+2e-5

k.x.ffgp <- function(dis.x,tau.x) matern(sqrt(dis.x),tau.x,5/2)
k.y.ffgp <- function(x1,x2,tau.y) exp(-abs(x1-x2)/tau.y)
t.y.ffgp <- function(tau.y){
  f <- function(x) (1-x^2*tau.y^2)*sin(x)+2*tau.y*x*cos(x)
  u = vector()
  for(i in 1:k){
    result1 <- uniroot(f, interval = c(pi*(i-1)+0.01, pi*(i)-0.01))
    u[i] = result1$root
  }
  beta = sort(as.matrix(2*tau.y/(1+u^2*tau.y^2)),decreasing=T)
  u = u[order(as.matrix(2*tau.y/(1+u^2*tau.y^2)),decreasing=T)]
  
  w1<-function(t0) as.matrix(u*tau.y*cos(u*t0)+sin(u*t0))
  c1=sqrt(adaptIntegrate(function(t0) (w1(t0))^2,lower.t,upper.t,fDim=k)$integral)
  w<-function(t) w1(t)/c1
  
  result = list(beta=beta,w=w)
  return(result)
}

like.ffgp <- function(dis.obs,g.obs,n,th){
  mu.prior = mean(g.obs)
  sigma2=th[1]; tau.x=th[2]; tau.y=th[3]; lambda=th[4]
  c = 2*integrate(function(u)(1-u)*exp(-u/tau.y),lower=lower.t,upper=upper.t)$value
  
  R1 = matrix(sapply(dis.obs,function(dis.x) k.x.ffgp(dis.x,tau.x)),n,n)
  Sigma.g = sigma2*c*R1+lambda*diag(n)
  
  i11 = determinant(Sigma.g,logarithm=TRUE)$modulus
  i22 = t(g.obs-mu.prior)%*%solve(Sigma.g)%*%(g.obs-mu.prior)
  
  likeli = i11+i22
  return(likeli)
}

## predict FFGP
pre.f.ffgp <- function(dis.obs,dis.pre,y.obs,n,n.test,th){
  mu.prior <- function(t0) mean(y.obs(t0))
  sigma2=th[1]; tau.x=th[2]; tau.y=th[3]; lambda=th[4]
  
  re.t.y = t.y.ffgp(tau.y)
  beta = re.t.y$beta; ev.t = sum(beta); w = re.t.y$w
  
  w.y <- function(t0) (y.obs(t0)-mu.prior(t0))%*%t(w(t0))
  w.1.y.sum = lapply(t.test, w.y)
  w.1.y = Reduce(`+`, w.1.y.sum)/length(t.test)
  
  R1 = matrix(sapply(dis.obs,function(dis.x) k.x.ffgp(dis.x,tau.x)),n,n)
  eigen.x <- eigen(R1); ev.x <- eigen.x$values; evec.x <- eigen.x$vectors
  
  theta = kronecker(ev.x,beta)
  
  w.1.1 = t(evec.x)%*%w.1.y; w.2 = matrix(1/(theta+lambda/sigma2),n,k,byrow=T)
  w.u = as.matrix(w.1.1*w.2)
  it11 <- function(t0) evec.x%*%w.u%*%as.matrix(beta*w(t0))
  
  R.pre = matrix(sapply(dis.pre,function(dis.x) k.x.ffgp(dis.x,tau.x)),n,n.test)
  y.pre.test <- function(t0) t(R.pre)%*%it11(t0)+mu.prior(t0)
  
  result = list(y.pre.test=y.pre.test)
  return(result)
}

pre.g.ffgp <- function(dis.obs,dis.pre,g.obs,n,n.test,th){
  mu.prior = mean(g.obs)
  sigma2=th[1]; tau.x=th[2]; tau.y=th[3]; lambda=th[4]
  c = 2*tau.y*((upper.t-lower.t)-tau.y*(1-exp(-(upper.t-lower.t)/tau.y)))
  
  R1 = matrix(sapply(dis.obs,function(dis.x) k.x.ffgp(dis.x,tau.x)),n,n)
  R.pre = matrix(sapply(dis.pre,function(dis.x) k.x.ffgp(dis.x,tau.x)),n,n.test)
  
  Sigma.g = solve(R1+lambda/(sigma2*c)*diag(n))
  
  mu.g = mu.prior+t(R.pre)%*%Sigma.g%*%(g.obs-mu.prior)
  cov.g = (sigma2*c)*(1-t(R.pre)%*%Sigma.g%*%R.pre)
  
  result = list(mu.g=mu.g, cov.g=cov.g)
  return(result)
}


############## FFBO ##################
J.ffgp = 10
y.for.ffgp = x.for.ffgp = g.for.ffgp = list()
hyper.for.ffgp = lapply(1:J.ffgp, function(x) list())
beta.ffgp = ucb.for.ffgp = lapply(1:J.ffgp, function(x) list())
time.for.ffgp = list()

###
n0=10
b.ini <- ini.data$b.ini 
x.ini <- ini.data$x.ini 
y.ini.for <- ini.data$y.ini.for
g.ini.for <- ini.data$g.ini.for


delta=0.05
for(j.ffgp in 1:J.ffgp){
  i.y = j.ffgp
  ptm = proc.time()
  n=10; m=20
  
  x.for.ffgp[[j.ffgp]] = apply(s.sp,1,function(s0) unlist(x.ini(s0)))
  g.for.ffgp[[j.ffgp]] = g.ini.for[[j.ffgp]]
  y.for.ffgp[[j.ffgp]] = y.ini.for[[j.ffgp]]
  
  x.da.ffgp <- function(s0) apply(x.for.ffgp[[j.ffgp]],1,function(y.row) x.sp(y.row,s0))
  
  opar<-par(no.readonly = TRUE)
  par(mfrow = c(2, 2))
  par(mar = c(4.5,4.5,2,2),xpd=TRUE)
  y.ini.test = apply(t.test,1,function(t0) unlist(y.for.ffgp[[j.ffgp]](t0)))
  for(i in 1:n){
    plot(s.sp, x.for.ffgp[[j.ffgp]][i,],type="l", lwd=3, ylab="x")
    plot(t.test, y.ini.test[i,],type="l", lwd=3, ylab="f")
  }
  
  dis.obs = dis.ffgp(x.da.ffgp,x.da.ffgp,n,n)
  hyper.for.ffgp.old = directL(function(th) like.ffgp(dis.obs,
                                                      g.for.ffgp[[j.ffgp]],n,th),lower.th,upper.th,control=list(xtol_rel=1e-8, maxeval=1000))$par
  hyper.for.ffgp[[j.ffgp]][[1]] = bobyqa(hyper.for.ffgp.old, function(th) like.ffgp(dis.obs,
                                                                                    g.for.ffgp[[j.ffgp]],n,th),lower=lower.th,upper=upper.th)$par
  hyper.ffgp.new = hyper.for.ffgp[[j.ffgp]][[1]]
  
  for(i.ffgp in 1:m){
    beta.ffgp[[j.ffgp]][[i.ffgp]] = 2*log(i.ffgp^2*pi^2/3)*sqrt(log(2/delta))
    
    ## construct ucb criteria
    ucb.ffgp <- function(dis.pre){
      mu.prior = mean(g.for.ffgp[[j.ffgp]])
      sigma2=hyper.ffgp.new[1]; tau.x=hyper.ffgp.new[2]
      tau.y=hyper.ffgp.new[3]; lambda=hyper.ffgp.new[4]
      c = 2*tau.y*((upper.t-lower.t)-tau.y*(1-exp(-(upper.t-lower.t)/tau.y)))
      
      R1 = matrix(sapply(dis.obs,function(dis.x) k.x.ffgp(dis.x,tau.x)),n,n)
      R.pre = matrix(sapply(dis.pre,function(dis.x) k.x.ffgp(dis.x,tau.x)),n,1)
      
      sol.Sigma.g = solve(R1+lambda/(sigma2*c)*diag(n))
      al.1 = sol.Sigma.g%*%(g.for.ffgp[[j.ffgp]]-mu.prior)
      
      mu.g = mu.prior+t(R.pre)%*%al.1
      cov.g = abs((sigma2*c)*(1-t(R.pre)%*%sol.Sigma.g%*%R.pre)+1e-8)
      
      ucb = mu.g+sqrt(beta.ffgp[[j.ffgp]][[i.ffgp]]*cov.g/n^2)
      d.ucb.k = al.1-as.numeric(sqrt(beta.ffgp[[j.ffgp]][[i.ffgp]]/(n^2*cov.g)))*sol.Sigma.g%*%R.pre
      
      result = list(ucb=ucb, d.ucb.k=d.ucb.k)
      return(result)
    }
    
    Ta.ffgp = 50; eps.ffgp=1e-8; ucb.val.ffgp = vector()
    x.old.ffgp <- function(s0) x.da.ffgp(s0)[[which.max(g.for.ffgp[[j.ffgp]])]]-1e-2
    
    dis.pre.0 = dis.ffgp(x.da.ffgp,x.old.ffgp,n,1)
    
    re.ucb.val.ffgp = ucb.ffgp(dis.pre.0)
    ucb.val.ffgp[1] = re.ucb.val.ffgp$ucb
    
    d.k.dis.ffgp = jacobian(function(dis) k.x.ffgp(dis,hyper.ffgp.new[2]), dis.pre.0)
    d.ucb.dis.ffgp = d.k.dis.ffgp%*%re.ucb.val.ffgp$d.ucb.k
    
    d.dis.x.ffgp <- function(x,s0) 2*(x(s0)-x.da.ffgp(s0))
    d.ucb.x.ffgp <- function(x,s0) t(d.ucb.dis.ffgp)%*%as.matrix(d.dis.x.ffgp(x,s0))

    #########FGA Algorithm
    for(l.ffgp in 1:Ta.ffgp){
      grad.old = apply(s.sp,1,function(s0) d.ucb.x.ffgp(x.old.ffgp,s0)^2)
      sr.ffgp = 0.01/sqrt(l.ffgp)
      
      x.new.ffgp <- function(s0) x.old.ffgp(s0)+sr.ffgp*d.ucb.x.ffgp(x.old.ffgp,s0)
      
      ##sparsify
      x.sp.ffgp = apply(s.sp,1,x.new.ffgp)
      x.new.hat.ffgp <- function(s0) x.sp(x.sp.ffgp,s0)
      
      plot(s.sp,x.sp.ffgp,type="l",lwd=2,col=1)
      lines(s.sp,apply(s.sp,1,x.new.hat.ffgp),type="l",lwd=2,col=2)
      lines(s.sp,apply(s.sp,1,x.star),type="l",lwd=2,col=3)
      
      dis.pre.0 = dis.ffgp(x.da.ffgp,x.new.hat.ffgp,n,1)
      
      re.ucb.val.ffgp = ucb.ffgp(dis.pre.0)
      ucb.val.ffgp[l.ffgp+1] = re.ucb.val.ffgp$ucb
      d.k.dis.ffgp = jacobian(function(dis) k.x.ffgp(dis,hyper.ffgp.new[2]), dis.pre.0)
      d.ucb.dis.ffgp = d.k.dis.ffgp%*%re.ucb.val.ffgp$d.ucb.k
      print(paste(l.ffgp,ucb.val.ffgp[l.ffgp+1]))
      
      x.old.ffgp <- x.new.hat.ffgp
      if((ucb.val.ffgp[l.ffgp+1]-ucb.val.ffgp[l.ffgp])^2<eps.ffgp) break
    }
    plot(ucb.val.ffgp,lwd=2,type="l",col=4,ylab="UCB",pch=3)
    
    ## Update design
    x.for.ffgp[[j.ffgp]] = rbind(x.for.ffgp[[j.ffgp]],x.sp.ffgp)
    
    y.new.sp.ffgp <- function(t0) f.obj(x.new.hat.ffgp,t0)+err(t0*j.ffgp)
    g.new.sp.ffgp = g.obj(y.new.sp.ffgp)
    
    g.for.ffgp[[j.ffgp]] = c(g.for.ffgp[[j.ffgp]],g.new.sp.ffgp)
    
    n=n+1
    dis.obs = dis.ffgp(x.da.ffgp,x.da.ffgp,n,n)
    # if(i.ffgp %% 10 == 0){
    #   hyper.for.ffgp.new = bobyqa(hyper.for.ffgp.new, function(the) like.ffgp(dis.obs,
    #             g.da.ffgp,n,th),lower=lower.th,upper=upper.th)$par
    # }else{
    #   hyper.for.ffgp.new = hyper.for.ffgp.new
    # }
    hyper.ffgp.new = bobyqa(hyper.ffgp.new, function(th) like.ffgp(dis.obs,
                     g.for.ffgp[[j.ffgp]],n,th),lower=lower.th,upper=upper.th)$par
    hyper.for.ffgp[[j.ffgp]][[i.ffgp+1]] = hyper.ffgp.new
    
    print(paste("rep:", j.ffgp, ", sequential:", i.ffgp))
  }
  time.ffgp = (proc.time()-ptm)[1]
  time.for.ffgp[[j.ffgp]] = time.ffgp
  
}

re.ffgp = list(x.for.ffgp=x.for.ffgp, y.for.ffgp=y.for.ffgp, g.for.ffgp=g.for.ffgp,
               hyper.for.ffgp=hyper.for.ffgp, beta.ffgp=beta.ffgp, ucb.for.ffgp=ucb.for.ffgp,
               time.for.ffgp=time.for.ffgp) 
save(re.ffgp, file="s1-re-ffgp.RData")  






