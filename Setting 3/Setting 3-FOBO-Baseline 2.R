getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("lib.R")
library(fda) # fpca
library(mlegp) # gaussian process
library(randtoolbox) # sobol point
library(statmod) # quadrature
library(support) # support point
library(CompQuadForm) # generalized chi-square
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
  return(-20 * exp(-5 * dist) - 10 * sin(3 * pi * t0) * asym)
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
y.star.dis = apply(t.dis,1,function(t0) f.obj(x.star,t0)) 
g.obj.dis <- function(y) mean(y)

## FPCA
basis.breaks <- seq(0, 1, length.out = 21)
basisobj <- create.bspline.basis(rangeval = c(0, 1), breaks = seq(0, 1, length.out = 21))
fdParobj <- fdPar(basisobj, lambda = 2)


############## FFBO ##################
J.fogp = 10; n0=10
y.for.fogp = x.for.fogp = g.for.fogp = list()
f.for.fogp = lapply(1:J.fogp, function(x) list())
hyper.for.fogp = lapply(1:J.fogp, function(x) list())
beta.fogp = ucb.for.fogp = lapply(1:J.fogp, function(x) list())
time.for.fogp = list()

x.ini.dis = ini.data$x.ini.dis
y.ini.for.dis = ini.data$y.ini.for.dis
g.ini.for.dis = ini.data$g.ini.for.dis



for(j.fogp in 1:J.fogp){
  ptm = proc.time()
  n=10; m=20
  
  nodes <- seq(1, 2*n0 - 1, by = 2) / (2 * n0)  # midpoints of intervals
  weights <- rep(1/n0, n0)
  quad.node <- list(nodes = nodes, weights = weights)
  gchisq.sp <- NULL
  
  x.for.fogp[[j.fogp]] = x.ini.dis
  g.for.fogp[[j.fogp]] = g.ini.for.dis[[j.fogp]]
  y.for.fogp[[j.fogp]] = y.ini.for.dis[[j.fogp]]
  
  opar<-par(no.readonly = TRUE)
  par(mfrow = c(2, 2))
  par(mar = c(4.5,4.5,2,2),xpd=TRUE)
  for(i in 1:n){
    plot(s.dis, x.for.fogp[[j.fogp]][i,],type="l", lwd=3, ylab="x")
    plot(t.dis, y.for.fogp[[j.fogp]][i,],type="l", lwd=3, ylab="f")
  }
  
  mse.best.ind.fogp <- which.min(apply(y.for.fogp[[j.fogp]],1,function(x) mean(x^2)))
  mse.best.fogp = mean((y.for.fogp[[j.fogp]][mse.best.ind.fogp,])^2)
  x.new.fogp = x.for.fogp[[j.fogp]][mse.best.ind.fogp,]
  
  for(i.fogp in 1:m){
    basis <- smooth.basis(t.dis, t(y.for.fogp[[j.fogp]]), fdParobj)
    fpca <- pca.fd(basis$fd, nharm = basisobj$nbasis)
    nharm <- which(cumsum(fpca$varprop) >= 0.99)[1] 
    fpca <- pca.fd(basis$fd, nharm = nharm)
    fpca.eig.value <- fpca$values[1:nharm]
    
    fpca.std.score <- t(t(fpca$scores) / sqrt(fpca.eig.value))
    
    fpca.mean <- eval.fd(t.dis, fpca$meanfd)
    residual <- y.star.dis - fpca.mean
    target.basis <- smooth.basis(t.dis, residual, fdParobj)
    target.score <- inprod(target.basis$fd, fpca$harmonics)
    target.std.score <- target.score / sqrt(fpca.eig.value)
    
    
    models <- list()
    for (i in 1:nharm){
      invisible(capture.output(
        models[[i]] <- mlegp(x.for.fogp[[j.fogp]], fpca.std.score[,i])
      ))
    } 
    
    if (is.null(gchisq.sp)){
      gchisq.sp <- sp(100,nharm,dist.str=rep("normal",nharm))$sp
    } else{
      if (ncol(gchisq.sp)!=nharm){
        gchisq.sp <- sp(100,nharm,dist.str=rep("normal",nharm))$sp
      }
    }
    
    mse.best.ind.fogp <- which.min(apply(y.for.fogp[[j.fogp]],1,function(x) mean(x^2)))
    mse.best.fogp = mean((y.for.fogp[[j.fogp]][mse.best.ind.fogp,])^2)
    
    EI <- function(x){
      score.fit <- rep(0, nharm)
      score.se <- rep(0, nharm)
      for (i in 1:nharm){
        score <- predict(models[[i]], x, se.fit = T)
        score.fit[i] <- score$fit
        score.se[i] <- score$se.fit
      }
      if (any(score.se < 1e-9)) return (0)
      # compute the expected improvement
      delta <- (score.fit - target.std.score)^2 / (score.se)^2
      lambda <- fpca.eig.value * (score.se)^2 / mse.best.fogp
      
      # integral approximation by quadrature
      quad.node.y <- rep(NA, length(quad.node$nodes))
      for (i in 1:length(quad.node$nodes)){
        quad.node.y[i] <- pgchisq(q=quad.node$nodes[i], lambda=lambda, delta=delta, D=gchisq.sp)
        # quad.node.y[i] <- 1 - imhof(q=quad.node$nodes[i], lambda=lambda, delta=delta)$Qq
      }
      val <- sum(quad.node.y * quad.node$weights) * mse.best.fogp
      return (-val) # negative for minimization
    }
    
    x.new.fogp = optim(x.new.fogp, EI, method = "Nelder-Mead")$par
    fit.fogp <- smooth.spline(s.dis,x.new.fogp)
    x.sp.fogp <- function(s0) predict(fit.fogp,s0)$y
    y.new.fogp = apply(t.dis,1,function(t0) f.obj(x.sp.fogp,t0))
    g.new.fogp = g.obj.dis(y.new.fogp)
    
    plot(s.dis,x.new.fogp,type="l",lwd=2,col=1)
    lines(s.sp,apply(s.sp,1,x.star),type="l",lwd=2,col=3)
    
    x.for.fogp[[j.fogp]] = rbind(x.for.fogp[[j.fogp]],x.new.fogp)
    y.for.fogp[[j.fogp]] = rbind(y.for.fogp[[j.fogp]],y.new.fogp)
    g.for.fogp[[j.fogp]] = c(g.for.fogp[[j.fogp]],g.new.fogp)
    
    print(paste("rep:", j.fogp, ", sequential:", i.fogp))
  }
}


re.fogp = list(x.for.fogp=x.for.fogp, y.for.fogp=y.for.fogp, g.for.fogp=g.for.fogp,
               hyper.for.fogp=hyper.for.fogp, beta.fogp=beta.fogp, ucb.for.fogp=ucb.for.fogp,
               time.for.fogp=time.for.fogp) 
save(re.fogp, file="s3-re-fogp.RData")  













