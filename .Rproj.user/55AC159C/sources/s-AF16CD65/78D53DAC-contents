#skrypt do basic stochastic volatility
library(pomp)
library(quantmod)
library(beepr)
library(doParallel)
library(doSNOW)
#rm(list = ls())
setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane")

#ladowanie zwrotow z pliku
dane<-read.csv.zoo('C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane/wig_zwroty.csv', header = T,
                   sep=',')
dane=as.xts(dane)


############################################
#poczatek okresu badania 
data.poczatkowa='1996-01'
#koniec okresu badania
data.koncowa='2016-12'
############################################


zwroty=dane[paste(data.poczatkowa,"/",data.koncowa,sep="")]



Breto.SVL_statenames <- c("H","G","Y_state")
Breto.SVL_rp_names <- c("sigma_nu","mu_h","phi","sigma_eta")
Breto.SVL_ivp_names <- c("G_0")
Breto.SVL_paramnames <- c(Breto.SVL_rp_names,Breto.SVL_ivp_names)
Breto.SVL_covarnames <- "covaryt"


rproc1 <- "
double beta,omega,nu;
omega = rnorm(0,sigma_eta  * sqrt(1-tanh(G)*tanh(G)));
nu = rnorm(0, sigma_nu);
G += nu;
beta = Y_state * sigma_eta;
H = mu_h*(1 - phi) + phi*H + beta * tanh( G ) * exp(-H/2) + omega;
"
rproc2.sim <- "
Y_state = rnorm( 0,exp(H/2) );
"

rproc2.filt <- "
Y_state = covaryt;
"
Breto.SVL_rproc.sim <- paste(rproc1,rproc2.sim)
Breto.SVL_rproc.filt <- paste(rproc1,rproc2.filt)


Breto.SVL_initializer <- "
G = G_0;
H =rnorm(mu_h,sigma_eta/sqrt((1-phi*phi)));
Y_state = rnorm( 0,exp(H/2) );
"
Breto.SVL_rmeasure <- "
y=Y_state;
"

Breto.SVL_dmeasure <- "
lik=dnorm(y,0,exp(H/2),give_log);
"

Breto.SVL_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tsigma_nu = log(sigma_nu);
Tphi = logit(phi);
"

Breto.SVL_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tsigma_nu = exp(sigma_nu);
Tphi = expit(phi);
"

Breto.SVL.filt <- pomp(data=data.frame(y=as.vector(zwroty),
                                       time=1:length(zwroty)),
                       statenames=Breto.SVL_statenames,
                       paramnames=Breto.SVL_paramnames,
                       covarnames=Breto.SVL_covarnames,
                       times="time",
                       t0=0,
                       covar=data.frame(covaryt=c(0,as.vector(zwroty)),
                                        time=0:length(zwroty)),
                       tcovar="time",
                       rmeasure=Csnippet(Breto.SVL_rmeasure),
                       dmeasure=Csnippet(Breto.SVL_dmeasure),
                       rprocess=discrete.time.sim(step.fun=Csnippet(Breto.SVL_rproc.filt),delta.t=1),
                       initializer=Csnippet(Breto.SVL_initializer),
                       toEstimationScale=Csnippet(Breto.SVL_toEstimationScale), 
                       fromEstimationScale=Csnippet(Breto.SVL_fromEstimationScale)
)




run_level <- 3 
Breto.SVL_Np <-          c(100,1e3,1e3)
Breto.SVL_Nmif <-        c(10, 100,100)
Breto.SVL_Nreps_eval <-  c(4,  10, 5)
Breto.SVL_Nreps_local <- c(10, 20, 5)
Breto.SVL_Nreps_global <-c(10, 20, 5)

Breto.SVL.list<-list(Breto.SVL_Np ,Breto.SVL_Nmif,Breto.SVL_Nreps_eval,
              Breto.SVL_Nreps_local,Breto.SVL_Nreps_global )

Breto.SVL_rw.sd_rp <- 0.02
Breto.SVL_rw.sd_ivp <- 0.1
Breto.SVL_cooling.fraction.50 <- 0.5


##testowa wersja parametr?w
params_test <- c(
  mu_h = -0.21,       
  phi = .98,     
  sigma_nu = exp(-4.5),
  sigma_eta = .1550,
  G_0 = 0
)

pf1 <- pfilter(Breto.SVL.filt,params=params_test,
               Np=1000,filter.traj=T)
par(mfrow=c(3,1))
plot(zwroty,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
sigma=as.xts(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),order.by=index(zwroty))
plot(sigma,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
rho=as.xts(tanh(pf1@filter.traj[2,1,2:(dim(pf1@filter.traj)[3])]),order.by=index(zwroty))
plot(rho, minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
par(mfrow=c(1,1))


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)




t.if.Breto.SVL <- system.time({
  if.Breto.SVL <- foreach(i=1:Breto.SVL.list[[5]][run_level] ,
                          .packages='pomp', .combine=c,.export = "Breto.SVL.list", 
                          .options.multicore=list(set.seed=TRUE)) %dopar% try(
                            pomp::mif2(Breto.SVL.filt,start=params_test,Np=Breto.SVL.list[[1]][run_level] , Nmif=Breto.SVL.list[[2]][run_level] ,cooling.type="geometric",
                                       cooling.fraction.50=Breto.SVL_cooling.fraction.50,
                                       transform=TRUE,
                                       rw.sd = rw.sd(
                                         mu_h      = Breto.SVL_rw.sd_rp,
                                         phi       = Breto.SVL_rw.sd_rp,
                                         sigma_eta = Breto.SVL_rw.sd_rp,
                                         sigma_nu = Breto.SVL_rw.sd_rp,
                                         G_0       = ivp(Breto.SVL_rw.sd_ivp)
                                       )
                            )
                            
                          )
  
  L.if.Breto.SVL <- foreach(i=1:Breto.SVL.list[[5]][run_level] ,.packages='pomp', .export = "Breto.SVL.list", 
                            .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                              logmeanexp(
                                replicate(Breto.SVL.list[[3]][run_level] ,
                                          logLik(pfilter(Breto.SVL.filt,params=coef(if.Breto.SVL[[i]]),Np=Breto.SVL.list[[1]][run_level]  ))
                                ),
                                se=TRUE)
                            )
  
  H.if.Breto.SVL<- foreach(i=1:Breto.SVL.list[[5]][run_level] ,.packages='pomp', .export = "Breto.SVL.list", 
                           .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                             exp(pfilter(Breto.SVL.filt,params=coef(if.Breto.SVL[[i]]),Np=Breto.SVL.list[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                           )
})


stopCluster(cl)
beep(2)
plot(if.Breto.SVL)


#save(if.Breto.SVL,L.if.Breto.SVL, H.if.Breto.SVL, file="Breto.SVL_if_eval.rda")

r.if.Breto.SVL <- data.frame(logLik=L.if.Breto.SVL[,1],logLik_se=L.if.Breto.SVL[,2],t(sapply(if.Breto.SVL,coef)))
summary(r.if.Breto.SVL$logLik,digits=5)
r.if.Breto.SVL[which.max(r.if.Breto.SVL$logLik),]
pairs(~logLik+mu_h+phi+sigma_eta,data=r.if.Breto.SVL)



params_nowe<- c(
  mu_h        = as.numeric(coef(if.Breto.SVL[which.max(r.if.Breto.SVL$logLik)])[1,'mu_h']),    
  phi         = as.numeric(coef(if.Breto.SVL[which.max(r.if.Breto.SVL$logLik)])[1,'phi']),    
  sigma_nu    = as.numeric(coef(if.Breto.SVL[which.max(r.if.Breto.SVL$logLik)])[1,'sigma_nu']), 
  sigma_eta   = as.numeric(coef(if.Breto.SVL[which.max(r.if.Breto.SVL$logLik)])[1,'sigma_eta']), 
  G_0         = as.numeric(coef(if.Breto.SVL[which.max(r.if.Breto.SVL$logLik)])[1,'G_0'])
)

pf1 <- pfilter(Breto.SVL.filt,params=params_nowe,
               Np=Breto.SVL.list[[1]][run_level],filter.traj=T)

par(mfrow=c(3,1))
plot(zwroty,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
sigma=as.xts(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),order.by=index(zwroty))
plot(sigma,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
rho=as.xts(tanh(pf1@filter.traj[2,1,2:(dim(pf1@filter.traj)[3])]),order.by=index(zwroty))
plot(rho, minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
par(mfrow=c(1,1))

##--------Likelihood maximization using randomized starting values--------

Breto.SVL_box <- rbind(
  sigma_eta=c(0.1,1),
  sigma_nu=c(0.1,.1),
  phi    = c(0.9,0.95),
  mu_h = c(-1,1),
  G_0=c(-0.9,0.9)
)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


t.Breto.SVL.box <- system.time({
  if.Breto.SVL.box <- foreach(i=1:Breto.SVL.list[[5]][run_level],.packages='pomp', .export = "Breto.SVL.list",.combine=c,
                              .options.multicore=list(set.seed=TRUE)) %dopar%  
   
     pomp::mif2(
      Breto.SVL.filt,start=apply(Breto.SVL_box,1,function(x) runif(1,x[1],x[2])),Np=Breto.SVL.list[[1]][run_level] , Nmif=Breto.SVL.list[[2]][run_level] ,cooling.type="geometric",
      cooling.fraction.50=Breto.SVL_cooling.fraction.50,
      transform=TRUE,
      rw.sd = rw.sd(
        mu_h      = Breto.SVL_rw.sd_rp,
        phi       = Breto.SVL_rw.sd_rp,
        sigma_eta = Breto.SVL_rw.sd_rp,
        sigma_nu = Breto.SVL_rw.sd_rp,
        G_0       = ivp(Breto.SVL_rw.sd_ivp)
      )
    )
  
  L.Breto.SVL.box <- foreach(i=1:Breto.SVL.list[[5]][run_level] ,.packages='pomp', .export = "Breto.SVL.list",.combine=rbind,
                             .options.multicore=list(set.seed=TRUE)) %dopar% {
                               set.seed(87932+i)
                               logmeanexp(
                                 replicate(Breto.SVL.list[[3]][run_level] ,
                                           logLik(pfilter(Breto.SVL.filt,params=coef(if.Breto.SVL.box[[i]]),Np=Breto.SVL.list[[1]][run_level] ))
                                 ), 
                                 se=TRUE)
                             }
  
  H.Breto.SVL.box<- foreach(i=1:Breto.SVL.list[[5]][run_level] ,.packages='pomp', .export = "Breto.SVL.list", 
                            .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                              exp(pfilter(Breto.SVL.filt,params=coef(if.Breto.SVL.box[[i]]),Np=Breto.SVL.list[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                            )
})

stopCluster(cl)
beep(2)
plot(if.Breto.SVL.box )

#save(if.Breto.SVL.box,L.Breto.SVL.box,H.Breto.SVL.box, file="Breto.SVL_box_eval.rda")
#load(file="Breto.SVL_box_eval.rda")

r.box <- data.frame(logLik=L.Breto.SVL.box [,1],logLik_se=L.Breto.SVL.box [,2],t(sapply(if.Breto.SVL.box,coef)))
summary(r.box$logLik,digits=5)
round(r.box [which.max(r.box $logLik),],5)


params_nowe2<- c(
  mu_h        = as.numeric(coef(if.Breto.SVL.box [which.max(r.box $logLik)])[1,'mu_h']),    
  phi         = as.numeric(coef(if.Breto.SVL.box [which.max(r.box $logLik)])[1,'phi']),    
  sigma_nu    = as.numeric(coef(if.Breto.SVL.box [which.max(r.box $logLik)])[1,'sigma_nu']), 
  sigma_eta   = as.numeric(coef(if.Breto.SVL.box [which.max(r.box $logLik)])[1,'sigma_eta']), 
  G_0         = as.numeric(coef(if.Breto.SVL.box [which.max(r.box $logLik)])[1,'G_0'])
)



pf2 <- pfilter(Breto.SVL.filt,params=params_nowe2,
               Np=Breto.SVL.list[[1]][run_level],filter.traj=T)

par(mfrow=c(3,1))
plot(zwroty,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
sigma=as.xts(exp(pf2@filter.traj[1,1,2:(dim(pf2@filter.traj)[3])]/2),order.by=index(zwroty))
plot(sigma,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
rho=as.xts(tanh(pf2@filter.traj[2,1,2:(dim(pf2@filter.traj)[3])]),order.by=index(zwroty))
plot(rho, minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
par(mfrow=c(1,1))

round(params_nowe2,4)
############################################################################
############################################################################
############################################################################
#profile funkcji wiarygodnosci


params_nowe2<- c(
  mu_h        =  0.3465,    
  phi         =  0.9778,    
  sigma_nu    =  0.1776, 
  sigma_eta   =  0.0003, 
  G_0         =-.22  
)

#parametr mu_h
xx1<-seq(from=params_nowe2['mu_h']-1,to=params_nowe2['mu_h']+1,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.Breto.SVL.log<- foreach(i=1:length(xx1) ,.packages='pomp', .export = "Breto.SVL.list",.combine=rbind,
                          .options.multicore=list(set.seed=TRUE)) %dopar% {
                            set.seed(87932)
                                        logLik(pfilter(Breto.SVL.filt,params=c(mu_h=xx1[i], 
                                                                               phi=as.numeric(params_nowe2['phi']),
                                                                               sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                                               sigma_nu=as.numeric(params_nowe2['sigma_nu']),
                                                                               G_0=as.numeric(params_nowe2['G_0'])),
                                                       Np=Breto.SVL.list[[1]][run_level] ))
                          }

stopCluster(cl)
beep(2)

plot(xx1, L.Breto.SVL.log, type='l',xlab=expression(mu[h]),ylab="logLik")
points(xx1, L.Breto.SVL.log)
#points(r.box[,'mu_h'], r.box[,'logLik'] ,col='red')
p=loess(L.Breto.SVL.log~xx1,span=0.5)
lines(xx1,p$fitted,col='blue',lwd=2)


#parametr phi
xx2<-seq(from=0.95,to=0.999,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.Breto.SVL.log2<- foreach(i=1:length(xx2) ,.packages='pomp', .export = "Breto.SVL.list",.combine=rbind,
                          .options.multicore=list(set.seed=TRUE)) %dopar% {
                            set.seed(87932)

                                        logLik(pfilter(Breto.SVL.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                                               phi=xx2[i],
                                                                               sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                                               sigma_nu=as.numeric(params_nowe2['sigma_nu']),
                                                                               G_0=as.numeric(params_nowe2['G_0'])),
                                                       Np=Breto.SVL.list[[1]][run_level] ))
}

stopCluster(cl)
beep(1)

plot(xx2, L.Breto.SVL.log2, type='l',xlab=expression(mu[h]),ylab="logLik")
points(xx2, L.Breto.SVL.log2)
points(r.box[,'phi'], r.box[,'logLik'] ,col='red')
p2=loess(L.Breto.SVL.log2~xx2, span=0.5)
lines(xx2,p2$fitted,col='blue',lwd=2)


#parametr sigma_eta
xx3<-seq(from=0.1,to=0.2,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.Breto.SVL.log3<- foreach(i=1:length(xx3) ,.packages='pomp', .export = "Breto.SVL.list",.combine=rbind,
                           .options.multicore=list(set.seed=TRUE)) %dopar% {
                             set.seed(87932)
                                         logLik(pfilter(Breto.SVL.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                                                phi=as.numeric(params_nowe2['phi']),
                                                                                sigma_eta=xx3[i],
                                                                                sigma_nu=as.numeric(params_nowe2['sigma_nu']),
                                                                                G_0=as.numeric(params_nowe2['G_0'])),
                                                        Np=Breto.SVL.list[[1]][run_level] ))
                           }

stopCluster(cl)
beep(1)

plot(xx3, L.Breto.SVL.log3, type='l',xlab=expression(sigma[eta]),ylab="logLik")
points(xx3, L.Breto.SVL.log3)
#points(r.box[,'sigma_eta'], r.box[,'logLik'] ,col='red')
p3=loess(L.Breto.SVL.log3~xx3,span=0.75)
lines(xx3,p3$fitted,col='blue',lwd=2)


#parametr sigma_nu
xx4<-seq(from=0.001,to=.1,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.Breto.SVL.log4<- foreach(i=1:length(xx4) ,.packages='pomp', .export = "Breto.SVL.list",.combine=rbind,
                           .options.multicore=list(set.seed=TRUE)) %dopar% {
                             set.seed(87932)
                                         logLik(pfilter(Breto.SVL.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                                                phi=as.numeric(params_nowe2['phi']),
                                                                                sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                                                sigma_nu=xx4[i],
                                                                                G_0=as.numeric(params_nowe2['G_0'])),
                                                        Np=Breto.SVL.list[[1]][run_level] ))
                           }

stopCluster(cl)
beep(1)

plot(xx4+0.01, L.Breto.SVL.log4, type='l',xlab=expression(sigma[nu]),ylab="logLik")
points(xx4+0.01, L.Breto.SVL.log4)
#points(r.box[,'sigma_nu'], r.box[,'logLik'] ,col='red')
p4=loess(L.Breto.SVL.log4~(xx4),span=0.75)
lines(xx4+0.01,p4$fitted,col='blue',lwd=2)

#parametr G_0
xx5<-seq(from=params_nowe2['G_0']-.1,to=params_nowe2['G_0']+.1,length.out = 10)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.Breto.SVL.log5<- foreach(i=1:length(xx4) ,.packages='pomp', .export = "Breto.SVL.list",.combine=rbind,
                           .options.multicore=list(set.seed=TRUE)) %dopar% {
                             set.seed(87932)
                             logmeanexp(
                               replicate(Breto.SVL.list[[3]][run_level] ,
                                         logLik(pfilter(Breto.SVL.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                                                phi=as.numeric(params_nowe2['phi']),
                                                                                sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                                                sigma_nu=as.numeric(params_nowe2['sigma_eta']),
                                                                                G_0=xx5[i]),
                                                        Np=Breto.SVL.list[[1]][run_level] ))
                               ), 
                               se=FALSE)
                           }

stopCluster(cl)
beep(2)

plot(xx5, L.Breto.SVL.log5, type='l',xlab=expression(sigma[nu]),ylab="logLik",mgp=c(3, 1, 0))
points(xx5, L.Breto.SVL.log5)
points(r.box[,'G_0'], r.box[,'logLik'] ,col='red')


par(mfrow = c(2,2),
    oma = c(1,4,0,0) + 0.1,
    mar = c(2.4,1,1,1) + 0.1)
plot(xx1, L.Breto.SVL.log, type='l',xlab=expression(mu[h]),ylab="logLik",mgp=c(1.7, 1, 0))
points(xx1, L.Breto.SVL.log)
#points(r.box[,'mu_h'], r.box[,'logLik'] ,col='red')
p=loess(L.Breto.SVL.log~xx1,span=0.5)
lines(xx1,p$fitted,col='blue',lwd=2)
abline(v=params_nowe2['mu_h'],lty=2)

plot(xx2, L.Breto.SVL.log2, type='l',xlab=expression(phi),ylab="logLik",mgp=c(1.7, 1, 0))
points(xx2, L.Breto.SVL.log2)
#points(r.box[,'phi'], r.box[,'logLik'] ,col='red')
p2=loess(L.Breto.SVL.log2~xx2, span=0.5)
lines(xx2,p2$fitted,col='blue',lwd=2)
abline(v=params_nowe2['phi'],lty=2)

plot(xx3, L.Breto.SVL.log3, type='l',xlab=expression(sigma[eta]),ylab="logLik",mgp=c(1.7, 1, 0))
points(xx3, L.Breto.SVL.log3)
#points(r.box[,'sigma_eta'], r.box[,'logLik'] ,col='red')
p3=loess(L.Breto.SVL.log3~xx3,span=0.75)
lines(xx3,p3$fitted,col='blue',lwd=2)
abline(v=params_nowe2['sigma_eta'],lty=2) 

plot(xx4+0.01, L.Breto.SVL.log4, type='l',xlab=expression(sigma[nu]),ylab="logLik",mgp=c(1.7, 1, 0))
points(xx4+0.01, L.Breto.SVL.log4)
#points(r.box[,'sigma_nu'], r.box[,'logLik'] ,col='red')
p4=loess(L.Breto.SVL.log4~(xx4),span=0.75)
lines(xx4+0.01,p4$fitted,col='blue',lwd=2)
abline(v=params_nowe2['sigma_nu']+0.02,lty=2) 

par(mfrow=c(1,1))


############################################################################
############################################################################
############################################################################
#PMCM
############################################################################
############################################################################
############################################################################

params_test<- c(
mu_h        =  0.3465,    
phi         =  0.9778,    
sigma_eta   =  0.1, 
sigma_nu    =  0.1776, 
G_0         =-.22  
)

hyperparams <- list(min = c(-1,0.9,0.1,0,-1), max = c(1,1,.5,.1,1) )
bsv.dprior <- function (params, ..., log) {
  f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max,
                 log = TRUE))
  if (log) f else exp(f)
}


pmcmc1 <-   pmcmc(pomp(Breto.SVL.filt, dprior = bsv.dprior), start = params_test,
                  Nmcmc = 1000, Np = 100, max.fail = Inf,
                  proposal = mvn.diag.rw(c(mu_h = 0.01, phi = 0.01, sigma_eta = 0.01, sigma_nu=0.001 ,G_0=0.01)))

continue( pmcmc1 ,Nmcmc=5000,proposal=mvn.rw(covmat( pmcmc1 ))) -> pmcmc1 
plot(  pmcmc1 )
pf1 <- pfilter(Breto.SVL.filt,params=params_test,
               Np=1000,filter.traj=T)
plot(pf1)
