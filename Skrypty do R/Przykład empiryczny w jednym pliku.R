library(pomp)
library(beepr)
library(doParallel)
library(doSNOW)
library(quantmod)
library(KFK)
#rm(list = ls())
setwd("C:/Users/user/Documents/github/Time-varying-beta-with-stochastic-volatility/Dane")


#sciaganie danych
baza.danych.zwroty=read.csv.zoo("baza_danych_zwroty.csv",sep=',')
baza.danych.zwroty=as.xts(baza.danych.zwroty)
wig.zwroty=read.csv.zoo("wig_zwroty_sub.csv",sep=',')
wig.zwroty=as.xts(wig.zwroty)


#wybór tickera
tiker='PKO'

y=baza.danych.zwroty[,tiker]
y=y[-1]
x=wig.zwroty

n=length(x)

par(mfrow=c(2,1))
plot(x,main='WIG',major.ticks = "years",grid.ticks.on = "years")
plot(y,main=tiker,major.ticks = "years",grid.ticks.on = "years")
par(mfrow=c(1,1))


####nazwy
beta_statenames <- c("H","R_i","Beta","Alpha")
beta_rp_names <- c("mu_h","phi","sigma_eta","sigma_ksi","sigma_gamma")
beta_ivp_names <- c("Beta_0","H_0","Alpha_0")
beta_paramnames <- c(beta_rp_names,beta_ivp_names)
beta_covarnames <- c("R_i_dane","R_M")

#pytanie czy Beta_0 powinno byc estymowane?

rproc1 <- "
double omega;
double ksi;
double gamma;
omega = rnorm(0,sigma_eta );
ksi= rnorm(0,sigma_ksi );
gamma= rnorm(0,sigma_gamma );
Beta=Beta+ksi;
Alpha=Alpha+gamma;
H = mu_h*(1 - phi) + phi*H + omega;
"

####rownanie procesu pomiaru
rproc2.sim <- "
R_i = rnorm( Alpha+Beta*R_M,exp(H/2) );
"
###do wypelniania danych
rproc2.filt <- "
R_i = R_i_dane;
"

###symulacja modelu SVL
beta_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
beta_rproc.filt <- paste(rproc1,rproc2.filt)


#Y_state = rnorm( 0,exp(H/2) );
######inicalizacja
beta_initializer <- "
Beta=Beta_0;
H = H_0 ;
Alpha=Alpha_0;
R_i = rnorm(Alpha+Beta*R_M,exp(H/2) );
"
###????
beta_rmeasure <- "
y=R_i;
"

####rozk?ad warunkowy zmiennej Y
beta_dmeasure <- "
lik=dnorm(y,Alpha+Beta*R_M,exp(H/2),give_log);
"


####przeskalowanie parametr?w 
beta_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tsigma_ksi = log(sigma_ksi);
Tsigma_gamma = log(sigma_gamma);
Tphi = logit(phi);
"

beta_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tsigma_ksi = exp(sigma_ksi);
Tsigma_gamma = exp(sigma_gamma);
Tphi = expit(phi);
"


####wypelnianie modelu danymi
beta.filt <- pomp(data=data.frame(y=as.vector(y),
                                  time=1:length(y)),
                  statenames=beta_statenames,
                  paramnames=beta_paramnames,
                  covarnames=beta_covarnames,
                  times="time",
                  t0=0,
                  covar=data.frame(R_i_dane=c(0,as.vector(y)) ,R_M=c(0,as.vector(x)),
                                   time=0:length(y)),
                  tcovar="time",
                  rmeasure=Csnippet(beta_rmeasure),
                  dmeasure=Csnippet(beta_dmeasure),
                  rprocess=discrete.time.sim(step.fun=Csnippet(beta_rproc.filt),delta.t=1),
                  initializer=Csnippet(beta_initializer),
                  toEstimationScale=Csnippet(beta_toEstimationScale), 
                  fromEstimationScale=Csnippet(beta_fromEstimationScale)
)



###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 2

#liczba czasteczek
beta_Np <-          c(1000,1e3,1e3)
beta_Nmif <-        c(10, 50, 200)
beta_Nreps_eval <-  c(4,  10,  20)
beta_Nreps_local <- c(4, 10, 10)
beta_Nreps_global <-c(4, 10, 20)

betalist<-list(beta_Np ,beta_Nmif,beta_Nreps_eval,
               beta_Nreps_local,beta_Nreps_global )

#parametry do metody mif2
beta_rw.sd_rp <- 0.02
beta_rw.sd_ivp <- 0.1
beta_cooling.fraction.50 <- 0.5





beta_box <- rbind(
  sigma_eta=c(0.1,1),
  phi    = c(0.9,0.99),
  mu_h = c(-1,1),
  sigma_ksi=c(0.1,1),
  sigma_gamma=c(0.1,1),
  Beta_0=c(0,1),
  Alpha_0=c(0,1),  
  H_0=c(-1,1)
)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


start_time <- Sys.time()
t.beta.box <- system.time({
  if.beta.box <- foreach(i=1:betalist[[5]][run_level],.packages='pomp', .export = "betalist",.combine=c,
                         .options.multicore=list(set.seed=TRUE)) %dopar%  
    pomp::mif2(beta.filt,Np=betalist[[1]][run_level] , Nmif=betalist[[2]][run_level] ,
               start=apply(beta_box,1,function(x) runif(1,x[1],x[2])),cooling.type="geometric", cooling.fraction.50=beta_cooling.fraction.50,
               transform=TRUE,
               rw.sd = rw.sd(
                 mu_h      = beta_rw.sd_rp,
                 phi       = beta_rw.sd_rp,
                 sigma_eta = beta_rw.sd_rp,
                 sigma_ksi = beta_rw.sd_rp,
                 sigma_gamma = beta_rw.sd_rp,
                 Beta_0       = ivp(beta_rw.sd_ivp),
                 Alpha_0      = ivp(beta_rw.sd_ivp),
                 H_0       = ivp(beta_rw.sd_ivp)
               )
    )
  L.beta.box <- foreach(i=1:betalist[[5]][run_level] ,.packages='pomp', .export = "betalist",.combine=rbind,
                        .options.multicore=list(set.seed=TRUE)) %dopar% {
                          set.seed(87932)
                          logmeanexp(
                            replicate(betalist[[3]][run_level] ,
                                      logLik(pfilter(beta.filt,params=coef(if.beta.box[[i]]),Np=betalist[[1]][run_level] ))
                            ), 
                            se=TRUE)
                        }
  
  H.beta.box<- foreach(i=1:betalist[[5]][run_level] ,.packages='pomp', .export = "betalist", 
                       .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                         exp(pfilter(beta.filt,params=coef(if.beta.box[[i]]),Np=betalist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                       )
})



stopCluster(cl)
end_time <- Sys.time()
end_time - start_time



r.box <- data.frame(logLik=L.beta.box [,1],logLik_se=L.beta.box [,2],t(sapply(if.beta.box,coef)))
round(apply(as.matrix(r.box),MARGIN=2, FUN=mean ),4)
round(r.box [which.max(r.box $logLik),],4)


params_nowe2<- c(
  mu_h        = as.numeric(r.box [which.max(r.box $logLik),'mu_h']),    
  phi         = as.numeric(r.box [which.max(r.box $logLik),'phi']),    
  sigma_eta   = as.numeric(r.box [which.max(r.box $logLik),'sigma_eta']), 
  sigma_ksi   = as.numeric(r.box [which.max(r.box $logLik),'sigma_ksi']), 
  sigma_gamma   = as.numeric(r.box [which.max(r.box $logLik),'sigma_gamma']),  
  Alpha_0        = as.numeric(r.box [which.max(r.box $logLik),'Alpha_0']),
  Beta_0       = as.numeric(r.box [which.max(r.box $logLik),'Beta_0']),
  H_0       = as.numeric(r.box [which.max(r.box $logLik),'H_0'])
)

pf1 <- pfilter(beta.filt,params=params_nowe2,
               Np=betalist[[1]][run_level],filter.traj=T)

wersja1=pf1@filter.traj[3,1,2:(dim(pf1@filter.traj)[3])]
wersja1=as.xts(wersja1,order.by = index(x))

beep(2)
################################################################################
################################################################################
################################################################################



####nazwy
beta_statenames <- c("H","R_i","Beta")
beta_rp_names <- c("mu_h","phi","sigma_eta","sigma_ksi","alpha")
beta_ivp_names <- c("Beta_0","H_0")
beta_paramnames <- c(beta_rp_names,beta_ivp_names)
beta_covarnames <- c("R_i_dane","R_M")

#pytanie czy Beta_0 powinno byc estymowane?

rproc1 <- "
double omega;
double ksi;
omega = rnorm(0,sigma_eta );
ksi= rnorm(0,sigma_ksi );
Beta=Beta+ksi;
H = mu_h*(1 - phi) + phi*H + omega;
"

####rownanie procesu pomiaru
rproc2.sim <- "
R_i = rnorm( alpha+Beta*R_M,exp(H/2) );
"
###do wypelniania danych
rproc2.filt <- "
R_i = R_i_dane;
"

###symulacja modelu SVL
beta_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
beta_rproc.filt <- paste(rproc1,rproc2.filt)


#Y_state = rnorm( 0,exp(H/2) );
######inicalizacja
beta_initializer <- "
Beta=Beta_0;
H = H_0 ;
R_i = rnorm(alpha+Beta*R_M,exp(H/2) );
"
###????
beta_rmeasure <- "
y=R_i;
"

####rozk?ad warunkowy zmiennej Y
beta_dmeasure <- "
lik=dnorm(y,alpha+Beta*R_M,exp(H/2),give_log);
"


####przeskalowanie parametr?w 
beta_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tsigma_ksi = log(sigma_ksi);
Tphi = logit(phi);
"

beta_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tsigma_ksi = exp(sigma_ksi);
Tphi = expit(phi);
"


####wypelnianie modelu danymi
beta.filt <- pomp(data=data.frame(y=as.vector(y),
                                  time=1:length(y)),
                  statenames=beta_statenames,
                  paramnames=beta_paramnames,
                  covarnames=beta_covarnames,
                  times="time",
                  t0=0,
                  covar=data.frame(R_i_dane=c(0,as.vector(y)) ,R_M=c(0,as.vector(x)),
                                   time=0:length(y)),
                  tcovar="time",
                  rmeasure=Csnippet(beta_rmeasure),
                  dmeasure=Csnippet(beta_dmeasure),
                  rprocess=discrete.time.sim(step.fun=Csnippet(beta_rproc.filt),delta.t=1),
                  initializer=Csnippet(beta_initializer),
                  toEstimationScale=Csnippet(beta_toEstimationScale), 
                  fromEstimationScale=Csnippet(beta_fromEstimationScale)
)


###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 2

#liczba czasteczek
beta_Np <-          c(1000,1e3,1e3)
beta_Nmif <-        c(10, 50, 200)
beta_Nreps_eval <-  c(4,  10,  20)
beta_Nreps_local <- c(4, 10, 10)
beta_Nreps_global <-c(4, 10, 20)

betalist<-list(beta_Np ,beta_Nmif,beta_Nreps_eval,
               beta_Nreps_local,beta_Nreps_global )

#parametry do metody mif2
beta_rw.sd_rp <- 0.02
beta_rw.sd_ivp <- 0.1
beta_cooling.fraction.50 <- 0.5





beta_box <- rbind(
  sigma_eta=c(0.1,1),
  phi    = c(0.9,0.99),
  mu_h = c(-1,1),
  sigma_ksi=c(0.1,1),
  alpha=c(-0.1,0.1),
  Beta_0=c(0,1),
  H_0=c(-1,1)
)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


start_time <- Sys.time()
t.beta.box <- system.time({
  if.beta.box <- foreach(i=1:betalist[[5]][run_level],.packages='pomp', .export = "betalist",.combine=c,
                         .options.multicore=list(set.seed=TRUE)) %dopar%  
    pomp::mif2(beta.filt,Np=betalist[[1]][run_level] , Nmif=betalist[[2]][run_level] ,
               start=apply(beta_box,1,function(x) runif(1,x[1],x[2])),cooling.type="geometric", cooling.fraction.50=beta_cooling.fraction.50,
               transform=TRUE,
               rw.sd = rw.sd(
                 mu_h      = beta_rw.sd_rp,
                 phi       = beta_rw.sd_rp,
                 sigma_eta = beta_rw.sd_rp,
                 sigma_ksi = beta_rw.sd_rp,
                 alpha=beta_rw.sd_rp,
                 Beta_0       = ivp(beta_rw.sd_ivp),
                 H_0       = ivp(beta_rw.sd_ivp)
               )
    )
  L.beta.box <- foreach(i=1:betalist[[5]][run_level] ,.packages='pomp', .export = "betalist",.combine=rbind,
                        .options.multicore=list(set.seed=TRUE)) %dopar% {
                          set.seed(87932)
                          logmeanexp(
                            replicate(betalist[[3]][run_level] ,
                                      logLik(pfilter(beta.filt,params=coef(if.beta.box[[i]]),Np=betalist[[1]][run_level] ))
                            ), 
                            se=TRUE)
                        }
  
  H.beta.box<- foreach(i=1:betalist[[5]][run_level] ,.packages='pomp', .export = "betalist", 
                       .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                         exp(pfilter(beta.filt,params=coef(if.beta.box[[i]]),Np=betalist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                       )
})



stopCluster(cl)
end_time <- Sys.time()
end_time - start_time


r.box <- data.frame(logLik=L.beta.box [,1],logLik_se=L.beta.box [,2],t(sapply(if.beta.box,coef)))
round(apply(as.matrix(r.box),MARGIN=2, FUN=mean ),4)
round(r.box [which.max(r.box $logLik),],4)


params_nowe2<- c(
  mu_h        = as.numeric(r.box [which.max(r.box $logLik),'mu_h']),    
  phi         = as.numeric(r.box [which.max(r.box $logLik),'phi']),    
  sigma_eta   = as.numeric(r.box [which.max(r.box $logLik),'sigma_eta']), 
  sigma_ksi   = as.numeric(r.box [which.max(r.box $logLik),'sigma_ksi']), 
  alpha       = as.numeric(r.box [which.max(r.box $logLik),'alpha']),
  Beta_0       = as.numeric(r.box [which.max(r.box $logLik),'Beta_0']),
  H_0       = as.numeric(r.box [which.max(r.box $logLik),'H_0'])
)

pf1 <- pfilter(beta.filt,params=params_nowe2,
               Np=betalist[[1]][run_level],filter.traj=T)

wersja2=pf1@filter.traj[3,1,2:(dim(pf1@filter.traj)[3])]
wersja2=as.xts(wersja2,order.by = index(x))


beep(2)

################################################################################
################################################################################
################################################################################




####nazwy
beta_statenames <- c("H","R_i","Beta")
beta_rp_names <- c("mu_h","phi","sigma_eta","sigma_ksi","alpha","mu_beta","nu")
beta_ivp_names <- c("Beta_0","H_0")
beta_paramnames <- c(beta_rp_names,beta_ivp_names)
beta_covarnames <- c("R_i_dane","R_M")

#pytanie czy Beta_0 powinno byc estymowane?

rproc1 <- "
double omega;
double ksi;
omega = rnorm(0,sigma_eta );
ksi= rnorm(0,sigma_ksi );
Beta=mu_beta*(1 - nu)+nu*Beta+ksi;
H = mu_h*(1 - phi) + phi*H + omega;
"

####rownanie procesu pomiaru
rproc2.sim <- "
R_i = rnorm( alpha+Beta*R_M,exp(H/2) );
"
###do wypelniania danych
rproc2.filt <- "
R_i = R_i_dane;
"

###symulacja modelu SVL
beta_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
beta_rproc.filt <- paste(rproc1,rproc2.filt)


#Y_state = rnorm( 0,exp(H/2) );
######inicalizacja
beta_initializer <- "
Beta=Beta_0;
H = H_0 ;
R_i = rnorm(alpha+Beta*R_M,exp(H/2) );
"
###????
beta_rmeasure <- "
y=R_i;
"

####rozk?ad warunkowy zmiennej Y
beta_dmeasure <- "
lik=dnorm(y,alpha+Beta*R_M,exp(H/2),give_log);
"


####przeskalowanie parametr?w 
beta_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tsigma_ksi = log(sigma_ksi);
Tphi = logit(phi);
Tnu = logit(nu);
"

beta_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tsigma_ksi = exp(sigma_ksi);
Tphi = expit(phi);
Tnu = expit(nu);
"


####wypelnianie modelu danymi
beta.filt <- pomp(data=data.frame(y=as.vector(y),
                                  time=1:length(y)),
                  statenames=beta_statenames,
                  paramnames=beta_paramnames,
                  covarnames=beta_covarnames,
                  times="time",
                  t0=0,
                  covar=data.frame(R_i_dane=c(0,as.vector(y)) ,R_M=c(0,as.vector(x)),
                                   time=0:length(y)),
                  tcovar="time",
                  rmeasure=Csnippet(beta_rmeasure),
                  dmeasure=Csnippet(beta_dmeasure),
                  rprocess=discrete.time.sim(step.fun=Csnippet(beta_rproc.filt),delta.t=1),
                  initializer=Csnippet(beta_initializer),
                  toEstimationScale=Csnippet(beta_toEstimationScale), 
                  fromEstimationScale=Csnippet(beta_fromEstimationScale)
)





###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 2

#liczba czasteczek
beta_Np <-          c(1000,1e3,1e3)
beta_Nmif <-        c(10, 50, 200)
beta_Nreps_eval <-  c(4,  10,  20)
beta_Nreps_local <- c(4, 10, 10)
beta_Nreps_global <-c(4, 10, 20)

betalist<-list(beta_Np ,beta_Nmif,beta_Nreps_eval,
               beta_Nreps_local,beta_Nreps_global )

#parametry do metody mif2
beta_rw.sd_rp <- 0.02
beta_rw.sd_ivp <- 0.1
beta_cooling.fraction.50 <- 0.5





beta_box <- rbind(
  sigma_eta=c(0.1,1),
  phi    = c(0.9,0.99),
  mu_h = c(-1,1),
  mu_beta = c(-1,1),
  sigma_ksi=c(0.1,1),
  alpha=c(-0.1,0.1),
  nu= c(0.9,0.99),
  Beta_0=c(0,1),
  H_0=c(-1,1)
)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


start_time <- Sys.time()
t.beta.box <- system.time({
  if.beta.box <- foreach(i=1:betalist[[5]][run_level],.packages='pomp', .export = "betalist",.combine=c,
                         .options.multicore=list(set.seed=TRUE)) %dopar%  
    pomp::mif2(beta.filt,Np=betalist[[1]][run_level] , Nmif=betalist[[2]][run_level] ,
               start=apply(beta_box,1,function(x) runif(1,x[1],x[2])),cooling.type="geometric", cooling.fraction.50=beta_cooling.fraction.50,
               transform=TRUE,
               rw.sd = rw.sd(
                 mu_h      = beta_rw.sd_rp,
                 phi       = beta_rw.sd_rp,
                 sigma_eta = beta_rw.sd_rp,
                 sigma_ksi = beta_rw.sd_rp,
                 alpha=beta_rw.sd_rp,
                 mu_beta      = beta_rw.sd_rp,
                 nu      = beta_rw.sd_rp,
                 Beta_0       = ivp(beta_rw.sd_ivp),
                 H_0       = ivp(beta_rw.sd_ivp)
               )
    )
  L.beta.box <- foreach(i=1:betalist[[5]][run_level] ,.packages='pomp', .export = "betalist",.combine=rbind,
                        .options.multicore=list(set.seed=TRUE)) %dopar% {
                          set.seed(87932)
                          logmeanexp(
                            replicate(betalist[[3]][run_level] ,
                                      logLik(pfilter(beta.filt,params=coef(if.beta.box[[i]]),Np=betalist[[1]][run_level] ))
                            ), 
                            se=TRUE)
                        }
  
  H.beta.box<- foreach(i=1:betalist[[5]][run_level] ,.packages='pomp', .export = "betalist", 
                       .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                         exp(pfilter(beta.filt,params=coef(if.beta.box[[i]]),Np=betalist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                       )
})



stopCluster(cl)
end_time <- Sys.time()
end_time - start_time



r.box <- data.frame(logLik=L.beta.box [,1],logLik_se=L.beta.box [,2],t(sapply(if.beta.box,coef)))
round(apply(as.matrix(r.box),MARGIN=2, FUN=mean ),4)
round(r.box [which.max(r.box $logLik),],4)


params_nowe2<- c(
  mu_h        = as.numeric(r.box [which.max(r.box $logLik),'mu_h']),    
  phi         = as.numeric(r.box [which.max(r.box $logLik),'phi']),   
  nu        = as.numeric(r.box [which.max(r.box $logLik),'nu']), 
  sigma_eta   = as.numeric(r.box [which.max(r.box $logLik),'sigma_eta']), 
  sigma_ksi   = as.numeric(r.box [which.max(r.box $logLik),'sigma_ksi']), 
  mu_beta      = as.numeric(r.box [which.max(r.box $logLik),'mu_beta']),
  alpha       = as.numeric(r.box [which.max(r.box $logLik),'alpha']),
  Beta_0       = as.numeric(r.box [which.max(r.box $logLik),'Beta_0']),
  H_0       = as.numeric(r.box [which.max(r.box $logLik),'H_0'])
)

pf1 <- pfilter(beta.filt,params=params_nowe2,
               Np=betalist[[1]][run_level],filter.traj=T)


wersja3=pf1@filter.traj[3,1,2:(dim(pf1@filter.traj)[3])]
wersja3=as.xts(wersja3,order.by = index(x))
beep(2)







################################################################################
################################################################################
################################################################################
kmnk<-lm(y~x)

#RW prametru beta i alpha nieskorelowane (odpowiada SV wersja 1)

OUss <- function(sigma.alfa, sigma.beta, epsilon){
  Tt <- diag(2)
  Zt <-array(0,dim=c(1,2,n))
  for(i in 1:n) Zt[,,i]=c(1,as.numeric(x[i]))
  ct <- matrix(c(0),ncol=1)
  dt <- matrix(c(0,0), nrow = 2)
  GGt<- matrix(data=c(epsilon^2),nrow = 1,ncol=1)
  HHt<- matrix(data=c(sigma.alfa^2,0,0,sigma.beta^2),nrow=2,ncol=2)
  a0 <-  kmnk$coefficients[1:2]
  P0 <- matrix(data=c(0,0,0,0),nrow=2,ncol=2)
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
              HHt = HHt))
  
}

KF1 <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(as.vector(y), nrow=1,ncol=n))
  return(ans$att[2,])
}

KF1.log <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(as.vector(y), nrow=1,ncol=n))
  return(-ans$logLik)
}
KF1.err <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3] )
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(as.vector(y), nrow=1,ncol=n))
  return(sum(ans$vt[1,]^2))
}



KF1.opt<-optim(c(.1,1,1.1),KF1.log, 
               method="L-BFGS-B",hessian = T )
parametry<-KF1.opt$par
hesjan=sqrt(diag(solve(KF1.opt$hessian)))
log.KF1<-KF1.opt$value
print("SP1:")
print(round(parametry,4))
print(round(hesjan,4))
print(round(-log.KF1,3))





KF1.RW=as.xts(KF1(parametry),order.by = index(x))



################################################################################
################################################################################
################################################################################


#RW prametru beta, stały parametr alpha (odpowiada SV wersja 2)

OUss2 <- function(alpha, sigma.beta, epsilon){
  Tt <- diag(1)
  Zt <-array(0,dim=c(1,1,n))
  for(i in 1:n) Zt[,,i]=as.numeric(x[i])
  ct <- matrix(c(alpha),ncol=1)
  dt <- matrix(c(0), nrow = 1)
  GGt<- matrix(data=c(epsilon^2),nrow = 1,ncol=1)
  HHt<- matrix(data=c(sigma.beta^2),nrow=1,ncol=1)
  a0 <- kmnk$coefficients[2]
  P0 <- matrix(data=c(0),nrow=1,ncol=1)
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
              HHt = HHt))
  
}

KF2 <- function(theta) {
  sp <- OUss2(theta[1], theta[2], theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(as.vector(y), nrow=1,ncol=n))
  return(ans$att[1,])
}

KF2.log <- function(theta) {
  sp <- OUss2(theta[1], theta[2], theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(as.vector(y), nrow=1,ncol=n))
  return(-ans$logLik)
}

KF2.err <- function(theta) {
  sp <- OUss2(theta[1], theta[2], theta[3] )
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt,yt =matrix(as.vector(y), nrow=1,ncol=n))
  return(sum(ans$vt[1,]^2))
}


KF2.opt<-optim(c(0.005,.1,1.1),KF2.log, 
               method="L-BFGS-B",hessian = T,lower = c(-Inf,0,0))
parametry2<-KF2.opt$par
hesjan2=sqrt(diag(solve(KF2.opt$hessian)))
log.KF2<-KF2.opt$value
print("SP2:")
print(round(parametry2,4))
print(round(hesjan2,4))  
print(round(-log.KF2,3))


KF2.RW=as.xts(KF2(parametry2),order.by = index(x))






################################################################################
################################################################################
################################################################################


OUss5 <- function(alpha,rho, beta, sigma.beta, epsilon){
  Tt <- rho*diag(1)
  Zt <-array(0,dim=c(1,1,n))
  for(i in 1:n) Zt[,,i]=as.numeric(x[i])
  ct <- matrix(c(alpha),ncol=1)
  dt <- matrix(c((1-rho)*beta), nrow = 1)
  GGt<- matrix(data=c(epsilon^2),nrow = 1,ncol=1)
  HHt<- matrix(data=c(sigma.beta^2),nrow=1,ncol=1)
  a0 <- kmnk$coefficients[2]
  P0 <- matrix(data=c(0),nrow=1,ncol=1)
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
              HHt = HHt))
  
}

KF5 <- function(theta) {
  sp <- OUss5(theta[1], theta[2], theta[3], theta[4], theta[5])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt,  yt =matrix(as.vector(y), nrow=1,ncol=n))
  return(ans$att[1,])
}

KF5.log <- function(theta) {
  sp <- OUss5(theta[1], theta[2], theta[3], theta[4], theta[5])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt,  yt =matrix(as.vector(y), nrow=1,ncol=n))
  return(-ans$logLik)
}


KF5.err <- function(theta) {
  sp <- OUss5(theta[1], theta[2], theta[3], theta[4], theta[5])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt,  yt =matrix(as.vector(y), nrow=1,ncol=n))
  return(sum(ans$vt[1,]^2))
}


KF5.opt<-optim(c(0.005,0.2,1.15,0.5,1.1),KF5.log, 
               method="L-BFGS-B",hessian = T,lower=c(-Inf,-Inf,-Inf,0,0))
parametry5<-KF5.opt$par
hesjan5=sqrt(diag(solve(KF5.opt$hessian)))
log.KF5<-KF5.opt$value
print("SP5:")
print(round(parametry5,4))
print(round(hesjan5,4))  
print(round(-log.KF5,3))



KF3.RW=as.xts(KF5(parametry5),order.by = index(x))

beep(2)
################################################################################
################################################################################
################################################################################

save(wersja1,wersja2,wersja3,KF1.RW,KF2.RW,KF3.RW,file='wersjePFKf')



plot(wersja1,main=expression(beta),major.ticks = "years",grid.ticks.on = "years")
plot(wersja2,main=expression(beta),major.ticks = "years",grid.ticks.on = "years")
plot(wersja3,main=expression(beta),major.ticks = "years",grid.ticks.on = "years")
plot(KF1.RW,main=expression(beta),major.ticks = "years",grid.ticks.on = "years")
plot(KF2.RW,main=expression(beta),major.ticks = "years",grid.ticks.on = "years")
plot(KF3.RW,main=expression(beta),major.ticks = "years",grid.ticks.on = "years")

plot(wersja1,main=expression(beta),major.ticks = "years",grid.ticks.on = "years")
lines(wersj2,col='red')
lines(wersja3,col='blue')
lines(KF1.RW,col='orange')
lines(KF2.RW,col='green')
lines(KF3.RW,col='brown')