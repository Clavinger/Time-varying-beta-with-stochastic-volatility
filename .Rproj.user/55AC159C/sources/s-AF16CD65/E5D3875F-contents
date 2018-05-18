library(pomp)
library(beepr)
library(doParallel)
library(doSNOW)


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

plot(beta.filt)

#parametry
mu = -0.5
phi = 0.98
sigma = 0.2

params_test <- c(
  mu_h = mu,       
  phi = phi,     
  sigma_eta =sigma,
  sigma_ksi=sigma,
  sigma_gamma=sigma,
  Beta_0=1.5,
  Alpha_0=0.5,  
  H_0=0
)


pf1 <- pfilter(beta.filt,params=params_test,
               Np=1000,filter.traj=T)


plot(pf1)




###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 3

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
beep(2)
plot(if.beta.box )

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
plot(pf1)
par(mfrow=c(1,1))
plot(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),type='l')
plot(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])],type='l',ylab="H")
plot(pf1@filter.traj[3,1,2:(dim(pf1@filter.traj)[3])],type='l',ylab="Beta")
plot(pf1@filter.traj[4,1,2:(dim(pf1@filter.traj)[3])],type='l',ylab="Alpha")
wersja2=pf1@filter.traj[3,1,2:(dim(pf1@filter.traj)[3])]