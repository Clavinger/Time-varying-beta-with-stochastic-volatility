library(pomp)
library(beepr)
library(doParallel)
library(doSNOW)
library(quantmod)

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
  alpha=0,
  Beta_0=0.5,
  nu=phi,
  mu_beta=mu,
  H_0=0
)


pf1 <- pfilter(beta.filt,params=params_test,
               Np=1000,filter.traj=T)


plot(pf1)



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
beep(2)
plot(if.beta.box )


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
plot(pf1)
par(mfrow=c(1,1))
plot(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),type='l')
plot(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])],type='l',ylab="H")
plot(pf1@filter.traj[3,1,2:(dim(pf1@filter.traj)[3])],type='l',ylab="Beta")

wersja3=pf1@filter.traj[3,1,2:(dim(pf1@filter.traj)[3])]
wersja3=as.xts(wersja3,order.by = index(x))



par(mfrow=c(3,1))
plot(x,main='WIG',major.ticks = "years",grid.ticks.on = "years")
plot(y,main=tiker,major.ticks = "years",grid.ticks.on = "years")
plot(wersja3,main=expression(beta),major.ticks = "years",grid.ticks.on = "years")
par(mfrow=c(1,1))




#parametr mu_h
xx1<-seq(from=params_nowe2['mu_h']-1,to=params_nowe2['mu_h']+1,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.beta.log<- foreach(i=1:length(xx1) ,.packages='pomp', .export = "betalist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logLik(pfilter(beta.filt,params=c(mu_h=xx1[i], 
                                                         phi=as.numeric(params_nowe2['phi']),
                                                         sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                         sigma_ksi =as.numeric(params_nowe2['sigma_ksi']),
                                                         nu =as.numeric(params_nowe2['nu']),
                                                         mu_beta=as.numeric(params_nowe2['mu_beta']),
                                                         Beta_0=as.numeric(params_nowe2['Beta_0']),
                                                         alpha=as.numeric(params_nowe2['alpha']),
                                                         H_0=as.numeric(params_nowe2['H_0'])),
                                      Np=betalist[[1]][run_level] ))
                     }

stopCluster(cl)



plot(xx1, L.beta.log, type='l',xlab=expression(mu[h]),ylab="logLik")
points(xx1, L.beta.log)
points(r.box[,'mu_h'], r.box[,'logLik'] ,col='red')
p=loess(L.beta.log~xx1,span=0.5)
lines(xx1,p$fitted,col='blue',lwd=2)
abline(v=params_nowe2['mu_h'],lty=2)
beep(2)


#parametr phi
xx2<-seq(from=0.95,to=0.999,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.beta.log2<- foreach(i=1:length(xx2) ,.packages='pomp', .export = "betalist",.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932)
                        logLik(pfilter(beta.filt,params=c(mu_h=as.numeric(params_nowe2['phi']),
                                                          phi=as.numeric(xx2[i]), 
                                                          sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                          sigma_ksi =as.numeric(params_nowe2['sigma_ksi']),
                                                          nu =as.numeric(params_nowe2['nu']),
                                                          mu_beta=as.numeric(params_nowe2['mu_beta']),
                                                          Beta_0=as.numeric(params_nowe2['Beta_0']),
                                                          alpha=as.numeric(params_nowe2['alpha']),
                                                          H_0=as.numeric(params_nowe2['H_0'])),
                                       Np=betalist[[1]][run_level] ))
                      }

stopCluster(cl)



plot(xx2, L.beta.log2, type='l',xlab=expression(phi),ylab="logLik")
points(xx2, L.beta.log2)
points(r.box[,'phi'], r.box[,'logLik'] ,col='red')
p2=loess(L.beta.log2~xx2,span=0.5)
lines(xx2,p2$fitted,col='blue',lwd=2)
abline(v=params_nowe2['phi'],lty=2)
beep(2)


#parametr sigma_eta
xx3<-seq(from=0.01,to=.5,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.beta.log3<- foreach(i=1:length(xx3) ,.packages='pomp', .export = "betalist",.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932)
                        logLik(pfilter(beta.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                          phi=as.numeric(params_nowe2['phi']),
                                                          sigma_eta=xx3[i],
                                                          sigma_ksi =as.numeric(params_nowe2['sigma_ksi']),
                                                          nu =as.numeric(params_nowe2['nu']),
                                                          mu_beta=as.numeric(params_nowe2['mu_beta']),
                                                          Beta_0=as.numeric(params_nowe2['Beta_0']),
                                                          alpha=as.numeric(params_nowe2['alpha']),
                                                          H_0=as.numeric(params_nowe2['H_0'])),
                                       Np=betalist[[1]][run_level] ))
                      }

stopCluster(cl)



plot(xx3, L.beta.log3, type='l',xlab=expression(sigma[eta]),ylab="logLik")
points(xx3, L.beta.log3)
points(r.box[,'sigma_eta'], r.box[,'logLik'] ,col='red')
p3=loess(L.beta.log3~xx3,span=0.25)
lines(xx3,p3$fitted,col='blue',lwd=2)
abline(v=params_nowe2['sigma_eta'],lty=2)
beep(2)



#parametr sigma_ksi
xx4<-seq(from=0.01,to=.025,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.beta.log4<- foreach(i=1:length(xx4) ,.packages='pomp', .export = "betalist",.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932)
                        logLik(pfilter(beta.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                          phi=as.numeric(params_nowe2['phi']),
                                                          sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                          sigma_ksi =xx4[i],
                                                          nu =as.numeric(params_nowe2['nu']),
                                                          mu_beta=as.numeric(params_nowe2['mu_beta']),
                                                          Beta_0=as.numeric(params_nowe2['Beta_0']),
                                                          alpha=as.numeric(params_nowe2['alpha']),
                                                          H_0=as.numeric(params_nowe2['H_0'])),
                                       Np=betalist[[1]][run_level] ))
                      }

stopCluster(cl)



plot(xx4, L.beta.log4, type='l',xlab=expression(sigma[xi]),ylab="logLik")
points(xx4, L.beta.log4)
points(r.box[,'sigma_ksi'], r.box[,'logLik'] ,col='red')
p4=loess(L.beta.log4~xx4,span=0.25)
lines(xx4,p4$fitted,col='blue',lwd=2)
abline(v=params_nowe2['sigma_ksi'],lty=2)
beep(2)



#parametr alpha
xx5<-seq(from=-.1,to=.1,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.beta.log5<- foreach(i=1:length(xx5) ,.packages='pomp', .export = "betalist",.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932)
                        logLik(pfilter(beta.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                          phi=as.numeric(params_nowe2['phi']),
                                                          sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                          sigma_ksi =as.numeric(params_nowe2['sigma_ksi']),
                                                          nu =as.numeric(params_nowe2['nu']),
                                                          mu_beta=as.numeric(params_nowe2['mu_beta']),
                                                          Beta_0=as.numeric(params_nowe2['Beta_0']),
                                                          alpha=xx5[i],
                                                          H_0=as.numeric(params_nowe2['H_0'])),
                                       Np=betalist[[1]][run_level] ))
                      }

stopCluster(cl)



plot(xx5, L.beta.log5, type='l',xlab=expression(alpha),ylab="logLik")
points(xx5, L.beta.log5)
points(r.box[,'alpha'], r.box[,'logLik'] ,col='red')
p5=loess(L.beta.log5~xx5,span=0.5)
lines(xx5,p5$fitted,col='blue',lwd=2)
abline(v=params_nowe2['alpha'],lty=2)
beep(2)





#parametr Beta_0
xx6<-seq(from=0.1,to=0.5,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.beta.log6<- foreach(i=1:length(xx6) ,.packages='pomp', .export = "betalist",.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932)
                        logLik(pfilter(beta.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                          phi=as.numeric(params_nowe2['phi']),
                                                          sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                          sigma_ksi =as.numeric(params_nowe2['sigma_ksi']),
                                                          nu =as.numeric(params_nowe2['nu']),
                                                          mu_beta=as.numeric(params_nowe2['mu_beta']),
                                                          Beta_0=xx6[i],
                                                          alpha=as.numeric(params_nowe2['alpha']),
                                                          H_0=as.numeric(params_nowe2['H_0'])),
                                       Np=betalist[[1]][run_level] ))
                      }

stopCluster(cl)



plot(xx6, L.beta.log6, type='l',xlab=expression(beta[0]),ylab="logLik")
points(xx6, L.beta.log6)
points(r.box[,'Beta_0'], r.box[,'logLik'] ,col='red')
p6=loess(L.beta.log6~xx6,span=0.5)
lines(xx6,p6$fitted,col='blue',lwd=2)
abline(v=params_nowe2['Beta_0'],lty=2)
beep(2)


#parametr H_0
xx7<-seq(from=-.5,to=-0.1,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.beta.log7<- foreach(i=1:length(xx7) ,.packages='pomp', .export = "betalist",.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932)
                        logLik(pfilter(beta.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                          phi=as.numeric(params_nowe2['phi']),
                                                          sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                          sigma_ksi =as.numeric(params_nowe2['sigma_ksi']),
                                                          nu =as.numeric(params_nowe2['nu']),
                                                          mu_beta=as.numeric(params_nowe2['mu_beta']),
                                                          Beta_0=as.numeric(params_nowe2['Beta_0']),
                                                          alpha=as.numeric(params_nowe2['alpha']),
                                                          H_0=xx7[i]),
                                       Np=betalist[[1]][run_level] ))
                      }

stopCluster(cl)



plot(xx7, L.beta.log7, type='l',xlab=expression(H[0]),ylab="logLik")
points(xx7, L.beta.log7)
points(r.box[,'H_0'], r.box[,'logLik'] ,col='red')
p7=loess(L.beta.log7~xx7,span=0.5)
lines(xx7,p7$fitted,col='blue',lwd=2)
abline(v=params_nowe2['H_0'],lty=2)
beep(2)




#parametr nu
xx9<-seq(from=.8,to=.99,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.beta.log9<- foreach(i=1:length(xx5) ,.packages='pomp', .export = "betalist",.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932)
                        logLik(pfilter(beta.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                          phi=as.numeric(params_nowe2['phi']),
                                                          sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                          sigma_ksi =as.numeric(params_nowe2['sigma_ksi']),
                                                          nu =as.numeric(params_nowe2['nu']),
                                                          mu_beta=xx9[i],
                                                          Beta_0=as.numeric(params_nowe2['Beta_0']),
                                                          alpha=as.numeric(params_nowe2['alpha']),
                                                          H_0=as.numeric(params_nowe2['H_0'])),
                                       Np=betalist[[1]][run_level] ))
                      }

stopCluster(cl)



plot(xx9, L.beta.log9, type='l',xlab=expression(nu),ylab="logLik")
points(xx9, L.beta.log9)
points(r.box[,'mu_beta'], r.box[,'logLik'] ,col='red')
p9=loess(L.beta.log9~xx9,span=0.5)
lines(xx9,p9$fitted,col='blue',lwd=2)
abline(v=params_nowe2['mu_beta'],lty=2)
beep(2)