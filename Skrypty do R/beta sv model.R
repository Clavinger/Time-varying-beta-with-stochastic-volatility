library(pomp)
library(beepr)
library(doParallel)
library(doSNOW)
#rm(list = ls())
n=147

#dane za okres 2013-04-22 do 2018-02-09
dane2<-read.csv('C:/Users/user/Dropbox/phd/Dane/wig_m.csv', header = T)
str(dane2)
dane3<-matrix(dane2$Zamkniecie,nrow=n,ncol=1)
str(dane3)
 



dane4<-read.csv('C:/Users/user/Dropbox/phd/Dane/bzw_m.csv', header = T)
str(dane4)
dane5<-matrix(dane4$Zamkniecie,nrow=n,ncol=1)
dane=matrix(0,nrow=n,ncol=2)
dane[,1]=dane3[,1]
dane[,2]=dane5[,1]
dim(dane)
lr.wbk=1:(n-1)
lr.wig=1: (n-1)
for(i in 2:(n)){
  lr.wbk[i-1]=(log(dane[i,2])-log(dane[i-1,2]))*100
  lr.wig[i-1]=(log(dane[i,1])-log(dane[i-1,1]))*100
}

y=lr.wbk
x=lr.wig

n=length(lr.wbk)

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
beta.filt <- pomp(data=data.frame(y=lr.wbk,
                                 time=1:length(lr.wbk)),
                 statenames=beta_statenames,
                 paramnames=beta_paramnames,
                 covarnames=beta_covarnames,
                 times="time",
                 t0=0,
                 covar=data.frame(R_i_dane=c(0,lr.wbk) ,R_M=c(0,lr.wig),
                                  time=0:length(lr.wbk)),
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
  H_0=0
)


pf1 <- pfilter(beta.filt,params=params_test,
               Np=1000,filter.traj=T)


plot(pf1)

plot(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),type='l',col='red',lty=2)
legend(x='topleft',legend = c('true','filtered'),
       col=c('black','red'),lty=c(1,2))


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
beep(2)
plot(if.beta.box )


#save(if.beta.box,L.beta.box,H.beta.box, file="beta_box_eval.rda")
#load( file="beta_box_eval.rda")

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
plot(pf1)
par(mfrow=c(1,1))
plot(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),type='l')
plot(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])],type='l',ylab="H")
plot(pf1@filter.traj[3,1,2:(dim(pf1@filter.traj)[3])],type='l',ylab="Beta")
wersja1=pf1@filter.traj[3,1,2:(dim(pf1@filter.traj)[3])]
