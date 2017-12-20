#evalor modificado

ev.mod<-function(dados,phi,r,peso,simu){
  n<-length(dados)
  tau<-peso^(1:r)
  x<-phi_barra(dados,phi,r)
  
  aux<-sum((n*x/sqrt(n+tau))^2)
  
  q<-rep(NA,simu)
  for (i in 1:simu){
    q[i]<-sum(((n*rnorm(r,sd=1/sqrt(n)))/sqrt(n+tau))^2)
  }
  
  prob<-sum(q<aux)/simu
  ev<-1-prob
  return(ev)
}


ev.nmod<-function(dados,phi,r,peso){
  n<-length(dados)
  tau<-peso^(1:r)
  x<-phi_barra(dados,phi,r)
  mu0<-rep(0,r)
  
  mu<-n/(n+tau)*x
  sigma<-diag(n+tau,r)
    
  c0<- t(mu0-mu)%*%sigma%*%(mu0-mu)
  ev<-1-pchisq(c0,df=r)
  return(ev)
}


#H_0 é falsa
evs<-rep(NA,1000)
nevs<-rep(NA,1000)
for (i in 1:1000){
  dados<-rnorm(10000,mean=5,sd=1)
  evs[i]<-ev.mod(dados,phi,50,5,1000)
  nevs[i]<-ev.nmod(dados,phi,50,5)
}
hist(evs)
hist(nevs)

#H_0 é verdadeira
evs<-rep(NA,1000)
nevs<-rep(NA,1000)
for (i in 1:1000){
  dados<-runif(10000)
  evs[i]<-ev.mod(dados,phi,50,5,1000)
  nevs[i]<-ev.nmod(dados,phi,50,5)
}
hist(nevs)
hist(evs)


#evalor nao modificado vs evalor modificado, quem ganha? sob H0 e H1.