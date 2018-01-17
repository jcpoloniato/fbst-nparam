setwd('/home/joao_polo/Dropbox/Disciplinas/Mestrado/Dissertação/fbst-nparam')
source('funcoesuteis.R') 


#evalor modificado

ev.mod<-function(dados,phi,r,peso,simu=1000){
  n<-length(dados)
  tau<-peso^(1:r)
  x<-phi_barra(dados,phi,r)
  aux<-evalue(x,n,tau)
  
  q<-rep(NA, simu)
  for (i in 1:simu){
    q[i] <- sum(((n*rnorm(r,sd=1/sqrt(n)))/sqrt(n+tau))^2)
  }
  
  1 - mean(q < aux)
}

ev.nmod<-function(dados,phi,r,peso){
  n<-length(dados)
  tau<-peso^(1:r)
  x<-phi_barra(dados,phi,r)
  aux<-evalue(x,n,tau)

  1 - pchisq(aux, df=r)
}

#comparacao entre evalor modificado e nao modificado 

#H_0 falsa
evm<-rep(NA,1000)
evnm<-rep(NA,1000)
for (i in 1:1000){
  dados<-rnorm(10000,mean=5,sd=1)
  evm[i]<-ev.mod(dados,phi,50,5)
  evnm[i]<-ev.nmod(dados,phi,50,5)
}
par(mfrow=c(1,2))
hist(evs)
hist(nevs)

#H_0 verdadeira
evm1<-rep(NA,1000)
evnm1<-rep(NA,1000)
for (i in 1:1000){
  dados<-runif(10000)
  evm1[i]<-ev.mod(dados,phi,50,5)
  evnm1[i]<-ev.nmod(dados,phi,50,5)
}
hist(evm1)
hist(evnm1)