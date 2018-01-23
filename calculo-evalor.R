setwd('/home/joao_polo/Dropbox/Disciplinas/Mestrado/Dissertação/fbst-nparam')
source('funcoes-uteis.R') 

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
hist(evm)
hist(evnm)

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

####################################
# gerar dados e aplicar uma F nos dados e 
# verificar se os dados modificado tem distribuicao uniforme

# funcao que vai calcular o e valor nos dados transformados

ev.mod.pf<-function(dados, phi, r, peso, pf){
  dados.inv <- inversa(dados, pf)
  ev.mod(dados.inv, phi, r, peso)
}

library(nortest)
# funcao que vai calcular o pvalor do ks nos dados transformados

pvalor.ks<-function(dados, pf){
  dados.inv <- inversa(dados, pf)
  Test <- ks.test(dados.inv, punif)
  round(Test$'p.value', 4)
}

# essa funcao compara o ks com o evalor, por isso dados deve ser uma matriz de dados
evalor.vs.ks<-function(dados, phi, r, peso, pf){
  
  ev<-apply(dados,2,ev.mod.pf,phi=phi,r= r,peso= peso,pf= pf)
  pv<-apply(dados,2,pvalor.ks,pf=pf)
  cbind(ev,pv)
}

dados <- matrix(rnorm(10000,0,1),ncol=100)
evalor.vs.ks(dados,phi,50,5,pnorm)
evalor.vs.ks(dados,phi,50,5,pexp)

dados <- matrix(rnorm(10000,5,1),ncol=100)
evalor.vs.ks(dados,phi,50,5,pnorm)

dados<- matrix(rgamma(10000,5,1),ncol=100)
evalor.vs.ks(dados,phi,50,5,pnorm)
evalor.vs.ks(dados,phi,50,5,pgamma)
