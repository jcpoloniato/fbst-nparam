source('funcoes uteis.R') 

#### ESSE DEU CERTO \o/ 

#######################################################
### METODO 2 - FAZENDO BOOTSTRAP PARA DEFINIR O CORTE #
#######################################################

##funcoes
empirica<-function(dados1, dados2){
  
  n<-length(dados1)
  m<-length(dados2)
  
  dado <- matrix(c(rep(1,n),rep(2,m),dados1,dados2),ncol=2,byrow = F)
  
  dado[,1]<-dado[,1][order(dado[,2])] #order pegaa posicao em que esta o elemento no vetor se o vetor for ordenado
  dado[,2]<-sort(dado[,2]) #ordenar a coluna em que estao os dados
  
  count=0
  aux<-rep(NA,n)
  j<-1
  for (i in 1:(n+m)){
    if(dado[i,1]==1){
      aux[j]<-count/m
      j<-j+1
    }
    else
      count=count+1
  }
  
  aux
}

calculo_ev_mod_boot<-function(dados1, dados2, phi, r=50, peso = function(x) x^2, simu=1000){
  n1<-length(dados1)
  n2<-length(dados2)  
  dados_transf<- empirica(dados1,dados2)
  n <- length(dados_transf)
  tau <- peso(1:r)
  x <- phi_barra(dados_transf, phi, r)
  aux <- testeFBST(x, n, tau)
  
  q <- rep(NA, simu)
  junto<-c(dados1,dados2)
  for (i in 1:simu){
    #resample <- sample(junto)
    dat1<-sample(junto, n1, replace =T)
    dat2<- sample(junto, n2, replace =T)
    AUX <- empirica(dat1,dat2)  
    X <- phi_barra(AUX, phi, r)  
    q[i] <- testeFBST(X, n, tau)
  }
  
  1 - mean(q < aux)
}

simulacao<-1000

##########
# NORMAL #
##########
medias<-seq(-2,2,by=0.1)

###2 pop Normal 50
poder <- matrix(NA, length(medias),3)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=3,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rnorm(50)
    data2<-rnorm(50, mean=medias[b])
    AUX <- empirica(data1,data2)
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- t.test(x=data1,y=data2, var.equal = T)$p.value
    resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias, medias),
                poder=c(poder[,1], poder[,2], poder[,3]), 
                teste=c(rep('ev', length(medias)), rep('pv.trv', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popnormal_tam50-50')

###2 pop Normal 100
poder <- matrix(NA, length(medias),3)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=3,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rnorm(100)
    data2<-rnorm(100, mean=medias[b])
    AUX <- empirica(data1,data2)
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- t.test(x=data1,y=data2, var.equal = T)$p.value
    resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias, medias),
                poder=c(poder[,1], poder[,2], poder[,3]), 
                teste=c(rep('ev', length(medias)), rep('pv.trv', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popnormal_tam100-100')

## 2 pop Normal 500
poder <- matrix(NA, length(medias),3)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=3,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rnorm(500)
    data2<-rnorm(500, mean=medias[b])
    AUX <- empirica(data1,data2)
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- t.test(x=data1,y=data2, var.equal = T)$p.value
    resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias, medias),
                poder=c(poder[,1], poder[,2], poder[,3]), 
                teste=c(rep('ev', length(medias)), rep('pv.trv', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popnormal_tam500-500')


###############
# EXPONENCIAL #
###############
medias<-seq(1,5,by=0.1)

###2 pop exp 50 
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rexp(50, rate = 2)
    data2<-rexp(50, rate =medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- ks.test(data1,data2)$'p.value'
    #resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poder[,1], poder[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popexp_tam50-50')

###2 pop exp 100
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rexp(100, rate = 2)
    data2<-rexp(100, rate =medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- ks.test(data1,data2)$'p.value'
    #resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poder[,1], poder[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popexp_tam100-100')

###2 pop exp 500 
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rexp(500, rate = 2)
    data2<-rexp(500, rate =medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- ks.test(data1,data2)$'p.value'
    #resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poder[,1], poder[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popexp_tam500-500')


##################
# TRIGONOMETRICA #
##################

densidade<-function(x,k){
  abs(sin(pi*k*x)*pi/2)
}
gerar<-function(n,k){
  gerado<-rep(NA,n)
  i<-1
  while(is.na(gerado[n])){
    aux1<- runif(1)
    aux2<- runif(1,0,pi/2)
    
    if(densidade(aux1,k) < aux2){
      gerado[i]<-aux1
      i<-i+1
    }
  }
  gerado
}
medias<-seq(1,30,by=1)

###2 pop trigo 50 
poder <- matrix(NA, length(medias),5)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=5,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-gerar(50, 5)
    data2<-gerar(50, medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^2)
    resul.teste[a,2] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^5)
    resul.teste[a,3] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 2^x)
    resul.teste[a,4] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 5^x)
    resul.teste[a,5] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=rep(1:30,5),
                poder=c(poder[,1], poder[,2], poder[,3], poder[,4], poder[,5]) ,
                teste=c(rep('x2', length(medias)), rep('x5', length(medias)), rep('2x', length(medias)), 
                        rep('5x', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2poptrigo_tam50-50')

###2 pop trigo 100 
poder <- matrix(NA, length(medias),5)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=5,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-gerar(50, 5)
    data2<-gerar(50, medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^2)
    resul.teste[a,2] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^5)
    resul.teste[a,3] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 2^x)
    resul.teste[a,4] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 5^x)
    resul.teste[a,5] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=rep(1:30,5),
                poder=c(poder[,1], poder[,2], poder[,3], poder[,4], poder[,5]) ,
                teste=c(rep('x2', length(medias)), rep('x5', length(medias)), rep('2x', length(medias)), 
                        rep('5x', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2poptrigo_tam100-100')

###2 pop trigo 500 
poder <- matrix(NA, length(medias),5)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=5,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-gerar(500, 5)
    data2<-gerar(500, medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^2)
    resul.teste[a,2] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^5)
    resul.teste[a,3] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 2^x)
    resul.teste[a,4] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 5^x)
    resul.teste[a,5] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=rep(1:30,5),
                poder=c(poder[,1], poder[,2], poder[,3], poder[,4], poder[,5]) ,
                teste=c(rep('x2', length(medias)), rep('x5', length(medias)), rep('2x', length(medias)), 
                        rep('5x', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2poptrigo_tam500-500')

#####################
# UNIFORME DE HAAR ##
#####################

unif.haar<- function(x, k) floor(k*x) %% 2
runif.haar<- function(n,k){
  gerado<- rep(NA,n)
  i<-1
  while(is.na(gerado[n])){
    aux<-runif(1)
    if(unif.haar(aux,k)){
      gerado[i]<-aux
      i <- i+1
    }
  }
  gerado
}
medias<-2^seq(1,30,by=1)

##2 pop haar 50
poder <- matrix(NA, length(medias),5)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=5,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-runif.haar(50, medias[4])
    data2<-runif.haar(50, medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^2)
    resul.teste[a,2] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^5)
    resul.teste[a,3] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 2^x)
    resul.teste[a,4] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 5^x)
    resul.teste[a,5] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=rep(1:30,5),
                poder=c(poder[,1], poder[,2], poder[,3], poder[,4], poder[,5]) ,
                teste=c(rep('x2', length(medias)), rep('x5', length(medias)), rep('2x', length(medias)), 
                        rep('5x', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2pophaar_tam50-50')

##2 pop haar 100
poder <- matrix(NA, length(medias),5)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=5,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-runif.haar(100, medias[4])
    data2<-runif.haar(100, medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^2)
    resul.teste[a,2] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^5)
    resul.teste[a,3] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 2^x)
    resul.teste[a,4] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 5^x)
    resul.teste[a,5] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=rep(1:30,5),
                poder=c(poder[,1], poder[,2], poder[,3], poder[,4], poder[,5]) ,
                teste=c(rep('x2', length(medias)), rep('x5', length(medias)), rep('2x', length(medias)), 
                        rep('5x', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2pophaar_tam100-100')

##2 pop haar 500
poder <- matrix(NA, length(medias),5)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=5,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-runif.haar(500, medias[4])
    data2<-runif.haar(500, medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^2)
    resul.teste[a,2] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) x^5)
    resul.teste[a,3] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 2^x)
    resul.teste[a,4] <- calculo_ev_mod_boot(data1,data2,phi=phi,peso=function(x) 5^x)
    resul.teste[a,5] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=rep(1:30,5),
                poder=c(poder[,1], poder[,2], poder[,3], poder[,4], poder[,5]) ,
                teste=c(rep('x2', length(medias)), rep('x5', length(medias)), rep('2x', length(medias)), 
                        rep('5x', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2pophaar_tam500-500')

########################
# MISTURA DE NORMAIS ###
########################

mist_norm <- function (n, medias, variancias, probs){
  k <- length(probs)
  vetor<-rep(NA,n)
  
  aux<-rmultinom(n, size = 1, prob = probs)
  delta<-rep(NA,n)
  j<-1
  while (j < (k+1)){
    delta[aux[j,]==1]<-j
    j <- j+1
  }
  
  j<-1
  while(j < (k+1)){
    vetor[delta==j]<-rnorm(sum(delta==j),medias[j],variancias[j])
    j<-j+1
  }
  
  vetor
}
#mistura de 2 normais
variancias <- c(1,1)
medias<-matrix(c(-seq(0,3,by=0.1),seq(0,3,by=0.1)),ncol=2)

#n = 50
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-mist_norm(50, medias=c(-1,1), variancias=c(1,1), probs=c(0.5,0.5))
    data2<-mist_norm(50, medias=medias[b,], variancias=c(1,1), probs=c(0.5,0.5))
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias[,2],medias[,2]),
                poder=c(poder[,1], poder[,2]) ,
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popmistnorm2_tam50-50')

#n = 100
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-mist_norm(100, medias=c(-1,1), variancias=c(1,1), probs=c(0.5,0.5))
    data2<-mist_norm(100, medias=medias[b,], variancias=c(1,1), probs=c(0.5,0.5))
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias[,2],medias[,2]),
                poder=c(poder[,1], poder[,2]) ,
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popmistnorm2_tam100-100')

#n = 500
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-mist_norm(500, medias=c(-1,1), variancias=c(1,1), probs=c(0.5,0.5))
    data2<-mist_norm(500, medias=medias[b,], variancias=c(1,1), probs=c(0.5,0.5))
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias[,2],medias[,2]),
                poder=c(poder[,1], poder[,2]) ,
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popmistnorm2_tam500-500')

#mistura de 3 normais
medias<-matrix(c(-seq(0,3,by=0.1),seq(0,3,by=0.1),rep(0,length(seq(0,3,by=0.1)))),ncol=3)

#n = 50
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-mist_norm(50, medias=c(-1,1,0), variancias=c(1,1,1), probs=c(1/3,1/3,1/3))
    data2<-mist_norm(50, medias=medias[b,], variancias=c(1,1,1), probs=c(1/3,1/3,1/3))
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias[,2],medias[,2]),
                poder=c(poder[,1], poder[,2]) ,
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popmistnorm3_tam50-50')

#n = 100
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-mist_norm(100, medias=c(-1,1,0), variancias=c(1,1,1), probs=c(1/3,1/3,1/3))
    data2<-mist_norm(100, medias=medias[b,], variancias=c(1,1,1), probs=c(1/3,1/3,1/3))
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias[,2],medias[,2]),
                poder=c(poder[,1], poder[,2]) ,
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popmistnorm3_tam100-100')

#n = 500
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-mist_norm(500, medias=c(-1,1,0), variancias=c(1,1,1), probs=c(1/3,1/3,1/3))
    data2<-mist_norm(500, medias=medias[b,], variancias=c(1,1,1), probs=c(1/3,1/3,1/3))
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2]<-ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias[,2],medias[,2]),
                poder=c(poder[,1], poder[,2]) ,
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popmistnorm3_tam500-500')


########################
### BETA   #############
########################

##beta 11
medias<-seq(0.5,5,by=0.1)

##n 50 
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rbeta(50, shape1 = 1, shape2=1)
    data2<-rbeta(50,shape1 = medias[b], shape2=medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poder[,1], poder[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popbeta11_tam50-50')

###n 100
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rbeta(100, shape1 = 1, shape2=1)
    data2<-rbeta(100,shape1 = medias[b], shape2=medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poder[,1], poder[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popbeta11_tam100-100')

###n 500 
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rbeta(500, shape1 = 1, shape2=1)
    data2<-rbeta(500,shape1 = medias[b], shape2=medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poder[,1], poder[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popbeta11_tam500-500')

##beta 1/2 1/2

##n 50 
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rbeta(50, shape1 = 1/2, shape2=1/2)
    data2<-rbeta(50,shape1 = medias[b], shape2=medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poder[,1], poder[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popbetameio_tam50-50')

##n 100 
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rbeta(100, shape1 = 1/2, shape2=1/2)
    data2<-rbeta(100,shape1 = medias[b], shape2=medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poder[,1], poder[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popbetameio_tam100-100')

##n 500
poder <- matrix(NA, length(medias),2)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=2,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rbeta(500, shape1 = 1/2, shape2=1/2)
    data2<-rbeta(500,shape1 = medias[b], shape2=medias[b])
    resul.teste[a,1] <- calculo_ev_mod_boot(data1,data2,phi=phi)
    resul.teste[a,2] <- ks.test(data1,data2)$'p.value'
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias),
                poder=c(poder[,1], poder[,2]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks', length(medias))))
saveRDS(dat, file='2popbetameio_tam500-500')