setwd('/home/joao_polo/Dropbox/Mestrado/Dissertacao/fbst-nparam')
source('funcoes-uteis.R') 


##funcao que aplica a empirica de uma amostra na outra outra amostra

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

calculo_ev_mod<-function(dados_transf, phi, r=50, peso = function(x) x^2, simu=1000){
  n <- length(dados_transf)
  tau <- peso(1:r)
  x <- phi_barra(dados_transf, phi, r)
  aux <- testeFBST(x, n, tau)
  
  q <- rep(NA, simu)
  for (i in 1:simu){
    q[i] <- sum(((n*rnorm(r,sd=1/sqrt(n)))/sqrt(n+tau))^2)
  }
  
  1 - mean(q < aux)
} 

ks<-function(vetor){
  ks.test(vetor, punif)$'p.value'
}

##gerando dados de duas pop normais e testando o poder sob h0
  simulacao<-500
  resul.teste<-matrix(NA,ncol=3,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rnorm(100)
    data2<-rnorm(1000)
    AUX <- empirica(data1,data2)
    resul.teste[a,1] <- calculo_ev_mod(AUX,phi=phi)
    resul.teste[a,2] <- ks(AUX)
    resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
    print(a)
  }
  colMeans(resul.teste<0.05)

##gerando dados de duas pop normais e tracando a funcao poder 
simulacao<-500
medias<-seq(-1,1,by=0.1)
poder <- matrix(NA, length(medias),3)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=3,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rnorm(100)
    data2<-rnorm(1000, mean=medias[b])
    AUX <- empirica(data1,data2)
    resul.teste[a,1] <- calculo_ev_mod(AUX,phi=phi)
    resul.teste[a,2] <- ks(AUX)
    resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
    #print(a)
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias, medias),
                poder=c(poder[,1], poder[,2], poder[,3]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks1', length(medias)), rep('pv.ks2', length(medias))))
saveRDS(dat, file='normal2pop')


plot(medias, poder[,3], type='l')
lines(medias, poder[,2], col = 'blue')
lines(medias, poder[,1], col= 'red')

##gerando dados de duas pop exp e tracando a funcao poder 
simulacao<-500
medias<-seq(1,5,by=0.2)
poder <- matrix(NA, length(medias),3)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=3,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-rexp(100, rate=2)
    data2<-rexp(1000, rate= medias[b])
    AUX <- empirica(data1,data2)
    resul.teste[a,1] <- calculo_ev_mod(AUX,phi=phi)
    resul.teste[a,2] <- ks(AUX)
    resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
    #print(a)
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias, medias),
                poder=c(poder[,1], poder[,2], poder[,3]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks1', length(medias)), rep('pv.ks2', length(medias))))
saveRDS(dat, file='exponencial2pop')


plot(medias, poder[,3], type='l')
lines(medias, poder[,2], col = 'blue')
lines(medias, poder[,1], col= 'red')

## trigonometrica

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

simulacao<-500
medias<-seq(1,10,by=1)
poder <- matrix(NA, length(medias),3)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=3,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-gerar(100, 3)
    data2<-gerar(100, medias[b])
    AUX <- empirica(data1,data2)
    resul.teste[a,1] <- calculo_ev_mod(AUX,phi=phi)
    resul.teste[a,2] <- ks(AUX)
    resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
    #print(a)
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias, medias),
                poder=c(poder[,1], poder[,2], poder[,3]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks1', length(medias)), rep('pv.ks2', length(medias))))
saveRDS(dat, file='trigo2popigual')

plot(medias, poder[,3], type='l')
lines(medias, poder[,2], col = 'blue')
lines(medias, poder[,1], col= 'red')

## uniforme descontinua

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

simulacao<-100
medias<-seq(1,10,by=1)
poder <- matrix(NA, length(medias),3)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=3,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-runif.haar(100, 3)
    data2<-runif.haar(1000, 1)
    AUX <- empirica(data1,data2)
    resul.teste[a,1] <- calculo_ev_mod(AUX,phi=phi)
    resul.teste[a,2] <- ks(AUX)
    resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
    print(a)
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias, medias),
                poder=c(poder[,1], poder[,2], poder[,3]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks1', length(medias)), rep('pv.ks2', length(medias))))
saveRDS(dat, file='haar2pop')

simulacao<-500
medias<-seq(1,10,by=1)
poder <- matrix(NA, length(medias),3)  
for(b in 1:length(medias)){
  resul.teste<-matrix(NA,ncol=3,nrow = simulacao)
  for(a in 1:simulacao){
    data1<-runif.haar(100, 3)
    data2<-runif.haar(100, medias[b])
    AUX <- empirica(data1,data2)
    resul.teste[a,1] <- calculo_ev_mod(AUX,phi=phi)
    resul.teste[a,2] <- ks(AUX)
    resul.teste[a,3]<-ks.test(data1,data2)$'p.value'
    #print(a)
  }
  poder[b,]<-colMeans(resul.teste<0.05)
  print(b)  
}

dat<-data.frame(medias=c(medias, medias, medias),
                poder=c(poder[,1], poder[,2], poder[,3]), 
                teste=c(rep('ev', length(medias)), rep('pv.ks1', length(medias)), rep('pv.ks2', length(medias))))
saveRDS(dat, file='haar2popigual')

plot(medias, poder[,3], type='l')
lines(medias, poder[,2], col = 'blue')
lines(medias, poder[,1], col= 'red')


##graficos

dat<-readRDS(file='normal2pop')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line() +
  scale_colour_discrete(labels=c("FBST", "KS uniforme", "KS entre 2"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "μ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dat<-readRDS(file='exponencial2pop')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line() +
  scale_colour_discrete(labels=c("FBST", "KS uniforme", "KS entre 2"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "μ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dat<-readRDS(file='trigo2pop')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line() +
  scale_colour_discrete(labels=c("FBST", "KS uniforme", "KS entre 2"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "μ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))

dat<-readRDS(file='trigo2popigual')
ggplot(data = dat, aes(x=medias, y=poder, group=teste , colour=teste)) +
  geom_line() +
  scale_colour_discrete(labels=c("FBST", "KS uniforme", "KS entre 2"))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col='green') +
  ggtitle("") +
  labs(y = "Poder do teste", x = "μ")+
  theme_bw() +
  theme(#axis.line = element_line(colour = "black"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(),
    #legend.position="top",
    legend.title=element_blank(),
    axis.text=element_text(size=15),
    legend.text = element_text(size=15),
    axis.title.x = element_text(size = 20, angle = 0),
    axis.title.y = element_text(size = 20, angle = 90))




