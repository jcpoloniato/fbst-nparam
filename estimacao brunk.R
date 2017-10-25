phi<-function(i,x){
  return(sqrt(2)*cos(i*pi*x))
}

phi2<-function(i,x){
  return(ifelse(i%%2,sqrt(2)*sin((i+1)*pi*x),sqrt(2)*cos((i*pi*x))))
}


#dados = amostra da funcao f
#phi = base ortonormal
#chute = uma lista com os membros:
# f0 = chute de densidade que se ajusta a f real
# pf0 = acumulada da funcao chute
#r = quantidade de funções ortogonais a serem usadas
#d = peso da base de ordem (d^r)
#enp retorna uma lista com os membros: 
# f_chapeu: a densidade estimada 
# beta_chapeu: os coeficientes estimados da base 
enp<-function(dados,phi,chute,r,d){
  
  pf0<-chute$pf0
  f0<-chute$f0
  
  n<-length(dados)
  z<-pf0(dados)
  beta_chapeu<-rep(NA,r)
  phi_barra<-rep(NA,r)
  
  for (i in 1:r){
    phi_barra[i]<-mean(phi(i,z))
    beta_chapeu[i]<-(n/(n-1+d^i))*phi_barra[i]
  }
  
  p_chapeu<-function(x){
    #out<-rep(1,length(x))
    out<-1
    for (i in 1:r){
      out<-out+beta_chapeu[i]*phi(i,x)
    }
    return(out)
  }
  
  f_chapeu<-function(x){
    out<-p_chapeu(pf0(x))*f0(x)
    return(out)
  }
  
  plot(density(dados))
  y<-seq(min(dados),max(dados),by=0.01)
  lines(y,f_chapeu(y),col='blue')
 
  return(list(f_chapeu=f_chapeu,
              beta_chapeu=beta_chapeu))
}

#
chute.norm<-function(dados){
  m0<-mean(dados)
  s0<-sd(dados)
  return(list(f0=function(x) dnorm(x,mean=m0,sd=s0),
              pf0=function(x) pnorm(x,mean=m0,sd=s0)
              ))
}
chute.gama<-function(dados){
  v0<-sd(dados)^2
  m0<-mean(dados)
  s0<-v0/m0
  a0<-m0/s0
  return(list(pf0=pf0<-function(x) pgamma(x,shape=a0,scale=s0),
              f0=f0<-function(x) dgamma(x,shape=a0,scale=s0)
  ))
}
chute.exp<-function(dados){
  rate<-1/mean(dados)
  return(list(  pf0=pf0<-function(x) pexp(x,rate=rate),
                f0=f0<-function(x) dexp(x,rate=rate)))
}
chute.beta<-function(dados){
v0<-sd(dados)^2
m0<-mean(dados)
a0<-(-m0^3+m0^2-m0*v0)/v0
b0<-(m0^3-2*m0^2+m0*v0+m0-v0)/v0
return(list(pf0=pf0<-function(x) pbeta(x,shape1=a0,shape2=b0),
            f0=f0<-function(x) dbeta(x,shape1=a0,shape2=b0)))
}
#EM DESENVOLVIMENTO
chute.mv<-function(dados,theta0,f,pf){
  log.lik<-function(theta) dados %>% f(theta) %>% log %>% sum 
  #sum(log(f(dados,theta=theta)))
  theta_chapeu<-optim(theta0,log.lik)
  return(list(pf0=function(x) pf(x,theta=theta_chapeu),
              f0=function(x) f(x,theta=theta_chapeu)))
}

dados<-rnorm(100,mean=5,sd=10)
est<-enp(dados,phi,chute.norm(dados),50,5)

dados<-rnorm(100000,mean=5,sd=10)
pnorm2<-function(x,theta) pnorm(x,mean=theta[1],sd=theta[2])
dnorm2<-function(x,theta) dnorm(x,mean=theta[1],sd=theta[2])
est<-enp(dados,phi,chute.mv(dados,c(0,1),dnorm2,pnorm2),50,5)

dados<-rnorm(100000,mean=0,sd=10)
est<-enp(dados,phi,chute.norm(dados),50,5)

dados<-rnorm(100000,mean=0,sd=10)
chute<-list(pf0=function(x) pnorm(x,mean=5,sd=1),
            f0=function(x) dnorm(x,mean=5,sd=1))
est<-enp(dados,phi,chute,5,5)

dados<-rexp(100000,5)
est<-enp(dados,phi,chute.norm(dados),50,5)

dados<-rexp(100000,5)
est<-enp(dados,phi,chute.exp(dados),1,5)

dados<-rnorm(100000,mean=2)
est<-enp(dados,phi,chute.exp(dados),50,5)

dados<-rgamma(100000,5)
f_chapeu<-enp(dados, phi,chute.gama(dados),50,5)

dados<-rgamma(10000,shape=5,scale=1/2)
f_chapeu<-enp(dados,phi,chute.gama(dados),50,5)