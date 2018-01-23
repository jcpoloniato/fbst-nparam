#Funcoes uteis

#funcao de x para calcular a i-esimo funcao da base ortonormal de cossenos
phi<-function(i,x){
  return(sqrt(2)*cos(i*pi*x))
}

phi_barra<-function(dados,phi,r){
  phi_barra<-rep(NA,r)
  for (i in 1:r){
    phi_barra[i]<-mean(phi(i,dados))
  } 
  return(phi_barra)
}

#funcao que recebe como parametros os dados ja tranformados (phi barra)
#e calcula a estatistica de teste do FBST
evalue<-function(x,n,tau){
  eval<-sum((n*x/sqrt(n+tau))^2)
  #eval <- n*x/sqrt(n+tau) %>% power(2) %>% sum()
  eval
}

# funcao que aplica F nos dados
# essa funcao recebe como argumentos os dados e uma distribuicao acumulada F
# se a F for a verdadeira distribuicao que gerou os dados eh esperado que os 
# dados modificados tenham districuicao uniforme
inversa<-function(dados, pf) {
  pf(dados)
}