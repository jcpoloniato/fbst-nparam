#Funcoes uteis

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
