######### Premier script R du projet de technique numérique. 
library(magrittr)
########### Premierre chose a faire : setup certians paramètres : 
M = 10000 # nombre de simulations de monte-carlo.
set.seed(100) # seed

#1 Une seul actif et une seule periode 

#1.1 Estimation de la PDD
#1.1.1. Methode naive.

# Fonction de simulation de PDD.
rPDD <- function(n=1,S0=100,Sa=120,delta_t = 1/365,alpha=0.2,r=0.015,sigma=0.45,raffinement=FALSE){
  # Probabilité de franchissement du seuil dans le cas de pont brownien : 
  p_k <- function(x,y){
    return(exp(-2*log(x/Sa)*log(y/Sa)/(sigma^2*delta_t)))
  }
  
  # Calcul d'une PDD : 
  rPDD_unitaire <- function(naif = !raffinement){
    # Commençons par calculer la trajectoire de S : 
    # S suit une dynamique de black-sholes, donc la solution explicite de BS nous donne :
      alea <- rnorm(1/delta_t)
      S <- S0 * cumprod(1+r * delta_t + sigma * sqrt(delta_t) * alea)
    
    
    # Gestion du raffinement par pont brownien : 
      p = S^0
      if(!naif){
        p <- c(1,(p_k(c(NA,S),c(S,NA))[2:(length(S)-1)]))
        p <- p > runif((1/delta_t)-1)
        p[1:(floor(1/delta_t)+1)] <- TRUE
      }
    
    # Maintenant qu'on a la trajectoire de S, calculons la PDD_0 correspondante : 
      # Approximation du sup sur 6mois : 
        sup_6_mois <- max(S[(floor(length(S)/2)+1):length(S)])
      # Calcul final de la PDD : 
        lambda <-prod(p)*max(Sa-S[length(S)],0) * (S[length(S)] < (1-alpha) * Sa) * (sup_6_mois < Sa) 
        return(lambda)
  }
  
  # On retourne le résultat :
  return(rep(NA,n) %>% 
              sapply(.,function(x){
                        return(rPDD_unitaire(naif=!raffinement))
                      }
  ))
}

# Petite fonction d'affichage : 
affichage <- . %>% 
  {data.frame(Moyenne = mean(.),ProbaPositif = mean(. > 0))} %>%
  {print(.)}

# Test : ( devrais retourner 23/24/25)
rPDD(M) %>% affichage



#1.1.2 Méthode par pont brownien : 
rPDD(M,raffinement=TRUE) %>% affichage



#1.2 Calcul de la sensibilit� par diff�rences finies


sensibilite <- function(S0=100,Sa=120,delta_t = 1/365,alpha=0.2,r=0.015,sigma=0.45,raffinement=FALSE, greek="Vega", delta_theta= 0.00001){
  if (greek == "Vega"){
    PDD_plus <-  mean(rPDD(n=10000,sigma = sigma + delta_theta, raffinement = TRUE))
    PDD_moins <- mean(rPDD(n=10000,sigma = sigma - delta_theta, raffinement = TRUE))
  }
  else if (greek == "Rho"){
    PDD_plus <-  mean(rPDD(n=10000,r = r + delta_theta, raffinement = TRUE))
    PDD_moins <- mean(rPDD(n=10000,r = r - delta_theta, raffinement = TRUE))
  }
  else if (greek == "Delta"){
    PDD_plus <-  mean(rPDD(n=10000,S0 = S0 + delta_theta, raffinement = TRUE))
    PDD_moins <- mean(rPDD(n=10000,S0 = S0 - delta_theta, raffinement = TRUE))
  }
  return((PDD_plus-PDD_moins)/delta_theta)
}
  
  
sensibilite()
sensibilite(greek = "Rho")
sensibilite(greek = "Delta")  
  
  
#greek doit �tre "Vaga, "Rho" ou "Delta"

  





























########### Code diephan ; 

PDDRaffinee<-function(r,sigma,T,alpha){
  J=365
  N=100
  dt=T/J
  Simu<-matrix(0,N,J*T+1)
  Simu[,1]=S0
  for(i in 1:N){
    for (k in 2:(J*T+1))
    {
      Simu[i,k]<-Simu[i,k-1]*(1+r*dt+sigma*sqrt(dt)*rnorm(1,0,1))
      if (Simu[i,k]<0) Simu[i,k]=0
    }
  }
  PDDanraf<-matrix(0,N,T)
  VecSraf<-matrix(0,N,T)
  VecSraf[,1]<-Sa
  PDDraf<-rep(0,T)
  PDDtotraf=0
  condition<-matrix(1,N,T*J+1)
  probas<-rep(0,N)
  
  for (i in 1:J*T+1){
    if(i%%J>(J/2)){
      for(j in 1:N){
        if (Simu[j,i]>Sa||Simu[j,(i+1)]>Sa) condition[j,i]=0
        else{
          
          
######### La formule est ici. 
          probas[j]=exp(-2*log((Simu[j,i]/Sa))*(log(Simu[j,(i+1)]/Sa))/(sigma^2*dt))
          if(runif(1,0,1)<probas[j]){
            condition[j,i]=0
          }
        }
      }
    }
  }
  for (i in 1:N){
    if (0%in%condition[i,1:(N+1)]==FALSE){
      PDDanraf[i,1]<-(Sa-Simu[i,(J+1)])*(Sa-Simu[i,(J+1)]>0)*((Simu[i,J+1]<(1-alpha)*Sa))*(max(Simu[i,(J/2+1):(J+1)])<Sa)
    }
    else PDDanraf[i,1]=0
    if(T>1){
      for (j in 1:(T-1)){
        VecSraf[i,(j+1)]<-VecSraf[i,j]-PDDanraf[i,j]
        if (0%in%condition[i,(j*J+1):((j+1)*J+1)]==FALSE){
          PDDanraf[i,(j+1)]<-(VecSraf[i,(j+1)]-Simu[i,(J*(j+1)+1)])*(VecSraf[i,(j+1)]-Simu[i,(J*(j+1)+1)]>0)*(Simu[i,(J*(j+1)+1)]<(1-alpha)*Sa)*(max(Simu[i,((j+0.5)*J+1):(J*(j+1)+1)])<S0)
        }
        else PDDanraf[i,(j+1)]=0
      }
    }
  }
  for (j in 1:T){
    PDDraf[j]=(sum(PDDanraf[,j])/N)
    PDDtotraf=PDDtotraf+exp(-r*j)*PDDraf[j]
  }
  return(PDDtotraf)
}
PDDRaffinee(0.015,0.45,1,0.20)






