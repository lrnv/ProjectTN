######### Premier script R du projet de technique numérique. 
library(magrittr)
set.seed(100) # On fixe la seed pour la reproductibilité

### On source les générateurs (utiliser "setwd()" si besoin pour se mettre dans le bon dossier.)
source("Generateurs.R")


#1 Une seul actif et une seule periode 

#1.1 Estimation de la PDD
#1.1.1. Methode naive.

# Petite fonction d'affichage : 
affichage <- . %>% 
  {data.frame(Moyenne = mean(.),ProbaPositif = mean(. > 0))} #%>%
  #{print(.)}

M = 100
list(
    # Test : ( devrais retourner 23/24/25)
    NaifT1 = rPDD(M) %>% affichage,
    
    
    
    #1.1.2 Méthode par pont brownien : 
    RafiT1 = rPDD(M,raffinement=TRUE) %>% affichage,
    
    #1.1.
    RafiT2 = rPDD(M,raffinement = TRUE,Temps = 2) %>% affichage,
    RafiT5 = rPDD(M,raffinement = TRUE,Temps = 5) %>% affichage,
    RafiT10 = rPDD(M,raffinement = TRUE,Temps = 10) %>% affichage,
    RafiT15 = rPDD(M,raffinement = TRUE,Temps = 15) %>% affichage,
    RafiT100 = rPDD(M,raffinement = TRUE,Temps = 100) %>% affichage
) %>% 
  do.call(rbind,.) %>% 
  data.frame



# Par exemple, plutot que de faire du MC, on peut juste sortir
# un plot de la densitée empirique des data générées. 
vector <- rPDD(n=1000,
     S0=100,
     Sa=120,
     delta_t = 1/2,
     alpha=0.2,
     r=0.015,
     sigma=0.45,
     raffinement=FALSE,
     Temps = 10
     ) %T>% {print(sum(.>120) + sum(. < 0))}
# Pour cetains paramètres, ça dépasse 120. Ce qui est bizare. 
# Donc la formule va pas. 
# Meme pour Temps = 1 ça épasse des fois; SheiBe


rPDD(1000) %>% density %>% plot












#1.2 Calcul de la sensibilite par differences finies

rGreek(1000,"Vega") %>% hist
rGreek(1000,"Vega") %>% hist
rGreek(1000,"Rho") %>% hist
rGreek(1000,"Delta") %>% hist
rGreek(10000,"Vega") %>% density %>% plot

















########### Code diephan ; 

PDDRaffinee<-function(r,sigma,T,alpha,S0 = 100,Sa = 120){
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
  #return(PDDtotraf)
  return(PDDanraf)
}
PDDRaffinee(0.015,0.45,1,0.20,100,120) %>% density %>% plot
PDDRaffinee(0.015,0.45,1,0.20,100,120) %>% mean





