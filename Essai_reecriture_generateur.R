# Essai de réécriture


rPDD <- function(n=1,S0=100,Sa=120,delta_t = 1/365,alpha=0.2,r=0.015,sigma=0.45,raffinement=FALSE, Temps = 1){
  
  # Probabilité de franchissement du seuil dans le cas de pont brownien : 
  p_k <- function(x,y){ return(exp(-2*log(x/Sa)*log(y/Sa)/(sigma^2*delta_t))) }
  
  # Calcul d'une PDD : 
  rPDD_unitaire <- function(naif,S0,Sa){
    
    # Chacun des lambdas dépends de la valeur de l'actif en début de periode et 
    # de la valeur du lambda précédent. Le rete n'est que simulation de l'actif. 
    
    lambda <- function(S0,lambdaPrecedent = 0){
      
      # Commençons par calculer la trajectoire de S : 
      # S suit une dynamique de black-sholes, donc la solution explicite de BS nous donne :
      alea <- rnorm(1/delta_t)
      S <- S0 * cumprod(1+r * delta_t + sigma * sqrt(delta_t) * alea)
      S <- S*(S > 0) # Histoire que le sous-jacent reste bien positif.
      
      # Gestion du raffinement par pont brownien : 
      p = S^0
      if(!naif){
        p <- c(1,(p_k(c(NA,S),c(S,NA))[2:(length(S)-1)]))
        p <- p > runif((1/delta_t)-1) # donc p = true quand la barrière désactivante n'est pas franchie. 
        p[1:(floor(1/delta_t)+1)] <- TRUE # Les 6 premiers mois sont toujours valides.
      }
      
      # Maintenant qu'on a la trajectoire de S, calculons la PDD_0 correspondante : 
      # Approximation du sup sur 6mois : 
      sup_6_mois <- max(S[(floor(length(S)/2)+1):length(S)])
      # Calcul final de la PDD : 
      lambda <-prod(p)*max(Sa-lambdaPrecedent-S[length(S)],0) * (S[length(S)] < (1-alpha) * Sa) * (sup_6_mois < Sa) 
      return(c(S[length(S)],lambda))
    }
    
    lambda_rez <- list()
    lambda_rez[[1]] <- lambda(S0,0)
    if(Temps > 1){
      for (i in 2:Temps){
        lambda_rez[[i]] <- lambda(lambda_rez[[i-1]][1],lambda_rez[[i-1]][2])
        # La formule est juste la même : 
        # On relance avec un nouveau point de départ S0 correspondant a la derière valeur obtenue au cout d'avent, 
        # et le lambda précédent qui est 0 au début.
      }
    }
    lambda_rez %>% 
      lapply(function(x) return(x[[2]])) %>%
      unlist %>%
      # tail(1) %>% 
      multiply_by(exp(-r * (1:Temps))) %>%
      tail(1) %>%
      return
  }
  
  # On retourne le résultat :
  return(rep(NA,n) %>% 
           sapply(.,function(x){
             return(rPDD_unitaire(naif=!raffinement,S0=S0,Sa=Sa))
           }
           ))
}

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