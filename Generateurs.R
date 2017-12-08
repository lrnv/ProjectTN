##################################################################
###################### Fonction .rPDD_unitaire : générateur de PDD à aléa paramètrable.
##################################################################
.rPDD_unitaire <- function(naif=TRUE, S0=100, Sa=120, delta_t = 1/365,
                           alpha=0.2, r=0.015, sigma=0.45, Temps = 1, alea){
      
      # Le paramètres "aléa" attend une rélisation d'une variable aléatoire N(0,1) de dimention Temps/delta_t
  
      # Probabilité de franchissement du seuil dans le cas de pont brownien : 
        p_k <- function(x,y){ 
          return(exp(-2*log(x/Sa)*log(y/Sa)/(sigma^2*delta_t))) 
        }
        
        
      # Chacun des lambdas dépends de la valeur de l'actif en début de periode et 
      # de la valeur du lambda précédent. Le rete n'est que simulation de l'actif. 
        lambda <- function(S0,lambdaPrecedent = 0,periode=1){
          
          # Commençons par calculer la trajectoire de S : 
          # S suit une dynamique de black-sholes, donc la solution explicite de BS nous donne :
          PetitAlea <- alea[(1 + (periode-1)/delta_t):(periode/delta_t)]
          S <- S0 * cumprod(1+r * delta_t + sigma * sqrt(delta_t) * PetitAlea)
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
      
      # Application de la fonction lambda récursivement : 
        lambda_rez <- list()
        lambda_rez[[1]] <- lambda(S0 = S0,lambdaPrecedent = 0,periode=1)
        if(Temps > 1){
          for (i in 2:Temps){
            lambda_rez[[i]] <- lambda(S0 = lambda_rez[[i-1]][1],lambdaPrecedent = lambda_rez[[i-1]][2], periode = i)
            # La formule est juste la même : 
            # On relance avec un nouveau point de départ S0 correspondant a la derière valeur obtenue au cout d'avent, 
            # et le lambda précédent qui est 0 au début.
          }
        }
      
      # Récupération des résultat, somme et actualisation pour récupérer la PDD : 
        lambda_rez %>% 
          lapply(function(x) return(x[[2]])) %>%
          unlist %>%
          multiply_by(exp(-r * (1:Temps))) %>% # Actualisation des différents lambdas.
          tail(1) %>%
          return
      # Version sans magrittr : 
      # return(tail(unlist(lapply(lambda_rez,function(x) return(x[[2]])))*exp(-r * (1:Temps)),1))
}

##################################################################
###################### Fonction rPDD : générateur de PDD.
##################################################################
rPDD <- function(n=1, S0=100, Sa=120, delta_t = 1/365, alpha=0.2, 
                 r=0.015, sigma=0.45, raffinement=FALSE, Temps = 1){
  
  # On applique juste plusieurs fois la fonction .rPDD_unitaire, avec à chaque fois un nouvel alléa.
  # FOnction a apply plusieurs fois : 
  to_apply <- function(x){
    return(
      .rPDD_unitaire(
        naif=!raffinement,
        S0=S0,
        Sa=Sa,
        delta_t = delta_t,
        alpha=alpha,
        r=r,
        sigma=sigma,
        Temps = Temps,
        alea=rnorm(Temps/delta_t)
      )
    )
  }
  rep(NA,n) %>% sapply(.,to_apply) %>% return
  # version sans maggritr : 
  # return(sapply(1:n,to_apply))
}


##################################################################
###################### Fonction rGreek : Générateur des greques de la PDD. (Vega, Rho, Delta, Saga)
###################### Le paramètre "greek" doit être "Vega","Rho","Delta" ou "Saga"
##################################################################
rGreek <- function(n=1,greek="Vega",S0=100,Sa=120,delta_t = 1/365,alpha=0.2,r=0.015,
                   sigma=0.45,raffinement=FALSE, Temps = 1, delta_theta=10^(-5)){
  
  # Chaque greeque est une dérivée selon un certain paramètre, avec un alea fixé. 
  # on a juste besoin de passer en paramètres de .rPDD_unitaire les mêmes paramètres, 
  # avec un +/- delta_theta pour le paramètre choisis, 
  # et le même aléa dans les deux fonctions.
  
  # On a plus qu'a retourner la différence finie pour la greque voulue, le nombre de fois voulue.
  
  
  to_apply <- function(x){
    alea <- rnorm(Temps/delta_t)
    return(
        (.rPDD_unitaire(
          S0 = S0              + delta_theta*(greek=="Delta"),
          Sa = Sa              + delta_theta*(greek=="Saga"),
          r = r                + delta_theta*(greek=="Rho"),
          sigma = sigma        + delta_theta*(greek=="Vega"),
          delta_t = delta_t,
          alpha = alpha,
          Temps = Temps,
          naif = !raffinement,
          alea = alea
        )
        - .rPDD_unitaire(
          S0 = S0              - delta_theta*(greek=="Delta"),
          Sa = Sa              - delta_theta*(greek=="Saga"),
          r = r                - delta_theta*(greek=="Rho"),
          sigma = sigma        - delta_theta*(greek=="Vega"),
          delta_t = delta_t,
          alpha = alpha,
          Temps = Temps,
          naif = !raffinement,
          alea = alea
        )
      )/(2*delta_theta)
    )
  }
  rep(NA,n) %>% sapply(.,to_apply) %>% return
  # version sans maggritr : 
  # return(sapply(1:n,to_apply))

}
