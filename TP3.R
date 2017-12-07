# l'idée est de calculer une PDD. 

S.a <- 120

# dépréciation si les deux evennements suivants e réalisent : 
S.t <= 80% de S.a
max(S.u, pour u dans ]t-1/2,t]) <= S.a

la perte de valeur est ans ce cas de S.a - S.t

dymanique de lactif : 
dS.t = st (rdt + sigma dWt) avec S.0 = 100

r = 1,5%
sigma = 45%
W un brownien standard. 

On se propose destimer la PDD sur 1an (T = 1) définie comme : 
  
PDD = expérance de exp(-rT) * (S.a - S.T) indic(S.T <= 80% S.a) indic(max(S.u, pour u dans ]t-1/2,t]) <= S.a)

Méthode naive : On co,sidère que lactif suit une dynamique simple (celle ci-dessut), 
et on va utiliser montécarlo pour estimer la PDD, en se basant sur la discrétisation d'EULER de l'eds :

Shéma d'euler de la' :
  
  S_t_{K+1} = S_t_K (1 + r\Delta_t + \sigma \Delta_t)