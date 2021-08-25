### ANOVA MODEL ###
library(nntcalc)

set.seed(1)

n = 200

N = 400

tr = c( rep("11", n/4), rep("12", n/4), rep("13", n/4) )
cr = c( rep("00", n/4) )

gr = c(tr, cr)

mu_c = 0
mu_t = c(0.3, 0.6, 0.9)

sig  = 1

y = numeric()

# fun_g   = function(x){ ifelse( x > 0, 1/x, Inf ) }

tau     = 0

list_1.2 = list(); list_1.3 = list(); list_1.4 = list()

for(j in 1:N){

  ### treatment arm ###
  
  for( i in 1:(n/4) ){
    y[i] = mu_t[1] + rnorm(1, 0, sig) }
  
  for( i in ( (n/4) + 1 ):(2 * n/4) ){
    y[i] = mu_t[2] + rnorm(1, 0, sig) }
  
  for( i in ( (2 * n/4) + 1 ):(3 * n/4) ){
    y[i] = mu_t[3] + rnorm(1, 0, sig) }
  
  ### control arm ###
  
  for( i in ( ( 3 * n/4) + 1 ):(4 * n/4) ){
    y[i] = mu_c + rnorm(1, 0, sig) }
  
  dat1 = as.data.frame( as.matrix( cbind( as.numeric(round(y, 2)), as.factor(gr)  ) ) )
  colnames( dat1 ) = c( "y", "gr" )
  
  list_1.2[[j]] = nnt_x( model    = "anova",
                         response = dat1$y,
                         x        = dat1$gr,
                         adj      = x01,
                         cutoff   = tau,
                         base     = 1,
                         decrease = FALSE,
                         data     = dat1 )
  
  print(j)
  
}

all_est = data.frame( matrix( NA, ncol = (27 + 9*2), nrow = N ))

for( i in 1:N){
  
  all_est[i,] = c( list_1.2[[i]][1,],
                   list_1.2[[i]][2,],
                   list_1.2[[i]][3,],
                   list_1.2[[i]][4,],
                   list_1.2[[i]][5,] )
  
}


colnames( all_est ) = c( "NNTL",    "NNTL_CI_TR_L",    "NNTL_CI_TR_U",    "NNTL_CI_DL_L",   "NNTL_CI_DL_U",   "NNTL_CI_NBS_L",   "NNTL_CI_NBS_U", "NNTL_CI_PBS_L", "NNTL_CI_PBS_U",
                         "NNTML",   "NNTML_CI_TR_L",   "NNTML_CI_TR_U",   "NNTML_CI_DL_L",  "NNTML_CI_DL_U",  "NNTML_CI_NBS_L",  "NNTML_CI_NBS_U", "NNTML_CI_PBS_L", "NNTML_CI_PBS_U",
                         "NNT1.2",  "NNT1.2_CI_TR_L",  "NNT1.2_CI_TR_U",  "NNT1.2_CI_DL_L", "NNT1.2_CI_DL_U", "NNT1.2_CI_NBS_L", "NNT1.2_CI_NBS_U", "NNT1.2_CI_PBS_L", "NNT1.2_CI_PBS_U",
                         "NNT1.3",  "NNT1.3_CI_TR_L",  "NNT1.3_CI_TR_U",  "NNT1.3_CI_DL_L", "NNT1.3_CI_DL_U", "NNT1.3_CI_NBS_L", "NNT1.3_CI_NBS_U", "NNT1.3_CI_PBS_L", "NNT1.3_CI_PBS_U",
                         "NNT1.4",  "NNT1.4_CI_TR_L",  "NNT1.4_CI_TR_U",  "NNT1.4_CI_DL_L", "NNT1.4_CI_DL_U", "NNT1.4_CI_NBS_L", "NNT1.4_CI_NBS_U", "NNT1.4_CI_PBS_L", "NNT1.4_CI_PBS_U" )


write.csv( all_est, "anova_nnt200_wide_setseed1.csv", row.names = F )


library(tidyr)

long_nnt <- all_est %>% gather()

write.csv(long_nnt, "anova_nnt200_long_setseed1.csv", row.names = F)

