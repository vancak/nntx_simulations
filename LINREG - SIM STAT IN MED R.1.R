### LIN REG ###
library(nntcalc)

# rm(list = ls())

set.seed(1)

n = 200

N = 400

tr = rep(1, n/2)
cr = rep(0, n/2)

gr = c(tr, cr)

sig  = 1
b0   = 1
bt   = 1
bc   = 1/2

xxx   = rnorm(n, 3, 1.5)

x_seq = seq(0.1, 8, length = n) 

yt    = numeric()
yc    = numeric()

fun_g = function(x){ ifelse( x > 0, 1/x, Inf ) }

tau   = 3

p_t_true = function(x){ 1 - pnorm( ( tau - b0 - bt * x ) / sig )  }

p_c_true = function(x){ 1 - pnorm( ( tau - b0 - (bt - 1/2) * x ) / sig )  }

av_ps       = mean( p_t_true(xxx) - p_c_true(xxx) )

NNT_UN_TRUE = fun_g( av_ps )
NNT_UN_TRUE

NNT_X = function(x){ fun_g( p_t_true(x) - p_c_true(x) ) }

NNT_X(1.2); NNT_X(1.3); NNT_X(1.4)

x01 = 1.2; x02 = 1.3; x03 = 1.4

list_1.2 = list(); list_1.3 = list(); list_1.4 = list()

for(j in 1:N){

  ### treatment arm ###
for( i in 1 : (n/2) ){ 
  
  yt[i] = b0 + bt * xxx[i] + rnorm(1, 0, 1) }

### control arm ###  
for( i in 1 : (n/2) ){ 
  
  yc[i] = b0 + (bt - 1/2) * xxx[(n/2)+i] + rnorm(1, 0, 1) 
  }

  y  = c( yt, yc )
  d  = c( tr, cr )
  dat1 = data.frame( cbind( y = y, gr = d, x_var = xxx ) )
  
  list_1.2[[j]] = nnt_x( model    = "linreg",
                         response = dat1$y,
                         x        = dat1$x_var,
                         group    = dat1$gr,
                         adj      = x01,
                         cutoff   = tau,
                         decrease = FALSE,
                         data     = dat1 )
  
  list_1.3[[j]] = nnt_x( model    = "linreg",
                         response = dat1$y,
                         x        = dat1$x_var,
                         group    = dat1$gr,
                         adj      = x02,
                         cutoff   = tau,
                         decrease = FALSE,
                         data     = dat1 )
  
  list_1.4[[j]] = nnt_x( model    = "linreg",
                         response = dat1$y,
                         x        = dat1$x_var,
                         group    = dat1$gr,
                         adj      = x03,
                         cutoff   = tau,
                         decrease = FALSE,
                         data     = dat1 )
  
  print(j)
  
}

all_est = data.frame( matrix( NA, ncol = (27 + 9*2), nrow = N ))

for( i in 1:N){
  
  all_est[i,] = c( list_1.2[[i]][1,],
                   list_1.2[[i]][2,],
                   list_1.2[[i]][3,],
                   list_1.3[[i]][3,],
                   list_1.4[[i]][3,] )
  
}


colnames( all_est ) = c( "NNTL",    "NNTL_CI_TR_L",    "NNTL_CI_TR_U",    "NNTL_CI_DL_L",   "NNTL_CI_DL_U",   "NNTL_CI_NBS_L",   "NNTL_CI_NBS_U", "NNTL_CI_PBS_L", "NNTL_CI_PBS_U",
                         "NNTML",   "NNTML_CI_TR_L",   "NNTML_CI_TR_U",   "NNTML_CI_DL_L",  "NNTML_CI_DL_U",  "NNTML_CI_NBS_L",  "NNTML_CI_NBS_U", "NNTML_CI_PBS_L", "NNTML_CI_PBS_U",
                         "NNT1.2",  "NNT1.2_CI_TR_L",  "NNT1.2_CI_TR_U",  "NNT1.2_CI_DL_L", "NNT1.2_CI_DL_U", "NNT1.2_CI_NBS_L", "NNT1.2_CI_NBS_U", "NNT1.2_CI_PBS_L", "NNT1.2_CI_PBS_U",
                         "NNT1.3",  "NNT1.3_CI_TR_L",  "NNT1.3_CI_TR_U",  "NNT1.3_CI_DL_L", "NNT1.3_CI_DL_U", "NNT1.3_CI_NBS_L", "NNT1.3_CI_NBS_U", "NNT1.3_CI_PBS_L", "NNT1.3_CI_PBS_U",
                         "NNT1.4",  "NNT1.4_CI_TR_L",  "NNT1.4_CI_TR_U",  "NNT1.4_CI_DL_L", "NNT1.4_CI_DL_U", "NNT1.4_CI_NBS_L", "NNT1.4_CI_NBS_U", "NNT1.4_CI_PBS_L", "NNT1.4_CI_PBS_U" )


write.csv( all_est, "linreg_nnt200_wide_setseed1.csv", row.names = F )


library(tidyr)

long_nnt <- all_est %>% gather()

write.csv(long_nnt, "linreg_nnt200_long_setseed1.csv", row.names = F)

