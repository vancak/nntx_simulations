###########################################
### LOGISTIC REGRESSION STAT IN MED R.1 ###
###########################################

# rm(list = ls())
library(nntcalc)

### DATA GENERATING PROCESS ###
set.seed(1)

n = 200

N = 400

gr = c(tr, cr)

tr = rep(1, n/2)
cr = rep(0, n/2)

b0c   = -2
b0t   = 1
bc    = 1
bt    = 1/2

p_t_true = function(x){ 1 / ( 1 + exp( - ( b0c + 2 + (bc - 0.5) * x ) ) ) }      
p_c_true = function(x){ 1 / ( 1 + exp( - ( b0c + bc * x ) ) ) }      

xxx   = rnorm(n, 2, 1)

x_seq = seq(0.1, 6, length = n) 

yt    = numeric()
yc    = numeric()

yt = c( rbinom( n/2, 1, p_t_true(xxx[1 : (n/2)]) ) )    

mean(yt)

### control arm ###  
yc = c( rbinom( n/2, 1, p_c_true(xxx[(n/2 + 1)  : n]) ) )

mean(yc)  

y   = c( yt, yc )
d   = c( tr, cr )


fun_g = function(x){ ifelse( x > 0, 1/x, Inf ) }

# NNT_K   = matrix( NA, ncol = 3,        nrow = N )
# NNT_L   = vector( length = 100 )
# NNT_UN  = vector( length = 100 )
# 
# NNTCI_LN   = matrix( NA, ncol = 19,          nrow = N )
# COV_M      = matrix( NA, ncol = 19*2,        nrow = N )

summary( glm( y ~ d + xxx + xxx * d, family = "binomial"))

av_ps      = mean( p_t_true(xxx) - p_c_true(xxx) )

NNT_UN_TRUE = fun_g( av_ps )
NNT_UN_TRUE

NNT_X = function(x){ fun_g( p_t_true(x) - p_c_true(x) ) }

NNT_X(1.5); NNT_X(2); NNT_X(2.5)

x01 = 1.5; x02 = 2; x03 = 2.5

list_2.7 = list(); list_3.1 = list(); list_3.5 = list()

for(j in 1:N){
  
  ### treatment arm ###
  yt = c( rbinom( n/2, 1, p_t_true(xxx[1 : (n/2)]) ) )    
  
  mean(yt)
  
  ### control arm ###  
  yc = c( rbinom( n/2, 1, p_c_true(xxx[(n/2 + 1)  : n]) ) )
  
  mean(yc)  
  
  y   = c( yt, yc )
  d   = c( tr, cr )
  # xtr = d * xxx
  # xco = ( 1 - d ) * xxx
  
  dat1 = data.frame( cbind( y = y, x_var = xxx, gr = d ) )
  
  list_2.7[[j]] = nnt_x( model    = "logreg",
                         response = dat1$y,
                         x        = dat1$x_var,
                         group    = dat1$gr,
                         adj      = x01,
                         data     = dat1 )

  list_3.1[[j]] = nnt_x( model    = "logreg",
                         response = dat1$y,
                         x        = dat1$x_var,
                         group    = dat1$gr,
                         adj      = x02,
                         data     = dat1 )

  list_3.5[[j]] = nnt_x( model    = "logreg",
                         response = dat1$y,
                         x        = dat1$x_var,
                         group    = dat1$gr,
                         adj      = x03,
                         data     = dat1 )
  
    print(j)
  
}

all_est = data.frame( matrix( NA, ncol = (27 + 9*2), nrow = N ))

for( i in 1:N){
  
  all_est[i,] = c( list_2.7[[i]][1,],
                   list_2.7[[i]][2,],
                   list_2.7[[i]][3,],
                   list_3.1[[i]][3,],
                   list_3.5[[i]][3,])
  
}


colnames( all_est ) = c( "NNTL",    "NNTL_CI_TR_L",    "NNTL_CI_TR_U",    "NNTL_CI_DL_L",   "NNTL_CI_DL_U",   "NNTL_CI_NBS_L",   "NNTL_CI_NBS_U", "NNTL_CI_PBS_L", "NNTL_CI_PBS_U",
                         "NNTML",   "NNTML_CI_TR_L",   "NNTML_CI_TR_U",   "NNTML_CI_DL_L",  "NNTML_CI_DL_U",  "NNTML_CI_NBS_L",  "NNTML_CI_NBS_U", "NNTML_CI_PBS_L", "NNTML_CI_PBS_U",
                         "NNT2.7",  "NNT2.7_CI_TR_L",  "NNT2.7_CI_TR_U",  "NNT2.7_CI_DL_L", "NNT2.7_CI_DL_U", "NNT2.7_CI_NBS_L", "NNT2.7_CI_NBS_U", "NNT2.7_CI_PBS_L", "NNT2.7_CI_PBS_U",
                         "NNT3.1",  "NNT3.1_CI_TR_L",  "NNT3.1_CI_TR_U",  "NNT3.1_CI_DL_L", "NNT3.1_CI_DL_U", "NNT3.1_CI_NBS_L", "NNT3.1_CI_NBS_U", "NNT3.1_CI_PBS_L", "NNT3.1_CI_PBS_U",
                         "NNT3.5",  "NNT3.5_CI_TR_L",  "NNT3.5_CI_TR_U",  "NNT3.5_CI_DL_L", "NNT3.5_CI_DL_U", "NNT3.5_CI_NBS_L", "NNT3.5_CI_NBS_U", "NNT3.5_CI_PBS_L", "NNT3.5_CI_PBS_U" )


write.csv( all_est, "logreg_nnt200_wide_setseed1.csv", row.names = F )


library(tidyr)

long_nnt <- all_est %>% gather()

write.csv(long_nnt, "logreg_nnt200_long_setseed1.csv", row.names = F)