set.seed(1)

library(nntcalc)
library(survival)
library(survsim)
library(boot)

n = 400

N = 400

fun_g = function(x){ ifelse( x > 0, 1/x, Inf)}

ps_xt    = function(t, x){ 
  
  diff_abs =  abs(t -  dat_both$time) 
  
  p_t_mle = dat_both$surv_base[which(diff_abs == min( diff_abs ) ) ] ^ ( exp( coef(coxph_both)[1] + coef(coxph_both)[2] * x ) )  
  
  p_c_mle = dat_both$surv_base[which(diff_abs == min( diff_abs ) )] ^ ( exp( ( coef(coxph_both)[2] ) * x ) )  
  
  return( p_t_mle - p_c_mle )
}

### TRUE NNT(y|x) ###
alph = 0.5

ps_tx_true = function(t, x, b0, b1){ 
  
  #  base_wei  = pweibull(25, shape = 1, scale = exp(1), lower.tail = F )
  
  p_tx_true = pweibull( t,
                        scale =   exp( b0 * alph + ( b1 * alph) * x )  ,
                        shape =   1, 
                        lower.tail = F)  
  
  p_cx_true = pweibull( t,
                        scale =  exp(  ( b1 * alph) * x )  ,
                        shape =  1,
                        lower.tail = F)
  
  #   p_tx_true = base_wei ^ ( exp(- b0 - b1 * x) )
  
  #   p_cx_true = base_wei ^ ( exp(- b1 * x) )
  
  return( p_tx_true - p_cx_true )
}

nnt_true = function(t, x, b0, b1) { fun_g( ps_tx_true(t, x, b0, b1) ) }

nnt_true(8, 0.6, 1.5, 2.5); nnt_true(8, 0.8, 1.5, 2.5); nnt_true(8, 1, 1.5, 2.5); 

xx = rnorm( 10 ^ 6, 1, 1 )

nnt_untrue = function(t, x, b0, b1) { fun_g( mean( ps_tx_true(t, x, b0, b1)  ) ) }

nnt_untrue(8, xx, 1.5, 2.5)


### POINT ESTIMATORS  GEN LOOP ###
NNT0.10.5 = numeric()
NNT0.11   = numeric()
NNT0.11.5 = numeric()
NNT_UN    = numeric()
NNT_KM    = numeric()

list_2.7 = list(); list_3.1 = list(); list_3.5 = list()

for( j in 1 : N ) { 
  
  dist.ev    = "weibull"
  anc.ev     = 0.5
  beta0.ev   = 1
  dist.cens  = "weibull"
  anc.cens   = 0.01
  beta0.cens = 12
  x          = list(c("bern", 0.5), c("normal", -1.5, 1))
  beta       = list(2, 1)
  
  simple.dat =  simple.surv.sim(n,
                                365,
                                dist.ev,
                                anc.ev,
                                beta0.ev,
                                dist.cens,
                                anc.cens,
                                beta0.cens,
                                z = NULL,
                                beta,
                                x)
  
  # summary(simple.dat)

  # mean number of events per arm
  # mean( simple.dat[simple.dat$x == 1, "status"] )
  # mean( simple.dat[simple.dat$x == 0, "status"] )
  #
  # mean( simple.dat$x )
  # mean( simple.dat$stop )
  # mean( simple.dat[simple.dat$x == 1, "stop"] )
  # mean( simple.dat[simple.dat$x == 0, "stop"] )
  # median( simple.dat[simple.dat$x == 1, "stop"] )
  # median( simple.dat[simple.dat$x == 0, "stop"] )
  #
  
  # View(simple.dat)
  
    list_2.7[[j]] = nnt_survreg( response = simple.dat$stop,
                                 status   = simple.dat$status,
                                 x        = simple.dat$x.1,
                                 group    = simple.dat$x,
                                 adj      = -2.5,
                                 time.point = 8,
                                 data = simple.dat )
    
    
    list_3.1[[j]] = nnt_survreg( response = simple.dat$stop,
                                 status   = simple.dat$status,
                                 x        = simple.dat$x.1,
                                 group    = simple.dat$x,
                                 adj      = -2,
                                 time.point = 8,
                                 data = simple.dat )
    
    
    list_3.5[[j]] = nnt_survreg( response = simple.dat$stop,
                                 status   = simple.dat$status,
                                 x        = simple.dat$x.1,
                                 group    = simple.dat$x,
                                 adj      = -1.5,
                                 time.point = 8,
                                 data = simple.dat )
  
    print(j)
    
  }
  
  all_est = data.frame( matrix( NA, ncol = (27 + 9*2), nrow = N ))
  
  for( i in 1:N){
    
    all_est[i,] = c( list_2.7[[i]][1,],
                     list_2.7[[i]][2,],
                     list_2.7[[i]][3,],
                     list_3.1[[i]][3,],
                     list_3.5[[i]][3,]
                  )
    
  }
  
  
  colnames( all_est ) = c( "NNTL",    "NNTL_CI_TR_L",    "NNTL_CI_TR_U",    "NNTL_CI_DL_L",   "NNTL_CI_DL_U",   "NNTL_CI_NBS_L",   "NNTL_CI_NBS_U", "NNTL_CI_PBS_L", "NNTL_CI_PBS_U",
                           "NNTML",   "NNTML_CI_TR_L",   "NNTML_CI_TR_U",   "NNTML_CI_DL_L",  "NNTML_CI_DL_U",  "NNTML_CI_NBS_L",  "NNTML_CI_NBS_U", "NNTML_CI_PBS_L", "NNTML_CI_PBS_U",
                           "NNT2.7",  "NNT2.7_CI_TR_L",  "NNT2.7_CI_TR_U",  "NNT2.7_CI_DL_L", "NNT2.7_CI_DL_U", "NNT2.7_CI_NBS_L", "NNT2.7_CI_NBS_U", "NNT2.7_CI_PBS_L", "NNT2.7_CI_PBS_U",
                           "NNT3.1",  "NNT3.1_CI_TR_L",  "NNT3.1_CI_TR_U",  "NNT3.1_CI_DL_L", "NNT3.1_CI_DL_U", "NNT3.1_CI_NBS_L", "NNT3.1_CI_NBS_U", "NNT3.1_CI_PBS_L", "NNT3.1_CI_PBS_U",
                           "NNT3.5",  "NNT3.5_CI_TR_L",  "NNT3.5_CI_TR_U",  "NNT3.5_CI_DL_L", "NNT3.5_CI_DL_U", "NNT3.5_CI_NBS_L", "NNT3.5_CI_NBS_U", "NNT3.5_CI_PBS_L", "NNT3.5_CI_PBS_U" )
  
  
  write.csv( all_est, "coxreg_nnt800_wide_setseed1_N400_highcens.csv", row.names = F )
  
  
  library(tidyr)
  
  long_nnt <- all_est %>% gather()
  
  write.csv(long_nnt, "coxreg_nnt800_long_setseed1_N400_highcens.csv", row.names = F)