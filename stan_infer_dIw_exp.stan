

data { 
  real <lower=0.0> dur_detect;
  int <lower=0> idx_T_UB;
  real <lower=0.0> T_UB_real;
  
  real <lower=0.0> x_0;
  
  int <lower=0> idx_dtVec_UB;
  vector <lower=0.0> [idx_dtVec_UB] dtVec;
  
  int <lower=0> T_dt_idxs_LB [idx_T_UB];
  int <lower=0> T_dt_idxs_UB [idx_T_UB];
  real <lower=0.0> dt;
  
  real <lower=0.0> popSize_wuhan;
  
  real <lower=0.0> W_BKK_day;
  real <lower=0.0> W_CNX_day;
  real <lower=0.0> W_KTM_day;
  real <lower=0.0> W_HAN_day;
  real <lower=0.0> W_CHI_day;
  real <lower=0.0> W_SEA_day;
  real <lower=0.0> W_SIN_day;
  real <lower=0.0> W_ICN_day;
  real <lower=0.0> W_HND_day;
  real <lower=0.0> W_TPE_day;
  real <lower=0.0> W_SYD_day;
  real <lower=0.0> W_global_day;
  
  int arr_time_BKK_1D [4];
  int arr_time_CNX; // Chiang Mai
  int arr_time_KTM;
  int arr_time_HAN;
  int arr_time_CHI;
  int arr_time_SEA;
  int arr_time_SIN;
  int arr_time_ICN_1D [2];
  int arr_time_HND_1D [2];
  int arr_time_TPE_1D [3];
  int arr_time_SYD_1D [2];  
}


transformed data {
  real growthRate_logLB = log( 0.01 );
  real growthRate_logUB = log( 0.5 );
  
  real n_seed_log_LB = log( 5.0 );
  real n_seed_log_UB = log( 20.0 );
  
  real outflux_rate_BKK = W_BKK_day / popSize_wuhan;
  real outflux_rate_CNX = W_CNX_day / popSize_wuhan;
  real outflux_rate_KTM = W_KTM_day / popSize_wuhan;
  real outflux_rate_HAN = W_HAN_day / popSize_wuhan;
  real outflux_rate_CHI = W_CHI_day / popSize_wuhan;
  real outflux_rate_SEA = W_SEA_day / popSize_wuhan;
  real outflux_rate_SIN = W_SIN_day / popSize_wuhan;
  real outflux_rate_ICN = W_ICN_day / popSize_wuhan;
  real outflux_rate_HND = W_HND_day / popSize_wuhan;
  real outflux_rate_TPE = W_TPE_day / popSize_wuhan;
  real outflux_rate_SYD = W_SYD_day / popSize_wuhan;
  real outflux_rate_global = W_global_day / popSize_wuhan;
}


parameters {
  real <lower=growthRate_logLB, upper=growthRate_logUB> growthRate_log;
  real <lower=n_seed_log_LB, upper=n_seed_log_UB> n_seed_log;
}

transformed parameters {
  real growthRate = exp( growthRate_log );
  real n_seed = exp( n_seed_log );
}


model {

  vector [idx_dtVec_UB] dIw_u_vec = n_seed * exp( growthRate * dtVec );
  
  vector [idx_T_UB] Iw_T_vec = rep_vector(0.0, idx_T_UB);
  int T_dt_LB = 0;
  int T_dt_UB = 0;
  
  // sampling statement
  growthRate_log ~ uniform(growthRate_logLB, growthRate_logUB);
  n_seed_log ~ uniform(n_seed_log_LB, n_seed_log_UB);
  
  // log-likelihood  
  for ( idx_T in 1:idx_T_UB ) {
    T_dt_LB = T_dt_idxs_LB[idx_T];
    T_dt_UB = T_dt_idxs_UB[idx_T];
    Iw_T_vec[idx_T] = sum( dIw_u_vec[ T_dt_LB : T_dt_UB ] ) * dt;
    }
  
  target += 
  sum( log( outflux_rate_BKK * Iw_T_vec[arr_time_BKK_1D] ) ) +
  sum( log( outflux_rate_ICN * Iw_T_vec[arr_time_ICN_1D] ) ) +
  sum( log( outflux_rate_HND * Iw_T_vec[arr_time_HND_1D] ) ) +
  sum( log( outflux_rate_TPE * Iw_T_vec[arr_time_TPE_1D] ) ) +
  sum( log( outflux_rate_SYD * Iw_T_vec[arr_time_SYD_1D] ) ) +
  log( outflux_rate_CNX * Iw_T_vec[arr_time_CNX] ) +
  log( outflux_rate_KTM * Iw_T_vec[arr_time_KTM] ) +
  log( outflux_rate_HAN * Iw_T_vec[arr_time_HAN] ) +
  log( outflux_rate_CHI * Iw_T_vec[arr_time_CHI] ) +
  log( outflux_rate_SEA * Iw_T_vec[arr_time_SEA] ) +
  log( outflux_rate_SIN * Iw_T_vec[arr_time_SIN] ) -
  outflux_rate_global * n_seed / growthRate / growthRate * (
    exp(growthRate * T_UB_real) - exp(x_0 * growthRate) + exp((x_0 - dur_detect) * growthRate) - exp(growthRate * (T_UB_real - dur_detect))
   );
    
}




