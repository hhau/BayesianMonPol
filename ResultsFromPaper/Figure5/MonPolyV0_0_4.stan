functions {
  vector convolve(vector a, vector b){ 
    int na;
    int nb;  
    int nc;
    vector[rows(b)] rev_b;
    vector[rows(a)+rows(b)-1] c;

    na = rows(a);
    nb = rows(b);
    nc = na+nb-1;
    for(i in 1:nb){
      rev_b[i] = b[nb-i+1];
    }
    for(j in 1:nc){
      int istart;
      int istop;
      istart = max(0, j-nb) + 1;
      istop = min(na, j);
      c[j] = sum(a[istart:istop] .* rev_b[(nb-j+istart):(nb-j+istop)]);
    }
    return c;
  } 

  vector horner(vector x, vector beta){ 
    int nb;
    vector[rows(x)] res;

    nb = rows(beta);
    res = rep_vector(beta[nb], rows(x));
    for(i in 1:(nb-1)){
      res = res .* x + beta[nb-i];
    }
    return res;
  } 
}

data { 
  # number of data points
  int <lower=1> N;

  # degree of polynomial we want to fit
  int <lower=1> q;

  # degree of sub polynomials, q = 2K or q = 2K + 1
  int <lower=0> K;

  # mode of operation, derived in the R code
  # 1 == whole real line
  # 2 == (-inf, b]
  # 3 == [a, Inf)
  # 4 == [a, b]
  int <lower = 1, upper = 4> operation_mode;

  # y values
  vector[N] y;

  # x values
  vector[N] x;

  # lower bound (could be -Infty / negative_infinity() )
  real a;

  # upper bound (could be Infty / positive_infinity() )
  real <lower=a> b;

  # alpha (polynomial is increasing or decreasing)
  int <lower=-1, upper=1> alpha;

  # number of new data points at which to predict
  int <lower=1> Nnew;
  vector[Nnew] xnew;
} 

transformed data { 
  int p1_length;
  int p2_length;
  int gamma_1_length;
  int gamma_2_length;
  int gamma_1_temp_length;
  int gamma_2_temp_length;
  vector[2] lower_bound_vector;
  vector[2] upper_bound_vector;

  gamma_1_temp_length = 0;
  gamma_2_temp_length = 0;
  if ((operation_mode == 1) && (q % 2 != 0)) {
    p1_length = K + 1;
    p2_length = K + 1;
  } else if (operation_mode == 2 || operation_mode == 3) {
      if (q % 2 == 0) {
        p1_length = K;
        p2_length = K;
      } else {
        p1_length = K + 1;
        p2_length = K;
      }
      gamma_2_temp_length = 1;
  } else {
      if (q % 2 == 0) {
        p1_length = K;
        p2_length = K;
        gamma_1_temp_length = 1;
        gamma_2_temp_length = 1;
      } else {
        p1_length = K + 1;
        p2_length = K;
        gamma_2_temp_length = 2;
      }
  }
  gamma_1_length = (2 * p1_length) - 1;
  gamma_2_length = (2 * p2_length) - 1;
  if( gamma_1_temp_length != 0) 
    gamma_1_temp_length =  gamma_1_temp_length + gamma_1_length;
  if( gamma_2_temp_length != 0) 
    gamma_2_temp_length =  gamma_2_temp_length + gamma_2_length;

  # This should be the vector (-a, 1) for (a - x) convolution
  lower_bound_vector[1] = -a;
  lower_bound_vector[2] = 1;

  # This should be the vector (b, -1) for (b - x) convolution
  upper_bound_vector[1] = b;
  upper_bound_vector[2] = -1;
} 

parameters { 
  vector[p1_length] beta_1;
  vector[p2_length] beta_2;

  real<lower=0> sd_y;
  real beta_zero;
} 

transformed parameters { 
  vector[gamma_1_length] gamma_1;
  vector[gamma_2_length] gamma_2;
  vector[gamma_1_temp_length] gamma_1_temp;
  vector[gamma_2_temp_length] gamma_2_temp;
  vector[q] gamma;
  vector[q + 1] beta_final;
  vector[N] mu;

  gamma_1 = rep_vector(0, gamma_1_length);
  gamma_2 = rep_vector(0, gamma_2_length);
  gamma_1_temp = rep_vector(0, gamma_1_temp_length);
  gamma_2_temp = rep_vector(0, gamma_2_temp_length);
  gamma = rep_vector(0, q);

  # self convolute betas to get gammas
  gamma_1 = convolve(beta_1, beta_1);
  gamma_2 = convolve(beta_2, beta_2);

  # convolute with lower and upper bounds.
  # combine together to get gamma
  # do this at the same time to avoid extra control flow steps
  if (operation_mode == 1) {
    gamma = gamma_1 + gamma_2;
  } else if (operation_mode == 2) {
      gamma_2_temp = convolve(upper_bound_vector, gamma_2);
      if(q % 2 == 0) {
        gamma[1:(gamma_1_length)] = gamma_1;
        gamma = gamma + gamma_2_temp;
      } else {
        gamma[1:(gamma_2_temp_length)] = gamma_2_temp;
        gamma = gamma + gamma_1;
      }
  } else if (operation_mode == 3) {
      gamma_2_temp = convolve(lower_bound_vector, gamma_2);
      if(q % 2 == 0) {
        gamma[1:(gamma_1_length)] = gamma_1;
        gamma = gamma + gamma_2_temp;
      } else {
        gamma[1:(gamma_2_temp_length)] = gamma_2_temp;
        gamma = gamma + gamma_1;
      }
  } else {
      if (q % 2 == 0) {
        gamma_1_temp = convolve(upper_bound_vector, gamma_1);
        gamma_2_temp = convolve(lower_bound_vector, gamma_2);
        gamma = gamma_1_temp + gamma_2_temp;
      } else {
        gamma_2_temp = convolve(upper_bound_vector, convolve(lower_bound_vector, gamma_2));
        gamma = gamma_1 + gamma_2_temp;
      }
  }

  beta_final[1] = beta_zero;
  for(i in 1:q) {
    beta_final[i+1] = alpha * gamma[i] / i;
  }

  # use horner() to evaluate and get mu
  mu = horner(x, beta_final);
} 

model { 
 y ~ normal(mu, sd_y);
} 

generated quantities { 
  vector[Nnew] mupred;
  vector[Nnew] ypred;

  mupred = horner(xnew, beta_final);
  ypred = multi_normal_cholesky_rng(mupred, diag_matrix(rep_vector(sd_y, Nnew)));
} 
