functions {
  # there is probably going to need to be a convolution function here
  vector convolve (vector a, vector b) {
    # length of the inital beta vector is needed to set length for convoluted beta
    int na;
    int nb;
    int nc;
    vector[(num_elements(a) * num_elements(b)) - 1] c;
    vector[num_elements(b)] rev_b;
    na <- num_elements(a);
    nb <- num_elements(b);
    nc <- (na * nb)  - 1;

    # just implement the sum as it is written in eqn 8 (in the new paper?)
    # am I missing something here?
    # not sure arrays / vectors are zero indexed like Berwin thinks? (documentation suggests no)
    for (i in 1:nb) {
      rev_b[i] <- b[nb - i + 1];
    }

    for (j in 1:nc) {
      int istart;
      int istop;
      int istart_2;
      int istop_2;
      istart <- max(0, j - nb) + 1;
      istop <- min(na, j);
      istart_2 <- nb - j + istart;
      istop_2 <- nb - j + istop;
      c[j] <- sum(a[istart:istop] .* rev_b[istart_2:istop_2]);
    }
    return(c);
  }

  vector horner(vector x_values, vector beta) {
    int n_beta;
    int n_x;
    vector[num_elements(x_values)] res;
    n_x <- num_elements(x_values);
    n_beta <- num_elements(beta);

    for (i in 1:num_elements(x_values)) {
      res[i] <- beta[n_beta];
    }


    #      for(q in 1:(n_beta - 1)) {
    #        res <- res .* x_values + beta[n_beta - q];      Vectorise ??
    #      }

    for (i in 1:n_x) {
      res[i] <- beta[n_beta];
      for (j in 1:(n_beta - 1)) {
        res[i] <- res[i] * x_values[i] + beta[n_beta - j];
      }
    }
    return(res);
  }
}

data {
  # number of data points
  int <lower=1> n_data;

  # degree of polynomial we want to fit
  int <lower=1> q;

  # degree of sub poly nomials and thangs
  int <lower=1> k;

  # sub poly order things
  int <lower=1> p1_length;
  int <lower=1> p2_length;
  int <lower=1> gamma_1_length;
  int <lower=1> gamma_2_length;
  int <lower=0> gamma_1_temp_length;
  int <lower=0> gamma_2_temp_length;
  int <lower=1> gamma_lenght;


  # mode of operation, derived in the r code
  # 1 = whole real line
  # 2 == (-inf, b]
  # 3 == [a, Inf)
  # 4 == [a, b]

  int <lower = 1, upper = 4> operation_mode;

  # y values
  vector[n_data] y_values;

  # x values
  vector[n_data] x_values;

  # lower bound (could be -Infty / negative_infinity() )
  vector[2] lower_bound_vector;
  # upper bound (could be Infty / positive_infinity() )
  vector[2] upper_bound_vector;

  # alpha (increasing or decreasing, want to do with (1/2, 1/2) priors on -1, 1 ?)
  int <lower=-1, upper=1> alpha;

  # number of new data points to predict with
  #int <lower=0> n_new_data_points;

  # probably also need to read in x.new as well? might be able to do in this once.
  #vector[n_new_data_points] x_new_values;

}

transformed data {

}


parameters {
  # underlying vectors, need to determine which more they are in though.
  vector[p1_length] beta_1;
  vector[p2_length] beta_2;

  real<lower=0> sd_y;
  real beta_zero;
}

transformed parameters {
  # declare gamma's, check mode and length for dimension
  vector[gamma_1_length] gamma_1;
  vector[gamma_2_length] gamma_2;
  vector[gamma_1_temp_length] gamma_1_temp;
  vector[gamma_2_temp_length] gamma_2_temp;
  vector[gamma_lenght] gamma
  vector[num_elements(gamma) + 1] beta_final;
  vector[n_data] mu;

  # self convolute beta to get gammas
  gamma_1 <- convolve(beta_1, beta_1);
  gamma_2 <- convolve(beta_2, beta_2);

  # convolute with lower and upper bounds.
  # combine together to get big gamma
  # do this at the same time to avoid extra control flow steps
  if (operation_mode == 1) {
    gamma <- gamma_1 + gamma_2;

  } else if (operation_mode == 2) {
    gamma_2_temp <- convolve(upper_bound_vector, gamma_2);

    if(q % 2 == 0) {
      gamma[1:(gamma_2_length + 1)] <- 0;
      gamma[1:(gamma_2_length)] <- gamma_1;
      gamma <- gamma + gamma_2_temp;

    } else {
      gamma[1:(gamma_2_length + 1)] <- 0;
      gamma[1:(gamma_2_length)] <- gamma_2_temp;
      gamma <- gamma + gamma_1;
    }

  } else if (operation_mode == 3) {
    gamma_2_temp <- convolve(lower_bound_vector, gamma_2);

    if(q % 2 == 0) {
      gamma[1:(gamma_2_length + 1)] <- 0;
      gamma[1:(gamma_2_length)] <- gamma_1;
      gamma <- gamma + gamma_2_temp;

    } else {
      gamma[1:(gamma_2_length + 1)] <- 0;
      gamma[1:(gamma_2_length)] <- gamma_2_temp;
      gamma <- gamma + gamma_1;
    }

  } else {
    if (q % 2 == 0) {
      gamma_1_temp <- convolve(upper_bound_vector, gamma_1);
      gamma_2_temp <- convolve(lower_bound_vector, gamma_2);
      gamma <- gamma_1_temp + gamma_2_temp;

    } else {
      gamma_2_temp <- convolve(upper_bound_vector, convolve(lower_bound_vector, gamma_2));
      gamma <- gamma_1 + gamma_2_temp;

    }

  }

  beta_final[1] <- beta_zero;
  for(i in 2:(num_elements(gamma) + 1)) {
    beta_final[i] <- alpha * gamma[i - 1] / (i - 1);
  }
  # use horner() to evaluate and get mu

  mu <- horner(x_values, beta_final);
  # doneski ? Feel like i'm forgetting something
}

model {
  y_values ~ normal(mu, sd_y);
}

generated quantities {

}