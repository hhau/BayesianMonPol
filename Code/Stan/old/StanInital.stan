functions {
  # there is probably going to need to be a convolution function here
  vector self_convolute (vector beta_vec) {
    # length of the inital beta vector is needed to set length for convoluted beta
    int q;
    q <- size(beta_vec);
    vector[(2*q) + 1] beta_sq;

    # just implement the sum as it is written in eqn 8 (in the new paper?)
    # am I missing something here?
    # not sure arrays / vectors are zero indexed like Berwin thinks? (documentation suggests no)
    for(j in 0:((2*q) + 1)) {
      int loop_lower;
      int loop_upper;
      loop_lower <- max(0, j - q) + 1;
      loop_upper <- min(q, j);
      for(i in loop_lower:loop_upper) {
        beta_sq[j] <- beta_vec[i] * beta_vec[j - i + 1];
      }
    }
   return(beta_sq);
  }

  vector horner(vector x_values, vector beta) {
      int n_x;
      int n_beta;

      n_x <- size(x_values);
      n_beta <- size(beta);

      vector[n_x] res;


      res <- beta[n_beta];

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
    int <lower=0> n_data_points;

  # degree of polynomial we want to fit
    int <lower=1> poly_order;

  # y values
    vector[n_data_points] y_values;

  # x values
    vector[n_data_points] x_values;

  # lower bound (could be -Infty / negative_infinity() )
    real lower_bound;

  # upper bound (could be Infty / positive_infinity() )
    real <lower=lower_bound> upper_bound;

  # alpha (increasing or decreasing, want to do with (1/2, 1/2) priors on -1, 1 ?)
    int <lower=-1, upper=1> alpha;

  # number of new data points to predict with
    #int <lower=0> n_new_data_points;

  # probably also need to read in x.new as well? might be able to do in this once.
    #vector[n_new_data_points] x_new_values;

}

transformed data {
  #(p_check degree)
  int <lower=0> poly_check_order;
  poly_check_order <- poly_order - 1;

  # this will depend on the upper and lower bounds as well
  int <lower=1> sub_poly_order_1; # this is effectivly k
  int <lower=1> sub_poly_order_2; # this is if we need k and k-1

  if (lower_bound == negative_infinity() && upper_bound == positive_infinity()) { # (-infty, Infty)
    if (poly_check_order % 2 == 0){
      # this is the okay example
      sub_poly_order_1 <- poly_check_order / 2;
      sub_poly_order_2 <- poly_check_order / 2;

    } else {
      reject("Invalid polynomial order for specified interval");

    }
  } else if (lower_bound == negative_infinity() && upper_bound < positive_infinity()) { # (-Infty, b)
    if (poly_check_order % 2 == 0) { # very unsure if this behaves the way I think it does?
                                      # I.e exactly the same as the (a, Infty) case.
      sub_poly_order_1 <- poly_check_order / 2;
      sub_poly_order_2 <- sub_poly_order_1 - 1;

    } else {
        sub_poly_order_1 <- (poly_check_order - 1) / 2;
        sub_poly_order_1 <- (poly_check_order - 1) / 2;
    }

  } else if (lower_bound < negative_infinity() && upper_bound == positive_infinity()) { # (a, Infty)
    if (poly_check_order % 2 == 0) { # note that we are using the order of poly_check, not q
      sub_poly_order_1 <- poly_check_order / 2;
      sub_poly_order_2 <- sub_poly_order_1 - 1;

    } else {
        sub_poly_order_1 <- (poly_check_order) / 2;
        sub_poly_order_2 <- (poly_check_order) / 2;
    }

  } else  { # (a,b)
    if (poly_check_order % 2 == 0) {
      sub_poly_order_1 <- ceil(poly_check_order / 2);
      sub_poly_order_2 <- sub_poly_order_1 - 1;
    } else {
        sub_poly_order_1 <- (poly_check_order - 1) / 2;
        sub_poly_order_2 <- (poly_check_order - 1) / 2;
    }
  }

}

parameters {
  # poly orders in here as well??

  # standard error for the data
  real <lower=0> sd_y;

  # actual betas

  # poly_1 / poly_2 betas
  vector[sub_poly_order_1] beta_p_1;
  vector[sub_poly_order_2] beta_p_2;


}

transformed parameters {
  #  post convolution parameters ??
  # these are almost certainly transformed parameters
  vector[poly_order] beta_final;

  # the dimension / order of this will need to be thought about carefully as it changes
  # in one specific case we will need to add a zero somewhere? or just be clever with the adding
  #
  int gamma_order;
  gamma_order <- 2 * max(sub_poly_order_1, sub_poly_order_2) + 1;

  # this next bit has about a zero % probablilty of working for th next little while
  vector[gamma_order] gamma_1;
  vector[gamma_order] gamma_2;
  gamma_1 <- self_convolute(beta_p_1);
  gamma_2 <- self_convolute(beta_p_2);

  vector[gamma_order] gamma_combined;
  gamma_combined <- gamma_1 + gamma_2;

  beta_final[1] <- gamma_combined[1]

  for (i in 1:(gamma_order - 1)) {
    beta_final[i + 1] <- alpha * gamma_combined[i] / i;
  }

  # this might be a matrix? Dimensions suggest it should be
  vector[] mu;

  #

  # write sd_y as diag_matrix , yes

}

model {
  # this just contains the top section in the bugs file, horner scheme stuff
    y_values ~ normal(mu, sd_y); # not sure this works dimension wise?


}

generated quantities {
  # simulated things, need x.new and things
  # probably the easiest, but not the quickest way to get intervals.
}
