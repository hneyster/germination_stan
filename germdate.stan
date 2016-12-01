data {
  int<lower=0> N; //number of obserations
  //int<lower=0> N2; // size of the new_X matrix 
  vector[N] y; //response 
  int<lower=0> K; //columns in the predictor matrix
  matrix[N,K] X; // the predictor matrix vector[N] origin;
  vector[N] temp;
  vector[N] strat;
  vector[N] origin;
}
transformed data {
vector[N] log_y;
vector[N] inter_ts;
vector[N] inter_to;
vector[N] inter_so;
vector[N] inter_tso;
matrix[N,2] m1;
matrix[N,3] m2;
matrix[N,4] inter;
matrix[N,8] X_int;

log_y=log(y);
inter_to= temp .* origin;
inter_so=strat .* origin;
inter_ts=strat .* temp;
inter_tso=origin .* strat .* temp; 
m1=append_col(inter_to, inter_so);
m2=append_col(m1, inter_ts);
inter=append_col(m2, inter_tso);
X_int = append_col(X, inter);
}
parameters {
  vector[8] beta; //the regression parameters
  real<lower=0> sigma; 
}
transformed parameters {
vector[N] mu;
mu = X_int*beta;
}
model {
  //priors 
  beta[1] ~ normal(0,10); //prior for the intercept
  for(i in 2:8)
   beta[i] ~ normal(0,5);//prior for the slopes 
  sigma ~ cauchy(0,5);
  // likelihood
  log_y ~ normal(mu, sigma);
}
