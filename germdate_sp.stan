data {
  int<lower=0> N; //number of obserations (total number of germinated seeds)
  vector[N] y; //response 
  //matrix[N,4] X; // the predictor matrix vector[N] origin;
//fixed effects
	vector[N] origin;
	vector[N] temp;
	vector[N] strat;
  
//random effects 
	//int<lower=0> nfamily; 
	//int<lower=0> nloc;  
	int<lower=0> nsp; //number of species
  
	int<lower=1, upper =nsp> sp[N]; //vector with species ID
	//int<lower=1, upper=nloc> loc[N];
	//int<lower=1, upper=nfamily> family[N];
}
transformed data {
vector[N] log_y; 		//logging the response 

vector[N] inter_ts;		//4 interaction terms 
vector[N] inter_to;
vector[N] inter_so;
vector[N] inter_tso;

log_y=log(y);
inter_to= temp .* origin;
inter_so=strat .* origin;
inter_ts=strat .* temp;
inter_tso=origin .* strat .* temp; 

//m1=append_col(inter_to, inter_so);
//m2=append_col(m1, inter_ts);
//inter=append_col(m2, inter_tso);
//X_int = append_col(X, inter);
}

parameters {
  vector[nsp] a_sp;
  //vector[nloc] a_loc;
  //vector[nfamily] a_family;
  vector[nsp] b_origin;
  vector[nsp] b_temp;
  vector[nsp] b_strat;
  vector[nsp] b_inter_to;
  vector[nsp] b_inter_so;
  vector[nsp] b_inter_ts;
  vector[nsp] b_inter_tso;
   
  real mu_b_origin;
  real mu_b_temp;
  real mu_b_strat;
  real mu_b_inter_to;
  real mu_b_inter_so;
  real mu_b_inter_ts;
  real mu_b_inter_tso;
  
  real<lower=0> sigma_b_origin;
  real<lower=0> sigma_b_temp;
  real<lower=0> sigma_b_strat;

  real<lower=0> sigma_b_inter_to;
  real<lower=0> sigma_b_inter_so;
  real<lower=0> sigma_b_inter_ts;
  real<lower=0> sigma_b_inter_tso;

real<lower=0> sigma_y;  
}

transformed parameters {
vector[N] y_hat;
		
	for(i in 1:N){
		y_hat[i] = a_sp[sp[i]] + 
		b_origin[sp[i]] * origin[i]+ 
		b_temp[sp[i]] * temp[i] + 
		b_strat[sp[i]] * strat[i] + 
	
		b_inter_to[sp[i]] * inter_to[i] +
		b_inter_so[sp[i]] * inter_so[i] +
		b_inter_ts[sp[i]] * inter_ts[i] +
		b_inter_tso[sp[i]] * inter_tso[i] 
		;
				
		}
	
}

model {
  // Priors
	mu_b_temp ~ normal(0, 4); // 
	mu_b_origin ~ normal(0, 4);
	mu_b_strat ~ normal(0, 4);


	mu_b_inter_to ~ normal(0, 4);
	mu_b_inter_ts ~ normal(0, 4);
	mu_b_inter_so ~ normal(0, 4);
	mu_b_inter_tso ~ normal(0, 4);	
	
	sigma_b_temp ~ normal(0, 4); 
	sigma_b_origin	~ normal(0, 4); 
	sigma_b_strat ~ normal(0, 4);

	sigma_b_inter_to ~ normal(0, 4);
	sigma_b_inter_ts ~ normal(0, 4);
	sigma_b_inter_so ~ normal(0, 4);
	sigma_b_inter_tso ~ normal(0, 4);	
	
	b_origin~ normal(mu_b_origin, sigma_b_origin);
	b_temp ~ normal(mu_b_temp, sigma_b_temp);
	b_strat ~ normal(mu_b_strat, sigma_b_strat);

	b_inter_to ~ normal(mu_b_inter_to, sigma_b_inter_to);
	b_inter_ts ~ normal(mu_b_inter_ts, sigma_b_inter_ts);
	b_inter_so ~ normal(mu_b_inter_so, sigma_b_inter_so);
	b_inter_tso ~ normal(mu_b_inter_tso, sigma_b_inter_tso);		
	
  
 log_y ~ normal(y_hat, sigma_y);
}
//generated quantities {

//vector[N2] y_pred;
  //y_pred = new_X * beta; //the y values predicted by the model
//}
