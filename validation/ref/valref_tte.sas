* Set working directory to current file;
data _null_;
	rc=dlgcdir();
	put rc=;
run;
data _null_;
	name="%sysget(SAS_EXECFILENAME).";
	path="%sysget(SAS_EXECFILEPATH).";
	path_len = length(path);
	name_len = length(name);
	pos = find(path, name, -path_len);
	   path1 = kupdate(path, pos ,name_len+1);
	rc = path1;
	put rc=;
run;
data _null_;
	rc=dlgcdir();
	put rc=;
run;
* VALIDATION;
proc iml;
* Set seed for generation of random examples;
CALL RANDSEED(22032023); 
* Number of examples;
n_examples = 10;
* Formula for the probability of moving to the next stage;
start pgo(x) global(d2, kappa);
   return(cdf('normal',
				(x - kappa)/
				sqrt(4/d2)
			  )
		 );
finish;
* Sample size formula;
start d3(x) global(z_alpha, z_beta);
   return( 4*(z_alpha + z_beta)**2
			/x**2);
finish;
* Probability density function;
start f(x) global(d2, theta);
	return(pdf('normal', x, theta, sqrt(4/d2)));
finish;
* Sample size formula times density;
start d3_times_f(x);
	return(d3(x)*f(x));
finish;

* Randomly generate input values;
HRgo_vec = randfun(n_examples, "Uniform", 0.5, 0.95);   
d2_vec = randfun(n_examples, "Integer", 10, 100);   
hr1_vec = randfun(n_examples, "Uniform", 0.5, 0.99);
alpha_vec = randfun(n_examples, "Uniform", 0.01, 0.2);
beta_vec = randfun(n_examples, "Uniform", 0.01, 0.3);

* Pre-allocate results vector;
pgo_sas = j(nrow(HRgo_vec),1,.);
d3_sas = j(nrow(HRgo_vec),1,.);
* Generate results for each parameter combination;
do i = 1 to nrow(d3_sas);
	* Generate derived variables;
	HRgo = HRgo_vec[i];
	kappa = -log(HRgo);
	d2 = d2_vec[i];
	theta = -log(hr1_vec[i]);
	* Quantiles of normal distribution;
	alpha = alpha_vec[i];
	z_alpha = quantile('NORMAL', 1 - alpha);
	beta = beta_vec[i];
	z_beta = quantile('NORMAL', 1 - beta);
	* Calculate probability to go to phase 3;
	pgo_sas[i] = pgo(theta);
	* Integrate for receiving expected sample size;
	limits = {0 .P};
	limits[1,1] = kappa;
	call quad(R, "d3_times_f", limits);
	* Write results;
	d3_sas[i] = R;
end;

* Write data set;
create dat_out var {"pgo_sas" "d3_sas" "HRgo_vec" "d2_vec" "hr1_vec" "alpha_vec" "beta_vec"};
append;
close dat_out;
* Save data set;
options replace;
libname out '.';
data out.valref_tte;
set dat_out;
run;
