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
start pgo(x) global(n2, kappa);
   return(cdf('normal',
				(x - kappa)/
				sqrt(4/n2)
			  )
		 );
finish;
* Sample size formula;
start n3(x) global(n2, p, Delta1, z_alpha, z_beta);
   return( 4*(z_alpha + z_beta)**2
			/x**2);
finish;
* Probability density function;
start f(x) global(n2, Delta1);
	return(pdf('normal', x, Delta1, sqrt(4/n2)));
finish;
* Sample size formula times density;
start n3_times_f(x);
	return(n3(x)*f(x));
finish;

* Randomly generate input values;
kappa_vec = randfun(n_examples, "Uniform", 0.02, 0.4);   
n2_vec = randfun(n_examples, "Integer", 10, 100);   
Delta1_vec = randfun(n_examples, "Uniform", 0.3, 0.7);
alpha_vec = randfun(n_examples, "Uniform", 0.01, 0.2);
beta_vec = randfun(n_examples, "Uniform", 0.01, 0.3);

* Pre-allocate results vector;
pgo_sas = j(nrow(kappa_vec),1,.);
n3_sas = j(nrow(kappa_vec),1,.);
* Generate results for each parameter combination;
do i = 1 to nrow(n3_sas);
	* Generate derived variables;
	kappa = kappa_vec[i];
	n2 = n2_vec[i];
	Delta1 = Delta1_vec[i];
	* Quantiles of normal distribution;
	alpha = alpha_vec[i];
	z_alpha = quantile('NORMAL', 1 - alpha);
	beta = beta_vec[i];
	z_beta = quantile('NORMAL', 1 - beta);
	* Calculate probability to go to phase 3;
	pgo_sas[i] = pgo(Delta1);
	* Integrate for receiving expected sample size;
	limits = {0 .P};
	limits[1,1] = kappa;
	call quad(R, "n3_times_f", limits);
	* Write results;
	n3_sas[i] = R;
end;

* Write data set;
create dat_out var {"pgo_sas" "n3_sas" "kappa_vec" "n2_vec" "Delta1_vec" "alpha_vec" "beta_vec"};
append;
close dat_out;
* Save data set;
options replace;
libname out '.';
data out.valref_normal;
set dat_out;
run;
