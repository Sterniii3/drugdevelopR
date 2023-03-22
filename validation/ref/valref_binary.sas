proc iml;
* Set seed for generation of random examples;
CALL RANDSEED(22032023); 
* Number of examples;
n_examples = 10;

* Sample size formula;
start n3(x) global(n2, p, p0, p11, z_alpha, z_beta);
   return( 2*(z_alpha*sqrt(2*(1-p)/p)
			  + z_beta*sqrt((1-p0)/p0
			  + (1-p11)/p11))**2
			/x**2);
finish;
* Probability density function;
start f(x) global(n2, p0, p11, rho);
	return(pdf('normal', x, rho, sqrt(2/n2*((1-p0)/p0 + (1-p11)/p11))));
finish;
* Sample size formula times density;
start n3_times_f(x);
	return(n3(x)*f(x));
finish;

* Randomly generate input values;
RRgo_vec = randfun(n_examples, "Uniform", 0.5, 0.95);   
n2_vec = randfun(n_examples, "Integer", 10, 100);   
p0_vec = randfun(n_examples, "Uniform", 0.5, 0.99);
p11_vec = randfun(n_examples, "Uniform", 0.1, 0.49);
alpha_vec = randfun(n_examples, "Uniform", 0.01, 0.2);
beta_vec = randfun(n_examples, "Uniform", 0.01, 0.3);

* Pre-allocate results vector;
res_sas = j(nrow(RRgo_vec),1,.);
* Generate results for each parameter combination;
do i = 1 to nrow(res_sas);
	* Generate derived variables;
	RRgo = RRgo_vec[i];
	kappa = -log(RRgo);
	n2 = n2_vec[i];
	p0 = p0_vec[i];
	p11 = p11_vec[i];
	rho = -log(p11/p0);
	p = (p0 + p11)/2;
	* Quantiles of normal distribution;
	alpha = alpha_vec[i];
	z_alpha = quantile('NORMAL', 1 - alpha);
	beta = beta_vec[i];
	z_beta = quantile('NORMAL', 1 - beta);
	* Integrate for receiving expected sample size;
	limits = {0 .P};
	limits[1,1] = kappa;
	call quad(R, "n3_times_f", limits);
	* Write results;
	res_sas[i] = R;
end;

* Write data set;
create dat_out var {"res_sas" "RRgo_vec" "n2_vec" "p0_vec" "p11_vec" "alpha_vec" "beta_vec"};
append;
close dat_out;
* Save data set;
options replace;
libname out '.';
data out.valref_binary;
set dat_out;
run;
