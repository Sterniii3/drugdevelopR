#### Values presented in table
* expected utility (u)
* assumed true hazard rate (HR)
* threshold value for the decision rule to go to phase III (HR<sub>go</sub>)
* total number of events for phase II (d<sub>2</sub>)
* total expected number of events for phase III (d<sub>3</sub>)
* total expected number of events in the program (d)
* total sample size for phase II (n<sub>2</sub>, rounded to the next equal natural number)
* total expected sample size for phase III (n<sub>3</sub>, rounded to the next equal natural number)
* total expected sample size in the program (n = n<sub>2</sub> + n<sub>3</sub>)
* maximal costs of the program (K)
* probability to go to phase III (p<sub>go</sub>)
* probability of a successful program (sProg)
* probability of a successful program with small treatment effect in Phase III (sProg1)
* probability of a successful program with medium treatment effect in Phase III (sProg2)
* probability of a successful program with large treatment effect in Phase III (sProg3)
* expected costs for phase II (K2)
* expected costs for phase III (K3)

and further input parameters.


#### References

Kirchner, M., Kieser, M., Goette, H., & Schueler, A. (2016). Utility-based optimization of phase II/III programs. <i>Statistics in Medicine</i>, 35(2), 305-316.

IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at https://www.iqwig.de/de/methoden/methodenpapier.3020.html, assessed last 15.05.19.

Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival distributions. <i>Biometrika</i>, 68(1), 316-319.

#### Additional Features
The true underlying treatment effect &theta; can also be modelled by a prior distribution (https://web.imbi.uni-heidelberg.de/prior/). Optimal planning can then be done with the function optimal_tte() of the R package drugdevelopR available at: https://github.com/Sterniii3/drugdevelopR.

#### Note

If the server is busy, you may need to double click the "Go"-button in order to see the updated plot.

#### Maintainer

Stella Preussler, Institute of Medical Biometry and Informatics, University of Heidelberg, email: Preussler@imbi.uni-heidelberg.de.

Beta version 0.5





