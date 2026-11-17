*******************************************************************************
* hmbivpois_examples.do
* Runs the examples of help hmbivpois on the simulated dataset
* hmbivpois_example.dta (1,000 individuals, presample + 6 periods, 3 states).
* Install the command and get the files first:
*   net install hmbivpois, from("https://raw.githubusercontent.com/plidonni/hmbivpois/main/") replace
*   net get hmbivpois, from("https://raw.githubusercontent.com/plidonni/hmbivpois/main/") replace
*******************************************************************************

version 16
clear all
set more off

capture log close hmbex
log using "hmbivpois_examples.log", replace text name(hmbex)

use hmbivpois_example, clear
xtset id time
describe
summarize y1 y2 x11 x12 x21 x22 z1
tabulate state if time > 0

* base model with three states
hmbivpois (y1 = x11 x12) (y2 = x21 x22), id(id) time(time) k(3) initcov(z1)

* postestimation
estat ic
estat transition
estat ame
predict double pp*, postprob
predict s_hat, viterbi
tabulate state s_hat if time > 0
estat lrtrans

* choice of the number of states
forvalues k = 1/4 {
    quietly hmbivpois (y1 = x11 x12) (y2 = x21 x22), id(id) time(time) k(`k') initcov(z1) nolog
    estat ic
}

* Frank copula, transitions with a covariate, log1p lags, slope of x11 by state
hmbivpois (y1 = x11 x12, statevary(x11)) (y2 = x21 x22), id(id) time(time) ///
    k(3) initcov(z1) dependence(frank) xtrans(z1) lagtransform(log1p)
estat stability

log close hmbex
