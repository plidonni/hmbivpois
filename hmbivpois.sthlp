{smcl}
{* *! version 1.0.1  18nov2026}{...}
{vieweralsosee "poisson" "help poisson"}{...}
{vieweralsosee "xtset" "help xtset"}{...}
{viewerjumpto "Syntax" "hmbivpois##syntax"}{...}
{viewerjumpto "Description" "hmbivpois##description"}{...}
{viewerjumpto "Model" "hmbivpois##model"}{...}
{viewerjumpto "Options" "hmbivpois##options"}{...}
{viewerjumpto "Postestimation" "hmbivpois##postestimation"}{...}
{viewerjumpto "Stored results" "hmbivpois##results"}{...}
{viewerjumpto "Examples" "hmbivpois##examples"}{...}
{viewerjumpto "References" "hmbivpois##references"}{...}
{viewerjumpto "Citation" "hmbivpois##citation"}{...}
{title:Title}

{phang}
{bf:hmbivpois} {hline 2} Hidden Markov bivariate Poisson regression for panel counts


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:hmbivpois}
{cmd:(}{it:depvar1} {cmd:=} [{it:indepvars1}] [{cmd:,} {it:eqopts}]{cmd:)}
{cmd:(}{it:depvar2} {cmd:=} [{it:indepvars2}] [{cmd:,} {it:eqopts}]{cmd:)}
{ifin}{cmd:,} {opth id(varname)} {opth time(varname)} {opt k(#)} [{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr:eqopts}
{synoptline}
{synopt:{opth statevary(varlist)}}regressors of the equation whose coefficients differ across states{p_end}
{synopt:{opt lagvary(depvars)}}lagged outcomes ({it:depvar1}, {it:depvar2}) whose coefficients differ across states{p_end}
{synoptline}

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opth id(varname)}}individual identifier{p_end}
{p2coldent:* {opth time(varname)}}integer time variable{p_end}
{p2coldent:* {opt k(#)}}number of latent states, 1 to 20{p_end}
{synopt:{opth initcov(varlist)}}covariates of the initial probabilities of the states{p_end}
{synopt:{opt dep:endence(form)}}dependence between the two counts: {opt famoye} (default), {opt frank} or {opt normal}{p_end}
{synopt:{opt lambdas:tate}}one dependence parameter per state{p_end}
{synopt:{opt lagt:ransform(g)}}transformation of the lagged counts: {opt linear} (default), {opt log1p}, {opt logc(#)}, {opt logzero}, {opt indicator}, {opt topcode(#)}{p_end}

{syntab:Transitions}
{synopt:{opth xtrans(varlist)}}covariates of the transition probabilities, measured at the period of arrival{p_end}
{synopt:{opt tvtrans}}transition probabilities that vary freely by period{p_end}
{synopt:{opt notr:ansition}}no transitions: the state of each individual is constant (latent class model){p_end}

{syntab:Estimation}
{synopt:{opt emiter:ate(#)}}maximum number of EM iterations; default {cmd:500}{p_end}
{synopt:{opt emtol:erance(#)}}tolerance of the EM on the relative change of the log likelihood; default {cmd:1e-8}{p_end}
{synopt:{opt emonly}}no direct maximization after the EM{p_end}
{synopt:{opt tech:nique(name)}}{opt bfgs} (default), {opt dfp} or {opt nr} for the direct maximization; specifying it forces the direct maximization{p_end}
{synopt:{opt iter:ate(#)}}maximum number of iterations of the direct maximization; default {cmd:300}{p_end}
{synopt:{opt from(matname)}}starting values, a row vector in the order of {cmd:e(b)}{p_end}
{synopt:{opt sort:by(how)}}labelling of the states: {opt eq1} (default), {opt eq2} or {opt none}{p_end}

{syntab:SE/Reporting}
{synopt:{opt vce(vcetype)}}{opt oim} (default) or {opt robust} (clustered on {cmd:id()}){p_end}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt nolo:g}}do not display the iteration log{p_end}
{synoptline}
{p 4 6 2}* {opt id()}, {opt time()} and {opt k()} are required.{p_end}
{p 4 6 2}{it:indepvars1}, {it:indepvars2}, {opt statevary()} and {opt xtrans()} may contain factor variables; see {help fvvarlist}.{p_end}
{p 4 6 2}{cmd:hmbivpois} is made of three files that must be in the same folder:
{cmd:hmbivpois.ado} (the whole command), {cmd:hmbivpois_p.ado} (for {cmd:predict}) and {cmd:hmbivpois_estat.ado} (for {cmd:estat}).{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:hmbivpois} fits the hidden Markov bivariate Poisson regression of Laudicella and Li Donni (2022) to
two counts observed on the same individuals over time, such as the visits to primary care and the
emergency admissions. Given a latent state, the two counts follow a bivariate distribution with Poisson
margins, whose means depend on regressors, on the lagged counts and on the state; the state follows a
first-order Markov chain.

{pstd}
The data are in long form, one observation per individual and period. The lagged counts are built by
the command: the first observed period of each individual is the presample period, which provides the
lags of the first estimation period and the presample counts used in the initial probabilities, and is
not an estimation period. The periods of an individual must be consecutive. Regressors that are missing
in the presample period are not a problem.

{pstd}
A technical guide with the model, the estimation method, the postestimation tools and Monte Carlo
evidence is in {browse "https://github.com/plidonni/hmbivpois/blob/main/doc/hmbivpois_guide.pdf":doc/hmbivpois_guide.pdf}
on the GitHub page of the command.


{marker model}{...}
{title:Model}

{pstd}
For individual {it:i}, period {it:t} and state {it:c} = 1, ..., {it:k}:

{p 8 8 2}
log mu_jitc = x_jit b_j + g_j1 G(y1_i,t-1) + g_j2 G(y2_i,t-1) + a_jc,  j = 1, 2

{pstd}
where G() is the transformation of {opt lagtransform()} and a_jc are the state intercepts (reported as
{cmd:state1}, ..., {cmd:state}{it:k} in each equation). With {opt statevary()} and {opt lagvary()} the
listed slopes are specific to the state and are reported in the equations {it:depvar}{cmd:_s1}, ...,
{it:depvar}{cmd:_s}{it:k}.

{pstd}
The joint probability of the two counts in a state has Poisson margins with means mu_1 and mu_2 and
dependence parameter lambda:

{p2colset 9 22 24 2}{...}
{p2col:{opt famoye}}p(y1,y2) = Po(y1;mu1) Po(y2;mu2) [1 + lambda (exp(-y1) - exp(-q mu1)) (exp(-y2) - exp(-q mu2))],
q = 1 - exp(-1) (Lakshminarayana, Pandit and Rao 1999; Famoye 2010). lambda must lie in an interval
that depends on the means; the command keeps it inside the interval and reports it in
{cmd:e(lambda_range)}. A lambda on a bound is held fixed in the variance estimate.{p_end}
{p2col:{opt frank}}p(y1,y2) is the probability that the Frank copula gives to the rectangle
(F1(y1-1), F1(y1)] x (F2(y2-1), F2(y2)], with F_j the Poisson distribution functions (Cameron et al. 2004).
The parameter {cmd:theta} is unrestricted; 0 means independence.{p_end}
{p2col:{opt normal}}the same with the normal copula. The parameter is {cmd:athrho} = atanh(rho); the
correlation rho is in {cmd:e(rho)}.{p_end}
{p2colreset}{...}

{pstd}
The initial probabilities of the states are a multinomial logit (base: state 1) on {opt initcov()},
the transformed presample counts and a constant (equations {cmd:init2}, ..., {cmd:init}{it:k}). The
transition probabilities are logits by row, with staying in the state as the base (equations
{cmd:tr}{it:c}{cmd:_}{it:d}): constant over time by default, functions of {opt xtrans()} and/or of the
period with {opt tvtrans}.

{pstd}
Estimation is by the EM algorithm, with scaled forward-backward recursions, followed by a direct
maximization of the log likelihood with its analytic score when the EM does not converge within
{opt emiterate()} iterations or when {opt technique()} is specified. The variance is the inverse of
the observed information, computed by differences of the analytic score; {cmd:vce(robust)} gives the
sandwich estimator clustered on the individuals. States are labelled in increasing order of the mean
fitted mu_1c ({cmd:sortby(eq1)}) or mu_2c ({cmd:sortby(eq2)}).


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opth initcov(varlist)} are covariates of the initial probabilities of the states, taken at the first
estimation period of each individual.

{phang}
{opt dependence(form)} sets the bivariate distribution in each state; see {help hmbivpois##model:Model}.

{phang}
{opt lambdastate} estimates one dependence parameter (lambda, theta or athrho) per state.

{phang}
{opt lagtransform(g)} applies G() to both lagged counts in both equations and to the presample counts in
the initial probabilities: {opt linear} G(y) = y; {opt log1p} ln(1+y); {opt logc(#)} ln(max(y,#)), 0 < # <= 1;
{opt logzero} ln(y) for y > 0, 0 otherwise, plus an indicator of y = 0; {opt indicator} 1(y > 0);
{opt topcode(#)} min(y,#). With linear lags the model can be explosive; see {cmd:estat stability}.

{dlgtab:Transitions}

{phang}
{opth xtrans(varlist)} adds covariates to the logits of the transitions into each period.

{phang}
{opt tvtrans} gives a separate transition matrix for each period of arrival.

{phang}
{opt notransition} fits the model without transitions between states; {cmd:estat lrtrans} tests it.

{dlgtab:Estimation}

{phang}
{opt emiterate(#)}, {opt emtolerance(#)}, {opt emonly}, {opt technique()} and {opt iterate(#)} control
the EM and the direct maximization. {cmd:emiterate(0) emonly} evaluates the model at {opt from()}.

{phang}
{opt from(matname)} gives starting values in the order and number of {cmd:e(b)}. By default the
command starts from Poisson regressions with the state intercepts spread around the constant and a
persistent chain.

{phang}
{opt sortby(how)} sets the labelling of the states; {opt none} keeps the order of the estimation.


{marker postestimation}{...}
{title:Postestimation}

{pstd}
The following commands use the estimates and the data in memory; the data must be those of the
estimation (the presample period of each individual included). Predictions are for the estimation sample.

{dlgtab:predict}

{p 8 16 2}
{cmd:predict} [{it:type}] {it:newvar} | {it:stub}{cmd:*} [{it:if}] [{it:in}] [{cmd:,} {it:statistic} {it:weights}]

{synoptset 20}{...}
{synopt:{opt mu1}}expected value of {it:depvar1}; the default{p_end}
{synopt:{opt mu2}}expected value of {it:depvar2}{p_end}
{synopt:{opt postprob}}probabilities of the states ({it:k} variables: {it:stub}{cmd:*} or {it:k} names){p_end}
{synopt:{opt state}}most likely state of each period{p_end}
{synopt:{opt viterbi}}most likely sequence of states of each individual (Viterbi algorithm){p_end}
{synopt:{opt rho}}correlation between the two counts given the state{p_end}
{synopt:{opt pr(a,b)}}probability that the two counts equal (a, b){p_end}

{pstd}
The states are weighted with their posterior probabilities given all the outcomes of the individual
(default), or with those given the outcomes up to the period ({opt filtered}), up to the period before
({opt onestep}) or given no outcome ({opt prior}). {opt instate(#)} gives the statistic in state #.
The sum over the sample of ln {cmd:pr(}{it:y1},{it:y2}{cmd:)} with {opt onestep} equals the log likelihood.

{dlgtab:estat}

{phang}
{cmd:estat ic} [{cmd:,} {opt n(#)}] reports the log likelihood, the number of parameters, AIC and BIC with
n = number of individuals, and the relative entropy of the classification into the states
(1 = states known with certainty).

{phang}
{cmd:estat transition} [{cmd:,} {opt period(#)} {opt level(#)}] reports the transition probabilities, the
expected duration of a stay in each state, 1/(1 - P(c -> c)), and the mean initial probabilities, with
delta-method standard errors; the confidence intervals of the probabilities are computed on the logit
scale. With heterogeneous transitions it reports the mean over the transition rows, or over those into
period # with {opt period()}.

{phang}
{cmd:estat lrtrans} [{cmd:,} {opt reps(#)} {opt seed(#)}] tests the model without transitions
({opt notransition}) against the model with homogeneous transitions in memory. The null hypothesis is
on the boundary of the parameter space: the LR statistic has a chi-bar-squared distribution, a mixture
of chi2(0), ..., chi2(k(k-1)), whose weights are estimated by {opt reps()} simulations (default 10000)
from the information of the off-diagonal transition probabilities at the restricted estimates
(Bartolucci 2006). The p-value of the chi2(k(k-1)) distribution, which is conservative, is also
reported. The restricted model is fitted from two sets of starting values and the larger log
likelihood is kept. Under the null hypothesis the maximum of the model with transitions is on the
boundary, where the EM converges slowly and may stop early, which makes LR too small and the test
conservative. No score test is reported: the individual scores of the transition probabilities are
sums of likelihood ratios between states, with very heavy tails, and a score statistic based on their
outer product is badly calibrated.

{phang}
{cmd:estat stability} [{cmd:,} {opt horizon(#)} {opt limit(#)} {opt nsim(#)} {opt seed(#)}] checks
whether the dynamics of the lagged counts can explode: the coefficients of the lags by state; their
spectral radius (log-type lags) or the index g_jj exp(a), which must not exceed exp(-1) (linear lags);
the share of deterministic paths, from every observation and in every state, whose expected count
exceeds {opt limit()} (default 1000) within {opt horizon()} periods (default 50); and the same share
for {opt nsim()} simulated paths per individual of the estimated model.

{phang}
{cmd:estat ame} [{cmd:,} {opt prior} {opt onestep} {opt filtered} {opt instate(#)} {opt level(#)}] reports
the average marginal effects on E({it:depvar1}) and E({it:depvar2}) with the probabilities of the states
held fixed: dy/dx for continuous regressors (variables coded 0/1 without {cmd:i.} are treated as
continuous), the discrete change from the base level for the levels of factor variables, and the
effect of one more count of each lagged outcome. Delta-method standard errors. The weights of the
states are those of {cmd:predict}.

{pstd}
{cmd:estat vce} and {cmd:estat summarize} work as after any estimation command. {cmd:margins} is not
supported; use {cmd:estat ame}.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:hmbivpois} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations (estimation periods){p_end}
{synopt:{cmd:e(N_g)}}number of individuals{p_end}
{synopt:{cmd:e(k)}}number of states{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(npar)}}number of parameters{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(aic)}, {cmd:e(bic)}}information criteria (BIC with n = individuals){p_end}
{synopt:{cmd:e(converged)}}1 if converged{p_end}
{synopt:{cmd:e(em_iter)}}number of EM iterations{p_end}
{synopt:{cmd:e(g_min)}, {cmd:e(g_avg)}, {cmd:e(g_max)}}periods per individual{p_end}
{synopt:{cmd:e(lambda_lo)}, {cmd:e(lambda_hi)}}admissible range of lambda (Famoye form){p_end}
{synopt:{cmd:e(n_bound)}}number of lambda parameters on a bound{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters ({cmd:vce(robust)}){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:hmbivpois}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}the two dependent variables{p_end}
{synopt:{cmd:e(ivar)}, {cmd:e(tvar)}}individual and time variables{p_end}
{synopt:{cmd:e(dependence)}}{cmd:famoye}, {cmd:frank} or {cmd:normal}{p_end}
{synopt:{cmd:e(deppar)}}{cmd:lambda}, {cmd:theta} or {cmd:athrho}{p_end}
{synopt:{cmd:e(transitions)}}{cmd:homogeneous}, {cmd:heterogeneous} or {cmd:none}{p_end}
{synopt:{cmd:e(lagtransform)}}transformation of the lags{p_end}
{synopt:{cmd:e(statevary1)}, {cmd:e(statevary2)}}terms with state-specific coefficients{p_end}
{synopt:{cmd:e(vce)}, {cmd:e(technique)}, {cmd:e(sortby)}}options used{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}coefficients and their variance{p_end}
{synopt:{cmd:e(Pi)}}transition matrix (mean over the transition rows if heterogeneous){p_end}
{synopt:{cmd:e(Pi_t)}}transition matrices by period of arrival (heterogeneous transitions){p_end}
{synopt:{cmd:e(pi0)}}mean initial probabilities{p_end}
{synopt:{cmd:e(omega)}, {cmd:e(omega_t)}}mean posterior probabilities of the states, overall and by period{p_end}
{synopt:{cmd:e(alpha)}}state intercepts{p_end}
{synopt:{cmd:e(lambda_range)}}admissible range of lambda by state (Famoye form){p_end}
{synopt:{cmd:e(rho)}}correlation of the normal copula{p_end}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}estimation sample{p_end}

{pstd}
{cmd:estat lrtrans} stores {cmd:r(lr)}, {cmd:r(p)}, {cmd:r(p_chi2)}, {cmd:r(df)}, {cmd:r(ll_1)},
{cmd:r(ll_0)}, {cmd:r(reps)} and {cmd:r(weights)}; {cmd:estat ic} {cmd:r(S)}, {cmd:r(entropy)} and
{cmd:r(N_obs)}; {cmd:estat transition} {cmd:r(P)}, {cmd:r(se_P)}, {cmd:r(pi0)}, {cmd:r(se_pi0)} and
{cmd:r(duration)}; {cmd:estat stability} {cmd:r(G)}, {cmd:r(radius)}, {cmd:r(index)}, {cmd:r(det)},
{cmd:r(sim)} and {cmd:r(q)}; {cmd:estat ame} {cmd:r(b)} and {cmd:r(V)}.


{marker examples}{...}
{title:Examples}

{pstd}
The examples use {cmd:hmbivpois_example.dta}, a simulated panel of 1,000 individuals observed in a
presample period and 6 estimation periods, generated from the model with three states, the Famoye
form (lambda = -0.5) and homogeneous transitions. The variable {cmd:state} holds the true latent state.
Download it into the current folder with{p_end}
{phang2}{cmd:. net get hmbivpois, from("https://raw.githubusercontent.com/plidonni/hmbivpois/main/")}{p_end}

{pstd}Setup{p_end}
{phang2}{cmd:. use hmbivpois_example, clear}{p_end}
{phang2}{cmd:. xtset id time}{p_end}

{pstd}Base model with three states{p_end}
{phang2}{cmd:. hmbivpois (y1 = x11 x12) (y2 = x21 x22), id(id) time(time) k(3) initcov(z1)}{p_end}

{pstd}Postestimation{p_end}
{phang2}{cmd:. estat ic}{p_end}
{phang2}{cmd:. estat transition}{p_end}
{phang2}{cmd:. estat ame}{p_end}
{phang2}{cmd:. predict double pp*, postprob}{p_end}
{phang2}{cmd:. predict s_hat, viterbi}{p_end}
{phang2}{cmd:. tabulate state s_hat if time > 0}{p_end}
{phang2}{cmd:. estat lrtrans}{p_end}

{pstd}Choice of the number of states{p_end}
{phang2}{cmd:. forvalues k = 1/4 {c -(}}{p_end}
{phang2}{cmd:.     quietly hmbivpois (y1 = x11 x12) (y2 = x21 x22), id(id) time(time) k(`k') initcov(z1) nolog}{p_end}
{phang2}{cmd:.     estat ic}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

{pstd}Frank copula, transitions with a covariate, log1p lags, slope of x11 by state{p_end}
{phang2}{cmd:. hmbivpois (y1 = x11 x12, statevary(x11)) (y2 = x21 x22), id(id) time(time) k(3) initcov(z1) dependence(frank) xtrans(z1) lagtransform(log1p)}{p_end}
{phang2}{cmd:. estat stability}{p_end}


{marker references}{...}
{title:References}

{phang}
Bartolucci, F. 2006. Likelihood inference for a class of latent Markov models under linear hypotheses
on the transition probabilities. {it:Journal of the Royal Statistical Society, Series B} 68: 155-178.

{phang}
Cameron, A. C., T. Li, P. K. Trivedi, and D. M. Zimmer. 2004. Modelling the differences in counted
outcomes using bivariate copula models with application to mismeasured counts.
{it:Econometrics Journal} 7: 566-584.

{phang}
Famoye, F. 2010. On the bivariate negative binomial regression model. {it:Journal of Applied Statistics}
37: 969-981.

{phang}
Lakshminarayana, J., S. N. N. Pandit, and K. Srinivasa Rao. 1999. On a bivariate Poisson distribution.
{it:Communications in Statistics - Theory and Methods} 28: 267-276.

{phang}
Laudicella, M., and P. Li Donni. 2022. The dynamic interdependence in the demand of primary and emergency
secondary care: A hidden Markov approach. {it:Journal of Applied Econometrics} 37: 521-536.


{marker citation}{...}
{title:Citation}

{pstd}
If you use {cmd:hmbivpois}, please cite Laudicella and Li Donni (2022) and the software:{p_end}
{phang2}Li Donni, P. 2026. hmbivpois: Stata module to fit hidden Markov bivariate Poisson regression
for panel counts. https://github.com/plidonni/hmbivpois.{p_end}


{title:Author}

{pstd}
Paolo Li Donni{break}
University of Palermo{break}
{browse "mailto:paolo.lidonni@unipa.it":paolo.lidonni@unipa.it}{break}
{browse "https://sites.google.com/site/plidonni/"}{break}
{browse "https://github.com/plidonni/hmbivpois"}{break}
Please report bugs and suggestions by email.{p_end}
