*! version 1.0.0  17nov2026
*! estat after hmbivpois. The program is in hmbivpois.ado (subprogram Estat),
*! which holds the whole command: this file only hands the request over.
*!
*!   estat lrtrans [, reps(#) seed(#)]
*!       likelihood-ratio test of no transitions between the states (latent
*!       class model) against homogeneous transitions, with the
*!       chi-bar-squared null distribution; reps() simulations for its weights
*!       (default 10000), seed() for the random numbers (the state of the
*!       generator is restored afterwards)
*!   estat ic [, n(#)]
*!       log likelihood, AIC, BIC with n = number of individuals (or n()),
*!       relative entropy of the classification into the states
*!   estat transition [, period(#) level(#)]
*!       transition probabilities, expected durations and mean initial
*!       probabilities, with delta-method standard errors; with heterogeneous
*!       transitions the mean over the transition rows (into period # with
*!       period())
*!   estat stability [, horizon(#) limit(#) nsim(#) seed(#)]
*!       diagnostics of explosive dynamics of the lagged counts: coefficients
*!       of the lags by state, spectral radius (log-type lags) or index of the
*!       univariate map (linear lags), deterministic and simulated paths over
*!       horizon() periods (default 50), share above limit() (default 1000)
*!   estat ame [, prior onestep filtered instate(#) level(#)]
*!       average marginal effects on E(y1) and E(y2), with the probabilities of
*!       the states held fixed (posterior by default): dy/dx for continuous
*!       regressors, discrete change from the base level for factor levels,
*!       one more count of each lagged outcome; delta-method standard errors
*!   estat vce, estat summarize: as after any estimation command

program define hmbivpois_estat
    version 16
    hmbivpois __estat `0'
end
