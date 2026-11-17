*! version 1.0.0  17nov2026
*! predict after hmbivpois. The program is in hmbivpois.ado (subprogram Predict),
*! which holds the whole command: this file only hands the request over.
*!
*!   predict [type] {newvar | stub* | newvarlist} [if] [in] [, statistic weights]
*!
*!   statistic   mu1        expected value of the first count (the default)
*!               mu2        expected value of the second count
*!               postprob   probabilities of the latent states (k variables)
*!               state      most likely state of each period
*!               viterbi    most likely sequence of states of each individual
*!               rho        correlation between the two counts given the state
*!               pr(a,b)    probability that the two counts equal (a, b)
*!   weights     (default)  posterior probabilities of the states, given all the
*!                          outcomes of the individual
*!               filtered   given the outcomes up to the period
*!               onestep    given the outcomes up to the period before
*!               prior      given the covariates and the presample outcomes only
*!               instate(#) the statistic in state #

program define hmbivpois_p
    version 16
    hmbivpois __predict `0'
end
