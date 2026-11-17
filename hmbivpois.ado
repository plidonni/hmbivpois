*! version 1.0.0  17nov2026
*! hmbivpois: hidden Markov bivariate Poisson regression for panel counts
*! Model: Laudicella & Li Donni (2022, Journal of Applied Econometrics 37: 521-536)
*! Author: Paolo Li Donni, https://github.com/plidonni/hmbivpois

/*
    Base model
    ----------
    y_it = (y1_it, y2_it) | U_it = c  ~  bivariate Poisson (Lakshminarayana et al. 1999; Famoye 2010)
        p(y1,y2) = Po(y1;mu1) Po(y2;mu2) [1 + lambda (exp(-y1)-exp(-q mu1)) (exp(-y2)-exp(-q mu2))]
        log mu_jitc = x_jit b_j + g_j1 y1_i,t-1 + g_j2 y2_i,t-1 + alpha_jc,   q = 1 - exp(-1)
    U_i1 ~ multinomial logit on (initcov, y1_i0, y2_i0), base = state 1
    U_it | U_i,t-1 ~ first-order Markov chain (logit rows, base = staying);
        homogeneous by default; xtrans(varlist) adds covariates measured at the
        arrival period, tvtrans adds period dummies (a free matrix per period)

    Lags are built by the command: the first observed period of each individual
    is the presample (y_i0) and is not used as an estimation period.
    lagtransform() applies g() to both lags, in both equations, and to y_i0 in the
    initial probabilities: linear (default), log1p = ln(1+y), logc(#) = ln(max(y,#)),
    logzero = ln(y) for y>0 plus a dummy for y=0, indicator = 1(y>0),
    topcode(#) = min(y,#).

    Dependence: lambda is common to the states; lambdastate estimates one lambda
    per state. Lambda must lie in the interval that keeps the pmf nonnegative at
    the fitted means (by state with lambdastate). The M-step projects it onto
    that interval; a lambda that ends on a bound is reported, and held fixed in
    the variance estimate.
    dependence(frank) and dependence(normal) replace the Famoye form with a copula
    applied to the Poisson distribution functions: p(y1,y2) is the probability
    that the copula gives to the rectangle (F1(y1-1), F1(y1)] x (F2(y2-1), F2(y2)].
    Frank: the parameter theta is unrestricted (theta = 0: independence).
    Normal: the parameter is athrho = atanh(rho), with rho the correlation of the
    copula (athrho = 0: independence); the probability of a cell is a
    one-dimensional integral, by 12-point Gauss-Legendre quadrature.
    With lambdastate the parameter is one per state. The M-step uses the outer
    product of the scores in place of the second derivatives.

    State-specific coefficients: inside each equation, statevary(varlist) names
    the regressors and lagvary(depvars) the lagged outcomes whose coefficients
    differ across states, (depvar = indepvars, statevary(...) lagvary(...)).
    They are reported in the equations <depvar>_s1, <depvar>_s2, ...
    States are labelled in increasing order of the mean over the sample of the
    fitted mu_1c (sortby(eq1), the default) or mu_2c (sortby(eq2)); with common
    slopes this is the order of the state intercepts.

    Estimation: EM (forward-backward with scaled recursions, vectorised over
    individuals by period). Direct maximisation of the marginal likelihood with
    the analytic score (Fisher identity) follows only if the EM does not converge
    within emiterate() or if technique() is specified; emonly never runs it. Variance: inverse of the observed
    information, from central differences of the analytic score; vce(robust)
    gives the sandwich clustered on individuals.

    Postestimation: predict (file hmbivpois_p.ado, which calls this program)
    rebuilds the data of the estimation from e(sample) and computes, for the
    estimation sample: the expected counts (mu1, mu2), the probabilities of the
    latent states (postprob), the most likely state of each period (state), the
    most likely sequence of states (viterbi), the correlation between the two
    counts given the state (rho) and the probability of a pair of counts
    (pr(a,b)). The states are weighted with their posterior probabilities given
    all the outcomes of the individual (the default), or with the filtered,
    one-step-ahead or prior probabilities (filtered, onestep, prior);
    instate(#) conditions on one state.

    notransition fits the model without transitions between states (the state
    of an individual does not change over time: latent class model).
    estat lrtrans (file hmbivpois_estat.ado, which calls this program) tests it
    against the model with homogeneous transitions: the LR statistic has a
    chi-bar-squared distribution, whose weights are estimated by simulation
    (Bartolucci 2006).
    estat ic: AIC, BIC (n = individuals) and relative entropy of the
    classification; estat transition: transition and initial probabilities
    with delta-method standard errors, expected durations; estat stability:
    diagnostics of explosive dynamics of the lagged counts.
    estat ame: average marginal effects on the expected counts, with the
    probabilities of the states held fixed, and delta-method standard errors.
*/

program define hmbivpois, eclass
    version 16
    if replay() {
        if "`e(cmd)'" != "hmbivpois" error 301
        Display `0'
        exit
    }
    // entry point of predict (hmbivpois_p.ado)
    gettoken sub rest : 0
    if `"`sub'"' == "__predict" {
        if "`e(cmd)'" != "hmbivpois" error 301
        Predict `rest'
        exit
    }
    // entry point of estat (hmbivpois_estat.ado)
    if `"`sub'"' == "__estat" {
        if "`e(cmd)'" != "hmbivpois" error 301
        Estat `rest'
        exit
    }
    Estimate `0'
end

program define Estimate, eclass sortpreserve
    version 16
    local cmdline `"hmbivpois `0'"'

    // ---------------------------------------------------------------- parsing
    gettoken eq1 0 : 0, match(paren)
    if "`paren'" != "(" {
        di as err "equation 1 must be enclosed in parentheses: (depvar1 = indepvars1)"
        exit 198
    }
    gettoken eq2 0 : 0, match(paren)
    if "`paren'" != "(" {
        di as err "equation 2 must be enclosed in parentheses: (depvar2 = indepvars2)"
        exit 198
    }
    ParseEq `eq1'
    local y1 `s(dep)'
    local x1 `s(indep)'
    local svspec1 `s(statevary)'
    local lvspec1 `s(lagvary)'
    ParseEq `eq2'
    local y2 `s(dep)'
    local x2 `s(indep)'
    local svspec2 `s(statevary)'
    local lvspec2 `s(lagvary)'
    if "`y1'" == "`y2'" {
        di as err "the two dependent variables must be different"
        exit 198
    }

    syntax [if] [in], ID(varname numeric) TIME(varname numeric)             ///
        [ K(integer 0) INITcov(varlist numeric) EMITERate(integer 500)     ///
          EMTOLerance(real 1e-8) ITERate(integer 300) TECHnique(string)    ///
          FROM(name) SORTby(string) VCE(string) Level(cilevel)             ///
          NOLOg EMONLY XTRANS(varlist numeric fv) TVTRANS             ///
          LAGTransform(string) LAMBDAState DEPendence(string)          ///
          NOTRansition ]

    if `k' == 0 {
        di as err "option k() required"
        exit 198
    }
    if `k' < 1 | `k' > 20 {
        di as err "k() must be between 1 and 20"
        exit 198
    }
    // direct maximisation: 0 = never, 1 = only if EM did not converge, 2 = always
    local doml = 1
    if "`technique'" != "" local doml = 2
    if "`emonly'" != ""    local doml = 0
    if "`technique'" == "" local technique "bfgs"
    if !inlist("`technique'", "bfgs", "dfp", "nr") {
        di as err "technique() must be bfgs, dfp or nr"
        exit 198
    }
    if "`sortby'" == "" local sortby "eq1"
    if !inlist("`sortby'", "eq1", "eq2", "none") {
        di as err "sortby() must be eq1, eq2 or none"
        exit 198
    }
    if "`vce'" == "" local vce "oim"
    if !inlist("`vce'", "oim", "robust") {
        di as err "vce() must be oim or robust"
        exit 198
    }
    if `emiterate' < 0 | `iterate' < 0 {
        di as err "emiterate() and iterate() must be nonnegative"
        exit 198
    }
    // form of the dependence between the two counts, and name of its parameter
    if "`dependence'" == "" local dependence "famoye"
    if !inlist("`dependence'", "famoye", "normal", "frank") {
        di as err "dependence() must be famoye, normal or frank"
        exit 198
    }
    local depcode = 0
    local depname "lambda"
    if "`dependence'" == "normal" {
        local depcode = 1
        local depname "athrho"
    }
    if "`dependence'" == "frank" {
        local depcode = 2
        local depname "theta"
    }
    // dependence parameter: common to the states, or one per state
    local lamst = ("`lambdastate'" != "")
    if `k' == 1 & `lamst' {
        di as txt "note: lambdastate ignored with k(1)"
        local lamst = 0
        local lambdastate
    }
    // state-specific coefficients: statevary() for regressors, lagvary() for lags
    if `k' == 1 & "`svspec1'`svspec2'`lvspec1'`lvspec2'" != "" {
        di as txt "note: statevary() and lagvary() ignored with k(1)"
        local svspec1
        local svspec2
        local lvspec1
        local lvspec2
    }
    forvalues j = 1/2 {
        // lagS`j': lags (1 = first dependent variable, 2 = second) whose
        // coefficients are state specific in equation `j'
        local lagS`j'
        foreach v of local lvspec`j' {
            capture unab vv : `v', max(1)
            if _rc local vv `v'
            if "`vv'" == "`y1'" {
                local lagS`j' `lagS`j'' 1
            }
            else if "`vv'" == "`y2'" {
                local lagS`j' `lagS`j'' 2
            }
            else {
                di as err "lagvary(): `v' is not one of the dependent variables"
                exit 198
            }
        }
        local lagS`j' : list uniq lagS`j'
    }
    local het = ("`xtrans'" != "" | "`tvtrans'" != "")
    if `k' == 1 & `het' {
        di as txt "note: xtrans() and tvtrans ignored with k(1)"
        local het = 0
        local xtrans
        local tvtrans
    }
    local tvt = ("`tvtrans'" != "")
    // no transitions: the state of each individual is constant over time
    local notr = ("`notransition'" != "")
    if `k' == 1 & `notr' {
        di as txt "note: notransition ignored with k(1)"
        local notr = 0
    }
    if `notr' & `het' {
        di as err "notransition may not be combined with xtrans() or tvtrans"
        exit 198
    }

    // transformation g() of the lagged outcomes (and of the presample values)
    local lagt : subinstr local lagtransform " " "", all
    if "`lagt'" == "" local lagt "linear"
    local lagfull "`lagt'"
    local lagarg
    if regexm("`lagt'", "^(logc|topcode)\(([0-9.]+)\)$") {
        local lagarg = regexs(2)
        local lagt   = regexs(1)
    }
    if !inlist("`lagt'", "linear", "log1p", "logc", "logzero", "indicator", "topcode") {
        di as err "lagtransform() must be linear, log1p, logc(#), logzero, indicator or topcode(#)"
        exit 198
    }
    if inlist("`lagt'", "logc", "topcode") & "`lagarg'" == "" {
        di as err "lagtransform(`lagt'(#)) requires a number"
        exit 198
    }
    if "`lagt'" == "logc" {
        if !(`lagarg' > 0 & `lagarg' <= 1) {
            di as err "lagtransform(logc(#)): # must be in (0, 1]"
            exit 198
        }
    }
    if "`lagt'" == "topcode" {
        if !(`lagarg' >= 1) {
            di as err "lagtransform(topcode(#)): # must be at least 1"
            exit 198
        }
    }
    if "`from'" != "" confirm matrix `from'
    local dolog = ("`nolog'" == "")

    // ---------------------------------------------------------------- sample
    // touse: rows usable as observed periods (outcomes, id, time)
    // xok  : rows where all regressors are also observed (estimation periods)
    marksample touse, novarlist
    markout `touse' `y1' `y2' `id' `time'

    tempvar xok first gap gap2 L1 L2 G1 G2 Z1 Z2 es Tn negT
    qui gen byte `xok' = `touse'
    if "`x1'`x2'" != "" {
        fvrevar `x1' `x2', list
        markout `xok' `r(varlist)'
    }
    if "`initcov'" != "" markout `xok' `initcov'
    if "`xtrans'" != "" {
        fvrevar `xtrans', list
        markout `xok' `r(varlist)'
    }

    capture assert `y1' >= 0 & `y1' == floor(`y1') & ///
                   `y2' >= 0 & `y2' == floor(`y2') if `touse'
    if _rc {
        di as err "dependent variables must be nonnegative integers"
        exit 459
    }
    capture assert `time' == floor(`time') if `touse'
    if _rc {
        di as err "`time' must be integer-valued"
        exit 459
    }

    sort `touse' `id' `time'
    qui by `touse' `id': gen byte `gap' = (`time' != `time'[_n-1] + 1) ///
        if `touse' & _n > 1
    qui count if `gap' == 1
    if r(N) > 0 {
        di as err "repeated or non-consecutive values of `time' within `id';"
        di as err "gaps in the observed periods are not allowed"
        exit 459
    }
    qui by `touse' `id': gen byte `first' = (_n == 1) if `touse'
    qui by `touse' `id': gen double `L1' = `y1'[_n-1] if `touse'
    qui by `touse' `id': gen double `L2' = `y2'[_n-1] if `touse'
    forvalues j = 1/2 {
        if "`lagt'" == "linear"    qui gen double `G`j'' = `L`j''
        if "`lagt'" == "log1p"     qui gen double `G`j'' = ln(1 + `L`j'')
        if "`lagt'" == "logc"      qui gen double `G`j'' = ln(max(`L`j'', `lagarg')) if `L`j'' < .
        if "`lagt'" == "indicator" qui gen double `G`j'' = (`L`j'' > 0) if `L`j'' < .
        if "`lagt'" == "topcode"   qui gen double `G`j'' = min(`L`j'', `lagarg') if `L`j'' < .
        if "`lagt'" == "logzero" {
            qui gen double `G`j'' = cond(`L`j'' > 0, ln(`L`j''), 0) if `L`j'' < .
            qui gen double `Z`j'' = (`L`j'' == 0) if `L`j'' < .
        }
    }
    local lagvars  `G1' `G2'
    local lagnames `=substr("lag_`y1'", 1, 32)' `=substr("lag_`y2'", 1, 32)'
    if "`lagt'" == "logzero" {
        local lagvars  `lagvars' `Z1' `Z2'
        local lagnames `lagnames' `=substr("lag0_`y1'", 1, 32)' `=substr("lag0_`y2'", 1, 32)'
    }
    qui gen byte `es' = `touse' & `first' == 0 & `xok'

    // estimation periods must be consecutive within individual
    sort `es' `id' `time'
    qui by `es' `id': gen byte `gap2' = (`time' != `time'[_n-1] + 1) ///
        if `es' & _n > 1
    qui count if `gap2' == 1
    if r(N) > 0 {
        di as err "missing regressors inside an individual's sequence of periods;"
        di as err "estimation periods must be consecutive"
        exit 459
    }
    qui by `es' `id': gen long `Tn' = _N if `es'
    // transition rows: estimation periods after the first one of each individual
    tempvar tr
    qui by `es' `id': gen byte `tr' = `es' & _n > 1
    if `het' {
        qui count if `tr'
        if r(N) == 0 {
            di as err "no transitions observed: xtrans() and tvtrans need T > 1"
            exit 459
        }
    }
    qui count if `es'
    local N = r(N)
    if `N' == 0 error 2000
    if inlist("`lagt'", "indicator", "logzero") {
        // the lag must take both zero and positive values in the sample
        forvalues j = 1/2 {
            qui count if `es' & `L`j'' == 0
            if r(N) == 0 | r(N) == `N' {
                di as err "lagtransform(`lagt'): lagged `y`j'' is never zero or always zero in the sample"
                exit 459
            }
        }
    }

    // ------------------------------------------------- regressors (fv aware)
    forvalues j = 1/2 {
        local xnames`j'
        local xvars`j'
        if "`x`j''" != "" {
            qui _rmcoll `x`j'' if `es', expand
            local full `r(varlist)'
            foreach t of local full {
                if !ustrregexm("`t'", "(^|#)[0-9]*[bo]\.") {
                    local xnames`j' `xnames`j'' `t'
                }
            }
            // one variable per term (TermVar: indicators of the levels of
            // factor variables are built explicitly)
            foreach t of local xnames`j' {
                if !ustrregexm("`t'", "[.#]") {
                    local xvars`j' `xvars`j'' `t'
                    continue
                }
                tempvar tv
                TermVar `tv' `t'
                local xvars`j' `xvars`j'' `tv'
            }
        }
        // regressors with state-specific coefficients: the terms of statevary(),
        // which must be among the regressors of the equation
        local svterm`j'
        if "`svspec`j''" != "" {
            fvexpand `svspec`j'' if `es'
            local svfull `r(varlist)'
            foreach t of local svfull {
                if !ustrregexm("`t'", "(^|#)[0-9]*[bo]\.") {
                    local isin : list t in xnames`j'
                    if !`isin' {
                        di as err "statevary(): `t' is not among the regressors of equation `j'"
                        exit 198
                    }
                    local svterm`j' `svterm`j'' `t'
                }
            }
        }
        // common terms (cnam, cvar) and state-specific terms (snam, svar):
        // regressors first, then lags
        local cnam`j'
        local cvar`j'
        local snam`j'
        local svar`j'
        local i = 0
        foreach t of local xnames`j' {
            local ++i
            local v : word `i' of `xvars`j''
            local isin : list t in svterm`j'
            if `isin' {
                local snam`j' `snam`j'' `t'
                local svar`j' `svar`j'' `v'
            }
            else {
                local cnam`j' `cnam`j'' `t'
                local cvar`j' `cvar`j'' `v'
            }
        }
        local i = 0
        foreach t of local lagnames {
            local ++i
            local v : word `i' of `lagvars'
            // terms 1 and 3 refer to the lag of the first count, 2 and 4 to the second
            local li = 2 - mod(`i', 2)
            local isin : list li in lagS`j'
            if `isin' {
                local snam`j' `snam`j'' `t'
                local svar`j' `svar`j'' `v'
            }
            else {
                local cnam`j' `cnam`j'' `t'
                local cvar`j' `cvar`j'' `v'
            }
        }
    }

    // covariates of the transition probabilities (measured at the arrival period)
    local xnamesT
    local xvarsT
    if "`xtrans'" != "" {
        qui _rmcoll `xtrans' if `tr', expand
        local full `r(varlist)'
        foreach t of local full {
            if !ustrregexm("`t'", "(^|#)[0-9]*[bo]\.") {
                local xnamesT `xnamesT' `t'
            }
        }
        foreach t of local xnamesT {
            if !ustrregexm("`t'", "[.#]") {
                local xvarsT `xvarsT' `t'
                continue
            }
            tempvar tv
            TermVar `tv' `t'
            local xvarsT `xvarsT' `tv'
        }
    }

    // ------------------------------------------------- starting values (k=1)
    tempname bs1 bs2
    if "`from'" == "" {
        // coefficients in the order: common terms, state-specific terms, _cons
        qui poisson `y1' `cvar1' `svar1' if `es'
        matrix `bs1' = e(b)
        qui poisson `y2' `cvar2' `svar2' if `es'
        matrix `bs2' = e(b)
    }

    // ------------------------------------------------------------ estimation
    qui gen long `negT' = -`Tn'
    sort `es' `negT' `id' `time'

    mata: hmbp_main("`y1'", "`y2'", "`cvar1'",                              ///
        "`cvar2'", "`initcov' `lagvars'", "`id'", "`es'", `k',              ///
        "`bs1'", "`bs2'", "`from'", `emiterate', `emtolerance', `iterate',  ///
        "`technique'", `dolog', "`sortby'", "`vce'", `doml',                ///
        "`xvarsT'", "`tr'", "`time'", `het', `tvt', "`svar1'", "`svar2'",   ///
        `lamst', `depcode', `notr')

    // --------------------------------------------------------------- results
    tempname b V Pi om omt al pi0 lrng
    matrix `b'   = __hmbp_b
    matrix `V'   = __hmbp_V
    matrix `Pi'  = __hmbp_Pi
    matrix `om'  = __hmbp_om
    matrix `omt' = __hmbp_omt
    matrix `al'  = __hmbp_al
    matrix `pi0' = __hmbp_pi0
    matrix `lrng' = __hmbp_lrng
    capture matrix drop __hmbp_b __hmbp_V __hmbp_Pi __hmbp_om __hmbp_omt ///
        __hmbp_al __hmbp_pi0 __hmbp_lrng

    // names: in each equation the common terms and the state intercepts; the
    // state-specific coefficients follow in the equations <depvar>_s1, _s2, ...
    local names
    forvalues j = 1/2 {
        foreach v of local cnam`j' {
            local names `names' `y`j'':`v'
        }
        forvalues c = 1/`k' {
            local names `names' `y`j'':state`c'
        }
        if "`snam`j''" != "" {
            forvalues c = 1/`k' {
                local eqs = substr("`y`j''", 1, 32 - strlen("_s`c'")) + "_s`c'"
                foreach v of local snam`j' {
                    local names `names' `eqs':`v'
                }
            }
        }
    }
    if `lamst' {
        forvalues c = 1/`k' {
            local names `names' `depname':state`c'
        }
    }
    else local names `names' `depname':_cons
    if `k' > 1 {
        forvalues c = 2/`k' {
            foreach v of local initcov {
                local names `names' init`c':`v'
            }
            foreach v of local lagnames {
                local names `names' init`c':`v'
            }
            local names `names' init`c':_cons
        }
        local trn `xnamesT'
        local tlab
        foreach tv of local hmbp_tvals {
            local tlab `tlab' t_`=subinstr("`tv'", "-", "m", .)'
        }
        if `tvt' local trn `trn' `tlab'
        else     local trn `trn' _cons
        if !`notr' {
            forvalues c = 1/`k' {
                forvalues d = 1/`k' {
                    if `c' != `d' {
                        foreach v of local trn {
                            local names `names' tr`c'_`d':`v'
                        }
                    }
                }
            }
        }
    }
    matrix colnames `b' = `names'
    matrix colnames `V' = `names'
    matrix rownames `V' = `names'

    local snames
    forvalues c = 1/`k' {
        local snames `snames' state`c'
    }
    matrix rownames `Pi'  = `snames'
    matrix colnames `Pi'  = `snames'
    matrix colnames `om'  = `snames'
    matrix colnames `omt' = `snames'
    matrix colnames `pi0' = `snames'
    matrix rownames `om'  = share
    matrix rownames `pi0' = prob
    matrix rownames `al'  = `snames'
    matrix colnames `al'  = `y1' `y2'
    matrix rownames `lrng' = `snames'
    matrix colnames `lrng' = lower upper
    local tnames
    forvalues t = 1/`=rowsof(`omt')' {
        local tnames `tnames' t`t'
    }
    matrix rownames `omt' = `tnames'

    if `het' {
        tempname Pit
        matrix `Pit' = __hmbp_Pit
        capture matrix drop __hmbp_Pit
        local pnames
        forvalues c = 1/`k' {
            forvalues d = 1/`k' {
                local pnames `pnames' p`c'_`d'
            }
        }
        matrix rownames `Pit' = `tlab'
        matrix colnames `Pit' = `pnames'
    }

    ereturn post `b' `V', esample(`es') obs(`N')
    ereturn scalar N_g       = __hmbp_n
    ereturn scalar k         = `k'
    ereturn scalar ll        = __hmbp_ll
    ereturn scalar npar      = __hmbp_npar
    ereturn scalar rank      = __hmbp_rank
    ereturn scalar bic       = -2*__hmbp_ll + __hmbp_npar*ln(__hmbp_n)
    ereturn scalar aic       = -2*__hmbp_ll + 2*__hmbp_npar
    ereturn scalar converged = __hmbp_conv
    ereturn scalar em_iter   = __hmbp_emit
    ereturn scalar eval_only = (`emiterate' == 0 & `doml' == 0)
    ereturn scalar g_min     = __hmbp_gmin
    ereturn scalar g_avg     = __hmbp_gavg
    ereturn scalar g_max     = __hmbp_gmax
    ereturn scalar lambda_lo = __hmbp_lamlo
    ereturn scalar lambda_hi = __hmbp_lamhi
    ereturn scalar n_bound   = __hmbp_nbnd
    if "`vce'" == "robust" {
        ereturn scalar N_clust = __hmbp_n
        ereturn local clustvar "`id'"
        ereturn local vcetype  "Robust"
    }
    capture scalar drop __hmbp_n __hmbp_ll __hmbp_npar __hmbp_rank   ///
        __hmbp_conv __hmbp_emit __hmbp_gmin __hmbp_gavg __hmbp_gmax  ///
        __hmbp_lamlo __hmbp_lamhi __hmbp_nbnd
    ereturn matrix Pi      = `Pi'
    if `het' ereturn matrix Pi_t = `Pit'
    ereturn matrix omega   = `om'
    ereturn matrix omega_t = `omt'
    ereturn matrix alpha   = `al'
    if "`dependence'" == "famoye" ereturn matrix lambda_range = `lrng'
    ereturn matrix pi0     = `pi0'
    ereturn local vce        "`vce'"
    ereturn local technique  "`technique'"
    ereturn local dependence "`dependence'"
    ereturn local deppar     "`depname'"
    if "`dependence'" == "normal" {
        // correlation of the normal copula, rho = tanh(athrho)
        tempname rho
        if `lamst' {
            matrix `rho' = J(1, `k', .)
            forvalues c = 1/`k' {
                matrix `rho'[1, `c'] = tanh(_b[athrho:state`c'])
            }
            matrix colnames `rho' = `snames'
        }
        else {
            matrix `rho' = J(1, 1, tanh(_b[athrho:_cons]))
            matrix colnames `rho' = rho
        }
        matrix rownames `rho' = rho
        ereturn matrix rho = `rho'
    }
    ereturn local lambdastate "`lambdastate'"
    ereturn local lagtransform "`lagfull'"
    ereturn local statevary1 "`snam1'"
    ereturn local statevary2 "`snam2'"
    if `het'       ereturn local transitions "heterogeneous"
    else if `notr' ereturn local transitions "none"
    else           ereturn local transitions "homogeneous"
    ereturn local xtrans     "`xtrans'"
    ereturn local tvtrans    "`tvtrans'"
    ereturn local sortby     "`sortby'"
    ereturn local initcov    "`initcov'"
    ereturn local tvar       "`time'"
    ereturn local ivar       "`id'"
    ereturn local depvar     "`y1' `y2'"
    ereturn local title      "Hidden Markov bivariate Poisson regression"
    // names of the terms, used by predict to rebuild the design matrices
    ereturn local common1      "`cnam1'"
    ereturn local common2      "`cnam2'"
    ereturn local lagnames     "`lagnames'"
    ereturn local xtrans_terms "`xnamesT'"
    ereturn local tvals        "`hmbp_tvals'"
    ereturn local marginsnotok "_ALL"
    ereturn local predict      "hmbivpois_p"
    ereturn local estat_cmd    "hmbivpois_estat"
    ereturn local cmdline    `"`cmdline'"'
    ereturn local cmd        "hmbivpois"

    Display, level(`level')
end

// Variable of one term in factor-variable notation, generated as `name':
// level.var gives the indicator of the level (missing where var is missing),
// c.var and var the variable, # the product of the parts. Other operators are
// left to fvrevar. The indicators are built here because fvrevar, given a list
// of levels without the base, may take the first level listed as the base.
program define TermVar
    version 16
    args name term
    mata: st_local("tex", hmbp_termexpr(`"`term'"'))
    if `"`tex'"' != "" {
        qui gen double `name' = `tex'
    }
    else {
        fvrevar `term'
        qui gen double `name' = `r(varlist)'
    }
end

program define ParseEq, sclass
    version 16
    sreturn clear
    local eq : subinstr local 0 "=" " ", all
    gettoken dep rest : eq
    capture unab dep : `dep', max(1)
    if _rc {
        di as err "`dep': dependent variable not found"
        exit 111
    }
    confirm numeric variable `dep'
    local 0 `rest'
    syntax [varlist(numeric fv default=none)]                              ///
        [, STATEVary(varlist numeric fv) LAGVary(namelist)]
    sreturn local dep       "`dep'"
    sreturn local indep     "`varlist'"
    sreturn local statevary "`statevary'"
    sreturn local lagvary   "`lagvary'"
end

program define Display
    version 16
    syntax [, Level(cilevel)]
    di _n as txt "`e(title)'" _col(49) "Number of obs" _col(67) "=" ///
        _col(69) as res %10.0fc e(N)
    local dep "Famoye (Sarmanov)"
    if "`e(dependence)'" == "frank"  local dep "Frank copula"
    if "`e(dependence)'" == "normal" local dep "normal copula"
    local dpar "`e(deppar)'"
    if "`dpar'" == "" local dpar "lambda"
    if "`e(lambdastate)'" != "" local dep "`dep', `dpar' by state"
    di as txt "Dependence: " as res "`dep'" _col(49) as txt ///
        "Number of groups" _col(67) "=" _col(69) as res %10.0fc e(N_g)
    di as txt "Group variable: " as res abbrev("`e(ivar)'", 16) ///
        _col(49) as txt "Obs per group:"
    di as txt "Lagged outcomes: " as res "`e(lagtransform)'"
    di as txt "Latent states:  " as res e(k) ///
        _col(63) as txt "min" _col(67) "=" _col(69) as res %10.0fc e(g_min)
    di _col(63) as txt "avg" _col(67) "=" _col(69) as res %10.1f e(g_avg)
    di _col(63) as txt "max" _col(67) "=" _col(69) as res %10.0fc e(g_max)
    di as txt "Log likelihood = " as res %12.0g e(ll) ///
        _col(49) as txt "BIC" _col(67) "=" _col(69) as res %10.0g e(bic)
    forvalues j = 1/2 {
        if "`e(statevary`j')'" != "" {
            local dv : word `j' of `e(depvar)'
            di as txt "State-specific coefficients in equation " as res "`dv'" ///
                as txt ": " as res "`e(statevary`j')'"
        }
    }
    di
    ereturn display, level(`level')
    if e(k) > 1 {
        if "`e(transitions)'" == "none" {
            di as txt _n "No transitions between the states (latent class model)"
        }
        else {
            if "`e(transitions)'" == "heterogeneous" {
                di as txt _n "Average transition probabilities (rows: state at t-1; columns: state at t)"
            }
            else {
                di as txt _n "Transition probabilities (rows: state at t-1; columns: state at t)"
            }
            matlist e(Pi), format(%9.4f) border(rows) rowtitle("from")
        }
        di as txt _n "Average posterior state shares"
        matlist e(omega), format(%9.4f) border(rows)
    }
    if "`e(dependence)'" == "normal" {
        di as txt _n "Correlation of the normal copula, rho = tanh(athrho)"
        matlist e(rho), format(%9.4f) border(rows)
    }
    // the parameter of a copula is unrestricted: the range concerns lambda only
    if !inlist("`e(dependence)'", "frank", "normal") {
        if "`e(lambdastate)'" != "" {
            di as txt _n "Admissible range of lambda at the fitted means, by state"
            matlist e(lambda_range), format(%9.4f) border(rows)
        }
        else {
            di as txt _n "Admissible range of lambda at the fitted means: [" ///
                as res %7.4f e(lambda_lo) as txt ", " as res %7.4f e(lambda_hi) as txt "]"
        }
    }
    if e(n_bound) > 0 & e(n_bound) < . {
        di as txt "note: " as res e(n_bound) as txt " lambda on a bound of its admissible range. It is held fixed in the"
        di as txt "      variance estimate (no standard error); the other standard errors are"
        di as txt "      conditional on it."
    }
    if e(eval_only) == 1 {
        di as txt "note: evaluation at the starting values, no iterations performed"
    }
    else if e(converged) == 0 {
        di as err "warning: convergence not achieved"
    }
end

// ---------------------------------------------------------------------------
// predict (called by hmbivpois_p.ado)
//   predict [type] {newvar | stub* | newvarlist} [if] [in] [, statistic weights]
//   statistic: mu1 (default), mu2, postprob (k variables), state, viterbi, rho,
//              pr(a,b)
//   weights:   posterior probabilities of the states (default), filtered,
//              onestep, prior, or instate(#)
// Predictions are computed for the estimation sample.
// ---------------------------------------------------------------------------
program define Predict
    version 16
    syntax anything(name=vlist id="newvarlist") [if] [in]                  ///
        [, MU1 MU2 POSTprob STate VITerbi RHO PR(string)                   ///
           PRIor ONEstep FILTered INSTate(string) ]

    local k  = e(k)
    local y1 : word 1 of `e(depvar)'
    local y2 : word 2 of `e(depvar)'

    // statistic
    local nstat = ("`mu1'" != "") + ("`mu2'" != "") + ("`postprob'" != "") + ///
        ("`state'" != "") + ("`viterbi'" != "") + ("`rho'" != "") + (`"`pr'"' != "")
    if `nstat' > 1 {
        di as err "only one statistic may be specified"
        exit 198
    }
    local stat "`mu1'`mu2'`postprob'`state'`viterbi'`rho'"
    local ya = 0
    local yb = 0
    if `"`pr'"' != "" {
        local stat "pr"
        local prs : subinstr local pr "," " ", all
        capture numlist "`prs'", min(2) max(2) integer range(>=0)
        if _rc {
            di as err "pr(a,b): a and b must be two nonnegative integers"
            exit 198
        }
        local prs `r(numlist)'
        local ya : word 1 of `prs'
        local yb : word 2 of `prs'
    }
    if `nstat' == 0 {
        local stat "mu1"
        di as txt "(option mu1 assumed; expected value of `y1')"
    }

    // probabilities of the states used as weights
    local nw = ("`prior'" != "") + ("`onestep'" != "") + ("`filtered'" != "")
    if `nw' > 1 {
        di as err "only one of prior, onestep and filtered may be specified"
        exit 198
    }
    local ptype = 0
    if "`filtered'" != "" local ptype = 1
    if "`onestep'"  != "" local ptype = 2
    if "`prior'"    != "" local ptype = 3
    local ist = 0
    if "`instate'" != "" {
        capture confirm integer number `instate'
        if _rc {
            di as err "instate() must be an integer between 1 and `k'"
            exit 198
        }
        local ist = `instate'
        if `ist' < 1 | `ist' > `k' {
            di as err "instate() must be an integer between 1 and `k'"
            exit 198
        }
        if inlist("`stat'", "postprob", "state", "viterbi") {
            di as err "instate() may not be combined with `stat'"
            exit 198
        }
        if `nw' > 0 {
            di as err "instate() may not be combined with prior, onestep or filtered"
            exit 198
        }
    }
    if "`stat'" == "viterbi" & `nw' > 0 {
        di as err "viterbi may not be combined with prior, onestep or filtered"
        exit 198
    }

    // new variables: k for postprob (stub* or k names), one otherwise
    local need = cond("`stat'" == "postprob", `k', 1)
    gettoken first rest : vlist
    local type "`c(type)'"
    if inlist("`first'", "byte", "int", "long", "float", "double") {
        local type "`first'"
        local vlist `rest'
    }
    local nv : word count `vlist'
    if `nv' == 1 {
        if substr("`vlist'", -1, 1) == "*" {
            local stub = substr("`vlist'", 1, strlen("`vlist'") - 1)
            local vlist
            forvalues c = 1/`need' {
                local vlist `vlist' `stub'`c'
            }
            local nv = `need'
        }
    }
    if `nv' != `need' {
        if `need' == 1 di as err "one new variable name required"
        else           di as err "postprob requires `need' new variable names, or stub*"
        exit 198
    }
    confirm new variable `vlist'

    // the results are computed in temporary variables on the estimation sample
    marksample touse, novarlist
    local rv
    forvalues i = 1/`need' {
        tempvar r`i'
        qui gen double `r`i'' = .
        local rv `rv' `r`i''
    }
    PostRun, task(predict) out(`rv') stat(`stat') ptype(`ptype')           ///
        instate(`ist') ya(`ya') yb(`yb')

    local wlab "posterior"
    if `ptype' == 1 local wlab "filtered"
    if `ptype' == 2 local wlab "one-step-ahead"
    if `ptype' == 3 local wlab "prior"
    local i = 0
    foreach v of local vlist {
        local ++i
        qui gen `type' `v' = `r`i'' if `touse'
        if "`stat'" == "postprob" label variable `v' "Pr(state `i'), `wlab'"
        if "`stat'" == "state"    label variable `v' "Most likely state, `wlab' probabilities"
        if "`stat'" == "viterbi"  label variable `v' "Most likely sequence of states (Viterbi)"
        if inlist("`stat'", "mu1", "mu2") {
            local dv = cond("`stat'" == "mu1", "`y1'", "`y2'")
            if `ist' > 0 label variable `v' "Expected `dv' in state `ist'"
            else         label variable `v' "Expected `dv', `wlab' probabilities"
        }
        if "`stat'" == "rho" {
            if `ist' > 0 label variable `v' "Correlation of the counts in state `ist'"
            else         label variable `v' "Correlation of the counts given the state, `wlab'"
        }
        if "`stat'" == "pr" {
            if `ist' > 0 label variable `v' "Pr(`y1'=`ya', `y2'=`yb') in state `ist'"
            else         label variable `v' "Pr(`y1'=`ya', `y2'=`yb'), `wlab' probabilities"
        }
    }
    local v1 : word 1 of `vlist'
    qui count if `v1' >= .
    if r(N) > 0 di as txt "(" r(N) " missing values generated)"
end

// ---------------------------------------------------------------------------
// Postestimation engine: rebuilds the data of the estimation from e() and the
// variables in memory, and runs a task in Mata. The estimation sample is
// e(sample); the lags of the first estimation period of each individual come
// from the observation of the period before, which must be identified by the
// individual and the period.
// ---------------------------------------------------------------------------
program define PostRun, sortpreserve
    version 16
    syntax , TASK(string) [ OUT(string) STAT(string) PTYPE(integer 0)      ///
        INSTATE(integer 0) YA(integer 0) YB(integer 0) REPS(integer 0)     ///
        SEED(integer 0) PER(string) HORizon(integer 50) LIMit(real 1000)   ///
        NSIM(integer 1) ]

    local y1 : word 1 of `e(depvar)'
    local y2 : word 2 of `e(depvar)'
    local id   "`e(ivar)'"
    local time "`e(tvar)'"
    local k = e(k)
    local initcov "`e(initcov)'"
    foreach v in `y1' `y2' `id' `time' `initcov' {
        capture confirm numeric variable `v'
        if _rc {
            di as err "variable `v' not found; the data must be those used by hmbivpois"
            exit 111
        }
    }
    local depcode = 0
    if "`e(dependence)'" == "normal" local depcode = 1
    if "`e(dependence)'" == "frank"  local depcode = 2
    local lamst = ("`e(lambdastate)'" != "")
    local het   = ("`e(transitions)'" == "heterogeneous")
    local tvt   = ("`e(tvtrans)'" != "")
    local notr  = ("`e(transitions)'" == "none")

    // ---------------------------------------------------------------- sample
    tempvar es t0 pre npre touse L1 L2 G1 G2 Z1 Z2 Tn negT tr
    qui gen byte `es' = e(sample)
    qui count if `es'
    if r(N) == 0 | r(N) != e(N) {
        di as err "estimation sample not available; the data must be those used by hmbivpois"
        exit 459
    }
    // the observation that precedes the first estimation period of each individual
    sort `id' `es' `time'
    qui by `id' `es': gen double `t0' = `time'[1] if `es'
    qui by `id': replace `t0' = `t0'[_N]
    qui gen byte `pre' = (`es' == 0 & `t0' < . & `time' == `t0' - 1 &      ///
        `y1' < . & `y2' < .)
    qui by `id': gen long `npre' = sum(`pre')
    qui by `id': replace `npre' = `npre'[_N]
    capture assert `npre' == 1 if `es'
    if _rc {
        di as err "cannot identify the period before the first estimation period of each"
        di as err "individual: `id' and `time' must identify the observations, and the data"
        di as err "must be those used by hmbivpois"
        exit 459
    }
    qui gen byte `touse' = (`es' | `pre')

    // lags and their transformation, as in the estimation
    local lagt "`e(lagtransform)'"
    local lagarg
    if regexm("`lagt'", "^(logc|topcode)\(([0-9.]+)\)$") {
        local lagarg = regexs(2)
        local lagt   = regexs(1)
    }
    sort `touse' `id' `time'
    qui by `touse' `id': gen double `L1' = `y1'[_n-1] if `touse'
    qui by `touse' `id': gen double `L2' = `y2'[_n-1] if `touse'
    forvalues j = 1/2 {
        if "`lagt'" == "linear"    qui gen double `G`j'' = `L`j''
        if "`lagt'" == "log1p"     qui gen double `G`j'' = ln(1 + `L`j'')
        if "`lagt'" == "logc"      qui gen double `G`j'' = ln(max(`L`j'', `lagarg')) if `L`j'' < .
        if "`lagt'" == "indicator" qui gen double `G`j'' = (`L`j'' > 0) if `L`j'' < .
        if "`lagt'" == "topcode"   qui gen double `G`j'' = min(`L`j'', `lagarg') if `L`j'' < .
        if "`lagt'" == "logzero" {
            qui gen double `G`j'' = cond(`L`j'' > 0, ln(`L`j''), 0) if `L`j'' < .
            qui gen double `Z`j'' = (`L`j'' == 0) if `L`j'' < .
        }
    }
    local lagvars `G1' `G2'
    if "`lagt'" == "logzero" local lagvars `lagvars' `Z1' `Z2'
    local lagnames `e(lagnames)'

    // design variables, from the names of the terms stored in e(): common
    // terms (cvar) and terms with state-specific coefficients (svar)
    forvalues j = 1/2 {
        local cvar`j'
        local svar`j'
        local ccode`j'
        local scode`j'
        foreach t in `e(common`j')' {
            local pos : list posof "`t'" in lagnames
            local ccode`j' `ccode`j'' `pos'
            if `pos' > 0 {
                local v : word `pos' of `lagvars'
            }
            else if !ustrregexm("`t'", "[.#]") local v `t'
            else {
                tempvar tv
                TermVar `tv' `t'
                local v `tv'
            }
            local cvar`j' `cvar`j'' `v'
        }
        foreach t in `e(statevary`j')' {
            local pos : list posof "`t'" in lagnames
            local scode`j' `scode`j'' `pos'
            if `pos' > 0 {
                local v : word `pos' of `lagvars'
            }
            else if !ustrregexm("`t'", "[.#]") local v `t'
            else {
                tempvar tv
                TermVar `tv' `t'
                local v `tv'
            }
            local svar`j' `svar`j'' `v'
        }
    }
    // lagged counts, read by estat stability
    local rawlags `L1' `L2'
    local xvarsT
    foreach t in `e(xtrans_terms)' {
        if !ustrregexm("`t'", "[.#]") {
            local xvarsT `xvarsT' `t'
            continue
        }
        tempvar tv
        TermVar `tv' `t'
        local xvarsT `xvarsT' `tv'
    }

    // estimation periods of each individual, sorted as in the estimation
    sort `es' `id' `time'
    qui by `es' `id': gen long `Tn' = _N if `es'
    qui by `es' `id': gen byte `tr' = (`es' & _n > 1)
    qui gen long `negT' = -`Tn'
    sort `es' `negT' `id' `time'

    // estat ame: counterfactual designs. For each effect, two versions (a and
    // b) of the four lists of design variables (equation 1 and 2, common and
    // state-specific terms): the variable shifted by +h and -h (continuous),
    // or set to a level and to the base level (factor). The lags are shifted
    // in Mata. The names, the equations where each effect enters and its type
    // are passed back to the caller.
    if "`task'" == "ame" {
        mata: hmbp_amevars("`e(common1)' `e(statevary1)' `e(common2)' `e(statevary2)'", "`lagnames'")
        local ne = 0
        local i  = 0
        foreach v of local amev {
            local ++i
            local isf : word `i' of `amef'
            if `isf' {
                qui levelsof `v' if `es', local(alllev)
                local lvs `amelv`i''
                local base
                foreach l of local alllev {
                    local inl : list l in lvs
                    if !`inl' & "`base'" == "" local base `l'
                }
                if "`base'" == "" {
                    di as err "estat ame: base level of `v' not found in the estimation sample"
                    exit 459
                }
                foreach l of local lvs {
                    local ++ne
                    tempvar sa`ne' sb`ne'
                    qui gen double `sa`ne'' = `l'
                    qui gen double `sb`ne'' = `base'
                    local ame_name`ne' "`l'.`v'"
                    local ame_t`ne' "x"
                    local ame_div`ne' = 1
                    local ame_v`ne' "`v'"
                }
            }
            else {
                qui summarize `v' if `es'
                local h = 1e-3*cond(r(sd) > 0 & r(sd) < ., r(sd), 1)
                local ++ne
                tempvar sa`ne' sb`ne'
                qui gen double `sa`ne'' = `v' + `h'
                qui gen double `sb`ne'' = `v' - `h'
                local ame_name`ne' "`v'"
                local ame_t`ne' "x"
                local ame_div`ne' = 2*`h'
                local ame_v`ne' "`v'"
            }
        }
        forvalues e = 1/`ne' {
            local v "`ame_v`e''"
            local ame_f1_`e' = 0
            local ame_f2_`e' = 0
            foreach sd in a b {
                local rep "`s`sd'`e''"
                forvalues j = 1/2 {
                    foreach kd in common statevary {
                        local base = cond("`kd'" == "common", "`cvar`j''", "`svar`j''")
                        local lst
                        local p = 0
                        foreach t in `e(`kd'`j')' {
                            local ++p
                            local w : word `p' of `base'
                            mata: st_local("tt", hmbp_subst(`"`t'"', "`v'", "`rep'"))
                            if `"`tt'"' != `"`t'"' {
                                tempvar tv
                                TermVar `tv' `tt'
                                local w `tv'
                                local ame_f`j'_`e' = 1
                            }
                            local lst `lst' `w'
                        }
                        local ame_`sd'`j'`kd'`e' `lst'
                    }
                }
            }
        }
        // one more count of each lagged outcome
        forvalues l = 1/2 {
            local ++ne
            local ame_name`ne' : word `l' of `lagnames'
            local ame_t`ne' "lag`l'"
            local ame_div`ne' = 1
            local ame_f1_`ne' = 1
            local ame_f2_`ne' = 1
        }
        local ame_n = `ne'
        c_local ame_n `ne'
        forvalues e = 1/`ne' {
            c_local ame_name`e' "`ame_name`e''"
            c_local ame_f1_`e' `ame_f1_`e''
            c_local ame_f2_`e' `ame_f2_`e''
        }
    }

    mata: hmbp_post("`task'", "`y1'", "`y2'", "`cvar1'", "`cvar2'",         ///
        "`initcov' `lagvars'", "`id'", "`es'", `k', "`xvarsT'", "`tr'",    ///
        "`time'", `het', `tvt', "`svar1'", "`svar2'", `lamst', `depcode',  ///
        "`out'", "`stat'", `ptype', `instate', `ya', `yb', `notr', `reps',  ///
        `seed')
end

// ---------------------------------------------------------------------------
// estat (called by hmbivpois_estat.ado)
//   estat lrtrans [, reps(#) seed(#)]   test of no transitions between states
//   estat ic [, n(#)]                   information criteria, entropy
//   estat transition [, period(#) level(#)]   transition probabilities
//   estat stability [, horizon(#) limit(#) nsim(#) seed(#)]   explosive lags
//   estat ame [, prior onestep filtered instate(#) level(#)]   marginal effects
//   other subcommands: those of estat for every estimation command (vce,
//   summarize)
// ---------------------------------------------------------------------------
program define Estat
    version 16
    gettoken sub rest : 0, parse(" ,")
    if `"`sub'"' == "lrtrans" {
        EstatLR `rest'
    }
    else if `"`sub'"' == "ic" {
        EstatIC `rest'
    }
    else if inlist(`"`sub'"', "transition", "trans") {
        EstatTrans `rest'
    }
    else if inlist(`"`sub'"', "stability", "stab") {
        EstatStab `rest'
    }
    else if `"`sub'"' == "ame" {
        EstatAME `rest'
    }
    else estat_default `0'
end

// Likelihood-ratio test of H0: no transitions between states (all the
// off-diagonal transition probabilities are zero; the restricted model is a
// latent class model). H0 is on the boundary of the parameter space, so the
// LR statistic has a chi-bar-squared null distribution: a mixture of chi2(j),
// j = 0..k(k-1), whose weights are the probabilities that the projection of
// a N(0, inv(I)) vector onto the nonnegative orthant has j positive
// components, with I the information of the off-diagonal transition
// probabilities net of the other parameters, at the restricted estimates
// (Bartolucci 2006). The information is estimated by the outer product of the
// individual scores; the weights by simulation.
// The restricted model is fitted twice, from the default starting values and
// from the estimates of the model with transitions, and the larger log
// likelihood is kept.
program define EstatLR, rclass
    version 16
    syntax [, REPS(integer 10000) SEED(integer 20261116) ]
    if e(k) < 2 {
        di as err "estat lrtrans requires a model with two or more states"
        exit 198
    }
    if "`e(transitions)'" == "none" {
        di as err "estat lrtrans: the model in memory has no transitions (notransition);"
        di as err "fit the model with transitions and run estat lrtrans after it"
        exit 198
    }
    if "`e(transitions)'" != "homogeneous" {
        di as err "estat lrtrans requires homogeneous transitions (no xtrans() or tvtrans)"
        exit 198
    }
    if `reps' < 100 {
        di as err "reps() must be at least 100"
        exit 198
    }
    local k  = e(k)
    local nr = `k'*(`k'-1)
    tempname ll1 ll0 lla llb b1 bfrom ba W full
    scalar `ll1' = e(ll)
    local N1 = e(N)
    matrix `b1' = e(b)
    local nf = colsof(`b1') - `nr'
    matrix `bfrom' = `b1'[1, 1..`nf']

    // command line of the restricted model: the same equations, sample and
    // options, without from(), emiterate(), emonly and sortby()
    local cmd `"`e(cmdline)'"'
    gettoken c0 cmd : cmd
    gettoken eq1 cmd : cmd, match(par)
    gettoken eq2 cmd : cmd, match(par)
    local 0 `"`cmd'"'
    syntax [if] [in] [, *]
    local ifin `"`if' `in'"'
    local 0 `", `options'"'
    syntax [, FROM(string) NOLOg EMITERate(string) EMONLY SORTby(string) *]
    local opts0 `"`options'"'
    local opts  `"`options'"'
    if "`emiterate'" != "" local opts `"`opts' emiterate(`emiterate')"'
    if "`emonly'" != ""    local opts `"`opts' emonly"'

    _estimates hold `full', copy restore
    di as txt "fitting the model without transitions..."
    capture quietly hmbivpois (`eq1') (`eq2') `ifin', `opts' notransition nolog
    local rca = _rc
    if `rca' == 1 exit 1
    scalar `lla' = .
    if `rca' == 0 {
        scalar `lla' = e(ll)
        matrix `ba' = e(b)
    }
    capture quietly hmbivpois (`eq1') (`eq2') `ifin', `opts' notransition nolog ///
        from(`bfrom')
    local rcb = _rc
    if `rcb' == 1 exit 1
    scalar `llb' = .
    if `rcb' == 0 scalar `llb' = e(ll)
    if `rca' & `rcb' {
        di as err "estat lrtrans: the model without transitions could not be fitted"
        exit 498
    }
    // the first fit is better: its estimates are made current again
    if `rca' == 0 & (`rcb' != 0 | `lla' > `llb' + 1e-8) {
        quietly hmbivpois (`eq1') (`eq2') `ifin', `opts0' notransition nolog ///
            from(`ba') emiterate(0) emonly sortby(none)
    }
    scalar `ll0' = e(ll)
    if e(N) != `N1' {
        di as err "estat lrtrans: the model without transitions has a different sample"
        exit 498
    }

    // weights of the chi-bar-squared distribution, at the restricted estimates
    PostRun, task(lrtrans) out(`W') reps(`reps') seed(`seed')
    local cn
    forvalues j = 0/`nr' {
        local cn `cn' chi2_`j'
    }
    matrix colnames `W' = `cn'
    matrix rownames `W' = weight

    tempname lr p pc
    scalar `lr' = 2*(`ll1' - `ll0')
    local neg = (`lr' < 0)
    if `neg' scalar `lr' = 0
    scalar `p' = 0
    if `lr' <= 0 scalar `p' = 1
    else {
        forvalues j = 1/`nr' {
            scalar `p' = `p' + `W'[1,`j'+1]*chi2tail(`j', `lr')
        }
    }
    scalar `pc' = chi2tail(`nr', `lr')

    di _n as txt "Likelihood-ratio test of no transitions between the states"
    di as txt "H0: the " as res `nr' as txt " transition probabilities off the diagonal are zero (latent class model)"
    di _n as txt "  log likelihood, model with transitions    = " as res %14.4f `ll1'
    di as txt "  log likelihood, model without transitions = " as res %14.4f `ll0'
    di as txt "  LR statistic                              = " as res %14.4f `lr'
    di _n as txt "  p-value, chi-bar-squared distribution     = " as res %10.4f `p'
    di as txt "  p-value, chi2(" as res `nr' as txt ") (conservative)" _col(45) "= " as res %10.4f `pc'
    di _n as txt "Weights of the chi-bar-squared distribution, mixture of chi2(j) (" ///
        as res `reps' as txt " simulations)"
    matlist `W', format(%7.4f) border(rows) names(columns)
    if `neg' {
        di as txt "note: the log likelihood of the model with transitions is below that of the"
        di as txt "      model without transitions; the model with transitions has not reached"
        di as txt "      its maximum, and LR is set to zero"
    }

    return scalar lr     = `lr'
    return scalar p      = `p'
    return scalar p_chi2 = `pc'
    return scalar df     = `nr'
    return scalar ll_1   = `ll1'
    return scalar ll_0   = `ll0'
    return scalar reps   = `reps'
    return matrix weights = `W'
end

// estat ic [, n(#)]: log likelihood, number of parameters, AIC and BIC with
// n = number of individuals (n() to change it), and the relative entropy of
// the classification of the observations into the states
program define EstatIC, rclass
    version 16
    syntax [, N(real 0) ]
    local nn = e(N_g)
    if `n' > 0 local nn = `n'
    tempname E S
    PostRun, task(ic) out(`E')
    local df = e(npar)
    local aic = -2*e(ll) + 2*`df'
    local bic = -2*e(ll) + `df'*ln(`nn')
    matrix `S' = (`nn', ., e(ll), `df', `aic', `bic')
    matrix colnames `S' = N ll0 ll df AIC BIC
    matrix rownames `S' = .

    di _n as txt "Akaike's and Bayesian information criteria"
    di as txt "{hline 13}{c TT}{hline 63}"
    di as txt "       Model {c |}  N indiv.    N obs.      ll(model)    df          AIC          BIC"
    di as txt "{hline 13}{c +}{hline 63}"
    di as txt "           . {c |}" as res %10.0fc `nn' %10.0fc e(N) "  " %13.3f e(ll) %6.0f `df' ///
        "  " %11.3f `aic' "  " %11.3f `bic'
    di as txt "{hline 13}{c BT}{hline 63}"
    di as txt "Note: BIC uses N = number of individuals" _c
    if `n' > 0 di as txt " (n() specified)."
    else       di as txt "."
    di as txt "Relative entropy of the classification (1 = states known with certainty): " ///
        as res %6.4f `E'

    return scalar entropy = `E'
    return scalar N_obs   = e(N)
    return matrix S       = `S'
end

// estat transition [, period(#) level(#)]: transition probabilities and mean
// initial probabilities, with delta-method standard errors; with
// heterogeneous transitions, the mean over the transition rows (of one
// period of arrival with period())
program define EstatTrans, rclass
    version 16
    syntax [, PERiod(string) Level(cilevel) ]
    local k = e(k)
    if `k' < 2 {
        di as err "estat transition requires a model with two or more states"
        exit 198
    }
    if "`e(transitions)'" == "none" {
        di as err "estat transition: the model has no transitions (notransition)"
        exit 198
    }
    local het = ("`e(transitions)'" == "heterogeneous")
    if "`period'" != "" {
        capture confirm number `period'
        if _rc {
            di as err "period() must be a number"
            exit 198
        }
        if !`het' {
            di as txt "note: period() ignored with homogeneous transitions"
            local period
        }
    }
    tempname R P SP P0 SP0
    PostRun, task(trans) out(`R') per(`period')
    local kk = `k'*`k'
    matrix `P'   = J(`k', `k', .)
    matrix `SP'  = J(`k', `k', .)
    local j0 = `kk' + 1
    local j1 = `kk' + `k'
    matrix `P0'  = `R'[1, `j0'..`j1']
    matrix `SP0' = `R'[2, `j0'..`j1']
    forvalues c = 1/`k' {
        forvalues d = 1/`k' {
            local j = (`c'-1)*`k' + `d'
            matrix `P'[`c', `d']  = `R'[1, `j']
            matrix `SP'[`c', `d'] = `R'[2, `j']
        }
    }
    local sn
    forvalues c = 1/`k' {
        local sn `sn' state`c'
    }
    matrix rownames `P'  = `sn'
    matrix colnames `P'  = `sn'
    matrix rownames `SP' = `sn'
    matrix colnames `SP' = `sn'
    matrix colnames `P0'  = `sn'
    matrix colnames `SP0' = `sn'
    matrix rownames `P0'  = prob
    matrix rownames `SP0' = se

    local z = invnormal(1 - (1 - `level'/100)/2)
    local ttl "Transition probabilities"
    if `het' & "`period'" == "" local ttl "Mean transition probabilities over the transition rows"
    if `het' & "`period'" != "" local ttl "Mean transition probabilities into period `period'"
    di _n as txt "`ttl' (delta-method standard errors;"
    di as txt "confidence intervals computed on the logit scale)"
    di as txt "{hline 13}{c TT}{hline 50}"
    di as txt "  from -> to {c |}      Prob.   Std. err.     [`level'% conf. interval]"
    di as txt "{hline 13}{c +}{hline 50}"
    forvalues c = 1/`k' {
        forvalues d = 1/`k' {
            local p  = `P'[`c', `d']
            local se = `SP'[`c', `d']
            // no interval for a probability at 0 or 1 (logit not defined)
            local lo = .
            local hi = .
            if `p' > 1e-8 & `p' < 1 - 1e-8 & `se' > 0 {
                local lo = invlogit(logit(`p') - `z'*`se'/(`p'*(1-`p')))
                local hi = invlogit(logit(`p') + `z'*`se'/(`p'*(1-`p')))
            }
            di as txt %12s "`c' -> `d'" " {c |}" as res %11.4f `p' %12.4f `se' ///
                %12.4f `lo' %11.4f `hi'
        }
        if `c' < `k' di as txt "{hline 13}{c +}{hline 50}"
    }
    di as txt "{hline 13}{c BT}{hline 50}"

    di _n as txt "Expected duration of a stay in each state, 1/(1 - P(c -> c)), in periods"
    di as txt "{hline 13}{c TT}{hline 24}"
    di as txt "       state {c |}   Duration   Std. err."
    di as txt "{hline 13}{c +}{hline 24}"
    tempname DU
    matrix `DU' = J(2, `k', .)
    forvalues c = 1/`k' {
        local p  = `P'[`c', `c']
        local du = 1/(1 - `p')
        local sd = `SP'[`c', `c']/(1 - `p')^2
        matrix `DU'[1, `c'] = `du'
        matrix `DU'[2, `c'] = `sd'
        di as txt %12s "`c'" " {c |}" as res %11.3f `du' %12.3f `sd'
    }
    di as txt "{hline 13}{c BT}{hline 24}"
    matrix colnames `DU' = `sn'
    matrix rownames `DU' = duration se

    di _n as txt "Initial probabilities, mean over the individuals"
    di as txt "{hline 13}{c TT}{hline 24}"
    di as txt "       state {c |}      Prob.   Std. err."
    di as txt "{hline 13}{c +}{hline 24}"
    forvalues c = 1/`k' {
        di as txt %12s "`c'" " {c |}" as res %11.4f `P0'[1, `c'] %12.4f `SP0'[1, `c']
    }
    di as txt "{hline 13}{c BT}{hline 24}"

    return matrix P      = `P'
    return matrix se_P   = `SP'
    return matrix pi0    = `P0'
    return matrix se_pi0 = `SP0'
    return matrix duration = `DU'
end

// estat stability [, horizon(#) limit(#) nsim(#) seed(#)]: diagnostics of the
// explosive behaviour of the lagged counts (section E6 of the document):
// coefficients of the lags by state; with log-type transformations the
// spectral radius of their 2 x 2 matrix, with linear lags the index of the
// univariate map y -> exp(a + g y); expected counts iterated for horizon()
// periods in each state from every observation; simulated paths of the
// model over horizon() periods from the last observed period
program define EstatStab, rclass
    version 16
    syntax [, HORizon(integer 50) LIMit(real 1000) NSIM(integer 1)        ///
        SEED(integer 20261116) ]
    if `horizon' < 1 | `nsim' < 1 | `limit' <= 0 {
        di as err "horizon() and nsim() must be positive integers, limit() a positive number"
        exit 198
    }
    local k = e(k)
    local y1 : word 1 of `e(depvar)'
    local y2 : word 2 of `e(depvar)'
    local lagt "`e(lagtransform)'"
    local lagt0 = regexr("`lagt'", "\(.*\)$", "")
    tempname R
    local stub = substr("`R'", 1, 25)
    PostRun, task(stability) out(`stub') horizon(`horizon') limit(`limit') ///
        nsim(`nsim') seed(`seed')
    tempname G RR U D S Q
    foreach m in G R U D S Q {
        local mm = cond("`m'" == "R", "RR", "`m'")
        matrix ``mm'' = `stub'`m'
        matrix drop `stub'`m'
    }
    local sn
    forvalues c = 1/`k' {
        local sn `sn' state`c'
    }
    matrix rownames `G'  = `sn'
    matrix colnames `G'  = `y1':lag_`y1' `y1':lag_`y2' `y2':lag_`y1' `y2':lag_`y2'
    matrix rownames `RR' = `sn'
    matrix colnames `RR' = radius norm
    matrix rownames `U'  = `sn'
    matrix colnames `U'  = `y1' `y2'
    matrix rownames `D'  = `sn'
    matrix colnames `D'  = share
    matrix colnames `S'  = share paths horizon limit
    matrix rownames `Q'  = sim_`y1' sim_`y2' fit_`y1' fit_`y2'
    matrix colnames `Q'  = p50 p90 p99 max

    di _n as txt "Stability of the dynamics of the counts (lagged outcomes: " as res "`lagt'" as txt ")"
    di _n as txt "1. Coefficients of the transformed lags, by state"
    matlist `G', format(%9.4f) border(rows) showcoleq(combined)
    if inlist("`lagt0'", "log1p", "logc", "logzero") {
        di _n as txt "   Spectral radius and norm of the 2 x 2 matrix of these coefficients. In logs"
        di as txt "   the dynamics is close to linear for large counts: a radius below 1 rules out"
        di as txt "   explosive paths, a norm below 1 is a sufficient condition"
        matlist `RR', format(%9.4f) border(rows)
    }
    if "`lagt0'" == "linear" {
        di _n as txt "   Index g_jj exp(a), with a the rest of the linear predictor at its 99th"
        di as txt "   percentile in the state: the map y -> exp(a + g_jj y) has a stable fixed point"
        di as txt "   only if the index is at most exp(-1) = 0.3679 (cross-lags ignored)"
        matlist `U', format(%9.4f) border(rows)
    }
    if inlist("`lagt0'", "indicator", "topcode") {
        di _n as txt "   The transformed lags are bounded: the expected counts cannot explode"
    }
    di _n as txt "2. Deterministic paths: from each observation, with its covariates, the counts"
    di as txt "   are replaced by their expected values for " as res `horizon' as txt " periods in each state."
    di as txt "   Share of paths whose expected count exceeds " as res `limit'
    matlist `D', format(%9.4f) border(rows)
    di _n as txt "3. Simulated paths of the model: " as res `S'[1,2] as txt " paths of " as res `horizon' as txt " periods, from"
    di as txt "   the last observed period of each individual (its covariates; state drawn"
    di as txt "   from its posterior probabilities)"
    di as txt "   Share of paths whose expected count exceeds " as res `limit' as txt ": " as res %6.4f `S'[1,1]
    di as txt "   Expected counts at the horizon (paths below the limit) and fitted in the sample"
    matlist `Q', format(%10.3f) border(rows)

    return matrix G      = `G'
    return matrix radius = `RR'
    return matrix index  = `U'
    return matrix det    = `D'
    return matrix sim    = `S'
    return matrix q      = `Q'
end

// estat ame [, prior onestep filtered instate(#) level(#)]: average marginal
// effects on the expected counts E(y1) and E(y2), holding the probabilities
// of the states fixed: dy/dx for continuous regressors, discrete change from
// the base level for the levels of factor variables, effect of one more count
// of each lagged outcome. The states are weighted with their posterior
// probabilities (default), or filtered, one-step-ahead or prior ones;
// instate(#) gives the effects in one state. Delta-method standard errors.
program define EstatAME, rclass
    version 16
    syntax [, PRIor ONEstep FILTered INSTate(integer 0) Level(cilevel) ]
    local k = e(k)
    local nw = ("`prior'" != "") + ("`onestep'" != "") + ("`filtered'" != "")
    if `nw' > 1 {
        di as err "only one of prior, onestep and filtered may be specified"
        exit 198
    }
    local ptype = 0
    if "`filtered'" != "" local ptype = 1
    if "`onestep'"  != "" local ptype = 2
    if "`prior'"    != "" local ptype = 3
    if `instate' != 0 {
        if `instate' < 1 | `instate' > `k' {
            di as err "instate() must be an integer between 1 and `k'"
            exit 198
        }
        if `nw' > 0 {
            di as err "instate() may not be combined with prior, onestep or filtered"
            exit 198
        }
    }
    local y1 : word 1 of `e(depvar)'
    local y2 : word 2 of `e(depvar)'

    tempname R b V
    local stub = substr("`R'", 1, 25)
    PostRun, task(ame) out(`stub') ptype(`ptype') instate(`instate')
    matrix `b' = `stub'
    matrix `V' = `stub'V
    matrix drop `stub' `stub'V
    // only the effects of the variables that enter each equation
    local names
    local keep
    forvalues j = 1/2 {
        forvalues e = 1/`ame_n' {
            if `ame_f`j'_`e'' {
                local names `names' `y`j'':`ame_name`e''
                local keep `keep' `=(`j'-1)*`ame_n' + `e''
            }
        }
    }
    mata: hmbp_msel("`b'", "`V'", "`keep'")
    matrix colnames `b' = `names'
    matrix colnames `V' = `names'
    matrix rownames `V' = `names'

    local wlab "posterior (given all the outcomes of the individual)"
    if `ptype' == 1 local wlab "filtered (given the outcomes up to the period)"
    if `ptype' == 2 local wlab "one-step-ahead (given the outcomes up to the period before)"
    if `ptype' == 3 local wlab "prior (given no outcome)"
    if `instate' > 0 local wlab "effects in state `instate'"
    di _n as txt "Average marginal effects on the expected counts, with the probabilities of"
    di as txt "the states held fixed. States: " as res "`wlab'"
    di as txt "dy/dx for continuous variables; discrete change from the base level for factor"
    local nl1 = `ame_n' - 1
    di as txt "levels; effect of one more count of the lagged outcome for " as res "`ame_name`nl1'' `ame_name`ame_n''"
    local z = invnormal(1 - (1 - `level'/100)/2)
    di as txt "{hline 13}{c TT}{hline 64}"
    di as txt "             {c |}            Delta-method"
    di as txt "             {c |}     Effect   std. err.      z    P>|z|     [`level'% conf. interval]"
    di as txt "{hline 13}{c +}{hline 64}"
    local c = 0
    forvalues j = 1/2 {
        di as res %-12s abbrev("E(`y`j'')", 12) as txt " {c |}"
        forvalues e = 1/`ame_n' {
            if `ame_f`j'_`e'' {
                local ++c
                local ef = `b'[1, `c']
                local se = sqrt(`V'[`c', `c'])
                local zz = `ef'/`se'
                di as txt %12s abbrev("`ame_name`e''", 12) " {c |}" as res "  " %9.0g `ef' "  " %9.0g `se' ///
                    %9.2f `zz' %8.3f 2*normal(-abs(`zz')) "    " %9.0g `ef' - `z'*`se' "   " %9.0g `ef' + `z'*`se'
            }
        }
        if `j' == 1 di as txt "{hline 13}{c +}{hline 64}"
    }
    di as txt "{hline 13}{c BT}{hline 64}"

    return matrix b = `b'
    return matrix V = `V'
end

// ===========================================================================
//                                   MATA
// ===========================================================================
capture mata: mata drop hmbp_*()
version 16
mata:

struct hmbp_dat {
    real scalar     n, N, k, p1, p2, s1, s2, nl, q, Tmax, m, nD, ntau, npar, het, qw, chk, dep, notr
    real colvector  y1, y2, ey1, ey2, lnfy, lfy1, lfy2, Ts, ns, trow, i00, i11
    real matrix     X1, X2, X1S, X2S, Z, info, W, Y1, Y2
    pointer(real colvector) colvector R, RT
}

// Famoye form: D, U and G = 1 + lambda U1 U2. Copulas: G is the probability of
// the cell, (L, H] the interval of each margin, EL and EH the derivatives of its
// endpoints with respect to ln(mu), Sg the sign that multiplies the parameter;
// normal copula: ZL and ZH are the standard normal quantiles of L and H
struct hmbp_mm {
    real matrix     lf, lf0, Mu1, Mu2, D1, D2, U1, U2, G
    real matrix     L1, H1, L2, H2, EL1, EH1, EL2, EH2, Sg
    real matrix     ZL1, ZH1, ZL2, ZH2
}

// ------------------------------------------------------------ parameters
real matrix hmbp_P(real rowvector tau, real scalar k)
{
    real matrix    P
    real rowvector e
    real scalar    c, cb, j

    P = J(k, k, 0)
    j = 0
    for (c=1; c<=k; c++) {
        e = J(1, k, 0)
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                j++
                e[cb] = tau[j]
            }
        }
        e = exp(e :- max(e))
        P[c,.] = e :/ sum(e)
    }
    return(P)
}

real rowvector hmbp_tau(real matrix P)
{
    real rowvector tau
    real scalar    k, c, cb, j

    k = rows(P)
    tau = J(1, k*(k-1), .)
    j = 0
    for (c=1; c<=k; c++) {
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                j++
                tau[j] = ln(P[c,cb]) - ln(P[c,c])
            }
        }
    }
    return(tau)
}

// Layout of the parameter vector b:
//   equation 1: common slopes (p1), state intercepts (k), state-specific slopes
//               (s1 for state 1, s1 for state 2, ...)
//   equation 2: the same, with p2 and s2
//   lambda (nl = 1 element, or nl = k: one per state); the measurement block
//               ends here, at position m
//   initial probabilities: q coefficients for each of the states 2..k
//   transitions: qw coefficients for each pair (from, to), from != to
// S1 and S2 are returned as s1 x k and s2 x k matrices (one column per state);
// lam is returned as a 1 x k row vector (a common lambda is repeated).
real rowvector hmbp_lam(struct hmbp_dat scalar M, real rowvector b)
{
    real rowvector lam

    lam = b[|M.m-M.nl+1 \ M.m|]
    if (M.nl == 1 & M.k > 1) lam = J(1, M.k, lam)
    return(lam)
}

void hmbp_mparts(struct hmbp_dat scalar M, real rowvector bm,
    real rowvector b1, real rowvector a1, real matrix S1,
    real rowvector b2, real rowvector a2, real matrix S2, real rowvector lam)
{
    real scalar o

    o = 0
    if (M.p1 > 0) {
        b1 = bm[|o+1 \ o+M.p1|]
        o  = o + M.p1
    }
    else b1 = J(1, 0, 0)
    a1 = bm[|o+1 \ o+M.k|]
    o  = o + M.k
    if (M.s1 > 0) {
        S1 = colshape(bm[|o+1 \ o+M.k*M.s1|], M.s1)'
        o  = o + M.k*M.s1
    }
    else S1 = J(0, M.k, 0)
    if (M.p2 > 0) {
        b2 = bm[|o+1 \ o+M.p2|]
        o  = o + M.p2
    }
    else b2 = J(1, 0, 0)
    a2 = bm[|o+1 \ o+M.k|]
    o  = o + M.k
    if (M.s2 > 0) {
        S2 = colshape(bm[|o+1 \ o+M.k*M.s2|], M.s2)'
        o  = o + M.k*M.s2
    }
    else S2 = J(0, M.k, 0)
    lam = hmbp_lam(M, bm)
}

void hmbp_unpack(struct hmbp_dat scalar M, real rowvector b,
    real rowvector b1, real rowvector a1, real matrix S1,
    real rowvector b2, real rowvector a2, real matrix S2,
    real rowvector lam, real matrix D, real rowvector tau)
{
    real scalar o

    hmbp_mparts(M, b, b1, a1, S1, b2, a2, S2, lam)
    o = M.m
    if (M.k > 1) {
        D = colshape(b[|o+1 \ o+M.nD|], M.q)'       // q x (k-1)
        o = o + M.nD
        if (M.ntau > 0) tau = b[|o+1 \ o+M.ntau|]
        else            tau = J(1, 0, .)
    }
    else {
        D   = J(M.q, 0, .)
        tau = J(1, 0, .)
    }
}

real matrix hmbp_pi0(real matrix Z, real matrix D, real scalar k)
{
    real matrix E

    if (k == 1) return(J(rows(Z), 1, 1))
    E = (J(rows(Z), 1, 0), Z*D)
    E = exp(E :- rowmax(E))
    return(E :/ rowsum(E))
}

// ------------------------------------------------------- measurement model
// means and everything else that does not depend on lambda
struct hmbp_mm scalar hmbp_means(struct hmbp_dat scalar M,
    real rowvector b1, real rowvector a1, real matrix S1,
    real rowvector b2, real rowvector a2, real matrix S2)
{
    struct hmbp_mm scalar mm
    real matrix  eta1, eta2, L, H, EL, EH, R1, R2
    real scalar  q

    q = 1 - exp(-1)
    eta1 = J(M.N, 1, a1)
    if (M.p1 > 0) eta1 = eta1 :+ (M.X1*b1')
    if (M.s1 > 0) eta1 = eta1 + M.X1S*S1
    eta2 = J(M.N, 1, a2)
    if (M.p2 > 0) eta2 = eta2 :+ (M.X2*b2')
    if (M.s2 > 0) eta2 = eta2 + M.X2S*S2
    mm.Mu1 = exp(eta1)
    mm.Mu2 = exp(eta2)
    if (M.dep == 0) {
        mm.D1  = exp(-q*mm.Mu1)
        mm.D2  = exp(-q*mm.Mu2)
        mm.U1  = M.ey1 :- mm.D1
        mm.U2  = M.ey2 :- mm.D2
        mm.lf0 = M.y1:*eta1 - mm.Mu1 + M.y2:*eta2 - mm.Mu2
    }
    else {
        hmbp_margin(mm.Mu1, M.Y1, M.lfy1, L, H, EL, EH, R1)
        mm.L1  = L
        mm.H1  = H
        mm.EL1 = EL
        mm.EH1 = EH
        hmbp_margin(mm.Mu2, M.Y2, M.lfy2, L, H, EL, EH, R2)
        mm.L2  = L
        mm.H2  = H
        mm.EL2 = EL
        mm.EH2 = EH
        // reflecting one margin only changes the sign of the parameter
        mm.Sg  = 1 :- 2*(R1 :!= R2)
        mm.lf0 = mm.H1 + mm.H2          // used only to detect missing values
        if (M.dep == 1) {
            mm.ZL1 = hmbp_zof(mm.L1)
            mm.ZH1 = hmbp_zof(mm.H1)
            mm.ZL2 = hmbp_zof(mm.L2)
            mm.ZH2 = hmbp_zof(mm.H2)
        }
    }
    return(mm)
}

// ------------------------------------------------------- copula dependence
// exp(x) - 1 and ln(1 + x), accurate for small |x|
real matrix hmbp_expm1(real matrix X)
{
    real matrix S

    S = abs(X) :< 1e-4
    return(S:*(X + X:^2/2 + X:^3/6 + X:^4/24) + (1:-S):*(exp(X) :- 1))
}

real matrix hmbp_ln1p(real matrix X)
{
    real matrix S

    S = abs(X) :< 1e-4
    return(S:*(X - X:^2/2 + X:^3/3 - X:^4/4) + (1:-S):*ln(1 :+ X))
}

// One margin of the copula models. The count y corresponds to an interval
// (L, H] of a uniform variable: in lower-tail form when y < mu, with L = F(y-1)
// and H = F(y); in upper-tail form when y >= mu, where the uniform is reflected,
// with L = P(Y > y) and H = P(Y >= y). In this way a cell with a small
// probability never has its endpoints close to 1, and its probability is
// computed without cancellation. EL and EH are the derivatives of L and H with
// respect to ln(mu); R marks the cells in upper-tail form.
void hmbp_margin(real matrix Mu, real matrix Y, real colvector lfy,
    real matrix L, real matrix H, real matrix EL, real matrix EH,
    real matrix R)
{
    real matrix    f, fm
    real colvector mu, y, r, l, idx
    real scalar    k

    k  = cols(Mu)
    f  = exp(Y:*ln(Mu) - Mu :- lfy)             // pmf at y
    fm = f:*Y:/Mu                                // pmf at y-1 (zero when y = 0)
    R  = (Y :>= Mu)
    mu = vec(Mu)
    y  = vec(Y)
    r  = vec(R)
    l  = J(rows(mu), 1, 0)
    idx = selectindex(r)
    if (rows(idx) > 0) l[idx] = poissontail(mu[idx], y[idx] :+ 1)
    idx = selectindex((r :== 0) :& (y :> 0))
    if (rows(idx) > 0) l[idx] = poisson(mu[idx], y[idx] :- 1)
    L  = rowshape(l, k)'
    H  = L + f
    EL = (R:*f - (1:-R):*fm):*Mu
    EH = (R:*fm - (1:-R):*f):*Mu
}

// Frank copula C(u,v) = -ln(1 + A B / Dn) / T, with A = exp(-T u) - 1,
// B = exp(-T v) - 1, Dn = exp(-T) - 1
real matrix hmbp_fC(real matrix A, real matrix B, real matrix Dn, real matrix T)
{
    return(-hmbp_ln1p(A:*B:/Dn):/T)
}

// the same for T close to zero: C = u v [1 + T (1-u)(1-v)/2]
real matrix hmbp_fCs(real matrix U, real matrix V, real matrix T)
{
    return(U:*V:*(1 :+ T:*(1:-U):*(1:-V)/2))
}

// derivatives of C with respect to u, v and T at the corner (u,v)
void hmbp_fd(real matrix U, real matrix A, real matrix V, real matrix B,
    real matrix Dn, real matrix T, real matrix C1, real matrix C2,
    real matrix Ct)
{
    real matrix g, den, gT

    g   = A:*B:/Dn
    den = Dn + A:*B
    C1  = (1:+A):*B:/den
    C2  = (1:+B):*A:/den
    gT  = (-U:*(1:+A):*B - V:*(1:+B):*A):/Dn + A:*B:*(1:+Dn):/Dn:^2
    Ct  = hmbp_ln1p(g):/T:^2 - gT:/(T:*(1:+g))
}

// the same for T close to zero
void hmbp_fds(real matrix U, real matrix V, real matrix T, real matrix C1,
    real matrix C2, real matrix Ct)
{
    C1 = V:*(1 :+ T:*(1:-V):*(1:-2*U)/2)
    C2 = U:*(1 :+ T:*(1:-U):*(1:-2*V)/2)
    Ct = U:*V:*(1:-U):*(1:-V)/2
}

// Frank copula: probability of the rectangle (L1,H1] x (L2,H2] in mm.G and its
// log in mm.lf. lam is the parameter by state (1 x k); T = lam x sign of the cell
void hmbp_frank(struct hmbp_mm scalar mm, real rowvector lam)
{
    real matrix    T, Ts, Ah, Al, Bh, Bl, Dn
    real rowvector sm

    T  = lam:*mm.Sg
    sm = abs(lam) :< 1e-6
    Ts = T :+ sm
    Ah = hmbp_expm1(-Ts:*mm.H1)
    Al = hmbp_expm1(-Ts:*mm.L1)
    Bh = hmbp_expm1(-Ts:*mm.H2)
    Bl = hmbp_expm1(-Ts:*mm.L2)
    Dn = hmbp_expm1(-Ts)
    mm.G = hmbp_fC(Ah, Bh, Dn, Ts) - hmbp_fC(Al, Bh, Dn, Ts) - hmbp_fC(Ah, Bl, Dn, Ts) + hmbp_fC(Al, Bl, Dn, Ts)
    if (any(sm)) {
        mm.G = (1:-sm):*mm.G + sm:*(hmbp_fCs(mm.H1, mm.H2, T) - hmbp_fCs(mm.L1, mm.H2, T) - hmbp_fCs(mm.H1, mm.L2, T) + hmbp_fCs(mm.L1, mm.L2, T))
    }
    // floor, so that the log is defined when the probability underflows
    mm.G  = mm.G:*(mm.G :> 1e-300) + 1e-300*(mm.G :<= 1e-300)
    mm.lf = ln(mm.G)
}

// Frank copula: scores of ln p with respect to eta1, eta2 and the parameter
void hmbp_frankd(struct hmbp_mm scalar mm, real rowvector lam, real matrix S1,
    real matrix S2, real matrix SL)
{
    real matrix    T, Ts, Ah, Al, Bh, Bl, Dn, ok
    real matrix    a1, a2, at, b1, b2, bt, c1, c2, ct, e1, e2, et
    real rowvector sm

    T  = lam:*mm.Sg
    sm = abs(lam) :< 1e-6
    Ts = T :+ sm
    Ah = hmbp_expm1(-Ts:*mm.H1)
    Al = hmbp_expm1(-Ts:*mm.L1)
    Bh = hmbp_expm1(-Ts:*mm.H2)
    Bl = hmbp_expm1(-Ts:*mm.L2)
    Dn = hmbp_expm1(-Ts)
    hmbp_fd(mm.H1, Ah, mm.H2, Bh, Dn, Ts, a1, a2, at)
    hmbp_fd(mm.L1, Al, mm.H2, Bh, Dn, Ts, b1, b2, bt)
    hmbp_fd(mm.H1, Ah, mm.L2, Bl, Dn, Ts, c1, c2, ct)
    hmbp_fd(mm.L1, Al, mm.L2, Bl, Dn, Ts, e1, e2, et)
    S1 = (a1 - c1):*mm.EH1 - (b1 - e1):*mm.EL1
    S2 = (a2 - b2):*mm.EH2 - (c2 - e2):*mm.EL2
    SL = mm.Sg:*(at - bt - ct + et)
    if (any(sm)) {
        hmbp_fds(mm.H1, mm.H2, T, a1, a2, at)
        hmbp_fds(mm.L1, mm.H2, T, b1, b2, bt)
        hmbp_fds(mm.H1, mm.L2, T, c1, c2, ct)
        hmbp_fds(mm.L1, mm.L2, T, e1, e2, et)
        S1 = (1:-sm):*S1 + sm:*((a1 - c1):*mm.EH1 - (b1 - e1):*mm.EL1)
        S2 = (1:-sm):*S2 + sm:*((a2 - b2):*mm.EH2 - (c2 - e2):*mm.EL2)
        SL = (1:-sm):*SL + sm:*(mm.Sg:*(at - bt - ct + et))
    }
    // a cell whose probability was floored carries no information
    ok = (mm.G :> 1e-300)
    S1 = ok:*S1:/mm.G
    S2 = ok:*S2:/mm.G
    SL = ok:*SL:/mm.G
}

// --------------------------------------------------------- normal copula
// standard normal quantile of a probability; -40 stands for minus infinity
// (U = 0), and a probability that rounds to 1 is taken just below 1
real matrix hmbp_zof(real matrix U)
{
    real matrix Us, Z

    Us = U + 0.5*(U :<= 0)
    Us = Us:*(Us :< 1) + (1 - 1.2e-16)*(Us :>= 1)
    Z  = (U :> 0):*invnormal(Us) - 40*(U :<= 0)
    return(Z:*(Z :> -40) - 40*(Z :<= -40))
}

// Phi(b) - Phi(a) for a < b, computed in the tail where the difference is accurate
real matrix hmbp_dnorm(real matrix a, real matrix b)
{
    real matrix neg

    neg = (a :> 0)
    return(normal(neg:*(-a) + (1:-neg):*b) - normal(neg:*(-b) + (1:-neg):*a))
}

// bivariate normal density with correlation r, s2 = 1 - r^2
real matrix hmbp_phi2(real matrix a, real matrix b, real matrix r,
    real matrix s2)
{
    return(exp(-(a:^2 - 2*r:*a:*b + b:^2):/(2*s2)):/(2*pi()*sqrt(s2)))
}

// nodes (column 1) and weights (column 2) of the 12-point Gauss-Legendre rule
real matrix hmbp_gl12()
{
    return((-0.98156063424671924, 0.047175336386511411 \
            -0.9041172563704748,  0.10693932599531907  \
            -0.76990267419430469, 0.16007832854334642  \
            -0.58731795428661748, 0.20316742672306573  \
            -0.36783149899818018, 0.23349253653835461  \
            -0.12523340851146891, 0.24914704581340269  \
             0.12523340851146891, 0.24914704581340269  \
             0.36783149899818018, 0.23349253653835461  \
             0.58731795428661748, 0.20316742672306573  \
             0.76990267419430469, 0.16007832854334642  \
             0.9041172563704748,  0.10693932599531907  \
             0.98156063424671924, 0.047175336386511411))
}

// Normal copula: probability of the rectangle (ZL1,ZH1] x (ZL2,ZH2] of a
// standard bivariate normal with correlation r = tanh(lam) x sign of the cell,
// in mm.G, and its log in mm.lf.
// Rows with y1 = y2 = 0 (both intervals start at minus infinity): binormal().
// Other rows: integral over the margin with the shorter interval (z) of
// phi(z) [Phi((zh - r z)/s) - Phi((zl - r z)/s)], where (zl, zh] is the
// interval of the other margin. The integrand is positive and smooth, so the
// probability is accurate in relative terms also when it is very small.
void hmbp_normalp(struct hmbp_dat scalar M, struct hmbp_mm scalar mm,
    real rowvector lam)
{
    real matrix R, p, zl1, zh1, zl2, zh2, r, s, o2, zlo, zho, zli, zhi
    real matrix c, h, z, q, gl
    real scalar m

    R = tanh(lam):*mm.Sg
    p = J(rows(R), cols(R), .)
    if (rows(M.i00) > 0) {
        p[M.i00, .] = binormal(mm.ZH1[M.i00, .], mm.ZH2[M.i00, .], R[M.i00, .])
    }
    if (rows(M.i11) > 0) {
        zl1 = mm.ZL1[M.i11, .]
        zh1 = mm.ZH1[M.i11, .]
        zl2 = mm.ZL2[M.i11, .]
        zh2 = mm.ZH2[M.i11, .]
        r   = R[M.i11, .]
        s   = sqrt(1 :- r:^2)
        o2  = (zh2 - zl2) :< (zh1 - zl1)
        zlo = o2:*zl2 + (1:-o2):*zl1
        zho = o2:*zh2 + (1:-o2):*zh1
        zli = o2:*zl1 + (1:-o2):*zl2
        zhi = o2:*zh1 + (1:-o2):*zh2
        c   = (zlo + zho)/2
        h   = (zho - zlo)/2
        q   = J(rows(c), cols(c), 0)
        gl  = hmbp_gl12()
        for (m=1; m<=12; m++) {
            z = c + gl[m,1]*h
            q = q + gl[m,2]*(normalden(z):*hmbp_dnorm((zli - r:*z):/s, (zhi - r:*z):/s))
        }
        p[M.i11, .] = h:*q
    }
    // floor, so that the log is defined when the probability underflows
    mm.G  = p:*(p :> 1e-300) + 1e-300*(p :<= 1e-300)
    mm.lf = ln(mm.G)
}

// Normal copula: scores of ln p with respect to eta1, eta2 and athrho
void hmbp_normald(struct hmbp_mm scalar mm, real rowvector lam, real matrix S1,
    real matrix S2, real matrix SL)
{
    real matrix    R, s2, s, ok
    real rowvector rho

    rho = tanh(lam)
    R   = rho:*mm.Sg
    s2  = 1 :- R:^2
    s   = sqrt(s2)
    // derivative of the probability with respect to an endpoint of a margin:
    // the conditional probability of the interval of the other margin
    S1 = hmbp_dnorm((mm.ZL2 - R:*mm.ZH1):/s, (mm.ZH2 - R:*mm.ZH1):/s):*mm.EH1 - hmbp_dnorm((mm.ZL2 - R:*mm.ZL1):/s, (mm.ZH2 - R:*mm.ZL1):/s):*mm.EL1
    S2 = hmbp_dnorm((mm.ZL1 - R:*mm.ZH2):/s, (mm.ZH1 - R:*mm.ZH2):/s):*mm.EH2 - hmbp_dnorm((mm.ZL1 - R:*mm.ZL2):/s, (mm.ZH1 - R:*mm.ZL2):/s):*mm.EL2
    // with respect to the correlation: the bivariate density at the four corners
    SL = hmbp_phi2(mm.ZH1, mm.ZH2, R, s2) - hmbp_phi2(mm.ZL1, mm.ZH2, R, s2) - hmbp_phi2(mm.ZH1, mm.ZL2, R, s2) + hmbp_phi2(mm.ZL1, mm.ZL2, R, s2)
    SL = SL:*mm.Sg:*(1 :- rho:^2)
    // a cell whose probability was floored carries no information
    ok = (mm.G :> 1e-300)
    S1 = ok:*S1:/mm.G
    S2 = ok:*S2:/mm.G
    SL = ok:*SL:/mm.G
}

// dependence term and log density given the means; lam is 1 x k
void hmbp_dep(struct hmbp_dat scalar M, struct hmbp_mm scalar mm,
    real rowvector lam)
{
    if (M.dep == 0) {
        mm.G  = 1 :+ lam:*(mm.U1:*mm.U2)
        mm.lf = mm.lf0 + ln(mm.G) :- M.lnfy
    }
    else if (M.dep == 1) hmbp_normalp(M, mm, lam)
    else hmbp_frank(mm, lam)
}

struct hmbp_mm scalar hmbp_meas(struct hmbp_dat scalar M,
    real rowvector b1, real rowvector a1, real matrix S1,
    real rowvector b2, real rowvector a2, real matrix S2, real rowvector lam)
{
    struct hmbp_mm scalar mm

    mm = hmbp_means(M, b1, a1, S1, b2, a2, S2)
    hmbp_dep(M, mm, lam)
    return(mm)
}

// admissible interval of lambda in each state (pmf nonnegative for every y at
// the current means): row 1 = lower bounds, row 2 = upper bounds, one column
// per state
real matrix hmbp_lambnd(struct hmbp_mm scalar mm)
{
    real matrix    A, B, X
    real rowvector lo, hi

    A  = mm.D1:*mm.D2
    B  = (1:-mm.D1):*(1:-mm.D2)
    X  = A:*(A:>=B) + B:*(A:<B)
    lo = colmax(-1:/X)
    A  = (1:-mm.D1):*mm.D2
    B  = mm.D1:*(1:-mm.D2)
    X  = A:*(A:>=B) + B:*(A:<B)
    hi = colmin(1:/X)
    return(lo \ hi)
}

// the same interval for the lambda parameters (2 x nl): by state, or the
// intersection over the states when lambda is common
real matrix hmbp_lamint(struct hmbp_dat scalar M, struct hmbp_mm scalar mm)
{
    real matrix bnd

    bnd = hmbp_lambnd(mm)
    if (M.nl == 1) return(max(bnd[1,.]) \ min(bnd[2,.]))
    return(bnd)
}

// ------------------------------------------------ forward / backward
// individuals are sorted by decreasing T_i: those active at period t are 1..ns[t]
real colvector hmbp_fwd(struct hmbp_dat scalar M, real matrix lf,
    real matrix Pi0, real matrix P, real matrix A, real matrix F)
{
    real colvector m, ll, s, r
    real matrix    a
    real scalar    t, nt

    m = rowmax(lf)
    F = exp(lf :- m)
    A = J(M.N, M.k, 0)
    r = *M.R[1]
    a = Pi0 :* F[r,.]
    s = rowsum(a)
    a = a :/ s
    A[r,.] = a
    ll = ln(s) + m[r]
    for (t=2; t<=M.Tmax; t++) {
        nt = M.ns[t]
        r  = *M.R[t]
        a  = (a[|1,1 \ nt,M.k|]*P) :* F[r,.]
        s  = rowsum(a)
        a  = a :/ s
        A[r,.] = a
        ll[|1 \ nt|] = ll[|1 \ nt|] + ln(s) + m[r]
    }
    return(ll)
}

void hmbp_bwd(struct hmbp_dat scalar M, real matrix F, real matrix A,
    real matrix P, real matrix W, real matrix Zs, real scalar wantxi,
    real matrix XiI)
{
    real matrix    B, G, Ap, bb
    real colvector r, rp, nrm
    real scalar    t, nt, k, c

    k  = M.k
    W  = A
    Zs = J(k, k, 0)
    if (wantxi) XiI = J(M.n, k*k, 0)
    if (M.Tmax >= 2) {
        B = J(M.ns[M.Tmax], k, 1)
        for (t=M.Tmax; t>=2; t--) {
            nt  = M.ns[t]
            r   = *M.R[t]
            rp  = (*M.R[t-1])[|1 \ nt|]
            G   = F[r,.] :* B
            Ap  = A[rp,.]
            nrm = rowsum((Ap*P) :* G)
            Ap  = Ap :/ nrm
            Zs  = Zs + P :* cross(Ap, G)
            if (wantxi) {
                for (c=1; c<=k; c++) {
                    XiI[|1,(c-1)*k+1 \ nt,c*k|] = XiI[|1,(c-1)*k+1 \ nt,c*k|] + (Ap[.,c] :* G) :* P[c,.]
                }
            }
            bb = G*P'
            bb = bb :/ rowsum(bb)
            B  = J(M.ns[t-1], k, 1)
            B[|1,1 \ nt,k|] = bb
            rp = *M.R[t-1]
            W[rp,.] = A[rp,.] :* B
        }
    }
    W = W :/ rowsum(W)
}

// log likelihood (and posteriors if post!=0); missing if b is not admissible
// transition probabilities by observation: one row per transition row,
// columns (c-1)*k+1..c*k = P(U_t = . | U_t-1 = c, w_it)
real matrix hmbp_PT(struct hmbp_dat scalar M, real rowvector tau)
{
    real matrix PT, TT, G, E
    real scalar c, cb, j, k

    k  = M.k
    TT = colshape(tau, M.qw)
    PT = J(rows(M.W), k*k, .)
    j  = 0
    for (c=1; c<=k; c++) {
        G = J(M.qw, k, 0)
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                j++
                G[.,cb] = TT[j,.]'
            }
        }
        E = M.W*G
        E = exp(E :- rowmax(E))
        PT[|1,(c-1)*k+1 \ .,c*k|] = E :/ rowsum(E)
    }
    return(PT)
}

// forward recursion with observation-specific transition probabilities
real colvector hmbp_fwdh(struct hmbp_dat scalar M, real matrix lf,
    real matrix Pi0, real matrix PT, real matrix A, real matrix F)
{
    real colvector m, ll, s, r, rt
    real matrix    a, ap
    real scalar    t, nt, c, k

    k = M.k
    m = rowmax(lf)
    F = exp(lf :- m)
    A = J(M.N, k, 0)
    r = *M.R[1]
    a = Pi0 :* F[r,.]
    s = rowsum(a)
    a = a :/ s
    A[r,.] = a
    ll = ln(s) + m[r]
    for (t=2; t<=M.Tmax; t++) {
        nt = M.ns[t]
        r  = *M.R[t]
        rt = *M.RT[t]
        ap = a[|1,1 \ nt,k|]
        a  = J(nt, k, 0)
        for (c=1; c<=k; c++) {
            a = a + ap[.,c] :* PT[rt, ((c-1)*k+1)..(c*k)]
        }
        a  = a :* F[r,.]
        s  = rowsum(a)
        a  = a :/ s
        A[r,.] = a
        ll[|1 \ nt|] = ll[|1 \ nt|] + ln(s) + m[r]
    }
    return(ll)
}

// backward recursion with observation-specific transition probabilities;
// XI (one row per transition row, k*k columns) holds the joint posteriors
void hmbp_bwdh(struct hmbp_dat scalar M, real matrix F, real matrix A,
    real matrix PT, real matrix W, real matrix XI)
{
    real matrix    B, G, Ap, bb
    real colvector r, rt, rp, nrm
    real rowvector cc
    real scalar    t, nt, k, c

    k  = M.k
    W  = A
    XI = J(rows(M.W), k*k, 0)
    if (M.Tmax >= 2) {
        B = J(M.ns[M.Tmax], k, 1)
        for (t=M.Tmax; t>=2; t--) {
            nt  = M.ns[t]
            r   = *M.R[t]
            rt  = *M.RT[t]
            rp  = (*M.R[t-1])[|1 \ nt|]
            G   = F[r,.] :* B
            Ap  = A[rp,.]
            bb  = J(nt, k, 0)
            nrm = J(nt, 1, 0)
            for (c=1; c<=k; c++) {
                cc = ((c-1)*k+1)..(c*k)
                bb[.,c] = rowsum(PT[rt, cc] :* G)
                nrm = nrm + Ap[.,c] :* bb[.,c]
            }
            for (c=1; c<=k; c++) {
                cc = ((c-1)*k+1)..(c*k)
                XI[rt, cc] = ((Ap[.,c] :/ nrm) :* PT[rt, cc]) :* G
            }
            bb = bb :/ rowsum(bb)
            B  = J(M.ns[t-1], k, 1)
            B[|1,1 \ nt,k|] = bb
            rp = *M.R[t-1]
            W[rp,.] = A[rp,.] :* B
        }
    }
    W = W :/ rowsum(W)
}

// log likelihood (and posteriors if post!=0); missing if b is not admissible.
// Homogeneous chain: P is the k x k transition matrix, Zs the expected
// transition counts and XiI (if wantxi) their individual contributions.
// Heterogeneous chain: P holds the probabilities by transition row and XiI
// the joint posteriors by transition row (both rows(M.W) x k*k).
real scalar hmbp_estep(struct hmbp_dat scalar M, real rowvector b,
    real matrix W, real matrix Zs, real matrix Pi0, real matrix P,
    struct hmbp_mm scalar mm, real scalar post, real scalar wantxi,
    real matrix XiI)
{
    real rowvector b1, a1, b2, a2, tau, lam
    real matrix    D, A, F, S1, S2, bnd
    real colvector ll

    hmbp_unpack(M, b, b1, a1, S1, b2, a2, S2, lam, D, tau)
    mm  = hmbp_meas(M, b1, a1, S1, b2, a2, S2, lam)
    if (M.chk & M.dep == 0) {
        // lambda must be inside its admissible interval
        bnd = hmbp_lambnd(mm)
        if (any(lam :<= bnd[1,.]) | any(lam :>= bnd[2,.])) return(.)
    }
    if (hasmissing(mm.lf)) return(.)
    Pi0 = hmbp_pi0(M.Z, D, M.k)
    if (M.het) {
        P  = hmbp_PT(M, tau)
        if (hasmissing(P)) return(.)
        ll = hmbp_fwdh(M, mm.lf, Pi0, P, A, F)
        if (hasmissing(ll)) return(.)
        if (post) {
            hmbp_bwdh(M, F, A, P, W, XiI)
            Zs = colshape(colsum(XiI), M.k)
        }
    }
    else {
        // without transitions the transition matrix is the identity
        if (M.notr) P = I(M.k)
        else        P = hmbp_P(tau, M.k)
        ll = hmbp_fwd(M, mm.lf, Pi0, P, A, F)
        if (hasmissing(ll)) return(.)
        if (post) hmbp_bwd(M, F, A, P, W, Zs, wantxi, XiI)
    }
    return(quadsum(ll))
}

// ------------------------------------------------- scores and derivatives
// observation-level scores of log f with respect to eta1, eta2, lambda
void hmbp_mscore(struct hmbp_dat scalar M, struct hmbp_mm scalar mm,
    real rowvector lam, real matrix S1, real matrix S2, real matrix SL)
{
    real scalar q

    if (M.dep == 1) {
        hmbp_normald(mm, lam, S1, S2, SL)
        return
    }
    if (M.dep == 2) {
        hmbp_frankd(mm, lam, S1, S2, SL)
        return
    }
    q  = 1 - exp(-1)
    S1 = M.y1 :- mm.Mu1 + lam:*(mm.U2:*(q*mm.Mu1:*mm.D1)):/mm.G
    S2 = M.y2 :- mm.Mu2 + lam:*(mm.U1:*(q*mm.Mu2:*mm.D2)):/mm.G
    SL = (mm.U1:*mm.U2):/mm.G
}

// analytic score of the marginal log likelihood (Fisher identity)
real rowvector hmbp_grad(struct hmbp_dat scalar M, real rowvector b,
    real matrix W, real matrix Zs, real matrix Pi0, real matrix P,
    struct hmbp_mm scalar mm, real matrix XiI)
{
    real rowvector g, gt, cc, lam
    real matrix    S1, S2, SL, wS1, wS2, W1, Xc, Gm
    real scalar    c, cb, j

    lam = hmbp_lam(M, b)
    hmbp_mscore(M, mm, lam, S1, S2, SL)
    wS1 = W:*S1
    wS2 = W:*S2
    g = hmbp_mgrad(M, wS1, wS2, W:*SL)
    if (M.k > 1) {
        W1 = W[M.info[.,1], .]
        g  = (g, vec(cross(M.Z, W1[|1,2 \ .,.|] - Pi0[|1,2 \ .,.|]))')
        if (M.het) {
            gt = J(1, 0, .)
            for (c=1; c<=M.k; c++) {
                cc = ((c-1)*M.k+1)..(c*M.k)
                Xc = XiI[., cc]
                Gm = cross(M.W, Xc - P[., cc] :* rowsum(Xc))
                for (cb=1; cb<=M.k; cb++) {
                    if (cb != c) gt = (gt, Gm[.,cb]')
                }
            }
        }
        else {
            gt = J(1, M.ntau, .)
            j  = 0
            if (M.ntau > 0) {
                for (c=1; c<=M.k; c++) {
                    for (cb=1; cb<=M.k; cb++) {
                        if (cb != c) {
                            j++
                            gt[j] = Zs[c,cb] - P[c,cb]*sum(Zs[c,.])
                        }
                    }
                }
            }
        }
        g = (g, gt)
    }
    return(g)
}

real rowvector hmbp_gradonly(struct hmbp_dat scalar M, real rowvector b)
{
    struct hmbp_mm scalar mm
    real matrix W, Zs, Pi0, P, XiI
    real scalar ll

    ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 1, 0, XiI)
    if (ll >= .) return(J(1, cols(b), .))
    return(hmbp_grad(M, b, W, Zs, Pi0, P, mm, XiI))
}

// individual contributions to the score (n x npar), for vce(robust)
real matrix hmbp_scores(struct hmbp_dat scalar M, real rowvector b)
{
    struct hmbp_mm scalar mm
    real matrix  W, Zs, Pi0, P, XiI, S1, S2, SL, wS1, wS2, E, St, Sc, Xc, Rc, Ef
    real colvector wprev
    real rowvector cc, lam
    real scalar  ll, c, cb, j, k

    k  = M.k
    ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 1, 1, XiI)
    lam = hmbp_lam(M, b)
    hmbp_mscore(M, mm, lam, S1, S2, SL)
    wS1 = W:*S1
    wS2 = W:*S2
    Sc = J(M.n, 0, .)
    if (M.p1 > 0) Sc = (Sc, panelsum(rowsum(wS1):*M.X1, M.info))
    Sc = (Sc, panelsum(wS1, M.info))
    if (M.s1 > 0) {
        for (c=1; c<=k; c++) Sc = (Sc, panelsum(wS1[.,c]:*M.X1S, M.info))
    }
    if (M.p2 > 0) Sc = (Sc, panelsum(rowsum(wS2):*M.X2, M.info))
    Sc = (Sc, panelsum(wS2, M.info))
    if (M.s2 > 0) {
        for (c=1; c<=k; c++) Sc = (Sc, panelsum(wS2[.,c]:*M.X2S, M.info))
    }
    if (M.nl == 1) Sc = (Sc, panelsum(rowsum(W:*SL), M.info))
    else           Sc = (Sc, panelsum(W:*SL, M.info))
    if (k > 1) {
        E = W[M.info[.,1], .] - Pi0
        for (c=2; c<=k; c++) Sc = (Sc, M.Z :* E[.,c])
        if (M.het) {
            for (c=1; c<=k; c++) {
                cc = ((c-1)*k+1)..(c*k)
                Xc = XiI[., cc]
                Rc = Xc - P[., cc] :* rowsum(Xc)
                for (cb=1; cb<=k; cb++) {
                    if (cb != c) {
                        Ef = J(M.N, M.qw, 0)
                        Ef[M.trow, .] = M.W :* Rc[.,cb]
                        Sc = (Sc, panelsum(Ef, M.info))
                    }
                }
            }
        }
        else {
            St = J(M.n, M.ntau, .)
            j  = 0
            if (M.ntau > 0) {
                for (c=1; c<=k; c++) {
                    wprev = rowsum(XiI[|1,(c-1)*k+1 \ .,c*k|])
                    for (cb=1; cb<=k; cb++) {
                        if (cb != c) {
                            j++
                            St[.,j] = XiI[.,(c-1)*k+cb] - P[c,cb]*wprev
                        }
                    }
                }
            }
            Sc = (Sc, St)
        }
    }
    return(Sc)
}

// score of the measurement block, in the layout of the parameter vector
real rowvector hmbp_mgrad(struct hmbp_dat scalar M, real matrix wS1,
    real matrix wS2, real matrix wSL)
{
    real rowvector g

    g = J(1, 0, .)
    if (M.p1 > 0) g = (g, cross(M.X1, rowsum(wS1))')
    g = (g, colsum(wS1))
    if (M.s1 > 0) g = (g, vec(cross(M.X1S, wS1))')
    if (M.p2 > 0) g = (g, cross(M.X2, rowsum(wS2))')
    g = (g, colsum(wS2))
    if (M.s2 > 0) g = (g, vec(cross(M.X2S, wS2))')
    if (M.nl == 1) return((g, quadsum(wSL)))
    return((g, quadcolsum(wSL)))
}

// H[r0.., c0..] = B (nothing to do when B is empty)
void hmbp_blk(real matrix H, real scalar r0, real scalar c0, real matrix B)
{
    if (rows(B) > 0 & cols(B) > 0) {
        H[|r0,c0 \ r0+rows(B)-1,c0+cols(B)-1|] = B
    }
}

// gradient and Hessian of the expected complete log likelihood, measurement block.
// Only the upper triangle of H is filled, block by block, and then reflected.
void hmbp_mderiv(struct hmbp_dat scalar M, real matrix W,
    struct hmbp_mm scalar mm, real rowvector lam, real rowvector g,
    real matrix H)
{
    real matrix S1, S2, SL, h1, h2, G2, r1, r2, wS1, wS2
    real matrix w11, w22, w12, w1L, w2L, wLL, A
    real scalar q, p1, p2, s1, s2, k, m, c
    real scalar o1, oa1, oS1, o2, oa2, oS2, oL, c1, c2, cL

    if (M.dep == 0) {
        q  = 1 - exp(-1)
        h1 = q*mm.Mu1:*mm.D1
        h2 = q*mm.Mu2:*mm.D2
        G2 = mm.G:^2
        r1 = lam:*(mm.U2:*h1):/mm.G
        r2 = lam:*(mm.U1:*h2):/mm.G
        S1 = M.y1 :- mm.Mu1 + r1
        S2 = M.y2 :- mm.Mu2 + r2
        SL = (mm.U1:*mm.U2):/mm.G
        // weighted second derivatives with respect to (eta1, eta2, lambda)
        w11 = W:*(-mm.Mu1 + r1:*(1:-q*mm.Mu1) - r1:^2)
        w22 = W:*(-mm.Mu2 + r2:*(1:-q*mm.Mu2) - r2:^2)
        w12 = W:*(((lam:*h1):*h2):/G2)
        w1L = W:*((mm.U2:*h1):/G2)
        w2L = W:*((mm.U1:*h2):/G2)
        wLL = W:*(SL:^2)
    }
    else {
        // copulas: minus the outer product of the scores takes the place of
        // the second derivatives
        hmbp_mscore(M, mm, lam, S1, S2, SL)
        w11 = -W:*S1:^2
        w22 = -W:*S2:^2
        w12 = -W:*S1:*S2
        w1L = -W:*S1:*SL
        w2L = -W:*S2:*SL
        wLL = W:*SL:^2
    }
    wS1 = W:*S1
    wS2 = W:*S2
    g = hmbp_mgrad(M, wS1, wS2, W:*SL)

    // first positions: common slopes, intercepts, state-specific slopes, lambda
    p1 = M.p1; p2 = M.p2; s1 = M.s1; s2 = M.s2; k = M.k
    o1  = 1
    oa1 = p1 + 1
    oS1 = p1 + k + 1
    o2  = p1 + k + k*s1 + 1
    oa2 = o2 + p2
    oS2 = o2 + p2 + k
    m   = M.m
    oL  = m - M.nl + 1
    H = J(m, m, 0)

    // equation 1
    if (p1 > 0) {
        hmbp_blk(H, o1, o1,  cross(M.X1, rowsum(w11), M.X1))
        hmbp_blk(H, o1, oa1, cross(M.X1, w11))
    }
    hmbp_blk(H, oa1, oa1, diag(colsum(w11)))
    // equation 2
    if (p2 > 0) {
        hmbp_blk(H, o2, o2,  cross(M.X2, rowsum(w22), M.X2))
        hmbp_blk(H, o2, oa2, cross(M.X2, w22))
    }
    hmbp_blk(H, oa2, oa2, diag(colsum(w22)))
    // equation 1 with equation 2
    if (p1 > 0 & p2 > 0) hmbp_blk(H, o1, o2, cross(M.X1, rowsum(w12), M.X2))
    if (p1 > 0) hmbp_blk(H, o1, oa2, cross(M.X1, w12))
    if (p2 > 0) hmbp_blk(H, oa1, o2, cross(w12, M.X2))
    hmbp_blk(H, oa1, oa2, diag(colsum(w12)))
    // lambda: one column when it is common, one per state otherwise
    if (M.nl == 1) {
        if (p1 > 0) hmbp_blk(H, o1, oL, cross(M.X1, rowsum(w1L)))
        hmbp_blk(H, oa1, oL, colsum(w1L)')
        if (p2 > 0) hmbp_blk(H, o2, oL, cross(M.X2, rowsum(w2L)))
        hmbp_blk(H, oa2, oL, colsum(w2L)')
        H[oL,oL] = -quadsum(wLL)
    }
    else {
        if (p1 > 0) hmbp_blk(H, o1, oL, cross(M.X1, w1L))
        hmbp_blk(H, oa1, oL, diag(colsum(w1L)))
        if (p2 > 0) hmbp_blk(H, o2, oL, cross(M.X2, w2L))
        hmbp_blk(H, oa2, oL, diag(colsum(w2L)))
        hmbp_blk(H, oL, oL, diag(-quadcolsum(wLL)))
    }
    // state-specific slopes: the block of state c starts at c1 (equation 1)
    // and at c2 (equation 2); cL is the column of the lambda of state c
    for (c=1; c<=k; c++) {
        c1 = oS1 + (c-1)*s1
        c2 = oS2 + (c-1)*s2
        cL = (M.nl == 1 ? oL : oL + c - 1)
        if (s1 > 0) {
            A = M.X1S :* w11[.,c]
            if (p1 > 0) hmbp_blk(H, o1, c1, cross(M.X1, A))
            hmbp_blk(H, oa1+c-1, c1, colsum(A))
            hmbp_blk(H, c1, c1, cross(M.X1S, A))
            A = M.X1S :* w12[.,c]
            if (p2 > 0) hmbp_blk(H, c1, o2, cross(A, M.X2))
            hmbp_blk(H, c1, oa2+c-1, colsum(A)')
            if (s2 > 0) hmbp_blk(H, c1, c2, cross(A, M.X2S))
            hmbp_blk(H, c1, cL, cross(M.X1S, w1L[.,c]))
        }
        if (s2 > 0) {
            A = M.X2S :* w22[.,c]
            if (p2 > 0) hmbp_blk(H, o2, c2, cross(M.X2, A))
            hmbp_blk(H, oa2+c-1, c2, colsum(A))
            hmbp_blk(H, c2, c2, cross(M.X2S, A))
            A = M.X2S :* w12[.,c]
            if (p1 > 0) hmbp_blk(H, o1, c2, cross(M.X1, A))
            hmbp_blk(H, oa1+c-1, c2, colsum(A))
            hmbp_blk(H, c2, cL, cross(M.X2S, w2L[.,c]))
        }
    }
    H = makesymmetric(H')
}

// ---------------------------------------------------------------- M-step
// Measurement block: Newton steps on the expected complete log likelihood,
// starting from bm0, where the measurement model is mm0 (from the E-step). The
// step is halved until the objective does not decrease.
// Famoye form: every trial point is made admissible by projecting lambda onto its
// admissible interval at the trial means. A lambda that sits on a bound of its
// interval with the gradient pointing outwards is left out of the Newton
// system, so that the other parameters keep moving (active set).
real rowvector hmbp_mmeas(struct hmbp_dat scalar M, real matrix W,
    real rowvector bm0, real scalar nsteps, struct hmbp_mm scalar mm0)
{
    struct hmbp_mm scalar mm, mn
    real rowvector bm, bn, g, step, s, b1, a1, b2, a2, lam, iL, lv, w, blk, fr
    real matrix    H, S1, S2, bnd
    real scalar    it, j, st, Q0, Q1, eps

    eps = 1e-6
    iL  = (M.m - M.nl + 1)..M.m
    bm  = bm0
    mm  = mm0
    for (it=1; it<=nsteps; it++) {
        lam = hmbp_lam(M, bm)
        Q0  = quadsum(W:*mm.lf)
        hmbp_mderiv(M, W, mm, lam, g, H)
        blk = J(1, M.nl, 0)
        if (M.dep == 0) {
            // lambdas on a bound, with the gradient pointing outwards
            bnd = hmbp_lamint(M, mm)
            w   = bnd[2,.] - bnd[1,.]
            lv  = bm[iL]
            blk = ((lv :>= bnd[2,.] - 2*eps*w) :& (g[iL] :> 0)) :| ((lv :<= bnd[1,.] + 2*eps*w) :& (g[iL] :< 0))
        }
        if (any(blk)) {
            fr = J(1, M.m, 1)
            fr[iL] = (blk :== 0)
            fr = selectindex(fr)
            s  = cholsolve(-H[fr, fr], g[fr]')'
            if (hasmissing(s)) s = (invsym(-H[fr, fr])*g[fr]')'
            if (hasmissing(s) | (g[fr]*s') <= 0) s = g[fr] :/ max(abs(diagonal(H)))
            step = J(1, M.m, 0)
            step[fr] = s
        }
        else {
            step = cholsolve(-H, g')'
            if (hasmissing(step)) step = (invsym(-H)*g')'
            if (hasmissing(step) | (g*step') <= 0) step = g :/ max(abs(diagonal(H)))
        }
        st = 1
        for (j=1; j<=40; j++) {
            bn = bm + st*step
            hmbp_mparts(M, bn, b1, a1, S1, b2, a2, S2, lam)
            mn = hmbp_means(M, b1, a1, S1, b2, a2, S2)
            if (!hasmissing(mn.lf0)) {
                if (M.dep == 0) {
                    // projection of lambda onto the admissible interval at the trial means
                    bnd = hmbp_lamint(M, mn)
                    w   = bnd[2,.] - bnd[1,.]
                    lv  = colmax((colmin((bn[iL] \ (bnd[2,.] - eps*w))) \ (bnd[1,.] + eps*w)))
                    if (lv != bn[iL]) {
                        bn[iL] = lv
                        lam = hmbp_lam(M, bn)
                    }
                }
                hmbp_dep(M, mn, lam)
                if (!hasmissing(mn.lf)) {
                    Q1 = quadsum(W:*mn.lf)
                    if (Q1 >= Q0 - 1e-12*abs(Q0)) break
                }
            }
            st = st/2
        }
        if (j > 40) break
        bm = bn
        mm = mn
        if (max(abs(st*step)) < 1e-10) break
    }
    return(bm)
}

// fractional multinomial logit for the initial probabilities (Newton-Raphson)
real matrix hmbp_minit(real matrix Z, real matrix W1, real matrix D0,
    real scalar k, real scalar nsteps)
{
    real matrix    D, Dn, Pn, Pi, Hm, blk
    real colvector g, step, wg
    real scalar    q, it, c, e, j, st, Q0, Q1

    q  = cols(Z)
    D  = D0
    Pi = hmbp_pi0(Z, D, k)
    Q0 = quadsum(W1:*ln(Pi))
    for (it=1; it<=nsteps; it++) {
        g  = vec(cross(Z, W1[|1,2 \ .,.|] - Pi[|1,2 \ .,.|]))
        Hm = J(q*(k-1), q*(k-1), 0)
        for (c=2; c<=k; c++) {
            for (e=c; e<=k; e++) {
                wg  = Pi[.,c] :* ((c==e) :- Pi[.,e])
                blk = -cross(Z, wg, Z)
                Hm[|(c-2)*q+1,(e-2)*q+1 \ (c-1)*q,(e-1)*q|] = blk
                Hm[|(e-2)*q+1,(c-2)*q+1 \ (e-1)*q,(c-1)*q|] = blk'
            }
        }
        step = cholsolve(-Hm, g)
        if (hasmissing(step)) step = invsym(-Hm)*g
        st = 1
        for (j=1; j<=30; j++) {
            Dn = D + st*colshape(step', q)'
            Pn = hmbp_pi0(Z, Dn, k)
            Q1 = quadsum(W1:*ln(Pn))
            if (Q1 >= Q0 - 1e-12*abs(Q0)) break
            st = st/2
        }
        if (j > 30) break
        D  = Dn
        Pi = Pn
        Q0 = Q1
        if (max(abs(st*step)) < 1e-10) break
    }
    return(D)
}

// ------------------------------------------------------------------ EM
// log of the row-wise softmax
real matrix hmbp_lnsoft(real matrix E)
{
    real colvector m

    m = rowmax(E)
    return(E :- m :- ln(rowsum(exp(E :- m))))
}

// weighted (fractional) multinomial logit by Newton-Raphson with step-halving.
// Y: n x k nonnegative responses (row sums are the weights); B0: q x k
// coefficients with a zero column for the base category. Returns B (q x k).
real matrix hmbp_fmlogit(real matrix X, real matrix Y, real matrix B0,
    real scalar base, real scalar nsteps)
{
    real matrix    B, Bn, lnP, lnPn, Pr, H, blk
    real colvector w, g, step, wg
    real rowvector free
    real scalar    k, q, it, a, e, j, st, Q0, Q1

    k    = cols(Y)
    q    = cols(X)
    w    = rowsum(Y)
    free = select(1..k, (1..k) :!= base)
    B    = B0
    lnP  = hmbp_lnsoft(X*B)
    Q0   = quadsum(Y :* lnP)
    for (it=1; it<=nsteps; it++) {
        Pr = exp(lnP)
        g  = vec(cross(X, Y[., free] - Pr[., free] :* w))
        H  = J(q*(k-1), q*(k-1), 0)
        for (a=1; a<=k-1; a++) {
            for (e=a; e<=k-1; e++) {
                wg  = w :* Pr[., free[a]] :* ((a==e) :- Pr[., free[e]])
                blk = -cross(X, wg, X)
                H[|(a-1)*q+1,(e-1)*q+1 \ a*q,e*q|] = blk
                H[|(e-1)*q+1,(a-1)*q+1 \ e*q,a*q|] = blk'
            }
        }
        step = cholsolve(-H, g)
        if (hasmissing(step)) step = invsym(-H)*g
        if (hasmissing(step)) break
        st = 1
        for (j=1; j<=30; j++) {
            Bn = B
            Bn[., free] = B[., free] + st*colshape(step', q)'
            lnPn = hmbp_lnsoft(X*Bn)
            Q1   = quadsum(Y :* lnPn)
            if (!hasmissing(lnPn) & Q1 >= Q0 - 1e-12*abs(Q0)) break
            st = st/2
        }
        if (j > 30) break
        B   = Bn
        lnP = lnPn
        Q0  = Q1
        if (max(abs(st*step)) < 1e-10) break
    }
    return(B)
}

// M-step for observation-specific transitions: one logit per origin state
real rowvector hmbp_mtrans(struct hmbp_dat scalar M, real rowvector tau,
    real matrix XI)
{
    real matrix    TT, B0, Bn
    real rowvector cc, idx
    real scalar    c, cb, j, i, k

    k  = M.k
    TT = colshape(tau, M.qw)
    j  = 0
    for (c=1; c<=k; c++) {
        cc  = ((c-1)*k+1)..(c*k)
        B0  = J(M.qw, k, 0)
        idx = J(1, k, 0)
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                j++
                B0[.,cb] = TT[j,.]'
                idx[cb]  = j
            }
        }
        Bn = hmbp_fmlogit(M.W, XI[., cc], B0, c, 5)
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                i = idx[cb]
                TT[i,.] = Bn[.,cb]'
            }
        }
    }
    return(vec(TT')')
}

real rowvector hmbp_em(struct hmbp_dat scalar M, real rowvector b0,
    real scalar maxit, real scalar tol, real scalar dolog,
    real scalar itout, real scalar llout)
{
    struct hmbp_mm scalar mm
    real rowvector b, bm, b1, a1, b2, a2, tau, lam
    real matrix    W, Zs, Pi0, P, XiI, D, Pn, S1, S2
    real colvector rs
    real scalar    it, ll, llold, c

    b = b0
    llold = .
    ll = .
    for (it=1; it<=maxit; it++) {
        ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 1, 0, XiI)
        if (ll >= .) {
            errprintf("EM: log likelihood not finite at iteration %5.0f\n", it)
            exit(430)
        }
        if (dolog) {
            printf("{txt}EM iteration %5.0f:{col 22}log likelihood = {res}%12.0g\n", it, ll)
            displayflush()
        }
        if (it > 1) {
            if (abs(ll - llold) <= tol*abs(llold)) break
        }
        llold = ll
        // measurement block: Newton steps with step-halving (two with the
        // Famoye form, one with the copulas, where an evaluation costs more)
        bm = hmbp_mmeas(M, W, b[|1 \ M.m|], (M.dep == 0 ? 2 : 1), mm)
        if (M.k > 1) {
            hmbp_unpack(M, b, b1, a1, S1, b2, a2, S2, lam, D, tau)
            if (M.het) {
                // transitions: one weighted multinomial logit per origin state
                tau = hmbp_mtrans(M, tau, XiI)
            }
            else if (!M.notr) {
                // transitions: closed form
                rs = rowsum(Zs)
                Pn = P
                for (c=1; c<=M.k; c++) {
                    if (rs[c] > 1e-300) Pn[c,.] = Zs[c,.] :/ rs[c]
                }
                Pn  = Pn :+ 1e-10
                Pn  = Pn :/ rowsum(Pn)
                tau = hmbp_tau(Pn)
            }
            // initial probabilities: fractional multinomial logit
            D = hmbp_minit(M.Z, W[M.info[.,1], .], D, M.k, 5)
            b = (bm, vec(D)', tau)
        }
        else b = bm
    }
    itout = (it > maxit ? maxit : it)
    llout = ll
    return(b)
}

// evaluator for optimize(), type d1
void hmbp_eval(real scalar todo, real rowvector b,
    pointer(struct hmbp_dat scalar) scalar pM,
    real scalar v, real rowvector g, real matrix H)
{
    struct hmbp_mm scalar mm
    real matrix W, Zs, Pi0, P, XiI

    v = hmbp_estep(*pM, b, W, Zs, Pi0, P, mm, (todo >= 1), 0, XiI)
    if (v >= .) return
    if (todo >= 1) g = hmbp_grad(*pM, b, W, Zs, Pi0, P, mm, XiI)
}

// Observed information: central differences of the analytic score with respect
// to the free parameters fr (a lambda on a bound is held fixed). The
// admissibility check of the E-step is switched off here: the log likelihood is
// smooth across the bounds, and the perturbed points next to a bound would
// otherwise be rejected.
real matrix hmbp_hess(struct hmbp_dat scalar M, real rowvector b,
    real rowvector fr)
{
    real matrix    H
    real rowvector bp, bq, gp, gq
    real scalar    i, j, h, nf

    nf = cols(fr)
    H  = J(nf, nf, .)
    M.chk = 0
    for (i=1; i<=nf; i++) {
        j  = fr[i]
        h  = 1e-5*(abs(b[j]) + 1)
        bp = b; bp[j] = b[j] + h
        bq = b; bq[j] = b[j] - h
        gp = hmbp_gradonly(M, bp)
        gq = hmbp_gradonly(M, bq)
        H[.,i] = (gp[fr] - gq[fr])' :/ (2*h)
    }
    M.chk = 1
    return((H + H')/2)
}

// relabel states in increasing order of the support points
// relabel the states in increasing order of the mean over the sample of mu_1c
// (how = "eq1") or of mu_2c (how = "eq2"). With common slopes this is the order
// of the state intercepts.
real rowvector hmbp_sort(struct hmbp_dat scalar M, real rowvector b,
    string scalar how)
{
    struct hmbp_mm scalar mm
    real rowvector b1, a1, b2, a2, tau, key, lam, lv
    real matrix    S1, S2, D, Df, TT, TN, pos
    real colvector perm
    real scalar    c, cb, j, k

    if (M.k == 1 | how == "none") return(b)
    k = M.k
    hmbp_unpack(M, b, b1, a1, S1, b2, a2, S2, lam, D, tau)
    mm   = hmbp_meas(M, b1, a1, S1, b2, a2, S2, lam)
    key  = (how == "eq2" ? mean(mm.Mu2) : mean(mm.Mu1))
    perm = order(key', 1)
    a1 = a1[perm']
    a2 = a2[perm']
    if (M.s1 > 0) S1 = S1[., perm']
    if (M.s2 > 0) S2 = S2[., perm']
    lv = b[|M.m-M.nl+1 \ M.m|]
    if (M.nl > 1) lv = lv[perm']
    Df = (J(M.q, 1, 0), D)[., perm']
    Df = Df :- Df[., 1]
    D  = Df[|1,2 \ .,.|]
    // transitions: the base of each row is "staying", which relabelling preserves
    pos = J(k, k, 0)
    j = 0
    for (c=1; c<=k; c++) {
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                j++
                pos[c,cb] = j
            }
        }
    }
    if (M.ntau == 0) return((b1, a1, vec(S1)', b2, a2, vec(S2)', lv, vec(D)'))
    TT = colshape(tau, M.qw)
    TN = TT
    for (c=1; c<=k; c++) {
        for (cb=1; cb<=k; cb++) {
            if (cb != c) TN[pos[c,cb], .] = TT[pos[perm[c],perm[cb]], .]
        }
    }
    return((b1, a1, vec(S1)', b2, a2, vec(S2)', lv, vec(D)', vec(TN')'))
}

// ------------------------------------------------------------------ driver
// data matrix for a list of variable names (N x 0 when the list is empty)
real matrix hmbp_getx(string scalar names, string scalar touse, real scalar N)
{
    if (strtrim(names) == "") return(J(N, 0, .))
    return(st_data(., tokens(names), touse))
}

// everything that depends only on the two counts
void hmbp_sety(struct hmbp_dat scalar M, real colvector y1, real colvector y2)
{
    M.y1   = y1
    M.y2   = y2
    M.ey1  = exp(-y1)
    M.ey2  = exp(-y2)
    M.lnfy = lngamma(y1 :+ 1) + lngamma(y2 :+ 1)
    if (M.dep > 0) {
        // copulas: counts as N x k matrices, log factorials by margin
        M.Y1   = J(1, M.k, y1)
        M.Y2   = J(1, M.k, y2)
        M.lfy1 = lngamma(y1 :+ 1)
        M.lfy2 = lngamma(y2 :+ 1)
        // rows where both counts are zero, and the other rows
        M.i00  = selectindex((y1 :== 0) :& (y2 :== 0))
        M.i11  = selectindex((y1 :> 0) :| (y2 :> 0))
    }
}

// Data and dimensions of the model, from the variables in memory: the rows
// selected by touse, sorted by decreasing T, individual and period. It is used
// by the estimation and by the postestimation tasks. tt (period of each
// transition row), u (distinct periods) and nwx (number of covariates of the
// transitions) are returned to the caller.
void hmbp_setup(struct hmbp_dat scalar M, string scalar y1n,
    string scalar y2n, string scalar x1n, string scalar x2n,
    string scalar zn, string scalar idn, string scalar touse,
    real scalar k, string scalar wn, string scalar trn,
    string scalar timen, real scalar het, real scalar tvt,
    string scalar xs1n, string scalar xs2n, real scalar lamst,
    real scalar dep, real scalar notr, real colvector tt, real colvector u,
    real scalar nwx)
{
    real colvector id, sel
    real matrix    zr, Wx
    real scalar    t

    M.k    = k
    M.dep  = dep
    M.notr = notr
    hmbp_sety(M, st_data(., y1n, touse), st_data(., y2n, touse))
    M.N    = rows(M.y1)
    M.X1   = hmbp_getx(x1n, touse, M.N)
    M.X2   = hmbp_getx(x2n, touse, M.N)
    M.X1S  = hmbp_getx(xs1n, touse, M.N)
    M.X2S  = hmbp_getx(xs2n, touse, M.N)
    id     = st_data(., idn, touse)
    zr     = st_data(., tokens(zn), touse)
    M.info = panelsetup(id, 1)
    M.n    = rows(M.info)
    M.Ts   = M.info[.,2] - M.info[.,1] :+ 1
    M.Tmax = max(M.Ts)
    M.Z    = (zr[M.info[.,1], .], J(M.n, 1, 1))
    M.p1   = cols(M.X1)
    M.p2   = cols(M.X2)
    M.s1   = cols(M.X1S)
    M.s2   = cols(M.X2S)
    M.q    = cols(M.Z)
    M.nl   = (lamst ? k : 1)
    M.chk  = 1
    M.m    = M.p1 + M.p2 + 2*k + k*(M.s1 + M.s2) + M.nl
    M.nD   = M.q*(k-1)
    M.het  = het
    M.qw   = 1
    nwx    = 0
    if (M.n > 1) {
        if (any(M.Ts[|2 \ M.n|] :> M.Ts[|1 \ M.n-1|])) {
            errprintf("hmbp_setup(): individuals must be sorted by decreasing T\n")
            exit(498)
        }
    }
    M.ns = J(M.Tmax, 1, .)
    M.R  = J(M.Tmax, 1, NULL)
    for (t=1; t<=M.Tmax; t++) {
        M.ns[t] = sum(M.Ts :>= t)
        M.R[t]  = &(M.info[|1,1 \ M.ns[t],1|] :+ (t-1))
    }

    // transition design: one row per transition row (estimation rows after the
    // first of each individual), in the same order as the data
    if (het) {
        tt = st_data(., timen, trn)
        if (strtrim(wn) != "") Wx = st_data(., tokens(wn), trn)
        else                   Wx = J(rows(tt), 0, .)
        nwx = cols(Wx)
        u   = uniqrows(tt)
        if (tvt) M.W = (Wx, (J(1, rows(u), tt) :== J(rows(tt), 1, u')))
        else     M.W = (Wx, J(rows(tt), 1, 1))
        M.qw = cols(M.W)
        if (rows(M.W) != M.N - M.n) {
            errprintf("hmbp_setup(): transition rows do not match the panel structure\n")
            exit(498)
        }
        sel = J(M.N, 1, 1)
        sel[M.info[.,1]] = J(M.n, 1, 0)
        M.trow = selectindex(sel)
        M.RT = J(M.Tmax, 1, NULL)
        for (t=2; t<=M.Tmax; t++) {
            M.RT[t] = &((*M.R[t]) - (1::M.ns[t]))
        }
    }
    M.ntau = (notr ? 0 : k*(k-1)*M.qw)
    M.npar = M.m + M.nD + M.ntau
}

void hmbp_main(string scalar y1n, string scalar y2n, string scalar x1n,
    string scalar x2n, string scalar zn, string scalar idn,
    string scalar touse, real scalar k, string scalar s1n,
    string scalar s2n, string scalar fromn, real scalar emit,
    real scalar emtol, real scalar maxit, string scalar tech,
    real scalar dolog, string scalar sortby, string scalar vce,
    real scalar doml, string scalar wn, string scalar trn,
    string scalar timen, real scalar het, real scalar tvt,
    string scalar xs1n, string scalar xs2n, real scalar lamst,
    real scalar dep, real scalar notr)
{
    struct hmbp_dat scalar M
    struct hmbp_mm scalar  mm
    real colvector tt, u
    real matrix    W, Zs, Pi0, P, XiI, H, Vo, V, Sc, D, omt, Pm, Pit, S1, S2
    real matrix    bnd, lint
    real rowvector b, bs1, bs2, off, a1, a2, b1, b2, tau, lam, iL, wl, atb, frv, fr
    real scalar    t, emits, ll, conv, rc, rank, nwx
    transmorphic   S

    // ---------------------------------------------------------- data
    hmbp_setup(M, y1n, y2n, x1n, x2n, zn, idn, touse, k, wn, trn, timen,
        het, tvt, xs1n, xs2n, lamst, dep, notr, tt, u, nwx)
    if (het) st_local("hmbp_tvals", invtokens(strofreal(u')))

    // ---------------------------------------------------------- start
    if (fromn != "") {
        b = st_matrix(fromn)
        if (cols(b) != M.npar) {
            errprintf("from(): vector has %5.0f elements; %5.0f required\n", cols(b), M.npar)
            exit(503)
        }
    }
    else {
        bs1 = st_matrix(s1n)
        bs2 = st_matrix(s2n)
        if (k > 1) off = 0.8*((1..k) :- (k+1)/2)/((k-1)/2)
        else       off = 0
        // bs1, bs2: Poisson regressions with (common terms, state-specific
        // terms, _cons); the state-specific slopes start equal across states
        b = J(1, 0, .)
        if (M.p1 > 0) b = (b, bs1[|1 \ M.p1|])
        b = (b, bs1[cols(bs1)] :+ off)
        if (M.s1 > 0) b = (b, J(1, k, 1) # bs1[|M.p1+1 \ M.p1+M.s1|])
        if (M.p2 > 0) b = (b, bs2[|1 \ M.p2|])
        b = (b, bs2[cols(bs2)] :+ off)
        if (M.s2 > 0) b = (b, J(1, k, 1) # bs2[|M.p2+1 \ M.p2+M.s2|])
        b = (b, J(1, M.nl, 0))
        // transitions start from a persistent chain (staying probability 0.9)
        if (k > 1) b = (b, J(1, M.nD, 0))
        if (M.ntau > 0) b = (b, J(1, k*(k-1), 1) #
                        (J(1, nwx, 0), J(1, M.qw - nwx, ln(0.1/(k-1)/0.9))))
    }
    ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 0, 0, XiI)
    if (ll >= .) {
        errprintf("log likelihood not finite at the starting values\n")
        exit(430)
    }

    // ---------------------------------------------------------- EM
    emits = 0
    if (emit > 0) {
        if (dolog) printf("{txt}\nEM algorithm\n")
        b = hmbp_em(M, b, emit, emtol, dolog, emits, ll)
    }
    conv = (emits < emit)

    // ---------------------------------------------------------- direct ML
    if (doml == 2 | (doml == 1 & conv == 0)) {
        if (dolog) printf("{txt}\nDirect maximisation of the marginal log likelihood\n")
        S = optimize_init()
        optimize_init_evaluator(S, &hmbp_eval())
        optimize_init_evaluatortype(S, "d1")
        optimize_init_argument(S, 1, &M)
        optimize_init_params(S, b)
        optimize_init_technique(S, tech)
        optimize_init_conv_maxiter(S, maxit)
        optimize_init_tracelevel(S, (dolog ? "value" : "none"))
        optimize_init_conv_warning(S, "off")
        optimize_init_singularHmethod(S, "hybrid")
        rc = _optimize(S)
        if (rc) {
            printf("{err}warning: direct maximisation failed (%s); EM estimates kept\n",
                optimize_result_errortext(S))
        }
        else {
            b    = optimize_result_params(S)
            conv = optimize_result_converged(S)
        }
    }

    // ---------------------------------------------------------- final
    b  = hmbp_sort(M, b, sortby)
    ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 1, 0, XiI)
    // a lambda on a bound of its admissible interval is held fixed in the
    // variance estimate: its rows and columns of V are zero
    iL   = (M.m - M.nl + 1)..M.m
    atb  = J(1, M.nl, 0)
    bnd  = J(2, k, .)
    if (M.dep == 0) {
        bnd  = hmbp_lambnd(mm)
        lint = hmbp_lamint(M, mm)
        wl   = lint[2,.] - lint[1,.]
        atb  = (b[iL] :>= lint[2,.] - 1e-4*wl) :| (b[iL] :<= lint[1,.] + 1e-4*wl)
    }
    frv  = J(1, M.npar, 1)
    frv[iL] = (atb :== 0)
    fr   = selectindex(frv)
    if (dolog) printf("{txt}\nComputing the observed information matrix...\n")
    displayflush()
    H  = hmbp_hess(M, b, fr)
    Vo = J(M.npar, M.npar, 0)
    Vo[fr, fr] = invsym(-H)
    rank = sum(diagonal(Vo) :!= 0)
    if (vce == "robust") {
        Sc = hmbp_scores(M, b)
        V  = Vo*((M.n/(M.n-1))*cross(Sc, Sc))*Vo
    }
    else V = Vo
    V = (V + V')/2

    omt = J(M.Tmax, k, .)
    for (t=1; t<=M.Tmax; t++) omt[t,.] = mean(W[*M.R[t], .])
    hmbp_unpack(M, b, b1, a1, S1, b2, a2, S2, lam, D, tau)
    if (M.het) {
        // average transition matrix, overall and by period of arrival
        Pm  = colshape(mean(P), k)
        Pit = J(rows(u), k*k, .)
        for (t=1; t<=rows(u); t++) {
            Pit[t,.] = mean(P[selectindex(tt :== u[t]), .])
        }
        st_matrix("__hmbp_Pit", Pit)
    }
    else Pm = P

    st_matrix("__hmbp_b",   b)
    st_matrix("__hmbp_V",   V)
    st_matrix("__hmbp_Pi",  Pm)
    st_matrix("__hmbp_om",  mean(panelsum(W, M.info) :/ M.Ts))
    st_matrix("__hmbp_omt", omt)
    st_matrix("__hmbp_al",  (a1', a2'))
    st_matrix("__hmbp_pi0", mean(Pi0))
    st_numscalar("__hmbp_n",     M.n)
    st_numscalar("__hmbp_ll",    ll)
    st_numscalar("__hmbp_npar",  M.npar)
    st_numscalar("__hmbp_rank",  rank)
    st_numscalar("__hmbp_conv",  conv)
    st_numscalar("__hmbp_emit",  emits)
    st_numscalar("__hmbp_gmin",  min(M.Ts))
    st_numscalar("__hmbp_gavg",  mean(M.Ts))
    st_numscalar("__hmbp_gmax",  M.Tmax)
    st_matrix("__hmbp_lrng", bnd')
    st_numscalar("__hmbp_lamlo", max(bnd[1,.]))
    st_numscalar("__hmbp_lamhi", min(bnd[2,.]))
    st_numscalar("__hmbp_nbnd",  sum(atb))
}

// ------------------------------------------------------------ postestimation
// index of the first maximum of each row
real colvector hmbp_argmax(real matrix X)
{
    real colvector best, bv, upd
    real scalar    c

    best = J(rows(X), 1, 1)
    bv   = X[.,1]
    for (c=2; c<=cols(X); c++) {
        upd  = X[.,c] :> bv
        best = best + upd:*(c :- best)
        bv   = bv + upd:*(X[.,c] - bv)
    }
    return(best)
}

// Probabilities of the states by observation: filtered (ptype = 1: given the
// outcomes up to the period), one-step-ahead (2: given the outcomes up to the
// period before) or prior (3: given no outcome). P is the transition matrix,
// or the probabilities by transition row when the chain is heterogeneous.
real matrix hmbp_sprob(struct hmbp_dat scalar M, real matrix lf,
    real matrix Pi0, real matrix P, real scalar ptype)
{
    real matrix    A, F, Q, a, ap
    real colvector ll, r, rp, rt
    real scalar    t, nt, c, k

    k = M.k
    if (ptype != 3) {
        if (M.het) ll = hmbp_fwdh(M, lf, Pi0, P, A, F)
        else       ll = hmbp_fwd(M, lf, Pi0, P, A, F)
        if (ptype == 1) return(A)
    }
    Q = J(M.N, k, 0)
    r = *M.R[1]
    Q[r,.] = Pi0
    a = Pi0
    for (t=2; t<=M.Tmax; t++) {
        nt = M.ns[t]
        r  = *M.R[t]
        // probabilities of the period before: filtered or prior
        if (ptype == 2) {
            rp = (*M.R[t-1])[|1 \ nt|]
            ap = A[rp, .]
        }
        else ap = a[|1,1 \ nt,k|]
        if (M.het) {
            rt = *M.RT[t]
            a  = J(nt, k, 0)
            for (c=1; c<=k; c++) {
                a = a + ap[.,c] :* P[rt, ((c-1)*k+1)..(c*k)]
            }
        }
        else a = ap*P
        Q[r,.] = a
    }
    return(Q)
}

// Most likely sequence of states of each individual (Viterbi algorithm, in
// logs). Dl holds the log probability of the best sequence that ends in each
// state, Psi the state of the period before along that sequence.
real colvector hmbp_viterbi(struct hmbp_dat scalar M, real matrix lf,
    real matrix Pi0, real matrix P)
{
    real matrix    lnP, Dl, Psi, d, dn, pb, cand
    real colvector r, r1, rt, S, cur, nxt, bv, arg, upd
    real scalar    t, nt, n1, k, c, e

    k   = M.k
    Dl  = J(M.N, k, .)
    Psi = J(M.N, k, 1)
    lnP = ln(P :+ 1e-300)
    r   = *M.R[1]
    d   = ln(Pi0 :+ 1e-300) + lf[r,.]
    Dl[r,.] = d
    for (t=2; t<=M.Tmax; t++) {
        nt = M.ns[t]
        r  = *M.R[t]
        d  = d[|1,1 \ nt,k|]
        dn = J(nt, k, .)
        pb = J(nt, k, 1)
        if (M.het) rt = *M.RT[t]
        for (e=1; e<=k; e++) {
            // best state of the period before, for arrival state e
            if (M.het) cand = d + lnP[rt, ((0..(k-1))*k) :+ e]
            else       cand = d :+ lnP[., e]'
            bv  = cand[.,1]
            arg = J(nt, 1, 1)
            for (c=2; c<=k; c++) {
                upd = cand[.,c] :> bv
                arg = arg + upd:*(c :- arg)
                bv  = bv + upd:*(cand[.,c] - bv)
            }
            dn[.,e] = bv
            pb[.,e] = arg
        }
        d = dn + lf[r,.]
        Dl[r,.]  = d
        Psi[r,.] = pb
    }
    // backtracking: individuals are sorted by decreasing T, so those observed
    // at t+1 are the first ns[t+1] of those observed at t
    S = J(M.N, 1, .)
    for (t=M.Tmax; t>=1; t--) {
        r   = *M.R[t]
        cur = hmbp_argmax(Dl[r,.])
        if (t < M.Tmax) {
            n1 = M.ns[t+1]
            r1 = *M.R[t+1]
            cur[|1 \ n1|] = rowsum(Psi[r1, .] :* (J(1, k, nxt) :== (1..k)))
        }
        S[r] = cur
        nxt  = cur
    }
    return(S)
}

// for each Poisson mean, a point a with P(Y > a) < 1e-10, where the sums over
// the support are stopped
real colvector hmbp_ptrunc(real colvector mu)
{
    real colvector a, todo, keep
    real scalar    j

    a    = floor(mu + 4*sqrt(mu))
    todo = selectindex(poissontail(mu, a :+ 1) :>= 1e-10)
    for (j=1; j<=5000; j++) {
        if (rows(todo) == 0) break
        a[todo] = a[todo] :+ 1
        keep    = selectindex(poissontail(mu[todo], a[todo] :+ 1) :>= 1e-10)
        if (rows(keep) == 0) break
        todo    = todo[keep]
    }
    return(a)
}

// copula C(u,v), element by element; par is theta (Frank) or the correlation
// (normal)
real colvector hmbp_cop(real scalar dep, real colvector u, real colvector v,
    real colvector par)
{
    real colvector sm, T

    if (dep == 1) return(binormal(hmbp_zof(u), hmbp_zof(v), par))
    sm = abs(par) :< 1e-6
    T  = par + sm
    return((1:-sm):*hmbp_fC(hmbp_expm1(-T:*u), hmbp_expm1(-T:*v), hmbp_expm1(-T), T) + sm:*hmbp_fCs(u, v, par))
}

// Copula models: correlation between the two counts in each state (N x k), from
//   Cov(Y1,Y2) = sum over a, b >= 0 of [C(F1(a),F2(b)) - F1(a) F2(b)]
// and the Poisson variances. Each cell (observation, state) is summed up to
// the points where the upper tails of its two margins are negligible; cells
// with a mean above 300 are returned as missing.
real matrix hmbp_cprho(struct hmbp_dat scalar M, real matrix Mu1,
    real matrix Mu2, real rowvector lam)
{
    real colvector mu1, mu2, par, ok, A, B, cov, f1, F1, idx, act, u, m2a
    real colvector pa, Ba, f2, F2, s, sub, rho
    real scalar    a, b, amax, bmax

    mu1 = vec(Mu1)
    mu2 = vec(Mu2)
    par = vec(J(M.N, 1, lam))
    if (M.dep == 1) par = tanh(par)
    ok  = (mu1 :<= 300) :& (mu2 :<= 300)
    A   = J(rows(mu1), 1, -1)
    B   = A
    idx = selectindex(ok)
    if (rows(idx) > 0) {
        A[idx] = hmbp_ptrunc(mu1[idx])
        B[idx] = hmbp_ptrunc(mu2[idx])
    }
    cov  = J(rows(mu1), 1, 0)
    f1   = exp(-mu1)
    F1   = f1
    amax = max(A)
    for (a=0; a<=amax; a++) {
        // cells whose sum over the first count has not ended
        act = selectindex(A :>= a)
        u   = F1[act]
        m2a = mu2[act]
        pa  = par[act]
        Ba  = B[act]
        f2  = exp(-m2a)
        F2  = f2
        s   = J(rows(act), 1, 0)
        bmax = max(Ba)
        for (b=0; b<=bmax; b++) {
            sub = selectindex(Ba :>= b)
            s[sub] = s[sub] + hmbp_cop(M.dep, u[sub], F2[sub], pa[sub]) - u[sub]:*F2[sub]
            f2 = f2:*m2a:/(b+1)
            F2 = F2 + f2
            F2 = F2:*(F2 :< 1) + (F2 :>= 1)
        }
        cov[act] = cov[act] + s
        f1 = f1:*mu1:/(a+1)
        F1 = F1 + f1
        F1 = F1:*(F1 :< 1) + (F1 :>= 1)
    }
    rho = (cov:/sqrt(mu1:*mu2)):/ok
    return(rowshape(rho, M.k)')
}

// correlation between the two counts in each state, by observation (N x k)
real matrix hmbp_rho(struct hmbp_dat scalar M, struct hmbp_mm scalar mm,
    real rowvector lam)
{
    real scalar q

    if (M.dep > 0) return(hmbp_cprho(M, mm.Mu1, mm.Mu2, lam))
    q = 1 - exp(-1)
    return(lam:*(q^2*sqrt(mm.Mu1:*mm.Mu2):*mm.D1:*mm.D2))
}

// Predictions for the estimation sample, one row per observation.
// stat: mu1, mu2, postprob, state, viterbi, rho, pr. The states are weighted
// with the probabilities ptype (0 posterior, 1 filtered, 2 one-step-ahead,
// 3 prior); ist > 0 conditions on state ist. pr: probability of the pair
// (ya, yb) at the regressors and lags of each observation.
real matrix hmbp_predict(struct hmbp_dat scalar M, real rowvector b,
    string scalar stat, real scalar ptype, real scalar ist,
    real scalar ya, real scalar yb)
{
    struct hmbp_mm scalar mm, m2
    real rowvector b1, a1, b2, a2, lam, tau
    real matrix    S1, S2, D, W, Zs, Pi0, P, XiI, Pr, f
    real scalar    ll

    ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 1, 0, XiI)
    if (ll >= .) {
        errprintf("the log likelihood cannot be evaluated at e(b) with the data in memory\n")
        exit(459)
    }
    if (reldif(ll, st_numscalar("e(ll)")) > 1e-8) {
        printf("{txt}note: the log likelihood at e(b) with the data in memory, {res}%12.0g{txt}, differs from e(ll);\n", ll)
        printf("{txt}      the data are not those used in the estimation\n")
    }
    if (stat == "viterbi") return(hmbp_viterbi(M, mm.lf, Pi0, P))

    // probabilities of the states that weight the statistic
    if (ist > 0) {
        Pr = J(M.N, M.k, 0)
        Pr[., ist] = J(M.N, 1, 1)
    }
    else if (ptype == 0) Pr = W
    else Pr = hmbp_sprob(M, mm.lf, Pi0, P, ptype)
    if (stat == "postprob") return(Pr)
    if (stat == "state")    return(hmbp_argmax(Pr))
    if (stat == "mu1")      return(rowsum(Pr:*mm.Mu1))
    if (stat == "mu2")      return(rowsum(Pr:*mm.Mu2))
    hmbp_unpack(M, b, b1, a1, S1, b2, a2, S2, lam, D, tau)
    if (stat == "rho")      return(rowsum(Pr:*hmbp_rho(M, mm, lam)))

    // pr: the measurement model is evaluated with the two counts set to (ya, yb)
    hmbp_sety(M, J(M.N, 1, ya), J(M.N, 1, yb))
    m2 = hmbp_meas(M, b1, a1, S1, b2, a2, S2, lam)
    if (M.dep == 0) {
        f = exp(m2.lf0 :- M.lnfy):*m2.G
        f = f:*(f :> 0)
    }
    else f = m2.G:*(m2.G :> 1e-300)
    return(rowsum(Pr:*f))
}

// driver of the postestimation tasks: the data are rebuilt by the Stata
// program PostRun, the parameters are e(b); the results are stored in the
// variables outn, on the rows of the estimation sample
void hmbp_post(string scalar task, string scalar y1n, string scalar y2n,
    string scalar x1n, string scalar x2n, string scalar zn,
    string scalar idn, string scalar touse, real scalar k,
    string scalar wn, string scalar trn, string scalar timen,
    real scalar het, real scalar tvt, string scalar xs1n,
    string scalar xs2n, real scalar lamst, real scalar dep,
    string scalar outn, string scalar stat, real scalar ptype,
    real scalar ist, real scalar ya, real scalar yb, real scalar notr,
    real scalar reps, real scalar seed)
{
    struct hmbp_dat scalar M
    real colvector tt, u
    real rowvector b
    real matrix    R
    real scalar    nwx

    hmbp_setup(M, y1n, y2n, x1n, x2n, zn, idn, touse, k, wn, trn, timen,
        het, tvt, xs1n, xs2n, lamst, dep, notr, tt, u, nwx)
    b = st_matrix("e(b)")
    if (cols(b) != M.npar) {
        errprintf("e(b) has %5.0f elements; the model rebuilt from the data has %5.0f\n", cols(b), M.npar)
        exit(459)
    }
    // the estimates are taken as they are, also with a lambda on a bound
    M.chk = 0
    if (task == "predict") {
        R = hmbp_predict(M, b, stat, ptype, ist, ya, yb)
        st_store(., tokens(outn), touse, R)
    }
    if (task == "lrtrans") st_matrix(outn, hmbp_lrw(M, b, reps, seed))
    if (task == "ic")      st_numscalar(outn, hmbp_entropy(M, b))
    if (task == "trans")   st_matrix(outn, hmbp_trse(M, b, tt))
    if (task == "stability") hmbp_stab(M, b, outn, seed, touse)
    if (task == "ame")     hmbp_ame(M, b, outn, ptype, ist, touse)
}

// ------------------------------------------------------------- estat lrtrans
// Individual derivatives of ln L_i with respect to the transition
// probabilities P_cd, at P = I (n x k*k, columns (c-1)*k+d): sum over t of
// A_t-1(c) F_t(d) B_t(d) / norm_t, from the scaled forward and backward
// variables of the model without transitions
real matrix hmbp_dpsc(struct hmbp_dat scalar M, real matrix lf,
    real matrix Pi0)
{
    real matrix    P, A, F, Dn, B, G, Ap, bb
    real colvector ll, r, rp, nrm
    real scalar    t, nt, k, c

    k  = M.k
    P  = I(k)
    ll = hmbp_fwd(M, lf, Pi0, P, A, F)
    Dn = J(M.n, k*k, 0)
    if (M.Tmax < 2) return(Dn)
    B = J(M.ns[M.Tmax], k, 1)
    for (t=M.Tmax; t>=2; t--) {
        nt  = M.ns[t]
        r   = *M.R[t]
        rp  = (*M.R[t-1])[|1 \ nt|]
        G   = F[r,.] :* B
        Ap  = A[rp,.]
        nrm = rowsum((Ap*P) :* G)
        Ap  = Ap :/ nrm
        for (c=1; c<=k; c++) {
            Dn[|1,(c-1)*k+1 \ nt,c*k|] = Dn[|1,(c-1)*k+1 \ nt,c*k|] + Ap[.,c] :* G
        }
        bb = G*P'
        bb = bb :/ rowsum(bb)
        B  = J(M.ns[t-1], k, 1)
        B[|1,1 \ nt,k|] = bb
    }
    return(Dn)
}

// min 0.5 th'Q th - c'th subject to th >= 0 (Lawson-Hanson active set on the
// normal equations); Q positive definite
real colvector hmbp_nnqp(real matrix Q, real colvector c)
{
    real colvector th, s, w, Ps, cand, idx, neg
    real scalar    r, it, inner, j, alpha

    r  = rows(c)
    th = J(r, 1, 0)
    Ps = J(r, 1, 0)
    w  = c
    for (it=1; it<=3*r+10; it++) {
        cand = selectindex((Ps :== 0) :& (w :> 1e-10))
        if (rows(cand) == 0) break
        j = cand[hmbp_argmax(w[cand]')]
        Ps[j] = 1
        for (inner=1; inner<=3*r+10; inner++) {
            idx = selectindex(Ps)
            s   = J(r, 1, 0)
            s[idx] = lusolve(Q[idx, idx], c[idx])
            if (all(s[idx] :> 1e-12)) {
                th = s
                break
            }
            neg   = idx[selectindex(s[idx] :<= 1e-12)]
            alpha = min(th[neg] :/ (th[neg] - s[neg]))
            if (alpha >= .) alpha = 0
            th = th + alpha*(s - th)
            Ps = Ps :& (th :> 1e-12)
            th = th :* Ps
        }
        w = c - Q*th
    }
    return(th)
}

// weights of the chi-bar-squared distribution of the LR test of no
// transitions (1 x k(k-1)+1, for j = 0, 1, ... positive components), at the
// estimates b of the model without transitions. The parameters are the
// off-diagonal transition probabilities, with P_cc = 1 - sum of the others;
// their information net of the other parameters is computed from the outer
// product of the individual scores. The state of the random-number generator
// is restored at the end.
real rowvector hmbp_lrw(struct hmbp_dat scalar M, real rowvector b,
    real scalar reps, real scalar seed)
{
    struct hmbp_mm scalar mm
    real matrix    W, Zs, Pi0, P, XiI, Sc, Dn, Sp, S, Jm, Q, L
    real colvector cv, th
    real rowvector w
    real scalar    ll, k, nr, c, e, j, i, np
    string scalar  rs

    k  = M.k
    ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 0, 0, XiI)
    if (ll >= .) {
        errprintf("the log likelihood of the model without transitions cannot be evaluated\n")
        exit(459)
    }
    Sc = hmbp_scores(M, b)
    Dn = hmbp_dpsc(M, mm.lf, Pi0)
    nr = k*(k-1)
    Sp = J(M.n, nr, .)
    j  = 0
    for (c=1; c<=k; c++) {
        for (e=1; e<=k; e++) {
            if (e != c) {
                j++
                Sp[.,j] = Dn[.,(c-1)*k+e] - Dn[.,(c-1)*k+c]
            }
        }
    }
    S  = (Sc, Sp)
    Jm = quadcross(S, S)
    np = cols(Sc)
    Q  = Jm[|np+1,np+1 \ .,.|] - Jm[|np+1,1 \ .,np|]*invsym(Jm[|1,1 \ np,np|])*Jm[|1,np+1 \ np,.|]
    Q  = (Q + Q')/2
    L  = cholesky(Q)
    if (hasmissing(L)) {
        errprintf("the information of the transition probabilities is not positive definite\n")
        exit(498)
    }
    rs = rngstate()
    rseed(seed)
    w  = J(1, nr+1, 0)
    for (i=1; i<=reps; i++) {
        cv = L*rnormal(nr, 1, 0, 1)
        th = hmbp_nnqp(Q, cv)
        j  = sum(th :> 1e-12)
        w[j+1] = w[j+1] + 1
    }
    rngstate(rs)
    return(w/reps)
}

// ------------------------------------------------------------------ estat ic
// relative entropy of the classification: 1 - sum(-w ln w)/(N ln k), with w
// the posterior probabilities of the states (1: states known with certainty)
real scalar hmbp_entropy(struct hmbp_dat scalar M, real rowvector b)
{
    struct hmbp_mm scalar mm
    real matrix W, Zs, Pi0, P, XiI
    real scalar ll

    if (M.k == 1) return(1)
    ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 1, 0, XiI)
    if (ll >= .) {
        errprintf("the log likelihood cannot be evaluated at e(b) with the data in memory\n")
        exit(459)
    }
    return(1 - quadsum(-W:*ln(W :+ (W :<= 0)))/(M.N*ln(M.k)))
}

// ---------------------------------------------------------- estat transition
// transition probabilities (by rows: element (c-1)*k+d = P(state d | state c))
// and mean initial probabilities. Heterogeneous chain: mean over the
// transition rows sel.
real rowvector hmbp_trf(struct hmbp_dat scalar M, real rowvector b,
    real colvector sel)
{
    real rowvector b1, a1, b2, a2, lam, tau, Pm
    real matrix    S1, S2, D, PT

    hmbp_unpack(M, b, b1, a1, S1, b2, a2, S2, lam, D, tau)
    if (M.het) {
        PT = hmbp_PT(M, tau)
        Pm = mean(PT[sel, .])
    }
    else if (M.notr) Pm = vec(I(M.k))'
    else             Pm = vec(hmbp_P(tau, M.k)')'
    return((Pm, mean(hmbp_pi0(M.Z, D, M.k))))
}

// the same quantities (row 1) and their delta-method standard errors (row 2),
// from e(V), with the Jacobian by central differences. per (local of the
// caller): period of arrival of the transition rows to average over
// (heterogeneous chain), all of them when empty
real matrix hmbp_trse(struct hmbp_dat scalar M, real rowvector b,
    real colvector tt)
{
    real colvector sel
    real rowvector f, fp, fq, bp, bq, idx
    real matrix    G, V
    real scalar    j, i, h, per

    sel = J(0, 1, .)
    if (M.het) {
        per = strtoreal(st_local("per"))
        if (per < .) sel = selectindex(tt :== per)
        else         sel = (1::rows(tt))
        if (rows(sel) == 0) {
            errprintf("period(): no transition into period %g in the estimation sample\n", per)
            exit(498)
        }
    }
    f   = hmbp_trf(M, b, sel)
    idx = (M.m+1)..M.npar
    G   = J(cols(f), cols(idx), 0)
    for (i=1; i<=cols(idx); i++) {
        j  = idx[i]
        h  = 1e-6*(abs(b[j]) + 1)
        bp = b; bp[j] = b[j] + h
        bq = b; bq[j] = b[j] - h
        fp = hmbp_trf(M, bp, sel)
        fq = hmbp_trf(M, bq, sel)
        G[., i] = (fp - fq)' :/ (2*h)
    }
    V = st_matrix("e(V)")
    V = V[idx, idx]
    return(f \ sqrt(diagonal(G*V*G'))')
}

// ----------------------------------------------------------- estat stability
// lag features (N x 4): g(y1), g(y2) and, with logzero, the indicators of
// y1 = 0 and y2 = 0 (zero otherwise). lagt and larg: transformation and its
// argument, as in lagtransform()
real matrix hmbp_lagf(real colvector z1, real colvector z2,
    string scalar lagt, real scalar larg)
{
    real matrix F

    F = J(rows(z1), 4, 0)
    if (lagt == "linear") F[., 1..2] = (z1, z2)
    else if (lagt == "log1p") F[., 1..2] = ln(1 :+ (z1, z2))
    else if (lagt == "logc") F[., 1..2] = (ln(rowmax((z1, J(rows(z1), 1, larg)))), ln(rowmax((z2, J(rows(z2), 1, larg)))))
    else if (lagt == "indicator") F[., 1..2] = ((z1 :> 0), (z2 :> 0))
    else if (lagt == "topcode") F[., 1..2] = (rowmin((z1, J(rows(z1), 1, larg))), rowmin((z2, J(rows(z2), 1, larg))))
    else if (lagt == "logzero") {
        F[., 1] = (z1 :> 0):*ln(z1 + (z1 :<= 0))
        F[., 2] = (z2 :> 0):*ln(z2 + (z2 :<= 0))
        F[., 3] = (z1 :<= 0)
        F[., 4] = (z2 :<= 0)
    }
    return(F)
}

// Pieces of the linear predictor of one equation: NL (N x k), the part that
// does not depend on the lags, intercepts included; G (4 x k), the
// coefficients of the four lag features in each state. cc and sc: codes of
// the common and of the state-specific terms (0 regressor, 1-4 lag feature)
void hmbp_lagparts(real matrix X, real matrix XS, real rowvector bc,
    real rowvector a, real matrix S, real rowvector cc, real rowvector sc,
    real scalar k, real matrix NL, real matrix G)
{
    real rowvector s0
    real scalar    m, c

    NL = J(rows(X), 1, a)
    if (cols(bc) > 0) {
        s0 = selectindex(cc :== 0)
        if (cols(s0) > 0) NL = NL :+ X[., s0]*bc[s0]'
    }
    if (rows(S) > 0) {
        s0 = selectindex(sc :== 0)
        if (cols(s0) > 0) NL = NL + XS[., s0]*S[s0, .]
    }
    G = J(4, k, 0)
    for (m=1; m<=4; m++) {
        if (cols(bc) > 0) G[m, .] = J(1, k, sum(bc :* (cc :== m)))
        if (rows(S) > 0) {
            for (c=1; c<=k; c++) G[m, c] = G[m, c] + sum(S[., c]' :* (sc :== m))
        }
    }
}

// smallest y with F(y; mu) >= u, element by element
real colvector hmbp_qpois(real colvector u0, real colvector mu)
{
    real colvector u, y, cum, pr, todo, hit
    real scalar    j

    u    = rowmin((u0, J(rows(u0), 1, 1 - 1e-12)))
    y    = J(rows(mu), 1, .)
    cum  = J(rows(mu), 1, 0)
    todo = J(rows(mu), 1, 1)
    pr   = exp(-mu)
    for (j=0; j<=20000; j++) {
        cum  = cum + pr
        hit  = todo :& (cum :>= u)
        if (any(hit)) y[selectindex(hit)] = J(sum(hit), 1, j)
        todo = todo :& (hit :== 0)
        if (!any(todo)) break
        pr   = pr:*mu:/(j+1)
    }
    if (any(todo)) y[selectindex(todo)] = J(sum(todo), 1, j)
    return(y)
}

// one draw of the pair of counts for each row, given the means and the
// dependence parameter of each row
void hmbp_draw(real scalar dep, real colvector mu1, real colvector mu2,
    real colvector lam, real colvector y1, real colvector y2)
{
    real colvector u1, u2, w, e1, e2, r, cum, todo, hit, pr, d1, uu, T, sm
    real scalar    q, y

    if (dep == 0) {
        // y2 is Poisson; y1 | y2 by inversion of its conditional pmf
        q    = 1 - exp(-1)
        y2   = rpoisson(1, 1, mu2)
        u2   = exp(-y2) - exp(-q*mu2)
        d1   = exp(-q*mu1)
        uu   = runiform(rows(mu1), 1)
        y1   = J(rows(mu1), 1, .)
        cum  = J(rows(mu1), 1, 0)
        todo = J(rows(mu1), 1, 1)
        pr   = exp(-mu1)
        for (y=0; y<=20000; y++) {
            cum  = cum + pr:*rowmax((1 :+ lam:*(exp(-y) :- d1):*u2, J(rows(mu1), 1, 0)))
            hit  = todo :& (cum :>= uu)
            if (any(hit)) y1[selectindex(hit)] = J(sum(hit), 1, y)
            todo = todo :& (hit :== 0)
            if (!any(todo)) break
            pr   = pr:*mu1:/(y+1)
        }
        if (any(todo)) y1[selectindex(todo)] = J(sum(todo), 1, y)
        return
    }
    if (dep == 1) {
        // normal copula: correlated standard normals
        r  = tanh(lam)
        e1 = rnormal(rows(mu1), 1, 0, 1)
        e2 = rnormal(rows(mu1), 1, 0, 1)
        u1 = normal(e1)
        u2 = normal(r:*e1 + sqrt(1 :- r:^2):*e2)
    }
    else {
        // Frank copula: u2 from the conditional distribution given u1
        u1 = runiform(rows(mu1), 1)
        w  = runiform(rows(mu1), 1)
        sm = abs(lam) :< 1e-6
        T  = lam + sm
        u2 = -hmbp_ln1p(w:*hmbp_expm1(-T):/(w + (1:-w):*exp(-T:*u1))):/T
        u2 = sm:*w + (1:-sm):*u2
    }
    y1 = hmbp_qpois(u1, mu1)
    y2 = hmbp_qpois(u2, mu2)
}

// p-th quantile of a column (the smallest value with at least a share p of
// the values below or at it)
real scalar hmbp_pctile(real colvector x, real scalar p)
{
    real colvector s

    s = sort(x, 1)
    return(s[max((1, ceil(p*rows(s))))])
}

// Stability of the dynamics. Results in Stata matrices whose names start with
// stub: stubG (k x 4, coefficients of the lags of y1 and y2 in the two
// equations), stubR (k x 2, spectral radius and norm of the 2 x 2 matrix of
// the coefficients of g(y1), g(y2)), stubU (k x 2, index gamma_jj
// exp(q99 of the rest of the linear predictor), equations 1 and 2), stubD
// (k x 1, share of deterministic paths that exceed the limit), stubS
// (simulated paths: share that exceeds the limit, and the quantiles of the
// expected counts at the horizon and in the sample).
// Options from the locals of the caller: horizon, limit, nsim, lagt, lagarg,
// rawlags (the lagged counts), ccode1, scode1, ccode2, scode2.
void hmbp_stab(struct hmbp_dat scalar M, real rowvector b, string scalar stub,
    real scalar seed, string scalar touse)
{
    struct hmbp_mm scalar mm
    real rowvector b1, a1, b2, a2, lam, tau, cc1, sc1, cc2, sc2, qq
    real matrix    S1, S2, D, W, Zs, Pi0, P, XiI, NL1, NL2, G1, G2, Lr, F, Gm
    real matrix    R, U, PT, PTi, Pm, Ri, Uk, Q
    real colvector z1, z2, m1, m2, dv, ok, last, s, nl1, nl2, y1, y2, ex
    real colvector lams, ind, trw, mx1, mx2, fit1, fit2, e
    real rowvector sv1, sv2
    real scalar    k, H, lim, nsim, larg, ll, c, h, i, n, tr, dt, chk
    string scalar  lagt, rs

    k    = M.k
    H    = strtoreal(st_local("horizon"))
    lim  = strtoreal(st_local("limit"))
    nsim = strtoreal(st_local("nsim"))
    lagt = st_local("lagt")
    larg = strtoreal(st_local("lagarg"))
    cc1  = strtoreal(tokens(st_local("ccode1")))
    sc1  = strtoreal(tokens(st_local("scode1")))
    cc2  = strtoreal(tokens(st_local("ccode2")))
    sc2  = strtoreal(tokens(st_local("scode2")))
    if (cols(cc1) == 0) cc1 = J(1, 0, .)
    if (cols(sc1) == 0) sc1 = J(1, 0, .)
    if (cols(cc2) == 0) cc2 = J(1, 0, .)
    if (cols(sc2) == 0) sc2 = J(1, 0, .)
    Lr   = st_data(., tokens(st_local("rawlags")), touse)

    ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 1, 0, XiI)
    if (ll >= .) {
        errprintf("the log likelihood cannot be evaluated at e(b) with the data in memory\n")
        exit(459)
    }
    hmbp_unpack(M, b, b1, a1, S1, b2, a2, S2, lam, D, tau)
    hmbp_lagparts(M.X1, M.X1S, b1, a1, S1, cc1, sc1, k, NL1, G1)
    hmbp_lagparts(M.X2, M.X2S, b2, a2, S2, cc2, sc2, k, NL2, G2)
    // the pieces must give back the fitted means at the observed lags
    F   = hmbp_lagf(Lr[., 1], Lr[., 2], lagt, larg)
    chk = max((mreldif(exp(NL1 + F*G1), mm.Mu1), mreldif(exp(NL2 + F*G2), mm.Mu2)))
    if (chk > 1e-8) {
        errprintf("hmbp_stab(): the linear predictor is not rebuilt (%g)\n", chk)
        exit(498)
    }

    // coefficients of the lags of the two counts, by state
    Gm = (G1[1, .]', G1[2, .]', G2[1, .]', G2[2, .]')
    st_matrix(stub + "G", Gm)

    // spectral radius and norm of [g11 g12 \ g21 g22] in each state
    R = J(k, 2, .)
    for (c=1; c<=k; c++) {
        Q  = (Gm[c, 1], Gm[c, 2] \ Gm[c, 3], Gm[c, 4])
        tr = Q[1,1] + Q[2,2]
        dt = Q[1,1]*Q[2,2] - Q[1,2]*Q[2,1]
        if (tr^2/4 - dt >= 0) R[c, 1] = max(abs((tr/2 + sqrt(tr^2/4 - dt), tr/2 - sqrt(tr^2/4 - dt))))
        else                  R[c, 1] = sqrt(dt)
        R[c, 2] = sqrt(max(symeigenvalues(Q'*Q)))
    }
    st_matrix(stub + "R", R)

    // index of the univariate map y -> exp(a + g y): explosive above exp(-1)
    U = J(k, 2, .)
    for (c=1; c<=k; c++) {
        U[c, 1] = Gm[c, 1]*exp(hmbp_pctile(NL1[., c], 0.99))
        U[c, 2] = Gm[c, 4]*exp(hmbp_pctile(NL2[., c], 0.99))
    }
    st_matrix(stub + "U", U)

    // deterministic paths: in each state, from the lags and with the
    // covariates of each observation, the counts are replaced by their
    // expected values for H periods
    D = J(k, 1, .)
    for (c=1; c<=k; c++) {
        z1 = Lr[., 1]
        z2 = Lr[., 2]
        dv = J(M.N, 1, 0)
        for (h=1; h<=H; h++) {
            F  = hmbp_lagf(z1, z2, lagt, larg)
            m1 = exp(NL1[., c] + F*G1[., c])
            m2 = exp(NL2[., c] + F*G2[., c])
            ok = (m1 :<= lim) :& (m2 :<= lim)
            dv = dv :| (ok :== 0)
            z1 = ok:*m1 + (1:-ok)*lim
            z2 = ok:*m2 + (1:-ok)*lim
        }
        D[c] = mean(dv)
    }
    st_matrix(stub + "D", D)

    // simulated paths: nsim per individual, from the last observed period,
    // with its covariates, the state drawn from its posterior probabilities
    // and the estimated chain
    rs   = rngstate()
    rseed(seed)
    n    = M.n
    last = M.info[., 2]
    ind  = J(nsim, 1, 1) # (1::n)
    last = last[ind]
    nl1  = NL1[last, .]
    nl2  = NL2[last, .]
    z1   = M.y1[last]
    z2   = M.y2[last]
    Uk   = uppertriangle(J(k, k, 1))
    s    = 1 :+ rowsum(runiform(rows(ind), 1) :> (W[last, .]*Uk))
    s    = rowmin((s, J(rows(ind), 1, k)))
    if (M.het) {
        // transition probabilities of the last transition row of each
        // individual; the mean over all the rows for those observed once
        PT  = P
        Pm  = mean(PT)
        trw = last :- ind
        PTi = J(rows(ind), k*k, .)
        for (i=1; i<=rows(ind); i++) {
            if (M.Ts[ind[i]] > 1) PTi[i, .] = PT[trw[i], .]
            else                  PTi[i, .] = Pm
        }
    }
    ex  = J(rows(ind), 1, 0)
    mx1 = J(rows(ind), 1, 0)
    mx2 = J(rows(ind), 1, 0)
    for (h=1; h<=H; h++) {
        F  = hmbp_lagf(z1, z2, lagt, larg)
        m1 = exp(rowsum(nl1 :* (J(1, k, s) :== (1..k))) + rowsum(F :* G1[., s]'))
        m2 = exp(rowsum(nl2 :* (J(1, k, s) :== (1..k))) + rowsum(F :* G2[., s]'))
        ok = (m1 :<= lim) :& (m2 :<= lim)
        ex = ex :| (ok :== 0)
        lams = (lam[s'])'
        // paths above the limit are not drawn (mean 1 in place of theirs)
        hmbp_draw(M.dep, ok:*m1 + (1:-ok), ok:*m2 + (1:-ok), lams, y1, y2)
        // paths that exceeded the limit stay at the limit
        z1 = ok:*y1 + (1:-ok)*lim
        z2 = ok:*y2 + (1:-ok)*lim
        if (h == H) {
            mx1 = m1
            mx2 = m2
        }
        // next state
        e = runiform(rows(ind), 1)
        if (M.het) {
            Ri = J(rows(ind), k, 0)
            for (c=1; c<=k; c++) Ri = Ri + (s :== c):*PTi[., ((c-1)*k+1)..(c*k)]
        }
        else Ri = P[s, .]
        s = 1 :+ rowsum(e :> (Ri*Uk))
        s = rowmin((s, J(rows(ind), 1, k)))
    }
    rngstate(rs)
    // quantiles of the expected counts at the horizon (paths below the limit)
    // and of the fitted expected counts in the sample (posterior weights)
    qq   = (0.5, 0.9, 0.99)
    ok   = selectindex(ex :== 0)
    sv1  = J(1, 4, .)
    sv2  = J(1, 4, .)
    if (rows(ok) > 0) {
        for (i=1; i<=3; i++) {
            sv1[i] = hmbp_pctile(mx1[ok], qq[i])
            sv2[i] = hmbp_pctile(mx2[ok], qq[i])
        }
        sv1[4] = max(mx1[ok])
        sv2[4] = max(mx2[ok])
    }
    fit1 = rowsum(W:*mm.Mu1)
    fit2 = rowsum(W:*mm.Mu2)
    st_matrix(stub + "S", (mean(ex), rows(ind), H, lim))
    st_matrix(stub + "Q", (sv1 \ sv2 \ (hmbp_pctile(fit1, 0.5), hmbp_pctile(fit1, 0.9), hmbp_pctile(fit1, 0.99), max(fit1)) \ (hmbp_pctile(fit2, 0.5), hmbp_pctile(fit2, 0.9), hmbp_pctile(fit2, 0.99), max(fit2))))
}

// ----------------------------------------------------------------- estat ame
// parts of an interaction term, split at #
string rowvector hmbp_split(string scalar t)
{
    string rowvector P
    string scalar    rest
    real scalar      pos

    P    = J(1, 0, "")
    rest = t
    pos  = strpos(rest, "#")
    while (pos > 0) {
        P    = (P, substr(rest, 1, pos-1))
        rest = substr(rest, pos+1, .)
        pos  = strpos(rest, "#")
    }
    return((P, rest))
}

// variables of a list of terms (lags excluded): st_local amev (names), amef
// (1 = factor variable) and amelv1, amelv2, ... (levels of each factor
// variable that appear in the terms)
void hmbp_amevars(string scalar terms, string scalar lagn)
{
    string rowvector T, L, P, vars, lev
    real rowvector   fac, w
    string scalar    p, pre, nm, lv
    real scalar      i, j, wi, isf

    T    = tokens(terms)
    L    = tokens(lagn)
    vars = J(1, 0, "")
    lev  = J(1, 0, "")
    fac  = J(1, 0, .)
    for (i=1; i<=cols(T); i++) {
        if (anyof(L, T[i])) continue
        P = hmbp_split(T[i])
        for (j=1; j<=cols(P); j++) {
            p = P[j]
            pre = ""
            nm  = p
            if (regexm(p, "^(.*)\.([^.]+)$")) {
                pre = regexs(1)
                nm  = regexs(2)
            }
            isf = regexm(pre, "^[0-9]+")
            lv  = (isf ? regexs(0) : "")
            w   = selectindex(vars :== nm)
            if (cols(w) == 0) {
                vars = (vars, nm)
                fac  = (fac, isf)
                lev  = (lev, "")
                wi   = cols(vars)
            }
            else wi = w[1]
            if (isf) {
                fac[wi] = 1
                lev[wi] = strtrim(lev[wi] + " " + lv)
            }
        }
    }
    st_local("amev", invtokens(vars))
    st_local("amef", invtokens(strofreal(fac)))
    for (i=1; i<=cols(vars); i++) {
        st_local("amelv" + strofreal(i), invtokens(uniqrows(tokens(lev[i])')'))
    }
}

// a term with the variable v replaced by r (t is a term as in e(common1):
// var, c.var, #.var, interactions joined by #)
string scalar hmbp_subst(string scalar t, string scalar v, string scalar r)
{
    string rowvector P
    string scalar    out, p
    real scalar      i

    P   = hmbp_split(t)
    out = ""
    for (i=1; i<=cols(P); i++) {
        p = P[i]
        if (regexm(p, "^(.*)\.([^.]+)$")) {
            if (regexs(2) == v) p = regexs(1) + "." + r
        }
        else if (p == v) p = r
        out = out + (i > 1 ? "#" : "") + p
    }
    return(out)
}

// linear predictor of one equation, N x k
real matrix hmbp_eta(real matrix X, real matrix XS, real rowvector bc,
    real rowvector a, real matrix S)
{
    real matrix eta

    eta = J(rows(X), 1, a)
    if (cols(bc) > 0) eta = eta :+ X*bc'
    if (rows(S) > 0)  eta = eta + XS*S
    return(eta)
}

// design matrix X with the columns of the lag features replaced by F
// (codes: 0 regressor, 1-4 column of F)
real matrix hmbp_lagrep(real matrix X, real rowvector cd, real matrix F)
{
    real matrix XN
    real scalar j

    XN = X
    for (j=1; j<=cols(cd); j++) {
        if (cd[j] > 0) XN[., j] = F[., cd[j]]
    }
    return(XN)
}

// columns o+1..o+w of X (N x 0 when w = 0)
real matrix hmbp_cols(real matrix X, real scalar o, real scalar w)
{
    if (w == 0) return(J(rows(X), 0, .))
    return(X[|1, o+1 \ ., o+w|])
}

// average marginal effects, equation 1 then equation 2 (1 x 2*nE): for each
// effect e, mean over the observations of sum_c Pr_c (mu_c(A) - mu_c(B)) / div,
// with the designs A and B side by side in XX: block e has the columns of
// eq. 1 common, eq. 1 by state, eq. 2 common, eq. 2 by state, for A then for B
real rowvector hmbp_amef(struct hmbp_dat scalar M, real rowvector b,
    real matrix XX, real rowvector div, real scalar ptype, real scalar ist)
{
    struct hmbp_mm scalar mm
    real rowvector b1, a1, b2, a2, lam, tau, wd
    real matrix    S1, S2, D, W, Zs, Pi0, P, XiI, Pr, out, mA, mB
    real scalar    ll, e, o, wb

    if (ist > 0) {
        Pr = J(M.N, M.k, 0)
        Pr[., ist] = J(M.N, 1, 1)
    }
    else {
        ll = hmbp_estep(M, b, W, Zs, Pi0, P, mm, 1, 0, XiI)
        if (ll >= .) return(J(1, 2*cols(div), .))
        if (ptype == 0) Pr = W
        else            Pr = hmbp_sprob(M, mm.lf, Pi0, P, ptype)
    }
    hmbp_unpack(M, b, b1, a1, S1, b2, a2, S2, lam, D, tau)
    wd  = (M.p1, M.s1, M.p2, M.s2)
    wb  = 2*sum(wd)
    out = J(2, cols(div), .)
    for (e=1; e<=cols(div); e++) {
        o  = wb*(e-1)
        mA = exp(hmbp_eta(hmbp_cols(XX, o, wd[1]), hmbp_cols(XX, o+wd[1], wd[2]), b1, a1, S1))
        mB = exp(hmbp_eta(hmbp_cols(XX, o+wb/2, wd[1]), hmbp_cols(XX, o+wb/2+wd[1], wd[2]), b1, a1, S1))
        out[1, e] = mean(rowsum(Pr:*(mA - mB)))/div[e]
        o  = o + wd[1] + wd[2]
        mA = exp(hmbp_eta(hmbp_cols(XX, o, wd[3]), hmbp_cols(XX, o+wd[3], wd[4]), b2, a2, S2))
        mB = exp(hmbp_eta(hmbp_cols(XX, o+wb/2, wd[3]), hmbp_cols(XX, o+wb/2+wd[3], wd[4]), b2, a2, S2))
        out[2, e] = mean(rowsum(Pr:*(mA - mB)))/div[e]
    }
    return((out[1, .], out[2, .]))
}

// estat ame: designs of the effects from the locals of the caller (ame_n,
// ame_t#, ame_div#, ame_[ab][12][common|statevary]#), estimates and their
// delta-method variance (Jacobian by central differences). Results: matrix
// stub (estimates) and stub + "V" (variance)
void hmbp_ame(struct hmbp_dat scalar M, real rowvector b, string scalar stub,
    real scalar ptype, real scalar ist, string scalar touse)
{
    real rowvector div, f, fp, fq, bp, bq, cc1, sc1, cc2, sc2, fr
    real matrix    Lr, FA, G, V, XX, XA, XB
    real scalar    nE, e, s, j, h, i
    string scalar  typ, sd, lagt
    string rowvector kd

    nE   = strtoreal(st_local("ame_n"))
    lagt = st_local("lagt")
    cc1  = strtoreal(tokens(st_local("ccode1")))
    sc1  = strtoreal(tokens(st_local("scode1")))
    cc2  = strtoreal(tokens(st_local("ccode2")))
    sc2  = strtoreal(tokens(st_local("scode2")))
    if (cols(cc1) == 0) cc1 = J(1, 0, .)
    if (cols(sc1) == 0) sc1 = J(1, 0, .)
    if (cols(cc2) == 0) cc2 = J(1, 0, .)
    if (cols(sc2) == 0) sc2 = J(1, 0, .)
    Lr   = st_data(., tokens(st_local("rawlags")), touse)
    div  = J(1, nE, .)
    kd   = ("1common", "1statevary", "2common", "2statevary")
    XX   = J(M.N, 0, .)
    for (e=1; e<=nE; e++) {
        typ    = st_local("ame_t" + strofreal(e))
        div[e] = strtoreal(st_local("ame_div" + strofreal(e)))
        if (typ == "x") {
            for (s=1; s<=2; s++) {
                sd = (s == 1 ? "a" : "b")
                for (j=1; j<=4; j++) {
                    XX = (XX, hmbp_getx(st_local("ame_" + sd + kd[j] + strofreal(e)), touse, M.N))
                }
            }
        }
        else {
            // one more count of the lagged outcome
            if (typ == "lag1") FA = hmbp_lagf(Lr[., 1] :+ 1, Lr[., 2], lagt, strtoreal(st_local("lagarg")))
            else               FA = hmbp_lagf(Lr[., 1], Lr[., 2] :+ 1, lagt, strtoreal(st_local("lagarg")))
            XA = (hmbp_lagrep(M.X1, cc1, FA), hmbp_lagrep(M.X1S, sc1, FA), hmbp_lagrep(M.X2, cc2, FA), hmbp_lagrep(M.X2S, sc2, FA))
            XB = (M.X1, M.X1S, M.X2, M.X2S)
            XX = (XX, XA, XB)
        }
    }
    f = hmbp_amef(M, b, XX, div, ptype, ist)
    if (hasmissing(f)) {
        errprintf("the marginal effects cannot be computed at e(b) with the data in memory\n")
        exit(459)
    }
    // Jacobian with respect to the parameters with nonzero variance
    V  = st_matrix("e(V)")
    fr = selectindex(diagonal(V)' :!= 0)
    G  = J(cols(f), cols(fr), 0)
    for (i=1; i<=cols(fr); i++) {
        j  = fr[i]
        h  = 1e-4*(abs(b[j]) + 1)
        bp = b; bp[j] = b[j] + h
        bq = b; bq[j] = b[j] - h
        fp = hmbp_amef(M, bp, XX, div, ptype, ist)
        fq = hmbp_amef(M, bq, XX, div, ptype, ist)
        G[., i] = (fp - fq)' :/ (2*h)
    }
    V = G*V[fr, fr]*G'
    st_matrix(stub, f)
    st_matrix(stub + "V", (V + V')/2)
}

// keeps the elements idx of the row vector bn and the rows and columns idx of
// the matrix vn (Stata matrices)
void hmbp_msel(string scalar bn, string scalar vn, string scalar idx)
{
    real rowvector k, bb
    real matrix    V

    k  = strtoreal(tokens(idx))
    bb = st_matrix(bn)
    V  = st_matrix(vn)
    st_matrix(bn, bb[k])
    st_matrix(vn, V[k, k])
}

// Stata expression of a term (see TermVar); "" for operators other than
// level indicators and c.
string scalar hmbp_termexpr(string scalar t)
{
    string rowvector P
    string scalar    ex, p, pre, nm, f
    real scalar      i

    P  = hmbp_split(t)
    ex = ""
    for (i=1; i<=cols(P); i++) {
        p = P[i]
        if (regexm(p, "^(.*)\.([^.]+)$")) {
            pre = regexs(1)
            nm  = regexs(2)
            if (pre == "c") f = nm
            else if (regexm(pre, "^[0-9]+$")) f = "cond(missing(" + nm + "),.," + nm + "==" + pre + ")"
            else return("")
        }
        else f = p
        ex = ex + (i > 1 ? "*" : "") + "(" + f + ")"
    }
    return(ex)
}

end
