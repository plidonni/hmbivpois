# hmbivpois

**Hidden Markov bivariate Poisson regression for panel counts in Stata**

`hmbivpois` fits the hidden Markov bivariate Poisson model of Laudicella and Li Donni (2022) to two counts observed on the same individuals over time, for example visits to primary care and emergency admissions. Given a latent state that follows a first-order Markov chain, the two counts have Poisson margins whose means depend on regressors, on the lagged counts and on the state.

## Installation

```stata
net install hmbivpois, from("https://raw.githubusercontent.com/plidonni/hmbivpois/main/") replace
```

To download the example dataset and the do-file with the examples into the current folder:

```stata
net get hmbivpois, from("https://raw.githubusercontent.com/plidonni/hmbivpois/main/") replace
```

Requires Stata 16 or later. Then type `help hmbivpois`.

## Quick start

```stata
use hmbivpois_example, clear
hmbivpois (y1 = x11 x12) (y2 = x21 x22), id(id) time(time) k(3) initcov(z1)
estat ic
estat transition
estat ame
predict double pp*, postprob
```

The data are in long form, one row per individual and period. The first observed period of each individual is the presample period: it provides the lagged counts of the first estimation period and the presample counts used in the initial probabilities of the states.

## Features

- **Dependence between the two counts:** Famoye form (default), Frank copula or normal copula, with one parameter or one per state (`lambdastate`).
- **Lagged counts:** linear, `log1p`, `logc(#)`, `logzero`, indicator or top-coded (`lagtransform()`).
- **State-specific slopes** for chosen regressors and lags (`statevary()`, `lagvary()`).
- **Latent chain:** homogeneous transitions, transitions depending on covariates (`xtrans()`), transitions free by period (`tvtrans`), or no transitions (latent class model, `notransition`).
- **Estimation:** EM algorithm with scaled forward-backward recursions, followed by direct maximization with the analytic score; observed-information or cluster-robust standard errors.
- **Postestimation:** `predict` (expected counts, posterior/filtered/one-step/prior state probabilities, most likely state, Viterbi path, within-state correlation, joint probabilities); `estat ic`, `estat transition`, `estat lrtrans` (chi-bar-squared test of no transitions), `estat stability` (explosive dynamics of the lags), `estat ame` (average marginal effects).

## Files

| File | Content |
|---|---|
| `hmbivpois.ado` | the command (estimation, `predict` and `estat` routines) |
| `hmbivpois_p.ado`, `hmbivpois_estat.ado` | entry points for `predict` and `estat` |
| `hmbivpois.sthlp` | help file |
| `hmbivpois_example.dta` | simulated example data |
| `hmbivpois_examples.do` | runs the examples of the help file |
| `stata.toc`, `hmbivpois.pkg` | installation files for `net install` / `net get` |

## Example data

`hmbivpois_example.dta` is simulated from the model with 1,000 individuals, a presample period and 6 estimation periods, three latent states, the Famoye form and homogeneous transitions. The variable `state` holds the true state. True values:

| Block | Parameters |
|---|---|
| Equation y1 | x11 = 0.4, x12 = −0.3, lag y1 = 0.03, lag y2 = −0.05, state intercepts (−0.3, 0.7, 1.5) |
| Equation y2 | x21 = 0.3, x22 = −0.4, lag y1 = −0.10, lag y2 = 0.10, state intercepts (−2.5, −1.5, −0.6) |
| Dependence | λ = −0.5 |
| Initial probabilities (states 2, 3 vs 1) | z1 = (0.5, 1.0), presample y1 = (0.15, 0.30), presample y2 = (0.30, 0.50), constant = (−0.5, −2.0) |
| Transition matrix | (0.85, 0.12, 0.03; 0.08, 0.84, 0.08; 0.03, 0.12, 0.85) |

## Citation

If you use `hmbivpois`, please cite:

- Laudicella, M., and P. Li Donni. 2022. The dynamic interdependence in the demand of primary and emergency secondary care: A hidden Markov approach. *Journal of Applied Econometrics* 37: 521–536.
- Li Donni, P. 2026. hmbivpois: Stata module to fit hidden Markov bivariate Poisson regression for panel counts. https://github.com/plidonni/hmbivpois.

## Development

The model and its extensions were specified by the author. Parts of the Stata and Mata code were written with the assistance of an AI coding tool (Claude, Anthropic); the code was tested by the author with Monte Carlo simulations.

## License

MIT License, see `LICENSE`.

## Author

Paolo Li Donni. Please report bugs and suggestions as [issues](https://github.com/plidonni/hmbivpois/issues).
