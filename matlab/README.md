# Inverse binomial sampling (IBS) for MATLAB

This folder contains MATLAB implementations and examples of IBS; see the [home page](https://github.com/lacerbi/ibs) for general information.

## Code

- `ibs_example.m` is an annotated tutorial with a full working example usage of IBS.
  - In order to run the tutorial, you will need to have installed Bayesian Adaptive Direct Search ([BADS](https://github.com/lacerbi/bads)). 
- `ibslike.m` is an advanced vectorized implementation of IBS, which supports several advanced features (still work in progress as we add more advanced features). Please read the documentation in the file for usage information. 
  - Note that `ibslike` returns the *negative* log-likelihood as it is meant to be used with an optimization method such as BADS.
  - Run `ibslike('test')` for a battery of unit tests.
  - If you want to run `ibslike` with [Variational Bayesian Monte Carlo](https://github.com/lacerbi/vbmc) (VBMC), note that you need to set the IBS options as
    - `OPTIONS.ReturnPositive = true` to return the *positive* log-likelihood;
    - `OPTIONS.ReturnStd = true` to return as second output the standard deviation (SD) of the estimate (as opposed to the variance).
- `ibs_basic.m` is a bare-bone implementation of IBS for didactic purposes.
- `psycho_gen.m` and `psycho_nll.m` are functions implementing, respectively, the generative model (simulator) and the negative log-likelihood function for the orientation discrimination model used as example in the tutorial (see also Section 5.2 in the IBS paper).
