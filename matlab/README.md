# Inverse binomial sampling (IBS) for MATLAB

This folder contains MATLAB implementations and examples of IBS; see the [home page](https://github.com/lacerbi/ibs) for general information.

## Code

- `ibs_example.m` is an annotated tutorial with a full working example usage of IBS.
  - In order to run the tutorial, you will need to have installed Bayesian Adaptive Direct Search ([BADS](https://github.com/lacerbi/bads)). 
- `ibslike.m` is an advanced vectorized implementation of IBS, which supports several advanced features (still work in progress). Read the documentation in the file for usage information. 
  - Note that `ibslike` returns the *negative* log-likelihood as it is meant to be used with an optimization method such as BADS.
- `ibs_basic.m` is a bare-bone implementation of IBS for didactic purposes.
