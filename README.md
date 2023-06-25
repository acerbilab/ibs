# Inverse binomial sampling (IBS)

This repository contains MATLAB implementations and examples of IBS. For Python implementations, see [PyIBS](https://github.com/acerbilab/pyibs).

## What is it?

Inverse binomial sampling (IBS) is a technique to obtain unbiased, efficient estimates of the log-likelihood of a model by simulation. [[1](#references)]

The typical scenario is the case in which you have a *simulator*, that is a model from which you can randomly draw synthetic observations (for a given parameter vector), but cannot evaluate the log-likelihood analytically or numerically. In other words, IBS affords likelihood-based inference for models without explicit likelihood functions (also known as implicit models).

## Quickstart


The main function is `ibslike.m`, which computes the estimate of the negative log-likelihood for a given simulator model and dataset via IBS.
We recommend to start with the tutorial in `ibs_example.m`, which contains a full walkthrough with working example usages of IBS.

IBS is commonly used as a part of an algorithm for maximum-likelihood estimation or Bayesian inference:
- For maximum-likelihood (or maximum-a-posteriori) estimation, we recommend to use IBS combined with [Bayesian Adaptive Direct Search (BADS)](https://github.com/acerbilab/bads). BADS is required to run the first part of the tutorial.
- For Bayesian inference of posterior distributions and model evidence, we recommend to use IBS with [Variational Bayesian Monte Carlo (VBMC)](https://github.com/acerbilab/vbmc). VBMC is required to run the second part of the tutorial.

For practical recommendations and any other question, check out the FAQ on the [IBS wiki](https://github.com/acerbilab/ibs/wiki).

## Code

We describe below the files in this repository:

- `ibs_basic.m` is a bare-bone implementation of IBS for didactic purposes.
- `ibslike.m` is an advanced vectorized implementation of IBS, which supports several advanced features (still work in progress as we add more features). Please read the documentation in the file for usage information, and refer to the [dedicated section on the FAQ](https://github.com/acerbilab/ibs/wiki#matlab-implementation-ibslike) for specific questions about this MATLAB implementation.
  - Note that `ibslike` returns the *negative* log-likelihood as it is meant to be used with an optimization method such as BADS.
  - Run `ibslike('test')` for a battery of unit tests.
  - If you want to run `ibslike` with [Variational Bayesian Monte Carlo](https://github.com/acerbilab/vbmc) (VBMC), note that you need to set the IBS options as
    - `OPTIONS.ReturnPositive = true` to return the *positive* log-likelihood;
    - `OPTIONS.ReturnStd = true` to return as second output the standard deviation (SD) of the estimate (as opposed to the variance).
- `ibs_example.m` is a tutorial with usage cases for IBS.
- `psycho_gen.m` and `psycho_nll.m` are functions implementing, respectively, the generative model (simulator) and the negative log-likelihood function for the orientation discrimination model used as example in the tutorial (see also Section 5.2 in the IBS paper).

The code used to produce results in the paper [[1](#references)] is available in the development repository [here](https://github.com/basvanopheusden/ibs-development).

## References

1. van Opheusden\*, B., Acerbi\*, L. & Ma, W.J. (2020). Unbiased and efficient log-likelihood estimation with inverse binomial sampling. *PLoS Computational Biology* 16(12): e1008483. (\* equal contribution) ([link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008483)) 

You can cite IBS in your work with something along the lines of

> We estimated the log-likelihood using inverse binomial sampling (IBS; van Opheusden, Acerbi & Ma, 2020), a technique that produces unbiased and efficient estimates of the log-likelihood via simulation. 

If you use IBS in conjunction with [Bayesian Adaptive Direct Search](https://github.com/acerbilab/bads), as recommended in the paper, you could add

> We obtained maximum-likelihood estimates of the model parameters via Bayesian Adaptive Direct Search (BADS; Acerbi & Ma, 2017), a hybrid Bayesian optimization algorithm which affords stochastic objective evaluations.

and cite the appropriate paper:

2. Acerbi, L. & Ma, W. J. (2017). Practical Bayesian optimization for model fitting with Bayesian Adaptive Direct Search. In *Advances in Neural Information Processing Systems* 30:1834-1844.

Similarly, if you use IBS in combination with [Variational Bayesian Monte Carlo](https://github.com/acerbilab/vbmc), you should cite these papers:

3. Acerbi, L. (2018). Variational Bayesian Monte Carlo. In *Advances in Neural Information Processing Systems* 31: 8222-8232.
4. Acerbi, L. (2020). Variational Bayesian Monte Carlo with Noisy Likelihoods. In *Advances in Neural Information Processing Systems* 33: 8211-8222.

Besides formal citations, you can demonstrate your appreciation for our work in the following ways:

- *Star* the IBS repository on GitHub;
- Follow us on Twitter ([Luigi](https://twitter.com/AcerbiLuigi), [Bas](https://twitter.com/basvanopheusden)) for updates about IBS and other projects we are involved with;
- Tell us about your model-fitting problem and your experience with IBS (positive or negative) in the [Discussions forum](https://github.com/orgs/acerbilab/discussions).

### BibTex

```
@article{vanOpheusden2020unbiased,
  title = {Unbiased and Efficient Log-Likelihood Estimation with Inverse Binomial Sampling},
  author = {van Opheusden, Bas and Acerbi, Luigi and Ma, Wei Ji},
  year = {2020},
  journal = {PLOS Computational Biology},
  volume = {16},
  number = {12},
  pages = {e1008483},
  publisher = {{Public Library of Science}},
  issn = {1553-7358},
  doi = {10.1371/journal.pcbi.1008483},
}
```

### License

The IBS code is released under the terms of the [MIT License](https://github.com/acerbilab/ibs/blob/master/LICENSE.txt).
