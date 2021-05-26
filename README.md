# Inverse binomial sampling (IBS)

## What is it

IBS is a technique to obtain unbiased, efficient estimates of the log-likelihood of a model by simulation. [[1](#references)]

The typical scenario is the case in which you have a *simulator*, that is a model from which you can randomly draw synthetic observations (for a given parameter vector), but cannot evaluate the log-likelihood analytically or numerically. In other words, IBS affords likelihood-based inference for models without explicit likelihood functions (also known as implicit models).

IBS is commonly used as a part of an algorithm for maximum-likelihood estimation or Bayesian inference.

- For maximum-likelihood (or maximum-a-posteriori) estimation, we recommend to use IBS combined with [Bayesian Adaptive Direct Search (BADS)](https://github.com/lacerbi/bads).
- For Bayesian inference of posterior distributions and model evidence, we recommend to use IBS with [Variational Bayesian Monte Carlo (VBMC)](https://github.com/lacerbi/vbmc).

## Code

This repository stores basic and advanced implementations and example usages of IBS in various programming languages for users of the method. At the moment, we have a MATLAB implementation and we plan to include other ones (e.g., Python).

- Go to the [MATLAB code page](https://github.com/lacerbi/ibs/tree/master/matlab).

The code used to produce results in the paper [[1](#references)] is available in the development repository [here](https://github.com/basvanopheusden/ibs-development).

For practical recommendations and any other question, check out the FAQ on the [IBS wiki](https://github.com/lacerbi/ibs/wiki).

## References

1. van Opheusden\*, B., Acerbi\*, L. & Ma, W.J. (2020). Unbiased and efficient log-likelihood estimation with inverse binomial sampling. *PLoS Computational Biology* 16(12): e1008483. (\* equal contribution) ([link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008483)) 

You can cite IBS in your work with something along the lines of

> We estimated the log-likelihood using inverse binomial sampling (IBS; van Opheusden, Acerbi & Ma, 2020), a technique that produces unbiased and efficient estimates of the log-likelihood via simulation. 

If you use IBS in conjunction with [Bayesian Adaptive Direct Search](https://github.com/lacerbi/bads), as recommended in the paper, you could add

> We obtained maximum-likelihood estimates of the model parameters via Bayesian Adaptive Direct Search (BADS; Acerbi & Ma, 2017), a hybrid Bayesian optimization algorithm which affords stochastic objective evaluations.

and cite the appropriate paper:

2. Acerbi, L. & Ma, W. J. (2017). Practical Bayesian optimization for model fitting with Bayesian Adaptive Direct Search. In *Advances in Neural Information Processing Systems* 30:1834-1844.

Similarly, if you use IBS in combination with [Variational Bayesian Monte Carlo](https://github.com/lacerbi/vbmc), you should cite these papers:

3. Acerbi, L. (2018). Variational Bayesian Monte Carlo. In *Advances in Neural Information Processing Systems* 31: 8222-8232.
4. Acerbi, L. (2020). Variational Bayesian Monte Carlo with Noisy Likelihoods. In *Advances in Neural Information Processing Systems* 33: 8211-8222.

Besides formal citations, you can demonstrate your appreciation for our work in the following ways:

- *Star* the IBS repository on GitHub;
- Follow us on Twitter ([Luigi](https://twitter.com/AcerbiLuigi), [Bas](https://twitter.com/basvanopheusden)) for updates about IBS and other projects we are involved with;
- Tell us about your model-fitting problem and your experience with IBS (positive or negative) at <luigi.acerbi@helsinki.fi> or <svo@princeton.edu> (putting 'IBS' in the subject of the email).

### License

The IBS code is released under the terms of the [MIT License](https://github.com/lacerbi/ibs/blob/master/LICENSE.txt).
