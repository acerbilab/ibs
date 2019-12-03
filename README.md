# Inverse binomial sampling (IBS)

## What is it

IBS is a technique to obtain unbiased, efficient estimates of the log-likelihood of a model by simulation. [[1](#references),[2](#references)]

The typical scenario is the case in which you have an *implicit* model from which you can randomly draw synthetic observations (for a given parameter vector), but cannot evaluate the log-likelihood analytically or numerically. In other words, IBS affords likelihood-based inference for likelihood-free models.

## References

1. de Groot, M. H. (1959). Unbiased sequential estimation for binomial populations. *The Annals of Mathematical Statistics*, 30(1):80–101.

2. van Opheusden\*, B., Acerbi\*, L. & Ma, W.J. (2019). Unbiased and efficient log-likelihood estimation with inverse binomial sampling. *arXiv preprint*. (\* equal contribution) ([preprint on arXiv]())

You can cite IBS in your work with something along the lines of

> We estimated the log-likelihood using inverse binomial sampling (IBS; de Groot, 1959; van Opheusden et al., 2019), a technique that produces unbiased and efficient estimates of the log-likelihood via simulation. 

If you use IBS in conjunction with [Bayesian Adaptive Direct Search](https://github.com/lacerbi/bads) (as recommended in the paper), you could add

> We obtained maximum-likelihood estimates of the model parameters via Bayesian Adaptive Direct Search (BADS; Acerbi & Ma, 2017), a hybrid Bayesian optimization algorithm which affords stochastic objective evaluations.

3. Acerbi, L. & Ma, W. J. (2017). Practical Bayesian optimization for model fitting with Bayesian Adaptive Direct Search. In *Advances in Neural Information Processing Systems* 30:1834-1844.

Besides formal citations, you can demonstrate your appreciation for our work in the following ways:

- *Star* the IBS repository on GitHub;
- Follow us on Twitter ([Luigi](https://twitter.com/AcerbiLuigi), [Bas](https://twitter.com/basvanopheusden)) for updates about IBS and other projects we are involved with;
- Tell us about your model-fitting problem and your experience with IBS (positive or negative) at <luigi.acerbi@gmail.com> or <svo213@nyu.edu> (putting 'IBS' in the subject of the email).

### License

IBS is released under the terms of the [MIT License](https://github.com/lacerbi/ibs/blob/master/LICENSE.txt).
