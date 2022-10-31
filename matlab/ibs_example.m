% IBS_EXAMPLE Example script for running inverse binomial sampling (IBS)

% In this example, we fit some simple orientation discrimination data using
% IBS and BADS. The model is described in the "Orientation discrimination" 
% section of the Results in the published paper.
% Note that this model is very simple and used only for didactic purposes;
% one should use the analytical log-likelihood whenever available.
%
%   Reference: 
%   van Opheusden*, B., Acerbi*, L. & Ma, W. J. (2020), "Unbiased and 
%   efficient log-likelihood estimation with inverse binomial sampling". 
%   (* equal contribution), PLoS Computational Biology 16(12): e1008483.
%   Link: https://doi.org/10.1371/journal.pcbi.1008483

%% Generate simulated dataset

% First, we simulate a synthetic dataset to fit.
% Normally, you would have your own (real or simulated) data here.

fprintf('Generating synthetic data set...\n');

Ntrials = 600;
eta = log(1);                   % Fake subject (log) sensory noise
bias = 0.2;                     % Fake subject bias
lapse = 0.03;                   % Fake subject lapse rate
theta_true = [eta,bias,lapse];  % Generating parameter vector

S = 3*randn(Ntrials,1);        % Generate stimulus orientation per trial
R = psycho_gen(theta_true,S);  % Generate fake subject responses

% We set the lower/upper bounds for the parameters (in particular, note that 
% we set a nonzero lower bound for the lapse rate)
LB = [log(0.1) -2 0.01];
UB = [log(10) 2 1];

% We also set the "plausible" lower/upper bounds representing where we
% would expect to find most parameters (this is not a prior, it just
% identifies a region e.g. for the starting points)
PLB = [log(0.2) -1 0.02];
PUB = [log(5) 1 0.2];


%% Maximum-likelihood estimation (MLE)

% We fit the data via maximum-likelihood estimation using Bayesian Adaptive 
% Direct Search (BADS), a particularly effective optimization algorithm for
% noisy target functions.
% If you do not have the BADS toolbox installed, you can freely download it 
% from here: https://github.com/acerbilab/bads
% Please ensure you have the latest version of BADS installed (v1.1.1 or
% more), as it incorporates tweaks that improve performance with IBS.

fprintf('Maximum-likelihood estimation with BADS using IBS...\n');

if isempty(which('bads.m'))                     % Check that you have BADS
    error('BADS not found. You can install it from here: https://github.com/acerbilab/bads');
end

old_bads = false;
try
    bads_version = bads('version');
    fprintf('BADS found (version %s).\n', bads_version);
    bads_version(bads_version == '.') = ' ';
    bads_version = str2num(bads_version);
    if bads_version(2) < 1
        old_bads = true;
    end
catch
    old_bads = true;
end

if old_bads
    error('You have installed an older version of BADS. You can find the latest version here: https://github.com/acerbilab/bads');        
end

fprintf('(press a key to continue)\n');
pause;

% We define the negative log-likelihood function via a call to IBSLIKE
% (IBSLIKE provides a vectorized implementation of IBS for MATLAB; check
% out IBS_BASIC for a bare bone implementation)

% We set 10 reps for the IBS estimator (see "Iterative multi-fidelity" 
% section in the Methods of the  paper). Note that this is also the default 
% for IBSLIKE
% We also tell IBSLIKE to return as a second output the standard deviation 
% of the estimate (as opposed to the variance). This is very important!

options_ibs = ibslike('defaults');
options_ibs.Nreps = 10;
options_ibs.ReturnStd  = true;

nllfun_ibs = @(theta) ibslike(@psycho_gen,theta,R,S,options_ibs);

% As a starting point for the optimization, we draw a sample inside the
% plausible box (in practice you should use multiple restarts!)
theta0 = rand(size(LB)).*(PUB-PLB) + PLB;

% We inform BADS that IBSLIKE returns noise estimate (SD) as second output
options = bads('defaults');
options.SpecifyTargetNoise = true;  

theta_ibs = bads(nllfun_ibs,theta0,LB,UB,PLB,PUB,[],options);

% Compare with MLE obtained using the analytical log-likelihood expression

fprintf('Maximum-likelihood estimation with BADS using the exact log-likelihood...\n');
fprintf('(press a key to continue)\n');
pause;

options = bads('defaults');
theta_exact = bads(@(x) psycho_nll(x,S,R),theta0,LB,UB,PLB,PUB,[],options);


%% Analysis of results

fprintf('Returned maximum-likelihood solutions with different methods:\n');

theta_ibs
theta_exact

fprintf('True data-generating parameter value:\n');
theta_true

% We compute now the log-likelihood at the returned solution with higher
% precision (100 repeats), and estimate its uncertainty

options_ibs = ibslike('defaults');
options_ibs.Nreps = 100;
options_ibs.ReturnStd  = true;

% IBSLIKE returns as second output the variance of the IBS estimate
% (please note that this might change in future versions)
[nll_ibs,nll_std] = ibslike(@psycho_gen,theta_ibs,R,S,options_ibs);
nll_exact = psycho_nll(theta_exact,S,R);

fprintf('Estimated log-likelihood at the IBS-found solution: %.2f +/- %.2f (exact value: %.2f).\n', ...
    -nll_ibs,nll_std,-psycho_nll(theta_ibs,S,R));

% We also evaluate the MLE found exactly
fprintf('Log-likelihood at the exact-found solution: %.2f.\n',-nll_exact);

% Note that even the exact MLE will differ from the true data-generating
% parameters due to finiteness of the dataset (we expect to recover the 
% true data generating parameters in the limit NTRIALS -> infinity)


%% Bayesian posterior estimation

close all;

% We now fit the same data via a method that computes the Bayesian posterior
% over model parameters, called Variational Bayesian Monte Carlo (VBMC).
% If you do not have the VBMC toolbox installed, you can freely download it 
% from here: https://github.com/acerbilab/vbmc

fprintf('Bayesian posterior estimation with Variational Bayesian Monte Carlo (VBMC) using IBS...\n');
fprintf('(press a key to continue)\n');
pause;

if isempty(which('vbmc.m'))                     % Check that you have VBMC
    error('VBMC not found. You can install it from here: https://github.com/acerbilab/vbmc');
end

folder = fileparts(which('vbmc.m'));
addpath([folder filesep 'utils']);  % Ensure that UTILS is on the path

% For VBMC, we need to return the POSITIVE log-likelihood as first output, 
% and the standard deviation (not variance!) of the estimate as second
% output.

% Options for IBSLIKE
options_ibs = ibslike('defaults');
options_ibs.Nreps = 100;            % Try and have a SD ~ 1 (and below 3)
options_ibs.ReturnPositive = true;  % Return *positive* log-likelihood
options_ibs.ReturnStd  = true;      % 2nd output is SD (not variance!)

llfun = @(theta) ibslike(@psycho_gen,theta,R,S,options_ibs);

% We add a trapezoidal or "tent" prior over the parameters (flat between 
% PLB and PUB, and linearly decreasing to 0 towards LB and UB). 
lpriorfun = @(x) log(mtrapezpdf(x,LB,PLB,PUB,UB));

% Since LLFUN has now two outputs (the log likelihood at X, and the
% estimated SD of the log likelihood at X), we cannot directly sum the log
% prior from LPRIORFUN and the log likelihood from LLFUN.
% Thus, we use an auxiliary function LPOSTFUN that does exactly this. 

fun = @(x) lpostfun(x,llfun,lpriorfun);     % Log joint

% FUN will have two outputs: the first output is the sum of LLFUN and 
% LPRIORFUN (the log joint), and the second output is the variability of
% LLFUN (that is, the second output of LLFUN). It is assumed that the prior
% can be computed analytically and adds no noise.

% Here we could use as starting point the result of a run of BADS
theta0 = 0.5*(PLB+PUB);
options = vbmc('defaults');
options.Plot = true;        % Plot iterations
options.SpecifyTargetNoise = true;  % Noisy function evaluations 

% Run VBMC
[vp,elbo,elbo_sd] = vbmc(fun,theta0,LB,UB,PLB,PUB,options);

close all;

Xs = vbmc_rnd(vp,3e5);  % Generate samples from the variational posterior

close all;
cornerplot(Xs,{'\eta (log noise)','bias','lapse'},[eta,bias,lapse]);

fprintf('  Have a look at the triangle-plot visualization of the approximate posterior.\n')
fprintf('  The black line represents the true generating parameters for the fake dataset.\n')
fprintf('  For more information, see tutorials and FAQ at <a href="https://github.com/acerbilab/vbmc">https://github.com/acerbilab/vbmc</a>.\n')


% Luigi Acerbi, 2021-2022
