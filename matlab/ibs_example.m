% IBS_EXAMPLE Example script for running inverse binomial sampling (IBS)

% In this example, we fit some simple orientation discrimination data using
% IBS and BADS. The model is described in Section 5.2 of the manuscript.
% Note that this model is very simple and used only for didactic purposes;
% one should use the analytical log-likelihood whenever available.

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

%% Maximum-likelihood estimation (MLE)

% We fit the data via maximum-likelihood estimation using Bayesian Adaptive 
% Direct Search (BADS), a particularly effective optimization algorithm.
% If you do not have the BADS toolbox installed, you can freely download it 
% from here: https://github.com/lacerbi/bads

% We set the lower/upper bounds for optimization (in particular, note that 
% we set a nonzero lower bound for the lapse rate)
LB = [log(0.1) -2 0.01];
UB = [log(10) 2 1];

% We also set the "plausible" lower/upper bounds used by BADS
PLB = [log(0.1) -1 0.01];
PUB = [log(5) 1 0.2];

% We define the negative log-likelihood function via a call to IBSLIKE
% (IBSLIKE provides a vectorized implementation of IBS for MATLAB; check
% out IBS_BASIC for a bare bone implementation)

% We set 10 reps for the IBS estimator (see Section 4.4 in the paper) -
% note that this is also the default for IBSLIKE
options_ibs.Nreps = 10;
nllfun_ibs = @(theta) ibslike(@psycho_gen,theta,R,S,options_ibs);

% As a starting point for the optimization, we draw a sample inside the
% plausible box (in practice you should use multiple restarts!)
theta0 = rand(size(LB)).*(PUB-PLB) + PLB;

fprintf('Maximum-likelihood estimation with BADS using IBS...\n');
fprintf('(press a key to continue)\n');
pause;

theta_ibs = bads(nllfun_ibs,theta0,LB,UB,PLB,PUB);

% Compare with MLE obtained using the analytical log-likelihood expression

fprintf('Maximum-likelihood estimation with BADS using the exact log-likelihood...\n');
fprintf('(press a key to continue)\n');
pause;

theta_exact = bads(@(x) psycho_nll(x,S,R),theta0,LB,UB,PLB,PUB);

%% Analysis of results

fprintf('Returned maximum-likelihood solutions with different methods:\n');

theta_ibs
theta_exact

fprintf('True data-generating parameter value:\n');
theta_true

% We compute now the log-likelihood at the returned solution with higher
% precision (100 repeats), and estimate its uncertainty
options_ibs.Nreps = 100;

% IBSLIKE returns as second output the variance of the IBS estimate
[nll_ibs,nll_var] = ibslike(@psycho_gen,theta_ibs,R,S,options_ibs);
nll_exact = psycho_nll(theta_exact,S,R);

fprintf('Estimated log-likelihood at the IBS-found solution: %.2f +/- %.2f (exact value: %.2f).\n', ...
    -nll_ibs,sqrt(nll_var),-psycho_nll(theta_ibs,S,R));

% We also evaluate the MLE found exactly
fprintf('Log-likelihood at the exact-found solution: %.2f.\n',-nll_exact);

% Note that even the exact MLE will differ from the true data-generating
% parameters due to finiteness of the dataset (we expect to recover the 
% true data generating parameters in the limit NTRIALS -> infinity)

% Luigi Acerbi, 2020