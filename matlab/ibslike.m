function [nlogL,nlogLvar,output] = ibslike(fun,params,respMat,designMat,options,varargin)
%IBSLIKE Unbiased negative log-likelihood via inverse binomial sampling.
%   NLOGL = IBLLIKE(FUN,PARAMS,RESPMAT,DESIGNMAT) returns unbiased estimate 
%   NLOGL of the negative of the log-likelihood for the simulated model 
%   and data, calculated using inverse binomial sampling (IBS). 
%   FUN is a function handle to a function that simulates the model's 
%   responses (see below). 
%   PARAMS is the parameter vector used to simulate the model's responses. 
%   RESPMAT is a "response" data matrix, where each row correspond to one 
%   observation or "trial" (e.g., a trial of a psychophysical experiment), 
%   and each column represents a different response feature (e.g., the 
%   subject's response and reported confidence level). Responses need to 
%   belong to a finite set.
%   DESIGNMAT is an optional experimental design matrix, where each row 
%   corresponds to one trial, and each column corresponds to a different 
%   trial feature (such as condition, stimulus value, etc.).
%
%   FUN takes as input a vector of parameters PARAMS and an experimental 
%   design matrix DMAT (one row per trial), and generates a matrix of 
%   simulated model responses (one row per trial, corresponding to rows of 
%   DMAT). DMAT is built by the algorithm out of rows of DESIGNMAT.
%
%   DESIGNMAT can be omitted or left empty, in which case FUN needs to 
%   accept a parameter vector PARAMS and an array of trial numbers T, and 
%   returns a matrix of simulated responses, where the i-th row contains 
%   the simulated response for trial T(i) (the indices in T may repeat).
%
%   NLOGL = IBSLIKE(FUN,PARAMS,RESPMAT,DESIGNMAT,OPTIONS) uses options in 
%   structure OPTIONS to replace default values. (To be explained...)
%
%   NLOGL = IBSLIKE(...,VARARGIN) additional arguments are passed to FUN.
%
%   [NLOGL,NLOGLVAR] = IBSLIKE(...) also returns an estimate NLOGVAR of the
%   variance of the log likelihood.
%
%   [NLOGL,NLOGLVAR,PC] = IBSLIKE(...) returns an array of probabilities of
%   correct response in each trial. These estimated probabilities are biased.
%
%   [NLOGL,NLOGLVAR,PC,OUTPUT] = IBSLIKE(...) returns a structure OUTPUT 
%   with additional information about the sampling.
%
%   OPTIONS = IBSLIKE('defaults') returns a basic default OPTIONS structure.
%
%   EXITFLAG = IBSLIKE('test') runs some tests. Here EXITFLAG is 0 if 
%   everything works correctly.
%   
%   Test code on binomial sampling:
%      p = 0.7; Ntrials = 100;                  % Define binomial probability
%      fun = @(x,dmat) rand(size(dmat)) < x;    % Simulating function
%      rmat = fun(p,NaN(Ntrials,1));            % Generate responses
%      [nlogL,nlogLvar,pc,output] = ibslike(fun,p,rmat);
%      nlogL_true = -log(p)*sum(rmat == 1) - log(1-p)*sum(rmat == 0);
%      fprintf('Ground truth: %.4g, Estimate: %.4g ± %.4g.\n',nlogL_true,nlogL,sqrt(nlogLvar));
%
%   Reference: 
%   van Opheusden*, B., Acerbi*, L. & Ma, W. J. (2020), "Unbiased and 
%   efficient log-likelihood estimation with inverse binomial sampling". 
%   (* equal contribution), arXiv preprint https://arxiv.org/abs/2001.03985
%
%   See also @.

%--------------------------------------------------------------------------
% IBS: Inverse Binomial Sampling for unbiased log-likelihood estimation
% To be used under the terms of the MIT License 
% (https://opensource.org/licenses/MIT).
%
%   Authors (copyright): Luigi Acerbi and Bas van Opheusden, 2020
%   e-mail: luigi.acerbi@{gmail.com,nyu.edu}, basvanopheusden@nyu.edu
%   URL: http://luigiacerbi.com
%   Version: 0.901
%   Release date: Jan 22, 2020
%   Code repository: https://github.com/lacerbi/ibs
%--------------------------------------------------------------------------

if nargin < 4; designMat = []; end
if nargin < 5; options = []; end

% Default options
% defopts.Display       = 'off';        % Level of display on screen
defopts.Nreps           = 10;           % # independent log-likelihood estimates per trial
defopts.NegLogLikeThreshold = Inf;      % Stop sampling if estimated nLL is above threshold (incompatible with vectorized sampling)
defopts.Vectorized      = 'auto';       % Use vectorized sampling algorithm with acceleration
defopts.Acceleration    = 1.5;          % Acceleration factor for vectorized sampling
defopts.NsamplesPerCall = 0;            % # starting samples per trial per function call (0 = choose automatically)
defopts.MaxIter         = 1e5;          % Maximum number of iterations (per trial and estimate)
defopts.ReturnPositive  = false;        % If true, the first returned output is the *positive* log-likelihood
defopts.ReturnStd       = false;        % If true, the second returned output is the standard deviation of the estimate

%% If called with no arguments or with 'defaults', return default options
if nargout <= 1 && (nargin == 0 || (nargin == 1 && ischar(fun) && strcmpi(fun,'defaults')))
    if nargin < 1
        fprintf('Basic default options returned (type "help ibslike" for help).\n');
    end
    nlogL = defopts;
    return;
end

%% If called with one argument which is 'test', run test
if nargout <= 1 && nargin == 1 && ischar(fun) && strcmpi(fun,'test')
    nlogL = runtest();
    return;
end

for f = fields(defopts)'
    if ~isfield(options,f{:}) || isempty(options.(f{:}))
        options.(f{:}) = defopts.(f{:});
    end
end

% Add hard-coded options
options.MaxSamples      = 250;          % Maximum number of samples per function call
options.AccelerationThreshold = 0.1;    % Keep accelerating until threshold is passed (in s)
options.VectorizedThreshold  = 0.1;     % Max threshold for using vectorized algorithm (in s)
options.MaxMem          = 1e6;          % Maximum number of samples for vectorized implementation

% NSAMPLESPERCALL should be a scalar integer
if ~isnumeric(options.NsamplesPerCall) || ~isscalar(options.NsamplesPerCall)
    error('ibslike:NsamplesPerCall','OPTIONS.NsamplesPerCall should be a scalar integer.');
end

% NEGLOGLIKETHRESHOLD should be a scalar greater than 0 (or Inf)
if ~isnumeric(options.Acceleration) || ~isscalar(options.Acceleration) || ...
        options.Acceleration < 1
    error('ibslike:Acceleration','OPTIONS.Acceleration should be a scalar equal or greater than one.');
end

% ACCELERATION should be a scalar equal or greater than 1
if ~isnumeric(options.NegLogLikeThreshold) || ~isscalar(options.NegLogLikeThreshold) || ...
        options.NegLogLikeThreshold <= 0
    error('ibslike:NegLogLikeThreshold','OPTIONS.NegLogLikeThreshold should be a positive scalar (including Inf).');
end

Ntrials = size(respMat,1);
Trials = (1:Ntrials)';
funcCount = 0;

simdata = []; elapsed_time = [];
if ischar(options.Vectorized) && options.Vectorized(1) == 'a'
    if isfinite(options.NegLogLikeThreshold)
        vectorized_flag = false;
    else
        % First full simulation to determine computation time
        fun_clock = tic;
        if isempty(designMat)    % Pass only trial indices
            simdata = fun(params,Trials(:),varargin{:});
        else                    % Pass full design matrix per trial
            simdata = fun(params,designMat(Trials(:),:),varargin{:});   
        end
        elapsed_time = toc(fun_clock);
        vectorized_flag = elapsed_time < options.VectorizedThreshold;
        funcCount = 1;
    end
else
    if options.Vectorized
        if isfinite(options.NegLogLikeThreshold)
            error('ibslike:UnsupportedVectorized','Vectorized sampling does not support negative log likelihood thresholding.');
        end
        vectorized_flag = true;
    else
        vectorized_flag = false;
    end
end

if vectorized_flag
    [nlogL,K,Nreps,Ns,fc] = vectorized_ibs_sampling(fun,params,respMat,designMat,simdata,elapsed_time,options,varargin{:});
else
    [nlogL,K,Nreps,Ns,fc] = loop_ibs_sampling(fun,params,respMat,designMat,simdata,elapsed_time,options,varargin{:});
end
funcCount = funcCount + fc;

% Variance of estimate per trial
if nargout > 1
    Ktab = -(psi(1,1:max(K(:)))' - psi(1,1));    
    LLvar = Ktab(K);   
    nlogLvar = sum(LLvar,2)./Nreps.^2;
end

% OUTPUT structure with additional information
if nargout > 2
    % output.Hits = hits;
    % output.TrialCount = sum(~isnan(hits(:)));
    % output.TargetHits = targetHits;
    output.funcCount = funcCount;
    output.NsamplesPerTrial = Ns/Ntrials;
    output.nlogL_trials = nlogL;
    output.nlogLvar_trials = nlogLvar;
end

% Return negative log-likelihood and variance summed over trials
nlogL = sum(nlogL);
if options.ReturnPositive; nlogL = -nlogL; end
if nargout > 1
    nlogLvar = sum(nlogLvar);
    if options.ReturnStd    % Return standard deviation instead of variance
        nlogLvar = sqrt(nlogLvar);
    end
end

end

%--------------------------------------------------------------------------
function [nlogL,K,Nreps,Ns,fc] = vectorized_ibs_sampling(fun,params,respMat,designMat,simdata0,elapsed_time0,options,varargin)

Ntrials = size(respMat,1);
Trials = (1:Ntrials)';
Ns = 0;
fc = 0;

BufferSize = 500;   % Initial size of hit matrix

% A "hit" is a sample draw from the model that matches the data
hits = NaN([BufferSize,Ntrials]);

% Set the first row to one -- these are not real hits, will be ignored 
% later, but it's necessary for how the run lengths K are computed 
hits(1,:) = 1;

offset = 1;
Nhits = zeros(Ntrials,1);
targetHits = options.Nreps(:).*ones(Ntrials,1);
lastsample = zeros(Ntrials,1);
MaxIter = options.MaxIter*max(options.Nreps);

% Starting samples
if options.NsamplesPerCall == 0
    samples_level = options.Nreps;
else
    samples_level = options.NsamplesPerCall;
end

for iter = 1:MaxIter
    % Pick trials that need more hits, sample multiple times
    T = Trials(Nhits < targetHits);
    if isempty(T); break; end
        
    % With accelerated sampling, might request multiple samples at once
    Nsamples = min(options.MaxSamples,max(1,round(samples_level)));
    MaxSamples = ceil(options.MaxMem / numel(T));
    Nsamples = min(Nsamples, MaxSamples);
    Tmat = repmat(T,[1,Nsamples]);
    
    % Simulate trials
    if iter == 1 && Nsamples == 1 && ~isempty(simdata0)
        simdata = simdata0;
        elapsed_time = elapsed_time0;
    else
        fun_clock = tic;
        if isempty(designMat)    % Pass only trial indices
            simdata = fun(params,Tmat(:),varargin{:});
            fc = fc + 1;
        else                    % Pass full design matrix per trial
            simdata = fun(params,designMat(Tmat(:),:),varargin{:});   
            fc = fc + 1;
        end
        elapsed_time = toc(fun_clock);        
    end
    Ns = Ns + numel(T);
    
    % Accelerated sampling
    if options.Acceleration > 0 && elapsed_time < options.AccelerationThreshold
        samples_level = samples_level*options.Acceleration;
    end
    
    % Add new hits to old ones
    hits_temp = all(respMat(Tmat(:),:) == simdata,2);    
    hits_new = reshape(hits_temp,size(Tmat));
    Nhits(T) = Nhits(T) + sum(hits_new,2);
    rows = [1,Nsamples]+offset;     % Rows in hits matrix
    if rows(2) > size(hits,1)       % Extend hits matrix if needed
        hits = [hits;NaN([max(rows(2)-Nsamples,BufferSize),Ntrials])];
    end
    hits(rows(1):rows(2),T) = hits_new';
    offset = offset + Nsamples;
    lastsample(T) = offset;
    
end

if ~isempty(T)
    error('ibslike:ConvergenceFail', ...
        'Maximum number of iterations reached and algorithm did not converge. Check FUN and DATA.');
end
    
% Log likelihood estimate per trial and run lengths K for each repetition
Nreps = min(targetHits,Nhits);
[nlogL,K] = estimateLog(hits(1:offset,:),Nreps);

end

%--------------------------------------------------------------------------
function [nlogL,K,Nreps,Ns,fc] = loop_ibs_sampling(fun,params,respMat,designMat,simdata0,elapsed_time0,options,varargin)

Ntrials = size(respMat,1);
Trials = (1:Ntrials)';
MaxIter = options.MaxIter;

nlogL = NaN(Ntrials,options.Nreps);
K = NaN(Ntrials,options.Nreps);
Ns = 0;
fc = 0;

for iRep = 1:options.Nreps
    
    offset = 1;
    nlogL_sum = 0;
    
    hits = zeros(Ntrials,1);
    for iter = 1:MaxIter
        % Pick trials that need more hits, sample multiple times
        T = Trials(hits < 1);
        if isempty(T); break; end        

        % Simulate trials
        if iter == 1 && iRep == 1 && ~isempty(simdata0)
            simdata = simdata0;
        elseif isempty(designMat)    % Pass only trial indices
            simdata = fun(params,T(:),varargin{:});
            fc = fc + 1;
        else                    % Pass full design matrix per trial
            simdata = fun(params,designMat(T(:),:),varargin{:});   
            fc = fc + 1;
        end
        Ns = Ns + numel(T); % Count samples
        hits_new = all(respMat(T(:),:) == simdata,2);    
        hits(T) = hits(T) + hits_new;
        
        K(T(hits_new),iRep) = offset; 
        nlogLs = (psi(K(T(hits_new),iRep)) - psi(1));
        
        nlogL(T(hits_new),iRep) = nlogLs;        
        nlogL_sum = nlogL_sum + sum(nlogLs);
        
        offset = offset + 1;
        
        % Terminate if negative log likelihood is above a given threshold
        if nlogL_sum > options.NegLogLikeThreshold
            T = Trials(hits < 1);
            if ~isempty(T)
                %K(T(hits_new),iRep) = offset; 
                %nlogL(T(hits_new),iRep) = psi(K(T(hits_new),iRep)) - psi(1);
                K(T,iRep) = offset;
                nlogL(T,iRep) = psi(K(T,iRep)) - psi(1);
            end
            T = [];
            break;
        end
    end    
end

if ~isempty(T)
    error('ibslike:ConvergenceFail', ...
        'Maximum number of iterations reached and algorithm did not converge. Check FUN and DATA.');
end
    
nlogL = mean(nlogL,2);
Nreps = options.Nreps;

end

%--------------------------------------------------------------------------
function [nlogL,K] = estimateLog(hitsmat,Nreps)
%ESTIMATELOG Estimate negative log-likelihood from matrix of hits.
%
%   [NLOGL,K] = ESTIMATELOG(HITSMAT,NREPS) estimates negative log-likelihood
%   via inverse binomial sampling using matrix of hits/misses HITSMAT.
%   HITSMAT is a binary matrix where each row is a trial, 0 stand for
%   misses and 1 for hits. NREPS is an array of requested number of
%   independent estimates per trial. NLOGL is a column array of negative 
%   log-likelihood estimates per trial. K is a matrix where each row is a 
%   trial and each column contains the number of samples needed to achieve 
%   a hit for that estimate. The columns of K are padded with NaNs if NREPS 
%   is unequal between trials.

[h,Ntrials] = size(hitsmat);

% This is the comprehensible version:
%
%   for iTrial = 1:Ntrials
%       index = find(hits(iTrial,:),targetHits(iTrial));
%       K = diff([0 index]);
%       logL(iTrial) = sum(Ktab(K))/numel(index);
%   end

% This is the efficient one, without for loops:
list = find(hitsmat(:) == 1)-1;
row = floor(list/h)+1;
col = mod(list,h)+1;
delta = diff([col;1]);
remidx = delta <= 0;
delta(remidx) = [];
row(remidx) = [];
indexcol = find(diff([0;row]));
col = 1 + (1:numel(row))' - indexcol(row);
K = ones(Ntrials,max(col));
K(row + (col-1)*size(K,1)) = delta;
Ktab = -(psi(1:max(K(:)))' - psi(1));   % Precompute digamma function   
mask = bsxfun(@le,repmat(1:size(K,2),[Ntrials,1]),Nreps);
K = max(mask.*K,1);
LL = Ktab(K);
nlogL = -sum(mask.*LL,2)./Nreps;

end

%--------------------------------------------------------------------------
function exitflag = runtest()

% Binomial probability model
Ntrials = 100;                  
p_true = 0.9*rand() + 0.05;             % True probability
p_model = 0.9*rand() + 0.05;            % Model probability
fun = @(x,dmat) rand(size(dmat)) < x;   % Simulating function
Nexps = 2e3;

zscores = zeros(1,Nexps);
for iter = 1:Nexps
    rmat = fun(p_true,NaN(Ntrials,1));            % Generate data
    [nlogL,nlogLvar] = ibslike(fun,p_model,rmat);
    nlogL_exact = -log(p_model)*sum(rmat == 1) - log(1-p_model)*sum(rmat == 0);
    zscores(iter) = (nlogL_exact - nlogL)/sqrt(nlogLvar);
end

edges = -4.75:0.5:4.75;
nz = histc(zscores,edges);
h(1) = bar(edges,nz,'histc');
hold on;
xx = linspace(-5,5,1e4);
h(2) = plot(xx,Nexps*exp(-xx.^2/2)/sqrt(2*pi)/2,'k-','LineWidth',2);

box off;
set(gca,'TickDir','out');
set(gcf,'Color','w');
xlabel('z-score');
ylabel('pdf')
xlim([-5 5]);
hl = legend(h,'z-scores histogram','expected pdf');
set(hl,'Location','NorthEast','Box','off');

fprintf('\n');
fprintf('Using IBS to compute the log-likelihood of a binomial distribution.\n');
fprintf('Parameters: p_true=%.2g, p_model=%.2g, %d trials per experiment.\n',p_true,p_model,Ntrials);
fprintf('The distribution of z-scores should approximate a standard normal distribution (mean 0, SD 1).\n');

exitflag = abs(mean(zscores)) > 0.15 || abs(std(zscores) - 1) > 0.1;

fprintf('\n');
if exitflag
    fprintf('Test FAILED. Something might be wrong.\n');    
else
    fprintf('Test PASSED. We verified that IBS is unbiased (~zero mean) and calibrated (SD ~1).\n');
end
fprintf('Distribution of z-scores (%d experiments). Mean: %.4g. Standard deviation: %.4g.\n',Nexps,mean(zscores),std(zscores));
fprintf('\n');


end

%   TODO:
%   - Fix help and documentation
%   - Optimal allocation of estimates?