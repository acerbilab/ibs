function L = ibs_basic(fun,theta,R,S)
%IBS_BASIC A barebone implementation of inverse binomial sampling (IBS).
%  IBS_BASIC returns an unbiased estimate L of the log-likelihood for the 
%  simulated model and data, calculated using inverse binomial sampling. 
%  FUN is a function handle to a function that simulates the model's 
%  responses (see below). 
%  THETA is the parameter vector used to simulate the model's responses. 
%  R is a "response" data matrix, where each row correspond to one 
%  observation or "trial" (e.g., a trial of a psychophysical experiment), 
%  and each column represents a different response feature (e.g., the 
%  subject's response and reported confidence level). Responses need to 
%  belong to a finite set.
%  S is an optional "stimulus" matrix, where each row corresponds to one 
%  trial, and each column corresponds to a different trial feature (such 
%  as condition, stimulus value, etc.).
%  FUN is the generative model or simulator, which takes as input a vector 
%  of parameters PARAMS and a row of S, and generates one row of the 
%  response matrix.
%
%  This is a slow, bare bone implementation of IBS which should be used only 
%  for didactic purposes. A proper vectorized implementation of IBS is 
%  offered in the function IBSLIKE.
%
%  See also IBSLIKE.

%  Luigi Acerbi 2020

N = size(R,1);
L = zeros(N,1);

for i = 1:N         % Loop over all trials (rows)
    K = 1;
    while any(fun(theta,S(i,:)) ~= R(i,:))
        K = K + 1;  % Sample until the generated response is a match
    end
    L(i) = -sum(1./(1:K-1));    % IBS estimator for the i-th trial
end

L = sum(L);     % Return summed log-likelihood

end