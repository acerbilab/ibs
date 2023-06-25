function L=psycho_nll(theta,S,R)
%PSYCHO_NLL Negative log-likelihood for psychometric function model.
%  L=PSYCHO_NLL(THETA,S,R) returns the negative log likelihood for a simple
%  orientation discrimination task; where THETA is a model parameter vector, 
%  with THETA(1) as eta=log(sigma), the log of the sensory noise; THETA(2) 
%  the bias term; THETA(3) is the lapse rate; S is a vector of stimulus 
%  orientations (in deg) for each trial; and R the vector of responses
%  per trial (1 for "rightwards" and -1 for "leftwards").

% Luigi Acerbi, 2020

sigma = exp(theta(1));
bias = theta(2);
lapse = theta(3);

% Likelihood per trial (analytical solution)
p_vec = lapse/2+(1-lapse)*((R==-1).*normcdf(-(S-bias)/sigma)+(R==1).*normcdf((S-bias)/sigma));

% Total negative log-likelihood
L = -sum(log(p_vec));

end


function z = normcdf(x)
%NORMCDF Normal cumulative distribution function (cdf)
z = 0.5 * erfc(-x ./ sqrt(2));
end