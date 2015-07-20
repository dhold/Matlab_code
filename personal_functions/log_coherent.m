function [wf,n_lims] = log_coherent(alpha,sd)
%generates the 
if nargin == 1
    sd = 5; %number of standard deviations to include
end

n0 = -ceil(abs(alpha)*sd + 10);
nend = ceil(abs(alpha)*sd + 10);

n_rng = ceil(abs(alpha).^2) + [n0:nend];
n_rng = n_rng(n_rng>0);
n_lims = [n_rng(1), n_rng(end)]; 

wf = -alpha.^2/2 - sum(log(2:n_rng(1)))/2+ log(alpha)*n_rng(1) +...
    cumsum([0,log(alpha)-log(n_rng(2:end))/2]);
end