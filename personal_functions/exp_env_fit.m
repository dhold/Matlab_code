function [aa_fit,first_guess,fsamp,ff] = ...
            exp_env_fit(tt,signal,est_nfl,over_est_w )
%fits an exponential envelope over a function

%example parameters
% fsamp = 10;
% tmax = 4*pi;
% tt=0:1/fsamp:tmax;
% LL = length(tt); NFFT = 2^nextpow2(LL); 
% ff = fsamp/2*linspace(0,1,NFFT/2+1);
% 
%  signal_damped_noise = (cos(1.4*tt) + cos(4.2*tt)+rand(1,LL)/10).*exp(-0.5*tt);
%  est_nfl = 0; 
%signal_undamped_noise = (cos(1.4*tt) + cos(4.2*tt)).*exp(-0.5*tt)  +randn(1,LL)/10;
%est_nfl = 1e-3; %some estimate of (non damped) noise level, 
%assumed to obey gaussian statistics with this the standard deviation
% sigfft = fft(signal,NFFT)/LL;

%shape to row vectors
signal = reshape(signal,1,length(signal));
tt = reshape(tt,1,length(tt));

 fsamp = 1/(tt(2)-tt(1));
 %tmax = tt(end);
 LL = length(tt); NFFT = 2^nextpow2(LL); 
 ff = fsamp/2*linspace(0,1,NFFT/2+1);
sigfft = fft(signal,NFFT)/LL;
sigred = 2*abs(sigfft(1:NFFT/2+1));
%make a guess based on the width of the maximum peak in the DFT
[a,b] = max(sigred); 

lg1 = sigred > a*3/4 & (1:length(ff) > b); lg{1} = find(lg1,1,'first');
lg1L = sigred > a*3/4 & (1:length(ff) > b); lg{2} = find(lg1L,1,'last');
lg2 = sigred > a*2/3 & (1:length(ff) > b); lg{3} = find(lg2,1,'first');
lg2L = sigred > a*2/3 & (1:length(ff) > b); lg{4} = find(lg2L,1,'last');
%use these four points to estimate widths and thus gamma in
% aGamma^2/(Gamma^2+(b'-b)^2) = k -> Gamma = sqrt(k/(1-k))*|b'-b|

G_est = zeros(4,1); est_mat = zeros(1,4);
w_vec = [3,3,2,2]; 
for k =1:4
    if ~isempty(lg{k})
     G_est(k) = sqrt(w_vec(k) * (ff(lg{k})-ff(b))^2);
    est_mat(k) = 1;
    end
end
if any(logical(est_mat))
first_guess = est_mat*G_est/sum(est_mat);
else
first_guess = 1; %this failed so I will have to take a generic guess
end
fit_fn = @(aa) (max(signal)+eps(max(signal))+7*est_nfl)*exp(-aa*tt); 
%function to fit to envelope
%cost function grows large when the fit function is smaller than the signal
%at any point and decreases very slowly as aa becomes larger
if nargin <4
over_est_w = 50/LL; %weighting value given to all signal over env fit
end
%adjust with noise level
shift = 2.5; 
step_fn = @(delta) heaviside(delta/est_nfl - shift).*...
                    tanh((delta/est_nfl).^2 - shift^2);
%tends to HSstep when no noise present, else switches on gradually
% shift factor in tanh allows for delta to drift to 2.5 s.d. of noise 
% before being included in the cost_fn

if est_nfl == 0  %e.g. theoretically generated data
cost_fn = @(aa) over_est_w*abs(aa)*sum((abs(signal)-fit_fn(aa))./fit_fn(aa)...
                .*heaviside(abs(signal) - fit_fn(aa))) - aa;
else
cost_fn = @(aa) over_est_w*abs(aa)*sum((abs(signal)-fit_fn(aa))...
                    ./max(fit_fn(aa),est_nfl)...
                .*step_fn(abs(signal)-fit_fn(aa))) - aa;
end

aa_fit = fminsearch(cost_fn,first_guess);

% figure
% plot(tt,abs(signal))
% hold on
% plot(tt,fit_fn(aa_fit),'--')

%%
if 1==0
fsamp = 1e5;
tmax = 0.1;
t=0:1/fsamp:tmax;
f = 12e3; %should be smaller than fsamp/2!
tau = 0.0765;
y=(sin(2 * pi * f * t) +rand(1,length(t))/20 ) .* exp(-t / tau);

%calculate running maximum
n = 20; %number of points to take max over
nblocks = floor(length(t) / n);
trun = mean(reshape(t(1:n*nblocks), n, nblocks), 1); %n-point mean
envelope = max(reshape(y(1:n*nblocks), n, nblocks), [], 1); %n-point max

%quick and dirty exponential fit, not the proper way in case of noise
p = polyfit(trun, log(envelope), 1);
tau_fit = -1/p(1);
k_fit = exp(p(2));
figure
plot(t, y, trun, envelope, 'or', t, k_fit * exp(-t / tau_fit), '-k')
title(sprintf('tau = %g', tau))

end