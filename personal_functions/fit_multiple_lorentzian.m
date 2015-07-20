function [peak_locations, peak_width,peak_amp] =fit_multiple_lorentzian(...
            data_input,x_rng, n_peaks0,n_peaks_max,tol);
%Unfinished

%find peaks in function 
[maxtab,~] = peakdet(data_input,min(real(data_input))+eps);

peak_amp = maxtab(:,2); peak_loc = maxtab(:,1); 
lg = peak_loc ~= 1 * peak_loc ~= length(x_rng) ;
peak_amp = peak_amp(lg); peak_loc= peak_loc(lg); %remove edge peaks
%select n_peaks0 if given

[peak_amp,b] = sort(peak_amp,'descend');

if ~isempty(n_peaks0)
    n_peaks0 = min(n_peaks0,length(peak_amp));
peak_loc = peak_loc(b(1:n_peaks0)); peak_amp = peak_amp(1:n_peaks0); 
else
    n_peaks0 = min(n_peaks_max,length(peak_amp));
peak_loc = peak_loc(b(1:n_peaks0)); peak_amp = peak_amp(1:n_peaks0); 
end

%calculate second derivative at the peaks to estimate width (bad for noisy
% or course grained data)
pp = spline(x_rng, data_input);
qq = ppdiff(pp,2); sec_der =ppval(qq,x_rng);

b_guess = (2*peak_amp/pi /sec_der(peak_loc))^(1/3);
init_guess = [peak_amp,b_guess,peak_loc];
%use these as Guesses for a multi lorentizan fit
'y=a1/((pi*b1) * (1+((x-c1)/b1)^2))+a2/((pi*b2) * (1+((x-c2)/b2)^2))')
1/((pi*G) * (1+((x-x0)/G)^2))
'y = a * sin (b * x)'
f = ezfit(data_input,fun_to_fit,init_guess)