function [tau,E,E_w,E_inc] = pulsefns(Hwhm_fs,type_of_pulse)

if nargin == 1 %gaussian_pulses 
    type_of_pulse = 'Gaussian';
end
if strcmp(type_of_pulse,'Gaussian')
tau = (Hwhm_fs/sqrt(2*log(2)))/1000*convfact; %pulse SD in inverse cm  
E = @(t) exp(-t.^2/tau^2/2) / sqrt(tau * sqrt(pi)); %init e field env
%intensity of which is normed to 1
E_w = @(om) exp(-tau^2.*om.^2/2)*sqrt(tau/sqrt(pi));
E_inc = @(om,t) exp(-tau^2.*om.^2/2)*(1+1i*erfi(tau^2*om-1i*t)/sqrt(2)/tau)/2;   
        
elseif strcmp(type_of_pulse,'Sech') %sech pulses
    
tau = Hwhm_fs/acosh(2)/2/1000*convfact; 
E = @(t) sech(t/2/tau) / sqrt(4*tau); %init e field env, 
%intensity of which is normed to 1
E_w = @(om) sech(pi*tau.*om)*sqrt(tau*pi/2);
E_u_inc = [];% don't know
else
    error('invalid pulse type specified')

end