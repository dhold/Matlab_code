function [convfact, beta,si_cnst,cnst_names] = inv_cm_unit_sys(T_kelvin)
%1 cm^-1 = 0.000123986 eV 
%hbar =  2*pi*c = 0.01m = 1 codification
 speed_unit = 2*pi*299792458; %units m s^{-1}, NOT cm 
 length_unit = 0.01; % units m
   hbar = 1.05457173*10^(-34) ;  %units Js
 h_plank = 6.62606957*10^(-34) ;  %units Js
 boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}

%in kelvin, conv fact is beta = 1/(0.695039*T (in K))
 beta = (hbar* speed_unit)/ ( T_kelvin* boltz_const* length_unit); %in cm^-1

%conversion factor for fs to to inv cm
convfact = speed_unit / length_unit  *10^(-12); %2 pi * "c" * 100m^-1 * 10^-12 
 %multiply time in ps by this to get it in inverse cm
 %t_in_ps * convfact == t_in_inverse_cm
 %t_in_inverse_cm/convfact = t_in_ps
 % freq_in_inverse_cm * convfact = freq_in_THz 

if nargout ==3 %output si unit constants

si_cnst = {speed_unit,length_unit,h_plank,boltz_const,hbar,};
end
if nargout == 4
cnst_names = {'speed_unit=c','length_unit=0.01m','h','boltz_const','hbar'};
end
end