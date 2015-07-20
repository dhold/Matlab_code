function [cc_mat,cc,cc2,vv_mat,vv,truncation_correction] = ...
    coeffients_from_brownian_alt2(Beta,Kappa,lam_dru,gam_dru)
warning('UNFINISHED')
%This function generates the coefficients for the exponential expansion of
%the time dependent correlation function using the PSD sum for the poles of
%the Bose function.  Drude only


if isempty(lam_dru) 
    cc_mat=[]; cc =[]; vv_mat=[]; vv=[];
    truncation_correction=0;  return
else
    cc = zeros(size(gam_dru)); vv =cc;  cc2 = cc; cc_mat = zeros(Kappa,1);
  for k = 1:length(gam_dru) %so overdamped treated with drude approx
    cc(cnt) = lam_dru(k)*gam_dru(k)*(cot(Beta*gam_dru(k)/2)-1i);       
     vv(cnt) =  gam_dru(k);            %frequency
    cc2(cnt) = conj(cc(cnt));
        cnt=cnt+1;

 truncation_correction = truncation_correction + lam_dru(k)*(...
         2-Beta*gam_dru(k)*cot(Beta*gam_dru(k)/2))/Beta/gam_dru(k);    
     if ~isempty(cc_mat)
             cc_mat = cc_mat - lam_dru(k).*gam_dru(k).*vv_mat./( gam_dru(k).^2-vv_mat.^2)/Beta;
     end
  end
end

 cc_mat = cc_mat*4; %pole residue, 
 %truncation_correction is  dot(cc1(j,Kappa+1:infinity),1./vv1(Kappa+1,2:infinity))

 
  truncation_correction = truncation_correction - dot(cc_mat,1./(vv_mat));