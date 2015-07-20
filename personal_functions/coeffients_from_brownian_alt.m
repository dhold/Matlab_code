function [cc_mat,cc,cc2,vv_mat,vv,truncation_correction] = ...
    coeffients_from_brownian_alt(lambda,gamma,om_0,Beta,Kappa,lam_dru,gam_dru)
%This function generates the coefficients for the exponential expansion of
%the time correlation function using Matsubara frequencies with over and
%under damped modes in the spectral density.

if isempty(lam_dru) && isempty(lambda)
    cc_mat=[]; cc =[]; vv_mat=[]; vv=[];
    truncation_correction=0;  return
end
% light_speed = 299792458; %units m s^{-1}, NOT cm 
%     length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
%     hbar = 1.05457173*10^(-34); %units Js
%     boltz_const = 1.3806488 * 10^(-23); 
% Beta = (2 * pi * hbar * light_speed * length_unit)/ ( Temp * boltz_const);

  %first calculate the poles from (1-exp(-beta *omega))^(-1)
 vv_mat = 2*pi*(1:Kappa)/Beta; %matsubara frequencies
  om_rng = -1i*vv_mat; %pole locations
 cc_mat = om_rng*0;   truncation_correction = 0;
    cnt =1; 
    
if ~isempty(om_0)    
lg = om_0./gamma/2 > 1; %under damped
lg2 = om_0./gamma/2==1; %critically damped
lg3 = om_0./gamma/2 < 1; %over damped


cc = zeros(1,2*sum(lg)+2*sum(lg3)+sum(lg2)+length(lam_dru)); 
vv = cc; cc2 = cc; %this is the cojugate of the residues, with the index switched for the BO

xi = zeros(size(lambda));


 for k = find(lg) %first underdamped most complicated due to the complex frequency
     
     xi(k) = sqrt(om_0(k)^2-gamma(k)^2/4);
     
     prefc = lambda(k)*om_0(k)^2/(2*xi(k));
     
     vv([cnt,cnt+1]) = gamma(k)/2 +[1i,-1i]*xi(k);
  cc([cnt,cnt+1]) = [1i,-1i].*prefc.*(cot(Beta*vv([cnt,cnt+1])/2)-1i);
  cc2([cnt,cnt+1]) = conj(cc([cnt+1,cnt])); %order switches
  
  cnt=cnt+2;
  
 end

  for k = find(lg2) %crit damped %not tested yet
              
    cc(cnt) = lambda(k)*om_0(k)*Beta/2/(sin(Beta*om_0(k)/2))^2;   %V^x coedfficient  
     vv(cnt) =  om_0(k);            %frequency
  cc2(cnt) = conj(cc(cnt));
        cnt=cnt+1;
              
  end
 
  for k = find(lg3) %overrrrdamped, poles imaginary %not tested yet
     
     xi(k) = sqrt(gamma(k)^2/4-om_0(k)^2);
     polepos1 = sqrt(om_0(k)^2 - gamma(k)^2/2 + gamma(k)*xi(k));  
     polepos2 = sqrt(om_0(k)^2 - gamma(k)^2/2 - gamma(k)*xi(k)); 
     
    vv(cnt) =  -1i*polepos1;
     tmp = lambda(k)*gamma(k)*om_0(k)^2/(2*om_0(k)^2-gamma(k)^2);
     %tmp = tmp*(1+1i*cot(Beta*polepos1/2)); 
    cc(cnt) = tmp*(1-1i*cot(Beta*polepos1/2));   
                %frequency
                    
      cnt=cnt+1;  
    cc(cnt) = tmp*(1-1i*cot(Beta*polepos2/2)); 
    cc2([cnt,cnt+1]) = conj(cc([cnt+1,cnt]));
     cnt=cnt+1;       
            
  end
 for k = 1:length(lambda)
 cc_mat = cc_mat - (lambda(k).*gamma(k).*om_0(k).^2 ).*(vv_mat./...
     ( (om_0(k).^2 + om_rng.^2).^2 - gamma(k).^2.*om_rng.^2));
      
 truncation_correction = truncation_correction + lambda(k)*(...
          (sin(Beta*gamma(k)/2) + gamma(k)*sinh(Beta*xi(k))/(2*xi(k)))/(...
          cos(Beta*gamma(k)/2)-cosh(Beta*xi(k))) + 2*gamma(k)/Beta/om_0(k)^2);

 end
  
end
if ~isempty(gam_dru)
  for k = 1:length(gam_dru) %so overdamped treated with drude approx
    cc(cnt) = lam_dru(k)*gam_dru(k)*(cot(Beta*gam_dru(k)/2)-1i)  ;     
     vv(cnt) =  gam_dru(k)  ;         %frequency
    cc2(cnt) = conj(cc(cnt));

 truncation_correction = truncation_correction + ...
                lam_dru(k)*(2/Beta/gam_dru(k)-1i) - cc(cnt)/vv(cnt);    
            
        cnt=cnt+1;
        
     if ~isempty(cc_mat)
             cc_mat = cc_mat + lam_dru(k).*gam_dru(k).*vv_mat./( vv_mat.^2-gam_dru(k).^2)/Beta;
     end
  end
end

 cc_mat = cc_mat*4; %pole residue, 
% cc_mat = -cc_mat*4;
 %truncation_correction is  dot(cc1(j,Kappa+1:infinity),1./vv1(Kappa+1,2:infinity))

 
  truncation_correction = truncation_correction - dot(cc_mat,1./(vv_mat));