function R_mod_red =  mod_redfield_calc3(B,gam_rng,om_0,...
                        lam_rng,gam_dru,lam_dru,c_nk,ex_frq,tol,tol2)
%calculates modified redfield rates, 
% g_sym is g_n(t) (assumed all the same in this program)
% g_deriv and g_sec_der are first and second derivatives 
%(same size in time), lam_tot is the total reorganisation energy
% c_nk is the participation of site n in exciton k
% ex_frq is the frequency of each exciton (relative to anything) t is time
% range
%
% R_kkk'k' = - 2 Re int_0^inf dt W'(w_kk',0,t)*{ (d^2/d_t^2 g_{k,k',k',k} - 
%    [d/d_t(g_{k',k,k',k'} - g_{k',k,k,k}+2it*lambda_{k',k,k',k'} )]^2  }
%               with g(t) the line broadening function
if nargin ==8
    tol2 = 1e-3;
    tol = 1e-6;
end
N = length(ex_frq); R_mod_red = zeros(N,N);
fctor = zeros(N);

%tol = 1e-9; 
g_broad =  @(t) line_broad_fn_full(B,gam_rng,om_0,lam_rng,gam_dru,lam_dru,t,tol);    
g_deriv  = @(t) line_broad_fn_deriv(B,gam_rng,om_0,lam_rng,gam_dru,lam_dru,t,tol);   
g_sec_der = @(t) line_broad_fn_sec_der(B,gam_rng,om_0,lam_rng,gam_dru,lam_dru,t,tol);   
lam_tot = sum(lam_rng) + sum(lam_dru);
for k =1:N %all sites same bath%all sites same bath
    for k2 = 1:N
        for k3= 1:N
            for k4=1:N
 fctor(k,k2,k3,k4)  = sum(c_nk(:,k).*c_nk(:,k2).*c_nk(:,k3).*c_nk(:,k4));
            end
        end
    end
end


for k = 1:N
    for j = 1:N %k'
        if k~=j
            
          dom = ex_frq(k)-ex_frq(j);
     exp_fct_imag = -dom + 2*(fctor(j,j,k,k) - fctor(j,j,j,j))*lam_tot ;
     exp_fct_real =  -fctor(k,k,k,k) - fctor(j,j,j,j) + 2*fctor(j,j,k,k);
           WW = @(t) exp(1i*exp_fct_imag.*t+exp_fct_real*g_broad(t));
        
      tmp  = @(t)  fctor(k,j,j,k)*g_sec_der(t) - ...
              ((fctor(j,k,j,j) -  fctor(j,k,k,k))*g_deriv(t) +...
              2i*fctor(j,k,j,j) *lam_tot).^2;
            
       R_tmp = @(t) WW(t).*tmp(t);    
        [Q] = integral(R_tmp,0,inf,'AbsTol',tol2);
        
       R_mod_red(k,j) = -2*real(Q);
        end
    end
end

R_mod_red = R_mod_red- diag(sum(R_mod_red,1));