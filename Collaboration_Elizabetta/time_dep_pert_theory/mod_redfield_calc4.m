function R_mod_red =  mod_redfield_calc4(c_nk,ex_frq,tol2)
    persistent   lam_tot                                                           
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
if nargin ==2
    tol2 = 1e-3; %absolute tolerance of elements
end

if nargin == 1 %passed g_broad g_deriv g_sec_der lam_tot as cell
    %calculate interpolating functions w/ pchip and pass these to the
    %functions values not interpolatable within the region can be
    %calculated later
    g_broad = c_nk{1} ; g_deriv = c_nk{2} ;
    g_sec_der = c_nk{3};  t_rng = c_nk{4}; 
    B  = c_nk{5}; gam = c_nk{6}; om_0 = c_nk{7};
    lambda= c_nk{8}; gam_dru= c_nk{9}; lam_dru= c_nk{10};
    lam_tot = sum(lambda) + sum(lam_dru);
    
    to_pass = {B ,gam,om_0,lambda,gam_dru,lam_dru};
   fun_0 = pchip(t_rng,g_broad); g_0(t_rng,fun_0,to_pass);
   fun_1 = pchip(t_rng,g_deriv); g_1(t_rng,fun_1,to_pass);
   fun_2 = pchip(t_rng,g_sec_der); g_2(t_rng,fun_2,to_pass);
   R_mod_red = []; return
end
     
N = length(ex_frq); R_mod_red = zeros(N,N);
fctor = zeros(N);

%tol = 1e-9; 
% g_broad =  @(t) line_broad_fn_full(B,gam_rng,om_0,lam_rng,gam_dru,lam_dru,t,tol);    
% g_deriv  = @(t) line_broad_fn_deriv(B,gam_rng,om_0,lam_rng,gam_dru,lam_dru,t,tol);   
% g_sec_der = @(t) line_broad_fn_sec_der(B,gam_rng,om_0,lam_rng,gam_dru,lam_dru,t,tol);   
% lam_tot = sum(lam_rng) + sum(lam_dru);


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
           WW = @(t) exp(1i*exp_fct_imag.*t+exp_fct_real*g_1(t));
        
      tmp  = @(t)  fctor(k,j,j,k)*g_2(t) - ...
              ((fctor(j,k,j,j) -  fctor(j,k,k,k))*g_0(t) +...
              2i*fctor(j,k,j,j) *lam_tot).^2;
            
       R_tmp = @(t) WW(t).*tmp(t);    
        [Q] = integral(R_tmp,0,inf,'AbsTol',tol2);
        
       R_mod_red(k,j) = -2*real(Q);
        end
    end
end

R_mod_red = R_mod_red- diag(sum(R_mod_red,1));
end
function out = g_0(t,passfun,passparams)
persistent fun_0 t_saved othervars
if nargin >1
   t_saved = t; fun_0 = passfun; othervars =passparams; out=[]; return
end
out = zeros(size(t)); lg = t<=max(t_saved) & t>=min(t_saved);
    out(lg) = ppval(fun_0,t(lg));
    if any(~lg)
        out(~lg) = line_broad_fn_full(othervars{1},othervars{2},...
            othervars{3},othervars{4},othervars{5},othervars{6},t(~lg), 1e-10);
    end
end
function out = g_1(t,passfun,passparams)
persistent fun_1 t_saved othervars
if nargin >1
   t_saved = t; fun_1 = passfun; out=[]; othervars =passparams;return
end
out = zeros(size(t)); lg = t<=max(t_saved) & t>=min(t_saved);
    out(lg) = ppval(fun_1,t(lg));
    if any(~lg)
        out(~lg) = line_broad_fn_deriv(othervars{1},othervars{2},...
            othervars{3},othervars{4},othervars{5},othervars{6},t(~lg), 1e-10);
    end
end
function out = g_2(t,passfun,passparams)
persistent fun_2 t_saved othervars
if nargin >1
   t_saved = t; fun_2 = passfun; out=[]; othervars =passparams; return
end
    out = zeros(size(t)); lg = t<=max(t_saved) & t>=min(t_saved);
    out(lg) = ppval(fun_2,t(lg));
    if any(~lg)
        out(~lg) = line_broad_fn_sec_der(othervars{1},othervars{2},...
            othervars{3},othervars{4},othervars{5},othervars{6},t(~lg), 1e-10);
    end
end
