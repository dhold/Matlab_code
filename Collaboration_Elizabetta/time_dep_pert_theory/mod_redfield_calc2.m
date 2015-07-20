function R_mod_red =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot...
            ,B,c_nk,ex_frq,t,working_tol)
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
N = length(ex_frq); R_mod_red = zeros(N,N);

t = reshape(t,1,size(g_broad,2));    
if sum(abs(diff(t,2)))<eps(length(t))
    use_simpsons = true;
else
    use_simpsons = false;
end

fctor = zeros(N);
for k =1:N %all sites same bath
    for k2 = 1:N
        for k3= 1:N
            for k4=1:N
 fctor(k,k2,k3,k4)  = sum(c_nk(:,k).*c_nk(:,k2).*c_nk(:,k3).*c_nk(:,k4));
            end
        end
    end
end

ex_frq = ex_frq-min(ex_frq); %scale to zero
for k = 1:N
    for j = k+1:N %k'
       % if k~=j
       dom = ex_frq(k)-ex_frq(j);     
       F_jk = fctor(k,k,k,k) + fctor(j,j,j,j) - 2*fctor(j,j,k,k);
       WW = exp(-1i*(dom+F_jk*lam_tot).*t - F_jk*g_broad); %based on Eds
       

      tmp  =   fctor(k,j,j,k)*g_sec_der - ...
              ((fctor(j,k,j,j) -  fctor(j,k,k,k))*g_deriv +...
              2i*fctor(j,k,j,j) *lam_tot).^2;
          
%           dom = ex_frq(k)-ex_frq(j);
%      exp_fct_imag = -dom + 2*(fctor(j,j,k,k) - fctor(j,j,j,j))*lam_tot ;
%      exp_fct_real =  -fctor(k,k,k,k) - fctor(j,j,j,j) + 2*fctor(j,j,k,k);
%            WW =  exp(1i*exp_fct_imag.*t+exp_fct_real*g_broad);
%                   
%        tmp  =   fctor(k,j,j,k)*g_sec_der - ...
%                 ((fctor(j,k,j,j) -  fctor(j,k,k,k))*g_deriv +...
%                 2i*fctor(j,k,j,j) *lam_tot).*...
%                 ((fctor(j,j,k,j) -  fctor(k,k,k,j))*g_deriv +...
%                 2i*fctor(j,j,k,j) *lam_tot);
       if nargin <9     
           if use_simpsons
       R_tmp = simpsons(WW.*tmp,t(1),t(end));        
           else
       R_tmp = trapz(t,WW.*tmp);    
           end
       else
       interp_fn = @(tt) ppval(pchip(t,WW.*tmp),tt);
        R_tmp = integral(interp_fn,min(t),max(t),'AbsTol',working_tol);  
       end
       R_mod_red(k,j) = -2*real(R_tmp);
       R_mod_red(j,k) = exp(B*dom)*R_mod_red(k,j); %detail balance cond
       % end
    end
end

R_mod_red = R_mod_red- diag(sum(R_mod_red,1));