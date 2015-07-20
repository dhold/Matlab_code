function R_mod_red =  mod_redfield_calc(g_sym,g_deriv,g_sec_der,lam_tot,c_nk,ex_frq,t)
%calculates modified redfield rates, doesn't currently seem to work
%g_sym is g_n(t) g_deriv and g_sec_der are first and second derivatives 
%(same size in time), lam_tot is the total reorganisation energy
%c_nk is the participation of site n in exciton k
%ex_frq is the frequency of each exciton (relative to anything) t is time
%range
N = length(ex_frq); R_mod_red = zeros(N,N);
if nargin~=7
   t = sym('t','real'); 
else
t = reshape(t,1,size(g_sym,2));    
end

for k = 1:N
    for kk = 1:N %k'
        if k~=kk
            
            if size(g_sym,1)==1 %all sites same bath
       tmp = sum(c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k)) * g_sec_der;
       
       tmp2 = sum(c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)-...
                c_nk(:,kk).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)) * g_deriv;
            
       tmp2 = tmp2 + 2i*sum(c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)) *lam_tot;

%        tmp3 = sum(c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk)...
%                 - c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,kk) ) * g_deriv;
%             
%        tmp3 = tmp3 + 2i*sum(c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk)) *lam_tot;
%        
%        tmp4 = tmp - tmp2.*tmp3;
        tmp4 = tmp - tmp2.^2;
       
       %now calculate W'(omega_{k,k'},0,t)
       deltaom = +ex_frq(k) - ex_frq(kk);
       %add (subtract) factor of 2i(lambda_{k'k'kk} - lambda_{k'k'k'k'})
      deltaom = deltaom - 2*sum(c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,k)-...
           c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk)) * lam_tot;
      %deltaom = deltaom + 2*sum(c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,k)-...
      %     c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk)) * lam_tot; %this isn't right, but doesn't give stupid results
       
       g_fct = sum(c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)+...
            c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk)-...
            2*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,k)) * g_sym;
       
       WW = exp(-1i*deltaom.*t-g_fct);    
            else
       tmp = (c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k)).' * g_sec_der;
       
       tmp2 = (c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)-...
                c_nk(:,kk).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)).' * g_deriv;
            
       tmp2 = tmp2 + 2i*(c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)).' *lam_tot;
       
       tmp3 = (c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk)...
                - c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,kk) ).' * g_deriv;
            
       tmp3 = tmp3 + 2i*(c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk)).' *lam_tot;
       
       tmp4 = tmp - tmp2.*tmp3;
       deltaom = ex_frq(k) - ex_frq(kk);
      deltaom = deltaom - 2*(c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,k)-...
           c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk)).' * lam_tot;
       
       WW = exp(-1i*deltaom.*t...
           -(c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)+...
            c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk)-...
            2*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,k)).' * g_sym);
            end
% 
%        if rand(1)>0.8
%                   figure
% plot(t,real(WW))
% hold on
% plot(t,imag(WW),'r')
%        end
       if nargin~=7
       R_tmp = double(int(WW*tmp4,t,0,inf));
       else
       R_tmp = trapz(t,WW.*tmp4);    
       end
       R_mod_red(k,kk) = -2*real(R_tmp);
        end
    end
end