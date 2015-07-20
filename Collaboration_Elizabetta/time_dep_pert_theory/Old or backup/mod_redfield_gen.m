function R_comp =  mod_redfield_gen(c_nk,ex_frq,site_dep)

%generate symbolic solution to modified Redfield tensor, also need to
%perform integral

N = length(c_nk);
R_comp = sym(zeros(N));

if site_dep
g_sym = sym('g',[N,1]);
g_deriv = sym('gg',[N,1]);
g_sec_der = sym('ggg',[N,1]);
lam_tot = sym('LL',[N,1]);
t = sym('t','real');


for k = 1:N
    for kk = 1:N %k'
        if k~=kk
       tmp = (c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k)).' * g_sec_der;
       
       tmp2 = (c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)-...
                c_nk(:,kk).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)).' * g_deriv;
            
       tmp2 = tmp2 + 2i*(c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)).' *lam_tot;
       
       tmp3 = (c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk)...
                - c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,kk) ).' * g_deriv;
            
       tmp3 = tmp3 + 2i*(c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk)).' *lam_tot;
       
       tmp4 = tmp - tmp2.*tmp3;
       deltaom = ex_frq(kk) - ex_frq(k);
       %include shift from the reorganisation energy differences
       deltaom = deltaom - (c_nk(:,k).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)-...
           c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)).' * lam_tot;
       
       
       WW = exp(-1i*deltaom.*t...
           -(c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)+...
            c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk)-...
            2*c_nk(:,k).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)).' * g_sym);
       
       R_comp(k,kk) = WW*tmp4;
        end
    end
end

else %all sites same bath
    
 g_sym = sym('g');
g_deriv = sym('gg');
g_sec_der = sym('ggg');
lam_tot = sym('LL');
t = sym('t','real');


for k = 1:N
    for kk = 1:N %k'
        if k~=kk
       tmp = sum(c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k)) * g_sec_der;
       
       tmp2 = sum(c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)-...
                c_nk(:,kk).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)) * g_deriv;
            
       tmp2 = tmp2 + 2i*sum(c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)) *lam_tot;
       
       tmp3 = sum(c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk)...
                - c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,kk) ) * g_deriv;
            
       tmp3 = tmp3 + 2i*sum(c_nk(:,kk).*c_nk(:,kk).*c_nk(:,k).*c_nk(:,kk)) *lam_tot;
       
       tmp4 = tmp - tmp2.*tmp3;
       deltaom = ex_frq(kk) - ex_frq(k);
      deltaom = deltaom - sum(c_nk(:,k).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)-...
           c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)) * lam_tot;
       
       WW = exp(-1i*deltaom.*t...
           -sum(c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)+...
            c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk).*c_nk(:,kk)-...
            2*c_nk(:,k).*c_nk(:,k).*c_nk(:,kk).*c_nk(:,kk)) * g_sym);
       
       R_comp(k,kk) = WW*tmp4;
        end
    end
end   
    
end