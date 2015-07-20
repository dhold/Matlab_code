function R_mod_red =  mod_redfield_calc_id_baths...
(g_broad,g_deriv,g_sec_der,lam_tot,B,QQ,ex_frq,t,working_tol)
%calculates modified redfield rates, 
% g_sym is g_n(t) (assumed all the same in this program)
% g_deriv and g_sec_der are first and second derivatives 
%(same size in time), lam_tot is the total reorganisation energy
% QQ  = sum_n <a | Q_n |b> <c|Q_n |d>
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


ex_frq = ex_frq-min(ex_frq); %scale to zero
for k = 1:N
    for j = k+1:N %k'
       % if k~=j
       dom = ex_frq(k)-ex_frq(j);     %frequency difference
       
       F_jk = QQ(k,k,k,k) + QQ(j,j,j,j) - 2*QQ(j,j,k,k); %will be real
       
       WW = exp(-1i*(dom + F_jk*lam_tot).*t - F_jk*g_broad); %based on Eds
       WWdecay = -F_jk*real(g_broad);
       WWphase = -F_jk*(lam_tot.*t+imag(g_broad));
       
      tmp  =   QQ(k,j,j,k)*g_sec_der - ...
              ((QQ(j,k,j,j) -  QQ(j,k,k,k))*g_deriv +2i*QQ(j,k,j,j)*lam_tot).^2;
          
       if nargin <9     
           if use_simpsons
       R_tmp = simpsons(WW.*tmp,t(1),t(end));        
           else
       R_tmp = trapz(t,WW.*tmp);    
           end
       else %interpolate each of the terms seperately as they will be smoother
           %also real valued, this might cause problems with 
           W_decay =  pchip(t,WWdecay); %real negative part of exp
           W_phase =  pchip(t,WWphase); %phase evolution excluding delta om
           W_prefct = pchip(t,tmp); %prefactor
           
       %interp_fn = @(tt) ppval(pchip(t,WW.*tmp),tt);
       interp_fn = @(tt) real(ppval(W_prefct,tt).*...
                    exp(ppval(W_decay,tt)+1i*(ppval(W_phase,tt)-dom.*tt)));
       
        R_tmp = integral(interp_fn,min(t),max(t),'AbsTol',working_tol);  
       end
       R_mod_red(k,j) = -2*R_tmp;
       R_mod_red(j,k) = exp(B*dom)*R_mod_red(k,j); %detail balance cond
       % end
    end
end

R_mod_red = R_mod_red- diag(sum(R_mod_red,1));