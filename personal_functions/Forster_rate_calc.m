function R_forster =  Forster_rate_calc(H_site,g_t,lam_tot,t)

N = length(H_site);

R_forster = zeros(size(H_site));

if numel(t) >1 && size(g_t,1) == 1
    same_bath_each_site = true;
else
    same_bath_each_site = false;
end

if length(t)~=1
   test = diff(t); 
   if sum(abs(test))/length(test) < eps
       use_simpson = true;
       if floor(length(t)/2)*2 ~= length(t)
           %extend time range to use
          t = [t,t(end)+test(end)] ;
          g_t = [g_t,zeros(size(g_t,1),1)]; %assume last point zero
       end
   else
       use_simpson = false;
   end

end

for k = 1:N
    for j = 1:N
        dom = H_site(k,k)-H_site(j,j);
        if same_bath_each_site 
            dom = dom + 2*lam_tot;
        else
            dom = dom + lam_tot(k)+lam_tot(j);
        end
        

       if length(t)~=1  
           
         if same_bath_each_site 
            g_fn = +2*g_t;
        else
            g_fn = +g_t(k,:)+g_t(j,:);
         end
        
        if use_simpson
       int_val = simpsons(exp(-g_fn-1i*dom*t),t(1),t(end),[]);         
        else %use trapz with unevenly spaced points
       int_val = trapz(t,exp(-g_fn-1i*dom*t));    
        end
        
       else %t is assumed to be a working tolerence and g_anon_fn
  
       integ_fn= @(tt) exp(-g_t(k,tt)-g_t(j,tt)-1i*dom*tt);
       int_val = integral(integ_fn,0,infinity,'AbsTol',t);  
       end        
        
        R_forster(k,j) = 2*abs(H_site(k,j))*real(int_val);        
        
    end
end