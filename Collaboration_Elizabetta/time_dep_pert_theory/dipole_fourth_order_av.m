function [alpha_av,CD_av,alpha_lin_pol] = ...
    dipole_fourth_order_av(mu,R,pol_L,pol_R,pol_linear,kk,simple)    
        if nargin <7
            simple = false; %take all transition combinations
        end
        
       N = size(mu,1); 
        
        if simple %take only averages of the form <mu_a,mu_a,mu_b,mu_b>
         alpha_L_4 = zeros(N); alpha_R_4 = zeros(N);
         alpha_L_5 = zeros(N); alpha_R_5 = zeros(N);              
        else
         alpha_L_4 = zeros(N,N,N,N); alpha_R_4 = zeros(N,N,N,N);
         alpha_L_5 = zeros(N,N,N,N); alpha_R_5 = zeros(N,N,N,N);   
         alpha_lin_pol = zeros(N,N,N,N);
        end
        
        if simple
for j1 = 1:N %pre save these
   for j2 = 1:N
       jj = [j1,j1,j2,j2];  
          alpha_lin_pol(j1,j2) = tensor_av(mu(jj,:),pol_linear); 
    alpha_L_4(j1,j2) = tensor_av(mu(jj,:),pol_L);
    alpha_R_4(j1,j2) = tensor_av(mu(jj,:),pol_R); 
    for j= 1:4
    alpha_L_5(j1,j2) = alpha_L_5(j1,j2) +...
        tensor_av([mu(jj,:);R(jj(j),:)],[pol_L;1i*kk(j,:)]);   
    alpha_R_5(j1,j2) = alpha_R_5(j1,j2) +...
        tensor_av([mu(jj,:);R(jj(j),:)],[pol_R;1i*kk(j,:)]);       
    end    
    
    
   end
end
        else
for j1 = 1:N %pre save these
   for j2 = 1:N
       for j3 = 1:N
           for j4=1:N
               
         jj = [j1,j2,j3,j4];             

        
    alpha_lin_pol(j1,j2,j3,j4) = tensor_av(mu(jj,:),pol_linear); 
    alpha_L_4(j1,j2,j3,j4) = tensor_av(mu(jj,:),pol_L);
    alpha_R_4(j1,j2,j3,j4) = tensor_av(mu(jj,:),pol_R);
    
    for j= 1:4
    alpha_L_5(j1,j2,j3,j4) = alpha_L_5(j1,j2,j3,j4) +...
        tensor_av([mu(jj,:);R(jj(j),:)],[pol_L;1i*kk(j,:)]);   
    alpha_R_5(j1,j2,j3,j4) = alpha_R_5(j1,j2,j3,j4) +...
        tensor_av([mu(jj,:);R(jj(j),:)],[pol_R;1i*kk(j,:)]);       
    end
           end
       end
   end
end
        end
    alpha_av = (alpha_L_4 + alpha_R_4)/2; 
    lg = alpha_L_5 + alpha_R_5; 
    if any(lg(:)~=0)
        warning('5th order contibution to absorbtion')
      alpha_av =  alpha_av+ lg /2;
    end


    CD_av = (alpha_L_5 - alpha_R_5); 
   lg2 = alpha_L_4 - alpha_R_4;  
   if any(lg2(:)~=0)
       warning('4th order contibution to CD')
    CD_av =  CD_av+ lg2;
   end