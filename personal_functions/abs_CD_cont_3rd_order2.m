 function [alpha_L_4,alpha_R_4,alpha_L_5, alpha_R_5] = abs_CD_cont_3rd_order2(mu,R,jj,pol,kk)
 %kk is set of wavevectors of light
 %
 mu = mu(jj,:);  R = R(jj,:); %reorder these terms to the order of interactions
 
 if length(pol) <  3
pol_L = [1,1i,0]/sqrt(2); pol_R = [1,-1i,0]/sqrt(2);
 else
pol_L = pol{3}; pol_R = conj(pol_L);     
 end
if size(kk,1) < 3
    kk = [kk;[0,0,1]]; %take probe along y with UNIT typical frequency
end
if length(pol) == 5 && length(kk)==4
pol_het_beam_L = pol{4}; pol_het_beam_R = pol{5}; %take another beam to be heterodyne sig
else
pol_het_beam_L = conj(pol_L); pol_het_beam_R = conj(pol_R); %self heterodyned
kk = [kk;-sum(kk)];  % conjugated by det scheme
end
alpha_L_4 = 0; alpha_R_4 = 0;  %4th order
alpha_L_5 = 0; alpha_R_5 = 0;   %5th order
%two required 4th order average functions
   xxxx_av = @(aa) (1/15)*(dot(aa(1,:),aa(2,:))*dot(aa(3,:),aa(4,:))+...
                dot(aa(1,:),aa(3,:))*dot(aa(2,:),aa(4,:))+...
                dot(aa(2,:),aa(3,:))*dot(aa(1,:),aa(4,:))); 

   yyxx_av = @(aa) (1/30)*(4*dot(aa(1,:),aa(2,:))*dot(aa(3,:),aa(4,:))-...
                dot(aa(1,:),aa(3,:))*dot(aa(2,:),aa(4,:))-...
                dot(aa(2,:),aa(3,:))*dot(aa(1,:),aa(4,:))); 
            
for x1 = 1:3 
    for x2 = 1:3
        for x3 = 1:3
            for x4 = 1:3
         pol_cmp_L = pol{1}(x1)*pol{2}(x2)*pol_L(x3)*pol_het_beam_L (x4);
         pol_cmp_R = pol{1}(x1)*pol{2}(x2)*pol_R(x3)*pol_het_beam_R (x4);
        if pol_cmp_L~=0
 co_set = [x1,x2,x3,x4]; 
[sorted_set,reorder] = sort(co_set);  
diffset = logical(diff(sorted_set)); 

if isequal(diffset,[0,0,0])%aaaa
term = xxxx_av(mu);  
elseif isequal(diffset,[0,1,0]) %aabb
    %reorder mu to the correct shape, note that it doesn't matter if I have
    %xxyy or yyxx so it doesn't matter which comes first 
 term = yyxx_av(mu(reorder,:));     
end
    alpha_L_4 = alpha_L_4 + term*pol_cmp_L;
    alpha_R_4 = alpha_R_4 + term*pol_cmp_R;
   
                for x5 = 1:3
               
co_set = [x1,x2,x3,x4,x5];
[sorted_set,reorder] = sort(co_set); 
diff_set = logical(diff(sorted_set)); inc_term = true; %logical whether non zero
if isequal(diff_set,[1,1,0,0]) %xxx
    %set already in the correct order %xxxyz
elseif isequal(diff_set,[0,1,1,0]) %yyy
    reorder = reorder([2:4,5,1]); %yyyzx
    %sorted_set = sorted_set(reorder);
elseif isequal(diff_set,[0,0,1,1]) %zzz
    reorder = reorder([3:5,1,2]); %zzzxy
    %sorted_set = sorted_set(reorder);
else
    inc_term = false;
end

if inc_term
    
                 for j  = 1:4
                    mu2 = [mu;R(j,:)];
    mu_reorder = mu2(reorder(1:4),:); %mu2(5,:) 
    mu_reorder(4,:) = cross(mu_reorder(4,:), mu(reorder(5),:))/2;
   term = xxxx_av(mu_reorder);     
       alpha_L_5 = alpha_L_5 + 1i*term*pol_cmp_L*kk(j,x5);                    
    alpha_R_5 = alpha_R_5 + 1i*term*pol_cmp_R*kk(j,x5); 
                 end
end
                end
        end
                 
            end
        end
    end
end

 end