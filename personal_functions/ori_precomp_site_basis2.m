function  [av_2,av_3,av2_3,av_4,av_5,av2_5,chk]=ori_precomp_site_basis2...
                                            (param_set1,param_set2,mu,R,mm);

% Precompute all the possible dipole expectation values 
% of the form <mu_4 dot e_out* ... mu_1 dot e_1  exp(i*k dot r_4..)>_{isotropic}
% these are computed as taylor expansions exp(1ix) ~ 1 + ix in terms of 4th
% and 5th order tensor averages
% param_set1 is a set of polarisations for the 
N = size(mu,1);
if nargin<5
    mm = [];
end
if nargin <4
    R = [];
end
if isempty(param_set1); param_set1 = zeros(4,3,0); end
if isempty(param_set2); param_set2 = zeros(7,3,0); end

av_4 = zeros(N,N,N,N,size(param_set1,3));  
av_2 = zeros(N,N,size(param_set1,3));  
if ~isempty(R)
av_5 = zeros(N,N,N,N,size(param_set2,3),3); %needs a factor of +i in front
av_3 = zeros(N,N,size(param_set2,3));  
else
    av_5 = [];
end
if ~isempty(mm)
av2_5 = zeros(N,N,N,N,size(param_set2,3),3); 
av2_3 = zeros(N,N,size(param_set2,3));  
else
    av2_5=[];
end

if nargout >4 
%check input polarizations and wavevectors constitute a valid set
chk = false(size(param_set2,3),3); 
%has value this scheme is valid for rephasing / nr / coh
              for lp2 = 1:size(param_set2,3)
      pol = param_set2(1:4,:,lp2);  k_set = param_set2(5:7,:,lp2);   
      lg = abs(dot(pol(1,:),k_set(1,:)))<=eps(2) ...
          & abs(dot(pol(2,:),k_set(2,:)))<=eps(2) ...
          &  abs(dot(pol(3,:),k_set(3,:)))<=eps(2) ;
      chk(lp2, 1) = lg & abs(dot(pol(4,:),[-1,1,1]*k_set))<=eps(2);
      chk(lp2, 2) = lg & abs(dot(pol(4,:),[1,-1,1]*k_set))<=eps(2);
      chk(lp2, 3) = lg & abs(dot(pol(4,:),[1,1,-1]*k_set))<=eps(2);
        if ~any(chk(lp2, :))
            lp2
warning('param_set2 contains beam configurations which are not possible') 
        end
              end
end
    
for k1 = 1:N %slow but only needs doing once
    for k2 = 1:N    
        
            for lp = 1:size(param_set1,3)

             pol = param_set1(1:2,:,lp); 
        
            mu_set = [mu(k1,:);mu(k2,:)];
        av_2(k1,k2,lp) = tensor_av(mu_set,pol(1:2,:));
            end
            
            for lp2 = 1:size(param_set2,3)
     %pol = param_set2(2:3,:,lp2) ;k_set = param_set2(1,:,lp2);   
      pol = param_set2(1:2,:,lp2) ;k_set = param_set2(3,:,lp2);    
            mu_set = [mu(k1,:);mu(k2,:)];
                     
                  if ~isempty(R)
        DeltaR = R(k1,:)-R(k2,:);
         av_3(k1,k2,lp) =  1i*(tensor_av([mu_set;DeltaR],[pol(1:2,:);k_set(1,:)]));
                  end
                  if ~isempty(mm)
         %now on site mag dipole moment        
         tmp1 = mu_set; tmp2 = mu_set; temp1 = pol(1:2,:);  temp2 = pol(1:2,:); 
        temp1(2,:)= -cross(k_set(1,:),pol(2,:));   temp2(1,:)= -cross(k_set(1,:),pol(1,:));  
         tmp1(2,:) = mm(k2,:); tmp2(1,:) = mm(k1,:); 
         av2_3(k1,k2,lp) = 1i*(tensor_av(tmp2,temp2)-tensor_av(tmp1,temp1));
                  end   
             end   
    end
end
    
if nargout >4 
for k1 = 1:N %slow but only needs doing once
    for k2 = 1:N
       for k3 = 1:N 
           for k4 = 1:N  
              lpset =[k1,k2,k3,k4];
              mu_set = [mu(k1,:);mu(k2,:);mu(k3,:);mu(k4,:)];
              for lp = 1:size(param_set1,3)

                pol = param_set1(1:4,:,lp); 
              av_4(k1,k2,k3,k4,lp) = tensor_av(mu_set,pol);
              %next calculate the additional averages
              end
              for lp2 = 1:size(param_set2,3)
     % pol = param_set2(4:7,:,lp2);  k_set = param_set2(1:3,:,lp2);     
      pol = param_set2(1:4,:,lp2);  k_set = param_set2(5:7,:,lp2);         
              for j = 1:3
                   k_j = k_set(j,:); 
                  if ~isempty(R)
        DeltaR = R(lpset(j),:)-R(lpset(4),:);
         av_5(k1,k2,k3,k4,lp2,j) =  1i*(tensor_av([mu_set;DeltaR],[pol;k_j]));
                  end
                  if ~isempty(mm)
         %now on site mag dipole moment        
         tmp1 = mu_set; tmp2 = mu_set; temp1 = pol;  temp2 = pol; 
        temp1(4,:)= -cross(k_j,pol(4,:));   temp2(j,:)= -cross(k_j,pol(j,:));  
         tmp1(4,:) = mm(lpset(4),:); tmp2(j,:) = mm(lpset(j),:); 
        
         av2_5(k1,k2,k3,k4,lp2,j) = 1i*(tensor_av(tmp2,temp2)-tensor_av(tmp1,temp1));
                  end
                  
                  %should also have a quadrupole moment but I have nothing
                  %good to calculate this yet
              end
              end
           end
       end
    end
end
end

% to get out the exciton basis averages from the coefficients C_ex do
% tmp = mtimesx(mtimesx(C_ex,'C',av_4),C_ex);
% tmp = mtimesx(mtimesx(C_ex,'C',squeeze(tmp)),C_ex);
% to get the 5th rank averages for a given set of wavevector amplitudes 
% tmp2 = mtimesx(mtimesx(C_ex,'C',av_5),C_ex);
% tmp2 = mtimesx(mtimesx(C_ex,'C',squeeze(tmp2)),C_ex);
% tmp2 = s(1)*tmp2(:,1).*norm(k_1) + s(2)*tmp2(:,2).*norm(k_2) + s(3)*tmp2(:,3).*norm(k_3)
% norms of these are proportional to the frequency
% integration over multiple k values may be required, but I hope not and
% that I can just take the carrier frequency.
 
%{
av_4 = zeros(N,N,N,N,size(param_set,3));  
%different averages will be different for rephasing, norephasing etc
av_5_rp = zeros(N,N,N,N,size(param_set,3),4);
av_5_nr = av_5_rp; av_5_coh = av_5_rp; %from interchromo-moment

%test to check both methods give the same answer
av2_5_rp = av_5_rp; av2_5_nr = av_5_rp; av2_5_coh = av_5_rp; %intrachromo mag moment

s_rf = [-1,+1,+1,-1]; s_nr = [+1,-1,+1,-1]; s_coh = [+1,+1,-1,-1];
%determines the order of conjugation etc

for lp = 1:size(param_set,3)

pol = param_set(1:4,:,lp); 
k_set = param_set(5:7,:,lp);
k_rp = s_rf(1).*k_set(1,:)+s_rf(2).*k_set(2,:)+s_rf(3).*k_set(3,:); 
k_nr = s_nr(1).*k_set(1,:)+s_nr(2).*k_set(2,:)+s_nr(3).*k_set(3,:); 
k_coh = s_coh(1).*k_set(1,:)+s_coh(2).*k_set(2,:)+s_coh(3).*k_set(3,:); 
for k1 = 1:N %slow but only needs doing once
    for k2 = 1:N
       for k3 = 1:N 
           for k4 = 1:N  
              lpset =[k1,k2,k3,k4];
              mu_set = [mu(k1,:);mu(k2,:);mu(k3,:);mu(k4,:)];
              av_4(k1,k2,k3,k4,lp) = tensor_av(mu_set,pol);
              %next calculate the additional averages
              for j = 1:3
                  k_j = k_set(j,:);
                  temp = tensor_av([mu_set;R(lpset(j),:)],[pol;k_j]);
           %calculate different contributions from RP non rephasing, coh
            av_5_rp(k1,k2,k3,k4,lp,j) =  s_rf(j)*temp;
            av_5_nr(k1,k2,k3,k4,lp,j) =  s_nr(j)*temp;
            av_5_coh(k1,k2,k3,k4,lp,j) = s_coh(j)*temp;
              end
av_5_rp(k1,k2,k3,k4,lp,4) =  - tensor_av([mu_set;R(lpset(4),:)],[pol;k_rp]);
av_5_nr(k1,k2,k3,k4,lp,4) =  - tensor_av([mu_set;R(lpset(4),:)],[pol;k_nr]);
av_5_coh(k1,k2,k3,k4,lp,4) = - tensor_av([mu_set;R(lpset(4),:)],[pol;k_coh]);           

        %also intrinsic mag moment
              for j = 1:3
                  k_j = k_set(j,:);  pol_tmp = pol; mm_set = mu_set;
                   pol_tmp(j,:) = cross(k_j,pol_tmp(j,:));
                   mm_set(j,:) = imag(mdip2(lpset(j),:));
                    temp = tensor_av(mm_set,pol_tmp);
           %calculate different contributions from RP non rephasing, coh
            av2_5_rp(k1,k2,k3,k4,lp,j) =  s_rf(j)*temp;
            av2_5_nr(k1,k2,k3,k4,lp,j) =  s_nr(j)*temp;
            av2_5_coh(k1,k2,k3,k4,lp,j) = s_coh(j)*temp;
              end
              pol_rp = pol; pol_rp(4,:) = cross(k_rp,pol_rp(j,:));
              pol_nr = pol; pol_nr(4,:) = cross(k_nr,pol_nr(j,:));
              pol_coh = pol; pol_coh(4,:) = cross(k_nr,pol_coh(j,:));
            mm_set = mu_set;  mm_set(4,:) = imag(mdip2(lpset(4),:)); 
            
            av2_5_rp(k1,k2,k3,k4,lp,j) =  -tensor_av(mm_set,pol_rp);
            av2_5_nr(k1,k2,k3,k4,lp,j) =  -tensor_av(mm_set,pol_nr);
            av2_5_coh(k1,k2,k3,k4,lp,j) = -tensor_av(mm_set,pol_coh);

           end
       end
    end
end
end
%}
