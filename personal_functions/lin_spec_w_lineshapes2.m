function  [OD,CD,FL,OD_ind,CD_ind,FL_ind] = lin_spec_w_lineshapes2...
            (w_ex,c_nk,g_t,lam,tau,d_n,R_n,P_k,t_rng,om_r_rng)
%Calculates absorption, Circular dichroism and Florescence linshapes
%this function takes only dipole moments d_n and positions R_n
%tau = 1/sum(-R_k'k'kk,{k'~=k}) = 1/R_kkkk
mid_freq = median(om_r_rng); %middle frequency
OD = zeros(size(om_r_rng)); CD = OD; FL = OD;
dt = t_rng(2)-t_rng(1);  tmax  = t_rng(end);
om_w_ft = 2*pi*(0:1/tmax:1/dt);
om_w_ft = om_w_ft -pi/dt + mid_freq;
N = length(w_ex);

%FFT stuff
%NFFT = 2^nextpow2(length(t_rng));
if nargout>3
   OD_ind = zeros(N,length(om_r_rng));
     CD_ind = OD_ind;  FL_ind = OD_ind;
end
for k = 1:N %exciton state loop

     if size(g_t,1)==1 %nosite dep baths
           g_ex = sum(c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k))*g_t;
           lam_ex = sum(c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k))*lam;
     else %site dependent baths
          g_ex = (c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)).'*g_t;
          lam_ex =(c_nk(:,k).*c_nk(:,k).*c_nk(:,k).*c_nk(:,k)).'*lam;
     end
    
            %ifft not fft as the integral is over exp(i om t)
 OD_tmp = real(fftshift(ifft...
            (exp(1i*(mid_freq-w_ex(k))*t_rng-g_ex -t_rng/2/tau(k)))));

FL_tmp = P_k(k)*real(fftshift(ifft(exp(1i*(mid_freq-w_ex(k))*t_rng...
            +2i*lam_ex*t_rng - conj(g_ex) - t_rng/2/tau(k)))));
           
%                if k==3
%                    figure
%             plot(om_w_ft,OD_tmp)
%             hold on 
%             plot(om_w_ft,FL_tmp,'r')
%                end
%      om_test_val = om_r_rng(round(end/2));
%      OD_test = real(trapz(t_rng,exp(1i*(om_test_val-w_ex(k))*t_rng-g_ex -t_rng/2/tau(k))));      
%     FL_test = P_k(k)*real(trapz(t_rng,exp(1i*(om_test_val-w_ex(k))*t_rng...
%                     +2i*t_rng*lam_ex - conj(g_ex) - t_rng/2/tau(k))));   
%                  (OD_tmp(round(end/2))-OD_test)/OD_test
%                  (FL_tmp(round(end/2))- FL_test)/FL_test

    OD_tmp = interp1(om_w_ft,OD_tmp,om_r_rng);
    FL_tmp = interp1(om_w_ft,N*FL_tmp,om_r_rng);
     
    dip_fct_1 = 0; CD_fct_1 = 0;
    for n1=1:N
        for n2 = 1:N
            dr = R_n(n1,:)-R_n(n2,:);
       dip_fct_1 = dip_fct_1 + c_nk(n1,k).*c_nk(n2,k).*dot(d_n(n1,:),d_n(n2,:));
       CD_fct_1 = CD_fct_1 + c_nk(n1,k).*c_nk(n2,k).*dot(d_n(n1,:),cross(d_n(n2,:), dr));
        end
    end
    dip_fct_1 = dip_fct_1/3; CD_fct_1 = CD_fct_1/3; 
    
     if nargout >3
         OD_ind(k,:) = om_r_rng.*OD_tmp*dip_fct_1;
         CD_ind(k,:) = om_r_rng.*OD_tmp*CD_fct_1;
         FL_ind(k,:) = om_r_rng.*FL_tmp*dip_fct_1;
         if k==N
          OD = sum(OD_ind,1);CD = sum(CD_ind,1);FL = sum(FL_ind,1);
         end
     else
    OD = OD + om_r_rng.*OD_tmp*dip_fct_1;
    FL = FL + om_r_rng.*FL_tmp*CD_fct_1;
    CD = CD + om_r_rng.*OD_tmp*dip_fct_1;
     end
end
      