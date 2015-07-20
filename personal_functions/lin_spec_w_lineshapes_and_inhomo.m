function  [OD,CD,FL,OD_ind,CD_ind,FL_ind] = lin_spec_w_lineshapes_and_inhomo...
    (w_ex,c_nk,g_t,lam,tau,d_kg, m_kg,P_k,t_rng,om_r_rng,sigma )
%Calculates absorption, Circular dichroism and Florescence linshapes
%tau = 1/sum(-R_k'k'kk,{k'~=k})
%sigma is the width of the Gaussian to integrate over to achieve the
%broadening
if nargin<11
   sigma = [];%no inhomobroadening  
end
if isempty(sigma)
   broad_fn = t_rng.^0; 
else %note if the broadening function is 
    broad_fn = exp(-t_rng.^2/sigma^2/2)/(sqrt(2*pi)*sigma);  
end
mid_freq = median(om_r_rng);
OD = zeros(size(om_r_rng)); CD = OD; FL = OD;
dt = t_rng(2)-t_rng(1);  tmax  = t_rng(end);
om_w_ft = 2*pi*(0:1/tmax:1/dt);
om_w_ft = fliplr(om_w_ft) -pi/dt + mid_freq;
N = length(w_ex);
d_kg = reshape(d_kg,N,3); m_kg= reshape(m_kg,N,3);
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
    
     if 1==1
   OD_tmp = (exp(1i*(mid_freq-w_ex(k))*t_rng-g_ex-t_rng/2/tau(k)).*broad_fn);
   OD_tmp = 2*real(fftshift(fft(OD_tmp )));            
   FL_tmp = P_k(k)*2*real(fftshift(fft(exp(1i*(mid_freq-w_ex(k))*t_rng...
                   +2i*lam_ex*t_rng - conj(g_ex) - t_rng/2/tau(k)).*broad_fn)));
               
    OD_tmp = interp1(om_w_ft,OD_tmp,om_r_rng);
    FL_tmp = interp1(om_w_ft,FL_tmp,om_r_rng);
     else
       for lp2 = 1:length(om_r_rng)
     OD_tmp(lp2) = real(trapz(t_rng,exp(1i*(om_r_rng(lp2)-w_ex(k))*t_rng-g_ex -t_rng/2/tau(k))));      
    FL_tmp(lp2) = P_k(k)*real(trapz(t_rng,exp(1i*(om_r_rng(lp2)-w_ex(k))*t_rng...
                    +2i*t_rng*lam_ex - conj(g_ex) - t_rng/2/tau(k))));           
       end
     end
     if nargout >3
         OD_ind(k,:) = om_r_rng.*OD_tmp*dot(d_kg(k,:),d_kg(k,:));
         CD_ind(k,:) = om_r_rng.*OD_tmp*dot(d_kg(k,:),m_kg(k,:));
         FL_ind(k,:) = om_r_rng.*FL_tmp*dot(d_kg(k,:),d_kg(k,:));
         if k==N
          OD = sum(OD_ind,1);CD = sum(CD_ind,1);FL = sum(FL_ind,1);
         end
     else
    OD = OD + om_r_rng.*OD_tmp*dot(d_kg(k,:),d_kg(k,:));
    FL = FL + om_r_rng.*FL_tmp*dot(d_kg(k,:),d_kg(k,:));
    CD = CD + om_r_rng.*OD_tmp*dot(d_kg(k,:),m_kg(k,:));
     end
end

      