%Used to calculate the third order response function for a dimer using
%Redfield theory and explicitly included, damped Harmonic oscillator modes
%get inverse cm units
Temp = 77; %temp in Kelvin
[convfact, beta,speed_unit]= inv_cm_unit_sys(Temp);
B = beta;
kpr = [0,0,1];  
%probe unit vectors

theta = atan(sqrt(2)); %angle between pump and probe
kpu = [sin(theta),0,cos(theta)];  
%pump unit vectors

%polarization vectors for pump beam
pol{1}=  [1,0,0];  pol{2} = pol{1} ; %others along x

simple_state_mix = true; %consider only direct mixings, none of the form 
% |g><e| ~> |g><e'| due to bath interactions etc

%either give half width at half max or sd
Hwhm_u_fs = 75; %half width at half maximum, pump
Hwhm_r_fs = 75; %half width at half maximum, probe

%I'm not sure this normalisation is a good idea, if I change the pulse
%width then I change the energy in the pulse... oh well
gaussian_pulses = true;
if gaussian_pulses
tau_u = (Hwhm_u_fs/sqrt(2*log(2)))/1000*convfact; %pumpp pulse SD in inverse cm
tau_r = (Hwhm_r_fs/sqrt(2*log(2)))/1000*convfact; %probe pulse SD in inverse cm    
E_u = @(t) exp(-t.^2/tau_u^2/2) / sqrt(tau_u * sqrt(pi)); %init e field env
%intensity of which is normed to 1
E_r = @(t) exp(-t.^2/tau_r^2/2) / sqrt(tau_r * sqrt(pi));
E_u_w = @(om) exp(-tau_u^2.*om.^2/2)*sqrt(tau_u/sqrt(pi));
E_u_inc = @(om,t) exp(-tau_u^2.*om.^2/2)*(1+1i*erfi(tau_u^2*om-1i*t)/sqrt(2)/tau_u)/2;
E_r_w = @(om) exp(-tau_r^2.*om.^2/2)*sqrt(tau_r/sqrt(pi));
else  %sech pulses
tau_u = Hwhm_u_fs/acosh(2)/2/1000*convfact; 
tau_r = Hwhm_r_fs/acosh(2)/2/1000*convfact;
E_u = @(t) sech(t/2/tau_u) / sqrt(4*tau_u); %init e field env, 
%intensity of which is normed to 1
E_r = @(t)  sech(t/2/tau_r) / sqrt(4*tau_r);  %init e field env
E_u_w = @(om) sech(pi*tau_u.*om)*sqrt(tau_u*pi/2);
%E_u_inc = @(om,t) don't know
E_r_w = @(om) sech(pi*tau_r.*om)*sqrt(tau_r*pi/2);   
    
end
%%
sys =4; %choose  system
[H_site,mu,R,lambda,gam,om_0,lam_dru, gam_dru,om_vib,displ,pdm,sd_shift] =sys_chooser(sys);
N = length(H_site); sz1 = N+1+N*(N-1)/2;
%lamtot = (cellfun(@sum ,lambda) + cellfun(@sum ,lam_dru));
om_r_rng = linspace(10^7./(700),10^7./(400),1000);
%set assumptions about how to solve, e.g. if modified Redfield used
use_markov = true; inc_double_ex = true; use_mod_red  =true;
site_dep_bath = false; %sites have same bath / diff baths
%write number of vibrational levels taken etc
%numvib = [4,4]; 
numvib=1
sz2 = prod(numvib);
points_to_save = 150000;
t1_range = linspace(0,15*tau_u,points_to_save);
t3_range = linspace(0,20*tau_r,points_to_save); g_cnst = zeros(N,1);
for n = 1:N
    g_cnst(n) = line_broad_fn_markov(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n});
end
Gamma_min = 1/min(real(g_cnst));
%construct time range to cover most of the interesting regions with similar
%weights
% t_int_rng = [linspace(0,B/(10*pi),2000),...
%              linspace(B/(10*pi),2*Gamma_min,10000),...
%             linspace(2*Gamma_min,10*Gamma_min,10000),...
%             linspace(10*Gamma_min,40*Gamma_min,5000)];
%  t_int_rng = t_int_rng([true,diff(t_int_rng)~=0]);
t_int_rng = linspace(0,20*Gamma_min,200000);
%calculate lin broadening functions in site basis
if site_dep_bath
g_t_1 = zeros(N,length(t1_range));g_t_3 = zeros(N,length(t3_range));
g_broad = zeros(N,length(t_int_rng)); 
g_deriv = g_broad;  g_sec_der = g_broad;  lam_tot = zeros(N,1);
t = sym('t','real');
for n =1:N

g_broad(n,:)   = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);    
g_deriv(n,:)  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);
g_sec_der(n,2:end) = line_broad_fn_sec_der(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng(2:end));
%just interpolate the last point, it diverges logarithmically in a way that
%is kind of pointless but screws up numerics
g_sec_der(n,1) = (3*g_sec_der(n,2)-g_sec_der(n,3))/2;
lam_tot(n) = sum(lambda{n}) + sum(lam_dru{n});
g_t_1(n,:)  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range);
g_t_3(n,:)  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t3_range);
end
else
    tic
 n=1; %bath at first site same as all others
 tol = 1e-13;
 g_broad   = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);    
g_deriv  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);
g_sec_der= line_broad_fn_sec_der(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);
%just interpolate the last point, it diverges logarithmically in a way that
%is kind of pointless but screws up numerics
g_sec_der(1) = (3*g_sec_der(2)-g_sec_der(3))/2;
lam_tot = sum(lambda{n}) + sum(lam_dru{n});
g_t_1  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range, tol);
g_t_3  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t3_range, tol); 
toc
end
[ex_basis,H0ex] = eig(H_site);
H0 = zeros(length(H_site)+1); H0(2:end,2:end) = H_site;
if site_dep_bath
H_site = H_site+diag(lam_tot);    
else
H_site = H_site+lam_tot*eye(N); %renormalise energies due to vib coupling
end
t_sep_range_ps = (0:10:4500)/1000;  Den_av = zeros(N,length(t_sep_range_ps));
%% Generate some Gaussian random noise for areas with static disorder

num_realisations = 6; %choose decent number of realisations. 

%If more are needed save old data and add more

%shifts of the site energies
site_shift = randn(num_realisations,N); %randn has variance 1, assumes indep
%site_shift = site_shift*0 %uncomment to remove static disorder
for j =  1:N
  site_shift(:,j) = site_shift(:,j)*sd_shift(j); %change variance
end
delta_w = sum(site_shift,2)/N; %mean shift to excited state manifold energy
site_shift_centred = site_shift - repmat(delta_w,1,N); %no centre broadening

%Could also put some static disorder in the relative dipole moments
%orientations, but I have no idea how and am pretty sure this wouldn't be
%independent of the site energies/interaction energies anyway
%Also if these can fluctuate they can be driven from equilibrium positions
%by interactions with a pump beam or interactions with the bath, leading to
%time dependent site-site interactions.  Leading to complex dynamics.

%om_u_rng = 1.0e+04 * [1.1910    1.2210    1.2410    1.2710];
om_u_rng = 1.0e+04*[1.9802];%,2.1];

points_to_save = 2^10+1; 

%Spp_alpha = zeros(length(om_r_rng),length(t_sep_rng),length(om_u_rng));
%Spp_CD = Spp_alpha; 

  
%%  Calculate linear spectra: 

num_linear = 2^16;
om_linear =  1.0e+04*linspace(1.3,2.5,num_linear);
tlin_range = linspace(0,20*Gamma_min,num_linear); tol = 1e-12;
if site_dep_bath
  g_lin_broad = zeros(N,length(tlin_range));  
    for n =1:N
 g_lin_broad(n,:) = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},tlin_range, tol); 
    end
 else
 g_lin_broad   = line_broad_fn_full(beta,gam{1},om_0{1},lambda{1},gam_dru{1},lam_dru{1},tlin_range, tol); 
end
mu_C_r = cross(mu,R);
  %need to convolve with Gaussian of appropriate width to deal with some
  %broadening due to static disorder
  centre_SD = mean(sd_shift)/sqrt(N);
  %shift_PD = @(w) exp(-w.^2/centre_SD^2/2)/sqrt(2*pi)/centre_SD;
  R_reduced_save = {} ; R_reduced_save{num_realisations} =[];
%% This part is effected by static disorder loops,
 OD_av =zeros(N,length(om_linear)); CD_av =OD_av ; FL_av = OD_av ;
 OD_av2 = OD_av; CD_av2 =OD_av ; FL_av2 = OD_av ; %test values
  for real_lp = 1:num_realisations %loop over realisations
      tic
  %shift energies with a centred shift (no mean drift, dealt with via a
  %convolution later
      H_site_shifted = H_site + diag(site_shift_centred(real_lp,:)); %shift site E
    [c_nk,E_ex] = eig(H_site_shifted); %diagonalise Hamiltonian
     E_ex = diag(E_ex);   c_nk = c_nk.'; %opposite definition for next bit
     E_ex = E_ex - ((cellfun(@sum ,lambda) + cellfun(@sum ,lam_dru))*c_nk.^4).';
if ~isempty(R_reduced_save{real_lp} )
     R_reduced = R_reduced_save{real_lp};
else     
     
     if use_mod_red %generator of dyanmics
    R_reduced =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,B,...
                      c_nk,E_ex,t_int_rng);
 %  R_reduced =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,B,...
 %                     c_nk,E_ex,t_int_rng,5e-3); %last argument is working tol
     else
     R_red = redfield_calc(H_site_shifted, fock_space_rep,beta,gam_dru,...
            lam_dru,gam,lambda,om_0,use_markov);   
         R_reduced = zeros(N,N);
        for k=1:N 
            for kk = 1:N              
  R_reduced(k,kk) = R_red(k,k,kk,kk);       
            end
        end
     end
   
       lg = R_reduced>0;lg2 = logical(eye(size(R_reduced)));
       if any(any(lg & ~lg2))

           warning('positive elements in Redfield off diag removed')
     R_reduced(R_reduced>0) = 0;   
       end    
    R_reduced  = R_reduced  - diag(sum(R_reduced,1)); %balance condition  
    
    R_reduced_save{real_lp} = R_reduced;
end 
 tau = 1./diag(R_reduced);

   eqstate = expm(-R_reduced*10*tlin_range(end)); %generic long timescale  
   eqstate = sum(eqstate,2)/N; %actual equilibrium state
    
mu_ex = zeros(3,N); %transition dipoles between exciton states
mag_ex = zeros(3,N); %mag transition dipoles between exciton states

mu_ex_test = zeros(3,N,1+N*(N-1)/2); 
for k = 1:N
    mu_ex(:,k) = mu.'*c_nk(:,k);
    mag_ex(:,k) = 1i*pi*E_ex(k)*mu_C_r.'*c_nk(:,k); %also factor of freq
end
 for k = 1:N
     mu_ex(:,k) = c_nk(k,:)*mu;
     mag_ex(:,k) = 1i*c_nk(k,:)*mu_C_r; %also factor of freq
 end  

[~,~,~,OD_ind,CD_ind,FL_ind,extras] = lin_spec_w_lineshapes(E_ex,c_nk,g_lin_broad,lam_tot,tau,mu_ex,...
                                   mag_ex,eqstate,tlin_range,om_linear,{mu,R});
OD_ind2 = extras{1};    CD_ind2 = extras{2};    FL_ind2 = extras{3};                          
% [~,~,~,OD_ind2,CD_ind2,FL_ind2] = lin_spec_w_lineshapes2(E_ex,c_nk,g_lin_broad,lam_tot,tau,...
%                             mu,R,eqstate,tlin_range,om_linear);                               
                             
                               %deal with broadening on entire signal
%[OD,CD,FL] = lin_spec_w_lineshapes_and_inhomo(E_ex,c_nk,g_lin_broad,lam_tot,tau,mu_ex,...
%                                   mag_ex,eqstate,tlin_range,om_linear,sigma);                               
 OD_av = OD_av + OD_ind; CD_av = CD_av + CD_ind; FL_av = FL_av + FL_ind;
  OD_av2 = OD_av2 + OD_ind2; CD_av2 = CD_av2 + CD_ind2; FL_av2 = FL_av2 + FL_ind2;
toc
%R_reduced_save{real_lp}=R_reduced;
  end
  CD_av2 = -CD_av2;
  CD_av  = -imag(CD_av );
  %% Plot if needed
 
[broaden_OD,renorm_fct1] = Gauss_broad_fn(centre_SD,om_linear,OD_av2,2);
[broaden_FL,renorm_fct2] = Gauss_broad_fn(centre_SD,om_linear,FL_av2,2);
[broaden_CD,renorm_fct3] = Gauss_broad_fn(centre_SD,om_linear,CD_av2,2);
  wl_nm = 10^7./om_linear;
  Alpha = sum(broaden_OD)/num_realisations;  
  CD = sum(broaden_CD)/num_realisations;
  
figure
%plotyy(wl_nm,[broaden_OD/num_realisations;Alpha],...
%    wl_nm,N/num_realisations*[broaden_FL;sum(broaden_FL)])
plotyy(wl_nm,Alpha,wl_nm,N/num_realisations*sum(broaden_FL))
xlabel('wavelength (nm)')
figure
%plot(wl_nm,real(broaden_CD)/num_realisations,'--')
hold on
plot(wl_nm,real(CD))
xlabel('wavelength (nm)')
ylabel('CD amplitude')
%% find local maxima and minima in spectra
peaks_abs = localMaximum(Alpha);
peaks_abs = peaks_abs(peaks_abs ~= 1 &  peaks_abs ~= length(wl_nm));
peaks = localMaximum(CD);
troughs = localMaximum(-CD);
peaks = peaks(peaks ~= 1 & peaks ~= length(wl_nm));
troughs = troughs(troughs ~= 1 & troughs ~= length(wl_nm));

min_CD =CD(troughs);  max_CD = CD(peaks);
max_abs = Alpha(peaks_abs);

  
%%
 n=1; %bath at first site same as all others
 tol = 1e-13;
% g_broad   = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);    
%g_deriv  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);
%g_sec_der = line_broad_fn_sec_der(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);
%lam_tot = sum(lambda{n}) + sum(lam_dru{n});
g_t_door  = line_broad_fn_full(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range, tol);
g_t  = line_broad_fn_full(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t3_range, tol); 
 t_sep_range = t_sep_range_ps*convfact;
D_mat_av = zeros(N,N,length(t_sep_range_ps)); %only for first pump freq
 
use_HEOM = false;
if use_HEOM
    Kappa=5; QQ_topass = zeros(N,2);
    for j = 1:N
     [cc1,cc2R,cc2I,vv1,vv2,QQ] = coeffients_from_brownian_new...
    (lambda{j},gam{j},om_0{j},Temp,Kappa,lam_dru{j},gam_dru{j});
  cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
 cc_acom{j}= [cc2I,cc1*0]; QQ_topass(j,:) = [QQ,0];
    end
    kap2 = 2;
    [ Qtrunc,H_prop_op,nn]=HEOM_propogator2(QQ_topass,cc_com,cc_acom,vv,kap2,1);
end

%%
for real_lp = 1: num_realisations
    %% calculate system parameters for this realisation of disorder
H_site_shifted = H_site + diag(site_shift(real_lp,:));
[H_ex_vib,fock_space_rep,~,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(H_site_shifted,om_vib,numvib,displ,mu,[]) ;
M_prj = indiv_op{5};

E = diag(H_exciton); 
E1 = E(2:N+1); E2 = E(N+2:end);

  c_nk = zeros(N); %as I have defined it c_nk = Proj(k,n)
  M_prj_sgl = M_prj(2:N+1,2:N+1);
  c_nm_f = zeros(N,N,N*(N-1)/2);
  M_prj_dbl = M_prj(N+2:end,N+2:end); fock_dbl = fock_space_rep(N+2:end,:);
      cnt=1;
      for m = 1:N
          c_nk(m,:) = M_prj_sgl(m,:);%M_prj_sgl(:,m);
      end
      for  f = 1:length(M_prj_dbl)
          tmp = find(fock_space_rep(N+1+f,:));
          n = tmp(1); m = tmp(2); %n always < m
              c_nm_f(n,m,:) =  M_prj_dbl(f,:); %M_prj_dbl(:,f);
      end 
      
 %shift energies back to bare exciton energies for use of redfield theory
%  E1 = E1 - ((cellfun(@sum ,lambda) + cellfun(@sum ,lam_dru))*c_nk.^4).';
% pickr = eye(N); 
%  for f = 1:length(M_prj_dbl)
%   coeff = zeros(N,1);   
%         for m=2:N
%             for n = 1:m-1
%                 for mm = 2:N
%                     for nn = 1:mm-1
%                         if n==nn || n==mm 
%                          coeff = coeff + (c_nm_f(n,m,f)*c_nm_f(nn,mm,f))^2*pickr(:,n);                       
%                         end
%                         if  m==mm || m==nn
%                          coeff = coeff + (c_nm_f(n,m,f)*c_nm_f(nn,mm,f))^2*pickr(:,m);                       
%                         end                                               
%                     end
%                 end       
%             end
%         end
% 
%   E2(f) = E2(f) - (cellfun(@sum ,lambda) + cellfun(@sum ,lam_dru))*coeff;       
%  end 
   
mu_C_r = cross(mu,R);  
mu_ex = zeros(3,N,1+N*(N-1)/2); %transition dipoles between exciton states
mu_cross_R_ex = zeros(3,N,1+N*(N-1)/2); %mag transition dipoles between exciton states

mu_ex_test = zeros(3,N,1+N*(N-1)/2); 
for k = 1:N
   
   % mu_ex_test(:,k,1) = c_nk(k,:)*mu;
    mu_ex(:,k,1) = mu.'*c_nk(:,k);
    mu_cross_R_ex(:,k,1) = H_exciton(k+1,k+1)*mu_C_r'*c_nk(:,k);
    
    for f= 1:N*(N-1)/2
        tmp = zeros(1,3); %tmp_test  = tmp; 
        tmp2 = tmp;
        for m = 2:N
            for n = 1:m-1      
        tmp = tmp + c_nm_f(n,m,f).*(c_nk(n,k)*mu(m,:) + c_nk(m,k)*mu(n,:));
        tmp2 = tmp2+ c_nm_f(n,m,f).*(c_nk(n,k)*mu_C_r(m,:) + c_nk(m,k)*mu_C_r(n,:));
            end
        end
        mu_ex(:,k,1+f) = tmp; %mu_ex_test(:,k,1+f) = tmp_test;
        mu_cross_R_ex (:,k,1+f) = tmp2;
    end
end

% [alpha_av,CD_av,alpha_lin_pol] = ...
%   dipole_fourth_order_av(mu,R,pol_L,pol_R,pol_linear,kk1,true);
pol_linear = [pol{1};pol{2};[1,0,0];[1,0,0]];
pol_L = [pol{1};pol{2};[1/sqrt(2),+1i/sqrt(2),0];[1/sqrt(2),-1i/sqrt(2),0]];
pol_R = [pol{1};pol{2};[1/sqrt(2),-1i/sqrt(2),0];[1/sqrt(2),+1i/sqrt(2),0]];
         alpha_ex_L = zeros(N,N,N*(N-1)/2+1); alpha_ex_R = zeros(N,N,N*(N-1)/2+1);
          xx_av_ex = zeros(N,1);  xy_av_ex =  xx_av_ex;
alpha_lin_pol_ex = alpha_ex_L;
for j1 = 1:N %pre save these
    xx_av_ex(j1) = dot(mu_ex(:,j1,1),mu_ex(:,j1,1))/3;
    xy_av_ex(j1) = dot(mu_ex(:,j1,1),mu_cross_R_ex(:,j1,1))/3;
   for j2 = 1:N
       for f2 = 1:N*(N-1)/2+1     
       mu_set=[mu_ex(:,j1,1),mu_ex(:,j1,1),mu_ex(:,j2,f2),mu_ex(:,j2,f2)].';
       
   alpha_lin_pol_ex(j1,j2,f2) = tensor_av(mu_set,pol_linear); 
    alpha_ex_L(j1,j2,f2) = tensor_av(mu_set,pol_L);
    alpha_ex_R(j1,j2,f2) = tensor_av(mu_set,pol_R);          
            
        end    
    
   end
end     
LL = length(H_exciton);
 manifold_chk = round(sum(double(fock_space_rep),2));
 %Q_j(1+a,1+b,n) = c_nk(a,n).*c_nk(b,n)
  Q_j = zeros(LL,LL,N); %mixes only within
 %the same excitation manifold
 %as I have defined it c_nk = ex_st(k,n)
 for lp = 1:N
     for lp2 = 2:LL
         occ_L = M_prj(:,lp2).*fock_space_rep(:,lp); %occupation of this state on bra
        for lp3 = 2:LL
            occ_R = M_prj(:,lp3).*fock_space_rep(:,lp); %occupation of this state on ket
            if manifold_chk(lp2) == manifold_chk(lp3) %else no mixing
                    Q_j(lp2,lp3,lp) =  occ_L'*occ_R;
            end
        end
     end
 end     
 tmp1 = reshape(Q_j,LL,LL,N,1,1,1);
 QQ = squeeze(mtimesx(permute(tmp1,[4,3,1,2,5,6]),permute(tmp1,[3,4,5,6,1,2])));
%this is the matrix QQ(a,b,c,d) := sum_n <a | Q_n |b> <c|Q_n |d>
% QQ(1+a,1+b,1+c,1+d) = sum_n c_nk(a,n)*c_nk(b,n)*c_nk(c,n)*c_nk(d,n), a<N 

%% calc spec params
use_fft = true;
[D_t,W_t,Wf_t,D0,Wgg,W1ee,W2ee] = spec_fun_simple(...
          t1_range,t3_range,om_r_rng,om_u_rng,tau_r,tau_u,E1,E2,...
          g_t,g_t_door,lam_tot,c_nk,c_nm_f,N,true,use_fft) ;

   
%% use modified Redfield theory
if ~use_HEOM
if ~exist('R_reduced_save','var') %presaved
   R_mod_red =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,B,...
                       M_prj(2:N+1,2:N+1)',diag(H_exciton(2:N+1,2:N+1)),t_int_rng,1e-2);
else
     R_mod_red = R_reduced_save{real_lp};
end
 %    R_mod_red =   mod_redfield_calc3(B,gam{1},om_0{1},lambda{1},gam_dru{1},...
 %                lam_dru{1},c_nk,E_ex,1e-9,5e-4);  
       lg = R_mod_red>0;lg2 = logical(eye(size(R_mod_red)));
       if any(any(lg &~lg2))
           warning('positive elements in Redfield off diag removed')
    R_mod_red(R_mod_red>0) = 0;
    R_mod_red  = R_mod_red  - diag(sum(R_mod_red,1)); %balance condition     
       end
end
 %% prop doorway function in delay time    
    poponly = true; ode_toler = 1e-12;
 
if use_HEOM
    for lp  =1:size(D0,1)
    rho_0 = diag(D0(1,:)); rho_0 = M_prj(2:N+1,2:N+1)*rho_0*M_prj(2:N+1,2:N+1)';
    tr_save = trace(rho_0); rho_0 = rho_0 / tr_save;
  [Time_units,rho_vec]=HEOM_exmanifold_solver(H_ex_vib(2:N+1,2:N+1),...
            Qtrunc,nn,H_prop_op,rho_0,length(t_sep_range),t_sep_range_ps(end),1);  
      rho_vec2 = reshape(rho_vec.',N,N,size(rho_vec,1));  
      rho_vec2 = mtimesx(M_prj(2:N+1,2:N+1),'C',rho_vec2);
      rho_vec2 = mtimesx(rho_vec2,M_prj(2:N+1,2:N+1));
      rho_vec2 = reshape(rho_vec2,N^2,size(rho_vec,1));
    end
else
       %pointless reshape to fit function I am using
    rho_neq_g = reshape(-D0,1,1,size(D0,1),size(D0,2));
    rho_neq_e = reshape(D0,1,1,size(D0,1),size(D0,2));
    
[rho_t_g,rho_t_e] = doorway_prop_basic(rho_neq_g,rho_neq_e,...
                t_sep_range,0,-R_mod_red,sz2,N,ode_toler,poponly ) ; 
            
 t_prop_mat = zeros(N,N,length(t_sep_range),length(om_u_rng)); %props in time
 D_t_ev = t_prop_mat; %not actually identical due to numerics and such like

for lp = 1:length(t_sep_range)
    t_prop_mat(:,:,lp) = expm(-R_mod_red*t_sep_range(lp));
    for lp2 = 1:length(om_u_rng)
         D_t_ev(:,:,lp,lp2) = t_prop_mat(:,:,lp)*diag(D0(lp2,:));
%          for lp3 = 1:N
%              tmp = 0*D0(lp2,:); tmp(lp3) = D0(lp2,lp3);
%              D_t_ev_test(:,lp3,lp,lp2) = t_prop_mat(:,:,lp)*tmp.';
%          end
    end
end                               
end


%calculate actual populations of ssecond order densisty matrix
D_actual = zeros(size(rho_t_e,1),size(rho_t_e,2),size(rho_t_e,3));

for j=1:N
D_actual(:,:,:) = D_actual(:,:,:) + xx_av_ex(j)*rho_t_e(:,:,:,1,j);
%D_actual(:,:,:) = D_actual(:,:,:) + xx_ex_av(j,j)*rho_t_e(:,:,:,1,j);
end
%comment out to keep in exciton basis
D_actual = mtimesx(M_prj(2:N+1,2:N+1),D_actual ); %proj to site basis
D_actual = mtimesx(D_actual,M_prj(2:N+1,2:N+1),'C');

D_mat_av = D_mat_av + D_actual / trace(D_actual(:,:,1));
%%  calculate transient spec            
   if 1==1 %old method
rho_t_ef_dip_av =  zeros(size(rho_t_e,3),size(rho_t_e,4),N,1+ N*(N-1)/2);
rho_t_g_dip_av =  zeros(size(rho_t_g,1),size(rho_t_g,2),N);
 %density matrix with weights coming from window prefactors as well
 %ground / double excited manifold transitions determine last index
for  lp = 1:N %section loop

    tmpe = zeros(size(rho_t_e,3),size(rho_t_e,4),1+ N*(N-1)/2,N);
    tmpg = zeros(size(rho_t_g,1),size(rho_t_g,2));
   for lp2 = 1:N %dipole element loop
%        tmpg = tmpg + rho_t_g(:,:,lp)*alpha_lin_pol_ex(lp,lp2,1);
% 
%        for f = 1:N*(N-1)/2+1
%            %excited state pop can transfer to ground (giving simulated 
%            %emission) or to the double ex-manifold -> absorption
%            %goes from state lp2 to state f, initial trans was to state lp
%            tmpe(:,:,f,lp2) =  rho_t_e(lp2,lp2,:,:,lp)*alpha_lin_pol_ex(lp,lp2,f);
%        %hole population can transfer to anything in exstate manifold and
%        %back
       tmpg = tmpg + rho_t_g(:,:,lp2)*alpha_lin_pol_ex(lp2,lp,1);

       for f = 1:N*(N-1)/2+1
           %excited state pop can transfer to ground (giving simulated 
           %emission) or to the double ex-manifold -> absorption
           %goes from lp
           tmpe(:,:,f,lp2) =  rho_t_e(lp,lp,:,:,lp2)*alpha_lin_pol_ex(lp2,lp,f);
           
        end
   end 
   tmpe = sum(tmpe,4);
    rho_t_ef_dip_av(:,:,lp,:)= tmpe;
    rho_t_g_dip_av(:,:,lp)= tmpg; 
end

GSB_sig2 = mtimesx(Wgg,permute(rho_t_g_dip_av,[3,1,2]));
SE_sig2 = -mtimesx(W1ee,permute(rho_t_ef_dip_av(:,:,:,1),[3,1,2]));
ESA_sig2 = zeros(size(SE_sig2));
for f = 2:N*(N-1)/2+1 %pathways featuring double ex states
    ESA_sig2  =  ESA_sig2  + mtimesx(W2ee(:,:,f-1),...
        permute(rho_t_ef_dip_av(:,:,:,1),[3,1,2]));
end
om_r_fct = reshape(om_r_rng,length(om_r_rng),1);
GSB_sig2  = GSB_sig2.*om_r_fct; 
om_r_fct = repmat(om_r_fct,1,length(t_sep_range),length(om_u_rng));
SE_sig2  = SE_sig2.*om_r_fct;
ESA_sig2  = ESA_sig2.*om_r_fct;

   end
%% Calculate the averaged density matrix
rho_t_g = squeeze(rho_t_g);  
D_dipole_av = zeros(N,N*(N-1)/2+1,length(t_sep_range),1);%length(om_u_rng));
%Dav_kk^(q)(tau,om_u) = sum_k' D_kk^(k')(tau,om_u)
%                       *<d_k'g.e1 d_k'g.e1 d_kq.e2 dkq.2>
rho_diag_only = zeros(N,size(rho_t_e,3),size(rho_t_e,4),size(rho_t_e,5));
for k = 1:N
    rho_diag_only(k,:,:,:) = rho_t_e(k,k,:,:,:);
end
for k =1:N %second exciton state
    for q = 1:N*(N-1)/2+1 %state 
        tmp= zeros(1,size(rho_diag_only,2),size(rho_diag_only,3));
        for kprime = 1:N %sum over initial exciton transfers
     tmp = tmp + rho_diag_only(k,:,:,kprime)*alpha_lin_pol_ex(kprime,k,q);
        end
      % tmp = squeeze(mtimesx(D_t_ev,alpha_lin_pol_ex(:,k,q)));
        D_dipole_av(k,q,:,:) = tmp;
       %D_dipole_av(:,q,:,:) = tmp;
    end
end

%% Calculate different signal contributions
%ground state bleaching
GSB_sig = -squeeze(mtimesx(Wgg,D_dipole_av(:,1,1,:))); %not time dep, just need init time
om_r_fct = reshape(om_r_rng,length(om_r_rng),1);
GSB_sig = GSB_sig.*repmat(om_r_fct,1,size(GSB_sig,2));

%Stimulated emission
SE_sig = -mtimesx(W1ee,D_dipole_av(:,1,:,:)); %1 is the ground state

%Excited state absorption
ESA_sig = zeros([size(SE_sig),N*(N-1)/2]);
for f = 2:N*(N-1)/2+1 %pathways featuring double ex states
    ESA_sig(:,:,:,f-1)  =  mtimesx(W2ee(:,:,f-1),D_dipole_av(:,f,:,:));
end

om_r_fct = repmat(om_r_fct,1,length(om_u_rng),length(t_sep_range));
SE_sig  = SE_sig.*om_r_fct;
om_r_fct = repmat(om_r_fct,1,1,1,N*(N-1)/2);
ESA_sig  = ESA_sig.*om_r_fct;
clear om_r_fct
SE_sig = squeeze(SE_sig); ESA_sig = squeeze(ESA_sig); %if om_u_rng = 1 pointless
%% Collect together
if real_lp==1
  GSB_sig_total = GSB_sig ;
SE_sig_total = SE_sig ;
 ESA_sig_total = sum(ESA_sig,ndims(ESA_sig)) ;  
else
GSB_sig_total = GSB_sig_total + GSB_sig ;
SE_sig_total = SE_sig_total + SE_sig ;
 ESA_sig_total =  ESA_sig_total+ sum(ESA_sig,ndims(ESA_sig));
end
if real_lp==1
  GSB_sig_total2 = GSB_sig2 ;
SE_sig_total2 = SE_sig2 ;
 ESA_sig_total2 =  ESA_sig2;  
else
GSB_sig_total2 = GSB_sig_total + GSB_sig2 ;
SE_sig_total2 = SE_sig_total + SE_sig2 ;
 ESA_sig_total2 =  ESA_sig_total+ ESA_sig2;
end
end

%%
totsig = (SE_sig_total + repmat(GSB_sig_total,1,length(t_sep_range_ps))...
    + ESA_sig_total)/num_realisations ;
% totsig = (SE_sig2 + repmat(GSB_sig2,1,length(t_sep_range_ps))...
%     + ESA_sig2)/num_realisations ;
broaden_sig = totsig;
%[broaden_sig,~] = Gauss_broad_fn(centre_SD,om_r_rng,totsig,1);
lam_rng = 10^7./om_r_rng;
lg = lam_rng<600 & lam_rng>500;
t_sep_ps = t_sep_range/convfact;
time_select  = [0.22,0.42,0.56,1.01,2.12,max(t_sep_ps)];
lg2 = time_select*0;
for k =1:length(time_select)
    lg2(k) = find(time_select(k)<=t_sep_ps,1,'first');
end

figure1 = figure;
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'all');
norm_fct = max(max(abs(broaden_sig(:,:,1))));
plot1 = plot(lam_rng(lg),broaden_sig(lg,lg2,1)/norm_fct,'Parent',axes1);
for k =1:length(time_select) 
set(plot1(k),'DisplayName',strcat(num2str(time_select(k)),'ps'));
end
%%
figure
plot(lam_rng(lg),SE_sig_total(lg,lg2,1)/norm_fct);
hold on
plot(lam_rng(lg),GSB_sig_total(lg)/norm_fct,'k');
plot(lam_rng(lg),ESA_sig_total(lg,lg2,1)/norm_fct,'--');
%%
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);

% figure
% pcolor(t_sep_range/convfact,om_r_rng,real(GSB_sig_total (:,:,1))) 
% shading flat
% set(gcf, 'renderer', 'zbuffer');
% xlabel('Time delay (fs)')
% ylabel('probe frequency \omega, cm^{-1}')  
%  colormap(CMRmap)
% colorbar
figure
plot(om_r_rng,real(GSB_sig_total (:,1,1))) 
xlabel('probe frequency \omega, cm^{-1}')
ylabel('GSB sig')  
hold on
plot(om_r_rng,real(SE_sig_total (:,1,1)),'r') 

figure
pcolor(t_sep_range/convfact,om_r_rng,real(SE_sig_total(:,:,1))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (ps)')
ylabel('probe frequency \omega, cm^{-1}')  
 colormap(CMRmap)
colorbar
figure
pcolor(t_sep_range/convfact,om_r_rng,real(ESA_sig_total(:,:,1))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (ps)')
ylabel('probe frequency \omega, cm^{-1}')  
 colormap(CMRmap)
colorbar

%%

tmp = reshape(D_mat_av,N^2,length(t_sep_range));
tmp = tmp(1:N+1:end,:);
tmp = tmp./sum(tmp(:,1));
figure
plot(t_sep_range_ps,tmp)
xlabel('\tau (ps)')
ylabel('Average site population (normalised)')


% D_mat_ex = mtimesx(M_prj(2:N+1,2:N+1),'C',D_mat_av); %proj to ex basis
% D_mat_ex = mtimesx(D_mat_ex,M_prj(2:N+1,2:N+1));
% figure
% tmp = reshape(D_mat_ex,N^2,length(t_sep_range));
% tmp = tmp(1:N+1:end,:);
% tmp = tmp./sum(tmp(:,1));
% plot(t_sep_range_ps,tmp)
% xlabel('\tau (ps)')
% ylabel('Average exciton population (normalised)')
% hold on
% tmp = eig(H_site); tmp = tmp-tmp(1);
% plot(t_sep_range_ps,repmat(exp(-beta*tmp)/sum(exp(-beta*tmp)),1,length(t_sep_range)),'--');
% 

