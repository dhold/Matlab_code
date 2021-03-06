%Used to calculate the third order response function for a dimer using
%Redfield theory and explicitly included, damped Harmonic oscillator modes
%get inverse cm units
Temp = 77; %temp in Kelvin
[convfact, beta,speed_unit]= inv_cm_unit_sys(Temp);
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
om_r_rng = linspace(10^7./(700),10^7./(400),1000);
%set assumptions about how to solve, e.g. if modified Redfield used
use_markov = true; inc_double_ex = true; use_mod_red  =true;
site_dep_bath = false; %sites have same bath / diff baths
%write number of vibrational levels taken etc
%numvib = [4,4]; 
numvib=1
sz2 = prod(numvib);
points_to_save = 10000;
t1_range = linspace(0,10*tau_u,points_to_save);
t3_range = linspace(0,20*tau_r,points_to_save); g_cnst = zeros(N,1);
for n = 1:N
    g_cnst(n) = line_broad_fn_markov(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n});
end
Gamma_min = 1/min(real(g_cnst));
t_int_rng = [linspace(eps,2*Gamma_min,5000),...
            linspace(2*Gamma_min,20*Gamma_min,5000),...
            linspace(20*Gamma_min,150*Gamma_min,5000)];
 t_int_rng = t_int_rng([true,diff(t_int_rng)~=0]);
 
 %t_int_rng = linspace(0,50*Gamma_min,50000);
%calculate lin broadening functions in site basis
if site_dep_bath
g_t_1 = zeros(N,length(t1_range));g_t_3 = zeros(N,length(t3_range));
g_broad = zeros(N,length(t_int_rng)); 
g_deriv = g_broad;  g_sec_der = g_broad;  lam_tot = zeros(N,1);
t = sym('t','real');
for n =1:N

g_broad(n,:)   = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);    
g_deriv(n,:)  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);
g_sec_der(n,:) = line_broad_fn_sec_der(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);
lam_tot(n) = sum(lambda{n}) + sum(lam_dru{n});
g_t_1(n,:)  = line_broad_fn_clderiv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range);
g_t_3(n,:)  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t3_range);
end
else
    tic
 n=1; %bath at first site same as all others
 tol = 1e-13;
 g_broad   = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);    
g_deriv  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);
g_sec_der = line_broad_fn_sec_der(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);
lam_tot = sum(lambda{n}) + sum(lam_dru{n});
g_t_1  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range, tol);
g_t_3  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t3_range, tol); 
toc
end
[ex_basis,H0ex] = eig(H_site);
H0 = zeros(length(H_site)+1); H0(2:end,2:end) = H_site;
H_site = H_site-lam_tot*eye(N); %renormalise energies

       t_sep_range_ps = (0:10:4500)/1000;  Den_av = zeros(N,length(t_sep_range_ps));
%% Generate some Gaussian random noise for areas with static disorder

num_realisations = 3; %choose decent number of realisations. 

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
om_u_rng = 1.0e+04*[1.9802,2.1];

points_to_save = 2^10+1; 

%Spp_alpha = zeros(length(om_r_rng),length(t_sep_rng),length(om_u_rng));
%Spp_CD = Spp_alpha; 

%%  Calculate dipole averages for the pump probe geometry
%Assuming no static disorder impacts the dipoles (apart from (random) 
%molecular alignment) we can precompute these values.  The actual values 
%will be averaged over exciton states though, which depend on site energies

 typ_op = 10^7/600; %typical frequency
k_u = 2*pi*abs(typ_op)*kpu; k_r = 2*pi*abs(typ_op)*kpr; %average amplitudes of these para
kk1 = [-k_u;k_u;k_r;-k_r]; kk2 = [k_u;-k_u;k_r;-k_r]; %order of interaction (left to right)
%take probe polarization wavevectors along x
%pol{3} = [1,0,0]; % pol_left = [1;-1i;0]/sqrt(2);  pol_right = [1;1i;0]/sqrt(2);


%pol{1} = [1/sqrt(2),1/sqrt(2),0]; pol{2} = pol{1} ; %45 degrees to probe apparently
pol_linear = [pol{1};pol{2};[1,0,0];[1,0,0]];
pol_L = [pol{1};pol{2};[1/sqrt(2),+1i/sqrt(2),0];[1/sqrt(2),-1i/sqrt(2),0]];
pol_R = [pol{1};pol{2};[1/sqrt(2),-1i/sqrt(2),0];[1/sqrt(2),+1i/sqrt(2),0]];

 [alpha_av,CD_av,alpha_lin_pol] = ...
   dipole_fourth_order_av(mu,R,pol_L,pol_R,pol_linear,kk1,false); 
   
   sz  = sz1*sz2; %total size of matrix
   
%% Two body averages, required for linear spec and Neq density populations

xx_av = zeros(N,N); yx_av = zeros(N,N);
for j1 = 1:N
    for j2 = 1:N
        xx_av(j1,j2) = dot(mu(j1,:),mu(j2,:))/3;
        %transition prop to k dot r, |k| omega = c-> |k| =
        %c/(omega_transition)
        yx_av(j1,j2) = 1i*pi*H_site(j2,j2)*dot(mu(j1,:),cross(mu(j2,:),R(j1,:)-R(j2,:)))/6;
    end
end
  
%%  Calculate linear spectra: 

num_linear = 2^15;
om_linear =  1.0e+04*linspace(1.3,2.5,num_linear);
tlin_range = linspace(0,10*Gamma_min,num_linear); tol = 1e-12;
if site_dep_bath
  g_lin_broad = zeros(N,length(tlin_range));  
    for n =1:N
 g_lin_broad(n,:)   = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},tlin_range, tol); 
    end
 else
 g_lin_broad   = line_broad_fn_full(beta,gam{1},om_0{1},lambda{1},gam_dru{1},lam_dru{1},tlin_range, tol); 
end
mu_C_r = cross(mu,R);
  %need to convolve with Gaussian of appropriate width to deal with some
  %broadening due to static disorder
  centre_SD = mean(sd_shift)/sqrt(N);
  %shift_PD = @(w) exp(-w.^2/centre_SD^2/2)/sqrt(2*pi)/centre_SD;
%% This part is effected by static disorder loops,
 OD_av =zeros(N,length(om_linear)); CD_av =OD_av ; FL_av = OD_av ;
  for real_lp = 1:num_realisations %loop over realisations
  %    tic
  %shift energies with a centred shift (no mean drift, dealt with via a
  %convolution later
      H_site_shifted = H_site + diag(site_shift_centred(real_lp,:)); %shift site E
    [c_nk,E_ex] = eig(H_site_shifted); %diagonalise Hamiltonian
     E_ex = diag(E_ex);  c_nk = c_nk.'; %opposite definition for next bit
     if use_mod_red %generator of dyanmics
   R_reduced =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,...
                      c_nk,E_ex,t_int_rng);
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
           
           R_reduced
           warning('positive elements in Redfield of diag removed')
     R_reduced(R_reduced>0) = 0;   
       end    
    R_reduced  = R_reduced  - diag(sum(R_reduced,1)); %balance condition         
 tau = zeros(N,1);
     for k = 1:N
     tau(k) = -sum(R_reduced([1:k-1,k+1:N],k)); %decay rates
     end    
   eqstate = expm(-R_reduced*10*tlin_range(end)); %generic long timescale  
   eqstate = sum(eqstate,2)/N; %actual equilibrium state
    
mu_ex = zeros(3,N); %transition dipoles between exciton states
mag_ex = zeros(3,N); %mag transition dipoles between exciton states

mu_ex_test = zeros(3,N,1+N*(N-1)/2); 
for k = 1:N
    mu_ex(:,k) = mu.'*c_nk(:,k);
    mag_ex(:,k) = 1i*pi*E_ex(k)*mu_C_r.'*c_nk(:,k); %also factor of freq
end
% for k = 1:N
%     mu_ex(:,k) = c_nk(k,:)*mu;
%     mag_ex(:,k) = 1i*c_nk(k,:)*mu_C_r; %also factor of freq
% end  

[~,~,~,OD_ind,CD_ind,FL_ind] = lin_spec_w_lineshapes(E_ex,c_nk,g_lin_broad,lam_tot,tau,mu_ex,...
                                   mag_ex,eqstate,tlin_range,om_linear);
                               %deal with broadening on entire signal
%[OD,CD,FL] = lin_spec_w_lineshapes_and_inhomo(E_ex,c_nk,g_lin_broad,lam_tot,tau,mu_ex,...
%                                   mag_ex,eqstate,tlin_range,om_linear,sigma);                               
 OD_av = OD_av + OD_ind; CD_av = CD_av + CD_ind; FL_av = FL_av + FL_ind;
% toc
  end
  %% Plot if needed
 
[broaden_OD,renorm_fct1] = Gauss_broad_fn(centre_SD,om_linear,OD_av,2);
[broaden_FL,renorm_fct2] = Gauss_broad_fn(centre_SD,om_linear,FL_av,2);
[broaden_CD,renorm_fct3] = Gauss_broad_fn(centre_SD,om_linear,CD_av,2);
  wl_nm = 10^7./om_linear;
  Alpha = sum(broaden_OD)/num_realisations;  
  CD = sum(imag(broaden_CD))/num_realisations;
  
figure
%plotyy(om_linear,OD_av/num_realisations,om_linear,N*FL_av/num_realisations)
plotyy(wl_nm,Alpha,wl_nm,N*sum(broaden_FL)/num_realisations)
figure
%plot(om_linear,imag(CD_av)/num_realisations)
plot(wl_nm,CD)
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
g_t_door  = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range, tol);
g_t  = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t3_range, tol); 
 t_sep_range = t_sep_range_ps*convfact;

for real_lp = 1:1% num_realisations
    %% calculate system parameters for this realisation of disorder
H_site_shifted = H_site + diag(site_shift(real_lp,:));
[H_ex_vib,fock_space_rep,~,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(H_site_shifted,om_vib,numvib,displ,mu,[]) ;

E = diag(H_exciton); E1 = E(2:N+1); E2 = E(N+2:end);

  c_nk = zeros(N); %as I have defined it c_nk = Proj(k,n)
  M_prj_sgl = M_prj(2:N+1,2:N+1);
  c_nm_f = zeros(N,N,N*(N-1)/2);
  M_prj_dbl = M_prj(N+2:end,N+2:end); fock_dbl = fock_space_rep(N+2:end,:);
      cnt=1;
      for m = 1:N
          c_nk(m,:) = M_prj_sgl(:,m);
      end
      for  f = 1:length(M_prj_dbl)
          tmp = find(fock_space_rep(N+1+f,:));
          n = tmp(1); m = tmp(2); %n always < m
              c_nm_f(n,m,:) =  M_prj_dbl(:,f);
      end 
%% calc spec params
[D_t,W_t,Wf_t,D0,Wgg,W1ee,W2ee] = spec_fun_simple(...
          t1_range,t3_range,om_r_rng,om_u_rng,tau_r,tau_u,E1,E2,...
          g_t,g_t_door,lam_tot,c_nk,c_nm_f,N,true,true) ;

   
%% use modified Redfield theory
   R_mod_red =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,...
                       M_prj(2:N+1,2:N+1)',diag(H_exciton(2:N+1,2:N+1)),t_int_rng);
     
       lg = any(any(R_mod_red>0));lg2 = logical(eye(size(R_mod_red)));
       if any(any(lg &~lg2))
           warning('positive elements in Redfield of diag removed')
    R_mod_red(R_mod_red>0) = 0;
    R_mod_red  = R_mod_red  - diag(sum(R_mod_red,1)); %balance condition     
       end
 %% prop doorway function in delay time    
    poponly = true; ode_toler = 1e-12;
    %pointless reshape to fit function I am using
    rho_neq_g = reshape(-D0,1,1,size(D0,1),size(D0,2));
    rho_neq_e = reshape(D0,1,1,size(D0,1),size(D0,2));
    
[rho_t_g,rho_t_e] = doorway_prop_basic(rho_neq_g,rho_neq_e,...
                t_sep_range,supop_g,-R_mod_red,sz2,N,ode_toler,poponly ) ;   

rho_t_g = squeeze(rho_t_g); 
%%  calculate transient spec            
   
rho_t_ef_dip_av =  zeros(size(rho_t_e,3),size(rho_t_e,4),N,1+ N*(N-1)/2);
rho_t_g_dip_av =  zeros(size(rho_t_e,3),size(rho_t_e,4),N);
 %density matrix with weights coming from window prefactors as well
 %ground / double excited manifold transitions determine last index
for  lp = 1:N %section loop

    tmpe = zeros(size(rho_t_e,3),size(rho_t_e,4),1+ N*(N-1)/2,N);
    tmpg = zeros(size(rho_t_g,3),size(rho_t_g,4));
   for lp2 = 1:N %dipole element loop
       tmpg = tmpg + rho_t_g(:,:,lp2)*alpha_lin_pol_ex(lp,lp2,1);
       for f = 1:N*(N-1)/2+1
           tmpe(:,:,f,lp2) = rho_t_e(lp,lp,:,:,lp2)*alpha_lin_pol_ex(lp,lp2,f);
       end
   end 
   tmpe = sum(tmpe,6);
    rho_t_ef_dip_av(:,:,:,:,lp,:)= tmpe;
    rho_t_g_dip_av(:,:,:,:,lp)= tmpg; 
end
end

GSB_sig = zeros(length(om_r_rng),length(t_sep_range),length(om_u_rng));
SE_sig = GSB_sig;  ESA_sig = GSB_sig;
om_r_fct = reshape(om_r_rng,length(om_r_rng),1);
om_r_fct = repmat(om_r_fct,1,length(t_sep_range),length(om_u_rng));

sz = size(rho_t_g_dip_av);
rho_t_g_rs = reshape(rho_t_g_dip_av,sz(1),sz(2),1,sz(3),sz(4),sz(5));
sz = size(V_gg_ft);
V_gg_ft_rs  = reshape(2*imag(V_gg_ft),sz(1),sz(2),sz(3),1,1,sz(4));

sz = size(rho_t_ef_dip_av);
rho_t_ef_rs = reshape(rho_t_ef_dip_av,sz(1),sz(2),1,sz(3),sz(4),sz(5),sz(6));

sz = size(V_ee_ft);
V_ee_ft_rs  = reshape(2*imag(V_ee_ft),sz(1),sz(2),sz(3),1,1,sz(4));
sz = size(V_eeff_ft);
V_eeff_ft_rs = reshape(2*imag(V_eeff_ft),sz(1),sz(2),sz(3),1,1,sz(4),sz(5));


for lp=1:N
   
       tmp = diagsum(mtimesx(V_gg_ft_rs(:,:,:,:,:,lp),...
             rho_t_g_rs(:,:,:,:,:,lp)),1,2);
    GSB_sig = GSB_sig - tmp; 

tmp2 = diagsum(mtimesx(V_ee_ft_rs(:,:,:,:,:,lp),rho_t_ef_rs(:,:,:,:,:,lp,1)),1,2);
 SE_sig  = SE_sig  - tmp2;
       for f = 2:N*(N-1)/2+1 %pathways featuring double ex states
           %average each one with the appropriate factor
         ESA_sig  = ESA_sig + diagsum(mtimesx(V_eeff_ft_rs(:,:,:,:,:,lp,f-1)...
                    ,rho_t_ef_rs(:,:,:,:,:,lp,f)),1,2);    
       end    
end
GSB_sig  = GSB_sig.*om_r_fct; SE_sig  = SE_sig.*om_r_fct;
ESA_sig  = ESA_sig.*om_r_fct;















% 
end





%% Calculate pump probe signal This part is effected by static disorder loops

for real_lp = 1: num_realisations
%% 

rho0 = zeros(size(H_site)); rho0(1)=1;  rho0 = reshape(rho0,numel(rho0),1);

H_site_shifted = H_site + diag(site_shift(real_lp,:));
%H0_renorm = H_site - diag(cellfun(@sum, lam_dru) - cellfun(@sum, lambda));

if ~inc_double_ex
[H_ex_vib,fock_space_rep,~,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(H_site_shifted,om_vib,numvib,displ,mu) ;    
else%include double excited states
[H_ex_vib,fock_space_rep,~,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(H_site_shifted,om_vib,numvib,displ,mu,[]) ;
end
%M_prj projects to the exciton basis

H_el = indiv_op{1};  H_vib = indiv_op{2};  
H_el_ex = indiv_op{4};  M_prj = indiv_op{5}; %projector to exciton basis
sz1 = length(indiv_op{1}) ; sz2 = length(indiv_op{2}); %H_vib
if sz2 == 1 %no vib levels
   H_exciton = H_el_ex; %take more accurate rep 
end
mu_ex_sym = indiv_op{6};

%have to generate this every loop as it is in the interaction basis, I
%could just project it each time but that probably isn't much faster
if~isempty(om_vib) %no explicit mode include
Lindblad_op = gen_lindblad(beta,om_vib,numvib,M_prj,gam,sz1,sz2);
end

%%
%calculate appropriate dipole averages
xx_ex_av = zeros(N,N); yx_ex_av = zeros(N,N);%exciton basis
for j1 = 1:N
    for j2 = 1:N
     for e1 = 1:N
         for e2=1:N
     xx_ex_av(j1,j2) = xx_ex_av(j1,j2) + ex_basis(j1,e1)*ex_basis(j2,e2)*xx_av(e1,e2) ;
     yx_ex_av(j1,j2) = yx_ex_av(j1,j2) + ex_basis(j1,e1)*ex_basis(j2,e2)*yx_av(e1,e2) ;
         end
     end     
    end
end
  c_nk = zeros(N); %as I have defined it c_nk = Proj(k,n)
  M_prj_sgl = M_prj(2:N+1,2:N+1);
  c_nm_f = zeros(N,N,N*(N-1)/2);
  M_prj_dbl = M_prj(N+2:end,N+2:end); fock_dbl = fock_space_rep(N+2:end,:);
      cnt=1;
      for m = 1:N
          c_nk(m,:) = M_prj_sgl(:,m);
      end
      for  f = 1:length(M_prj_dbl)
          tmp = find(fock_space_rep(N+1+f,:));
          n = tmp(1); m = tmp(2); %n always < m
              c_nm_f(n,m,:) =  M_prj_dbl(:,f);
      end 
  
mu_C_r = cross(mu,R);  
mu_ex = zeros(3,N,1+N*(N-1)/2); %transition dipoles between exciton states
mu_cross_R_ex = zeros(3,N,1+N*(N-1)/2); %mag transition dipoles between exciton states

mu_ex_test = zeros(3,N,1+N*(N-1)/2); 
for k = 1:N
   
   % mu_ex_test(:,k,1) = c_nk(k,:)*mu;
    mu_ex(:,k,1) = mu.'*c_nk(:,k);
    mu_cross_R_ex(:,k,1) = mu_C_r'*c_nk(:,k);
    
    for f= 1:N*(N-1)/2
        tmp = zeros(1,3); tmp_test  = tmp; 
        tmp2 = tmp;
        for m = 2:N
            for n = 1:m-1
         tmp_test  = tmp_test  + c_nm_f(n,m,f)*(c_nk(k,n)*mu(m,:) + c_nk(k,m)*mu(n,:));       
        tmp = tmp + c_nm_f(n,m,f)*(c_nk(n,k)*mu(m,:) + c_nk(m,k)*mu(n,:));
        tmp2 = tmp2+ c_nm_f(n,m,f)*(c_nk(n,k)*mu_C_r(m,:) + c_nk(m,k)*mu_C_r(n,:));
            end
        end
        mu_ex(:,k,1+f) = tmp; mu_ex_test(:,k,1+f) = tmp_test;
        mu_cross_R_ex (:,k,1+f) = tmp2;
    end
end

xx_av_ex = zeros(N,1+N*(N-1)/2); yx_av_ex = zeros(N,1+N*(N-1)/2);
for j1 = 1:N
    xx_av_ex(j1,1) = dot(mu_ex(:,j1,1),mu_ex(:,j1,1))/3;  
    yx_av_ex(j1,1) = 1i*typ_op*dot(mu_ex(:,j1,1),mu_cross_R_ex(:,j1,1))/3;
    for j2 = 2:N*(N-1)/2+1
        xx_av_ex(j1,j2) = dot(mu_ex(:,j1,j2),mu_ex(:,j1,j2))/3;
        yx_av_ex(j1,j2) = 1i*typ_op*dot(mu_ex(:,j1,j2),mu_cross_R_ex(:,j1,j2))/3;
    end
end


% [alpha_av,CD_av,alpha_lin_pol] = ...
%   dipole_fourth_order_av(mu,R,pol_L,pol_R,pol_linear,kk1,true); 
         alpha_ex_L = zeros(N,N,N*(N-1)/2+1); alpha_ex_R = zeros(N,N,N*(N-1)/2+1);
alpha_lin_pol_ex = alpha_ex_L;
for j1 = 1:N %pre save these
   for j2 = 1:N
       for f2 = 1:N*(N-1)/2+1
           
       
       mu_set=[mu_ex(:,j1,1),mu_ex(:,j1,1),mu_ex(:,j2,f2),mu_ex(:,j2,f2)].';
       
   alpha_lin_pol_ex(j1,j2,f2) = tensor_av(mu_set,pol_linear); 
    alpha_ex_L(j1,j2,f2) = tensor_av(mu_set,pol_L);
    alpha_ex_R(j1,j2,f2) = tensor_av(mu_set,pol_R);          
            
        end    
    
   end
end


%% no point in doing this anyway...
for j=1:N
    temp = zeros(sz1*sz2); 
    temp(1:sz2,(1+j*sz2):((j+1)*sz2)) = diag(ones(sz2,1));

V_ge{j} = temp  - diag(diag(temp));
V_eg{j} = temp' - diag(diag(temp));
end


%% Generate redfield prop op

%use_mod_red  =true;
if  ~use_mod_red
 %Calculate Markovian parameters from redfield
 %H0_renorm = H_site - diag(cellfun(@sum, lam_dru) + cellfun(@sum, lambda));
 tic
 if ~inc_double_ex
[R_red]= redfield_calc(H_el(2:N+1,2:N+1),beta,gam_dru,...
            lam_dru,{[],[]},{[],[]},{[],[]},use_markov);
 else %includes doubley excited
     %H_el is already renormalised as it were
H_single = H_el(2:N+1,2:N+1); 
H_double =  blkdiag(H_el(1,1),H_el(2+N:end,2+N:end));    %assumes no interactions 
%in general there are interactions between zero ex and double ex states
%unless one has already shifted to the ground state and rescaled dipole
%moments accordinly
%[R_red]= redfield_calc_2(H_single,H_double, fock_space_rep,beta,gam_dru,...
%            lam_dru,[],[],[],use_markov);     
[R_red]= redfield_calc_2(H_single,H_double, fock_space_rep,beta,gam_dru,...
            lam_dru,gam,lambda,om_0,use_markov);  

 end
  % check that the balance condition is satisfied and correct for errors
 tmp = zeros(N,1); tmp2 = tmp;
for k =1:length(R_red)
    tmp(k) = diagsum(abs(R_red(k,k,:,:)),3,4);
    tmp2(k) = trace(R_red(:,:,k,k));
    if abs(tmp2(k)) > 100*max(eps,eps(tmp(k)))
        warning('balance condition fails beyond what is expected')
        k
        tmp2(k)
    end
    R_red(k,k,k,k) =  R_red(k,k,k,k) - tmp2(k);
end    
%set balance condition for R_kkkk so this is matched
R_reduced = zeros(N);
for k = 1:N
    for kk = 1:N
        R_reduced (k,kk) = R_red(k+1,k+1,kk+1,kk+1);
    end
end

 %generate full op in Liouville space, sparse of courese
[R_red_op_full,ham_rs,R_red_sc]  = gen_redfield_Liouville(R_red,sz1,sz2,false);
%last parameter decides if I should bother rescaling my Hamiltonian with
%imaginary elements from the Redfield tensor
%toc
 H_exciton = H_exciton + diag(ham_rs);

%decoherence_op = Lindblad_op - R_red_op_full; 
tmp = sparse(1:length(H_exciton),1:length(H_exciton),ones(length(H_exciton),1));
L_op = -1i*(kron(tmp,sparse(H_exciton))-kron(sparse(H_exciton.'),tmp));    
if ~exist('Lindblad_op','var') %no vibrational states
supop = sparse(L_op - R_red_op_full);     
else
supop = sparse(L_op + Lindblad_op - R_red_op_full); 
end

else
    
   R_mod_red =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,...
                       M_prj(2:N+1,2:N+1)',diag(H_exciton(2:N+1,2:N+1)),t_int_rng);
     
       lg = any(any(R_mod_red>0));lg2 = logical(eye(size(R_mod_red)));
       if any(any(lg &~lg2))
           warning('positive elements in Redfield of diag removed')
       R_mod_red(R_mod_red>0) = 0;
    R_mod_red  = R_mod_red  - diag(sum(R_mod_red,1)); %balance condition     
       end
  if sz2 ~=1
      R_mod_red = kron(R_mod_red,ones(sz2));
  end
  
  
  tmp = sparse(1:length(H_exciton),1:length(H_exciton),ones(length(H_exciton),1));
  
if ~exist('Lindblad_op','var') %no explicit vibrational states in Ham
supop = -1i*(kron(tmp,sparse(H_exciton))-kron(sparse(H_exciton.'),tmp));       
else
supop = sparse(Lindblad_op ...
    -1i*(kron(tmp,sparse(H_exciton))-kron(sparse(H_exciton.'),tmp))); 
end


end


de_fn = @(t,v) supop*v;

%%  Get out specfic reduced operators
%reduced operators acting only on p_gg' elements
for kkkkkk=1:1
tmpg = zeros(sz1*sz2); tmpg(1:sz2,1:sz2)=1;
    [ag,bg,sg] = find(tmpg);
	tmpg = logical(reshape(tmpg,sz2^2*sz1^2,1));
    supop_g = supop(tmpg,tmpg);

%reduced operators acting only on 1st ex state manifold

    tmpe  = zeros(sz1*sz2); tmpe(sz2+1:sz2*(N+1),sz2+1:sz2*(N+1))=1;
    [ae,be,se] = find(tmpe);
	tmpe = logical(reshape(tmpe,sz2^2*sz1^2,1));

    supop_e = supop(tmpe,tmpe);   
    
%these matricies are always Hermitian so I only need to consider the upper
%right/lower left half of them to get the whole thing

   LI_ele_g = reshape(tril(true(1*sz2)) , 1,(1*sz2)^2);
   LI_ele_e = reshape(tril(true(N*sz2)) , 1,(N*sz2)^2);
   supop_g_1 = supop_g(LI_ele_g,LI_ele_g);  supop_e_1 = supop_e(LI_ele_e,LI_ele_e);
            %op_1 doesn't include anything that multiples the conjugates      
   supop_g_2 = supop_g(LI_ele_g,~LI_ele_g); supop_e_2 = supop_e(LI_ele_e,~LI_ele_e);
            %element that map the conjugate terms to non conjugate terms
            
 % This was written for my HEOM code and may be over general for what I am
 % doing here.
            
            count1g =1;  count2g = 0; count1e = 1; count2e= 0;
            testg = zeros(sz2); teste=zeros(sz2*N);
     for k = 1:(sz2-1)
     count2g = count2g + sz2-k;
     testg(1+k:sz2,k) =count1g:count2g;
     count1g= count2g+1; 
     end
     for k = 1:(sz2*N-1)
     count2e = count2e + sz2*N-k;
     teste(1+k:sz2*N,k) =count1e:count2e;
     count1e= count2e+1; 
     end
    test2g = reshape(testg',1,sz2^2);  test2e = reshape(teste',1,N^2*sz2^2);  
    [~,orderg] = sort(test2g(test2g~=0)); [~,ordere] = sort(test2e(test2e~=0));

            pickrg = reshape(tril(true(sz2),-1) , 1,sz2^2);
            pickre = reshape(tril(true(sz2*N),-1) , 1,N^2*sz2^2);            
            %select non linearly independent elements
            pickrg = pickrg(LI_ele_g); pickre = pickre(LI_ele_e);
            %elements which are on diagonal will be zero, not mapped from
            %but mapped to
%reorder to pick correct things
  supop_g_2 = supop_g_2(:,orderg);  supop_e_2=supop_e_2(:,ordere);              
clear  count2e count1e count2g count1g test2g testg test2e teste  
    
%reduced operators acting only on ground excited coherences, these don't
%mix to p_gg or p_ee'

    tmpge  = zeros(sz1*sz2); tmpge(1:sz2,sz2+1:sz2*(N+1))=1; %upper diag
    [aa,ba,sa] = find(tmpge);
	tmpge = logical(reshape(tmpge,sz2^2*sz1^2,1));
    %reshape(V_ge{1}(tmpge),sz2,sz2*(N))
    supop_ge = supop(tmpge,tmpge);
    
    tmpeg  = zeros(sz1*sz2);tmpeg(sz2+1:sz2*(N+1),1:sz2)=1; %lower diag
    [ab,bb,sb] = find(tmpeg);
	tmpeg = logical(reshape(tmpeg,sz2^2*sz1^2,1));
    %reshape(V_eg{1}(tmpeg),sz2*N,sz2)
    supop_eg = supop(tmpeg,tmpeg); %parts mixing between these elements only
    %note that if I add spont emmision this is not the case

    
%reduced operators acting only on 1st-2nd ex state manifold, 

    tmpef  = zeros(sz1*sz2); tmpef(sz2+1:sz2*(N+1),sz2*(N+1)+1:end)=1; %upper diag
    [aef,bef,sef] = find(tmpef);
	tmpef = logical(reshape(tmpef,sz2^2*sz1^2,1));
    % reshape(V_ef{1}(tmpef),sz2*(N),sz2*(N-1)*N/2)
    supop_ef = supop(tmpef,tmpef);
    
     tmpfe  = zeros(sz1*sz2); tmpfe(sz2*(N+1)+1:end,sz2+1:sz2*(N+1))=1; %lower diag
    [afe,bfe,sfe] = find(tmpfe);
	tmpfe = logical(reshape(tmpfe,sz2^2*sz1^2,1));
    % reshape(V_fe{1}(tmpfe),sz2*(N-1)*N/2,sz2*(N))
    supop_fe = supop(tmpfe,tmpfe);       
    
    %pretty sure I never need to consider states with population in double
    %excited state for 3rd order spec without any decay
    
    %important for first interaction propogation and window
    de_fn_ge = @(t,v) supop_ge*v; de_fn_eg = @(t,v) supop_eg*v;  
    
    %important for density matrix
    de_fn_gg = @(t,v) supop_g*v; de_fn_ee = @(t,v) supop_e*v;
    %only consider part of the density matrix
    if any(supop_g_2(:))
de_fn_gg_red = @(t,v) mtimesx(supop_g_1,v) + mtimesx(supop_g_2,v(pickrg),'G');
    else
de_fn_gg_red = @(t,v) mtimesx(supop_g_1,v);        
    end
de_fn_ee_red = @(t,v) mtimesx(supop_e_1,v) + mtimesx(supop_e_2,v(pickre),'G');
    
    %only important for third (window)
    de_fn_ef = @(t,v) supop_ef*v; de_fn_fe = @(t,v) supop_fe*v;
end%just put this loop so I can wrap it up
    
    
%% Calculate initial condition from the operator, should be sensible thermal dist
points_to_save = 40; t_end_ps = 30; t_end = t_end_ps*convfact;
%t_range = linspace(0,t_end,points_to_save);
rho_fl = zeros(length(supop),1); rho_fl(1) = 1;
%populate the vibrations thermally by allowing system to equilibrate

options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-16);

init_state = [1;zeros(length(supop_g)-1,1)];
init_state_red = init_state(LI_ele_g);

output_DE_fun(points_to_save,init_state_red,'notsavingnayway'); 
ode45(de_fn_gg_red,[0,t_end],init_state_red,options);
[tmp1,tmp22]  =  output_DE_fun(points_to_save,rho_fl,'get_data+clean');
tmp2 = 0*init_state; tmp2(~LI_ele_g) = conj(tmp22(end,pickrg));
tmp2(LI_ele_g) = tmp22(end,:); 

% 
% output_DE_fun(points_to_save,init_state,'notsavingnayway'); 
% ode45(de_fn_gg,[0,t_end],init_state,options);
% [tmp11,tmp22]  =  output_DE_fun(points_to_save,rho_fl,'get_data+clean');

rho_0 =    full(sparse(ag,bg,tmp2,sz1*sz2,sz1*sz2));
rho_fl = reshape(rho_0,numel(rho_0),1);

if abs(1-reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl) > eps(100)
    warning('eq density matrix evaluated to have trace neq 1')
    trace_discrep = 1-reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl
end
%renormalise to have trace 1, some eps level errors occur
    rho_fl = rho_fl/(reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl);
    rho_0 = rho_0/(reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl);
%% Calculate to second order the density matrix before the probe beam hits
%A density matrix associated to each path in Louiville space for
%averaging purposes
points_to_save = 2^12+1; 

if use_mod_red
    topass{1} = g_t_3; topass{2} = lam_tot ;
    topass{3} = c_nk;  topass{4} = supop ;  
[rho_neq_g,rho_neq_e] = doorway_fun_basic(rho_0,t1_range,tau_u...
          ,om_u_rng ,supop ,V_ge,V_eg,N,sz2) ;   
else
[rho_neq_g,rho_neq_e] = doorway_fun_basic(rho_0,t1_range,tau_u...
          ,pump_freq ,supop ,V_ge,V_eg,N,sz2) ;
end



%% Propogate in time     

      t_sep_range = t_sep_range_ps*convfact;
      poponly = true;
if use_mod_red
      tic
[rho_t_g,rho_t_e] = doorway_prop_basic(rho_neq_g,rho_neq_e,...
                t_sep_range,supop_g,-R_mod_red,sz2,N,1e-12,true ) ;   
      toc  
else

      tic
[rho_t_g,rho_t_e] = doorway_prop_basic(rho_neq_g,rho_neq_e,...
                t_sep_range,supop_g,supop_e,sz2,N,1e-12,poponly) ;   
      toc
end
D_actual = zeros(size(rho_t_e,1),size(rho_t_e,2),size(rho_t_e,3));


for j=1:N
D_actual(:,:,:) = D_actual(:,:,:) + xx_av_ex(j,1)*rho_t_e(:,:,:,1,j);
%D_actual(:,:,:) = D_actual(:,:,:) + xx_ex_av(j,j)*rho_t_e(:,:,:,1,j);
end
%comment out to keep in exciton basis
D_actual = mtimesx(M_prj(2:N+1,2:N+1),D_actual ); %proj to site basis
D_actual = mtimesx(D_actual,M_prj(2:N+1,2:N+1),'C');

test = D_actual / trace(D_actual(:,:,1));
tmp = zeros(N,length(t_sep_range));
for k = 1:N
    tmp(k,:) = test(k,k,:);
end
   Den_av = Den_av + tmp;
% test2 = diagsum(test,1,2);
 test=reshape(test,size(rho_t_e,1)^2,length(t_sep_range));
 if poponly
     tmp = eye(size(rho_t_e,1)); tmp = reshape(logical(tmp),size(rho_t_e,1)^2,1);
 test = test(tmp,:);    
 end
% 
% test=reshape(rho_t_e(:,:,:,1,1),size(rho_t_e,1)^2,length(t_sep_range));
% test = test/trace(rho_t_e(:,:,1,1,1));
% test2=reshape(rho_t_g(:,:,:,1,1),size(rho_t_g,1)^2,length(t_sep_range));
%  figure
%  plot(t_sep_range_ps,test)
% figure
% plot(t_sep_range_ps,test2-1,'--')
% figure
% plot(t_sep_range_ps,abs(test2))
%reduce size 

%%  Calculate the window function
% Assume the pulse is short enough that no time evolution occurs during the
% interaction.
clear topass
if use_mod_red
    if 1==0
    topass{1} = g_t_3; topass{2} = lam_tot ;
    topass{3} = c_nk;  topass{4} =  c_nm_f; topass{5} =supop ;

[V_gg_t,V_ee_t,V_eeff_t] = window_fun_basic(...
          t3_range,topass ,V_ge,V_eg,N,sz2,sz1,true); 
    else
   %   topass{1} =g_t_3;
     topass{1}  = line_broad_fn_full(beta,gam{1},om_0{1},lambda{1},gam_dru{1},lam_dru{1},t3_range, tol);  topass{6}=[];
     topass{2} = lam_tot ; 
    topass{3} = c_nk; topass{4} =  c_nm_f; topass{5} =supop ;       
[V_gg_t,V_ee_t,V_eeff_t] = window_fun_basic(...
          t3_range,topass ,V_ge,V_eg,N,sz2,sz1,true);     
    end
else
[V_gg_t,V_ee_t,V_eeff_t] = window_fun_basic(...
          t3_range,supop ,V_ge,V_eg,N,sz2,sz1,true);     
end

%% Calculate actual populations by integrating over the functions

E_comb = exp(-t3_range.^2/tau_r^2/4); %int_-inf^inf dt' E(t')E(t'+t)
E_comb = reshape(E_comb,1,1,length(t3_range));
%V_gg_om = zeros(sz2^2,length(om_r_rng),N);
%V_ee_om = zeros(sz2^2,length(om_r_rng),N);
%V_eeff_om =zeros(sz2^2,length(om_r_rng),N,N*(N-1)/2);

dt = t3_range(2)-t3_range(1);
tmax = t3_range(end);
om_rng = pi*(-1/dt:2/tmax:1/dt);


if gaussian_pulses

   om_typ = (H_exciton(sz2*(1+N),sz2*(1+N)-H_exciton(1,1)));
  % to_choose = om_rng >- 0.6e4 & om_rng < 0.8e4;  %~60000 wavenumbers either way
   om_rng = om_rng + om_typ;    
   
   om_r_fc = exp(-1i*om_typ*reshape(t3_range,1,1,length(t3_range)));
   V_gg_ft = fftshift(fft(V_gg_t.*repmat(E_comb.*om_r_fc,sz2,sz2,1,N),[],3),3)/length(t3_range);
   V_ee_ft =  fftshift(fft(V_ee_t.*repmat(E_comb.*om_r_fc,sz2,sz2,1,N),[],3),3)/length(t3_range);
   V_eeff_ft = fftshift(fft(V_eeff_t.*repmat(E_comb.*om_r_fc,sz2,sz2,1,N,N*(N-1)/2),[],3),3)/length(t3_range);   
else
    
end
%interpolate to desired frequency range
  V_gg_ft = permute(interp1(om_rng,permute(V_gg_ft,[3,2,1,4,5]),om_r_rng),[3,2,1,4,5]);
   V_ee_ft =  permute(interp1(om_rng,permute(V_ee_ft,[3,2,1,4,5]),om_r_rng),[3,2,1,4,5]);
   V_eeff_ft = permute(interp1(om_rng,permute(V_eeff_ft,[3,2,1,4,5,6]),om_r_rng),[3,2,1,4,5,6]);

%% Linear spec test
plot_lin_spec = false;
if plot_lin_spec
lin_resp =  reshape(diagsum(V_gg_ft(1:sz2,1:sz2,:,:)*rho0(1:sz2,1:sz2),1,2),length(om_rng),N);
    
alpha_sig = lin_resp*diag(xx_ex_av); 
beta_sig = lin_resp*diag(yx_ex_av);    
    
figure
plotyy(om_rng,imag(alpha_sig),om_rng,real(alpha_sig))
figure
plotyy(om_rng,real(beta_sig),om_rng,imag(beta_sig))


lin_resp = mtimesx(V_gg_t,rho_0(1:sz2,1:sz2));
lin_resp = diagsum(lin_resp,1,2); %trace

centre_wl = 800;%put centre wavelength in nm
centre_freq = 0*10^7/800; exp_shift = exp(-1i*centre_freq*t3_range.');
dt = t3_range(2)-t3_range(1);  Dt = t3_range(end)-t3_range(1);  
om_range = centre_freq + 2*pi*(-1/dt/2:1/Dt:1/dt/2);
lam_nm = 10^7./om_range;

alpha_sig = lin_resp*diag(xx_ex_av);  alpha_sig = alpha_sig.*exp_shift;
beta_sig = lin_resp*diag(yx_ex_av); beta_sig = beta_sig.*exp_shift;
alpha_sig_ft = fftshift(fft(alpha_sig)); beta_sig_ft = fftshift(fft(beta_sig)); 

Alpha = -imag(alpha_sig_ft);  ref_ind = real(alpha_sig_ft);
CD = real(beta_sig_ft);  OR = imag(beta_sig_ft);

lg = lam_nm > 450 & lam_nm < 600;

figure
plotyy(lam_nm(lg),Alpha(lg),lam_nm(lg),ref_ind(lg))
figure
plotyy(lam_nm(lg),CD(lg),lam_nm(lg),OR(lg))

%% find local maxima and minima in spectra
peaks_abs = localMaximum(Alpha(lg));
peaks_abs = peaks_abs(peaks_abs ~= 1 &  peaks_abs ~= length(lam_nm));
peaks = localMaximum(CD(lg));
troughs = localMaximum(-CD(lg));
peaks = peaks(peaks ~= 1 & peaks ~= length(lam_nm));
troughs = troughs(troughs ~= 1 & troughs ~= length(lam_nm));

max_CD = CD(lg); min_CD = max_CD(troughs);  max_CD = max_CD(peaks);
max_abs = Alpha(lg);  max_abs = max_abs(peaks_abs);
max_save{real_lp} = {max_CD,min_CD,max_abs};
alpha_CD_or_save{real_lp} = {Alpha,CD,OR};

end


%%  Calculate PP signal

if simple_state_mix
rho_t_ef_dip_av =  zeros(sz2,sz2,size(rho_t_e,3),size(rho_t_e,4),N,1+ N*(N-1)/2);
rho_t_g_dip_av =  zeros(sz2,sz2,size(rho_t_e,3),size(rho_t_e,4),N);
 %density matrix with weights coming from window prefactors as well
 %ground / double excited manifold transitions determine last index
for  lp = 1:N %section loop
    sec = (lp-1)*sz2+1:lp*sz2; 
    tmpe = zeros(sz2,sz2,size(rho_t_e,3),size(rho_t_e,4),1+ N*(N-1)/2,N);
    tmpg = zeros(sz2,sz2,size(rho_t_g,3),size(rho_t_g,4));
   for lp2 = 1:N %dipole element loop
       tmpg = tmpg + rho_t_g(:,:,:,:,lp2)*alpha_lin_pol_ex(lp,lp2,1);
       for f = 1:N*(N-1)/2+1
           tmpe(:,:,:,:,f,lp2) = rho_t_e(sec,sec,:,:,lp2)*alpha_lin_pol_ex(lp,lp2,f);
       end
   end 
   tmpe = sum(tmpe,6);
    rho_t_ef_dip_av(:,:,:,:,lp,:)= tmpe;
    rho_t_g_dip_av(:,:,:,:,lp)= tmpg; 
end
end

GSB_sig = zeros(length(om_r_rng),length(t_sep_range),length(om_u_rng));
SE_sig = zeros(length(om_r_rng),length(t_sep_range),length(om_u_rng));
ESA_sig = zeros(length(om_r_rng),length(t_sep_range),length(om_u_rng));
om_r_fct = reshape(om_r_rng,length(om_r_rng),1);
om_r_fct = repmat(om_r_fct,1,length(t_sep_range),length(om_u_rng));

sz = size(rho_t_g_dip_av);
rho_t_g_rs = reshape(rho_t_g_dip_av,sz(1),sz(2),1,sz(3),sz(4),sz(5));
sz = size(V_gg_ft);
V_gg_ft_rs  = reshape(2*imag(V_gg_ft),sz(1),sz(2),sz(3),1,1,sz(4));

%sz = size(rho_t_e_dip_av);
%rho_t_e_rs = reshape(rho_t_e_dip_av,sz(1),sz(2),1,sz(3),sz(4),sz(5));
sz = size(rho_t_ef_dip_av);
rho_t_ef_rs = reshape(rho_t_ef_dip_av,sz(1),sz(2),1,sz(3),sz(4),sz(5),sz(6));

sz = size(V_ee_ft);
V_ee_ft_rs  = reshape(2*imag(V_ee_ft),sz(1),sz(2),sz(3),1,1,sz(4));
sz = size(V_eeff_ft);
V_eeff_ft_rs = reshape(2*imag(V_eeff_ft),sz(1),sz(2),sz(3),1,1,sz(4),sz(5));


for e4=1:N
   if simple_state_mix  % e3 = e4;
      
     %  tic
%        for e1 = 1:N
%            e2 = e1;
%         dipole_fct = alpha_lin_pol_ex(e4,e1,1);
%         
%     tmp = diagsum(mtimesx(V_gg_ft_rs(:,:,:,:,:,e4),...
%     rho_t_g_rs(:,:,:,:,:,e1)),1,2);
%     GSB_sig = GSB_sig -  dipole_fct*om_r_fct .*tmp;       
       tmp = diagsum(mtimesx(V_gg_ft_rs(:,:,:,:,:,e4),...
             rho_t_g_rs(:,:,:,:,:,e4)),1,2);
    GSB_sig = GSB_sig - tmp; 

    %because elements can be mixed between the states in a time dependent
    %way we must sum over 
    
%red_rng = sz2*(e1-1)+1:sz2*(e1); %smaller range that can contribute
%        tmp2 = diagsum(mtimesx(V_ee_ft_rs(:,:,:,:,:,e4),...
%                rho_t_e_rs(red_rng,red_rng,:,:,:,e1)),1,2);

tmp2 = diagsum(mtimesx(V_ee_ft_rs(:,:,:,:,:,e4),rho_t_ef_rs(:,:,:,:,:,e4,1)),1,2);

 SE_sig  = SE_sig  - tmp2;

%         tmp = V_eeff_ft_rs(:,:,:,:,:,e4,1)*alpha_lin_pol_ex(e4,e1,2);
       for f = 2:N*(N-1)/2+1 %pathways featuring double ex states
           %average each one with the appropriate factor
         ESA_sig  = ESA_sig + diagsum(mtimesx(V_eeff_ft_rs(:,:,:,:,:,e4,f-1)...
                    ,rho_t_ef_rs(:,:,:,:,:,e4,f)),1,2);    
       end
%        for f = 3:N*(N-1)/2+1 %pathways featuring double ex states
%            %average each one with the appropriate factor
%          tmp=tmp + V_eeff_ft_rs(:,:,:,:,:,e4,f-1)*alpha_lin_pol_ex(e4,e1,f);
%        end
%             tmp2 = diagsum(mtimesx(tmp,...
%                 rho_t_e_rs(red_rng,red_rng,:,:,:,e1)),1,2);           

         %ESA_sig  = ESA_sig  +  dipole_fct*om_r_fct.*tmp2;
     %  end    
   %toc       
   else
              warning('nothing written yet')
%trace_ee = squeeze(mtimesx(permute(V_ee_ft,[2,1,3]),rho_neq_e));
%trace_ef = squeeze(mtimesx(permute(V_eeff_ft,[2,1,3,4]),rho_neq_e));  
%trace_gg = squeeze(mtimesx(permute(V_gg_ft,[2,1,3]),rho_neq_g));       
       

    end
end
GSB_sig  = GSB_sig.*om_r_fct; SE_sig  = SE_sig.*om_r_fct;
ESA_sig  = ESA_sig.*om_r_fct;
%% Collect together
if real_lp==1
  GSB_sig_total = GSB_sig ;
SE_sig_total = SE_sig ;
 ESA_sig_total =  ESA_sig;  
else
GSB_sig_total = GSB_sig_total + GSB_sig ;
SE_sig_total = SE_sig_total + SE_sig ;
 ESA_sig_total =  ESA_sig_total+ ESA_sig;
end
end

%%
totsig = (SE_sig_total + GSB_sig_total + ESA_sig_total)/num_realisations ;
[broaden_sig,~] = Gauss_broad_fn(centre_SD,om_r_rng,totsig,1);
lam_rng = 10^7./om_r_rng;
lg = lam_rng<690 & lam_rng>300;
t_sep_ps = t_sep_range/convfact;
time_select  = [0.22,0.42,0.56,1.01,2.12];
lg2 = time_select*0;
for k =1:length(time_select)
    lg2(k) = find(time_select(k)<t_sep_ps,1,'first');
end

figure
norm_fct = max(max(abs(totsig(:,:,1))));
plot(lam_rng(lg),totsig(lg,lg2,1)/norm_fct)

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


%%  Calculate the (Louiville space) inner product of the Window and Doorway


% 
% %final signal is x only in this scheme (assuming CD effects small etc)
% Spp_g = Spp_alpha*0; Spp_e = Spp_alpha*0; Spp_f = Spp_alpha*0; %sep conts
%  Spp_CDf = Spp_CD*0;Spp_CDe =Spp_CD*0; Spp_CDg = Spp_CD*0;
% 
%    
% for e1=1:N  %loop over all possible interactions
%     for e2 = 1:N
%         for e3 = 1:N
%             for e4 = 1:N
%                 c1 = zeros(1,N); c2 = c1; c3=c1;c4=c1;
%                 %work out average based on occupation of each different
%                 %site
%                     for lp = 1:N
%                         %calculate effective occupancy of sites
%                         %i.e. decomposition of ex dipole op into site 
%                         %somewhat general to be adapted for DE states
%      c1 = c1 + ex_basis(e1,lp)*double(fock_space_rep(1+lp,:));
%      c2 = c2 + ex_basis(e2,lp)*double(fock_space_rep(1+lp,:));
%      c3 = c3 + ex_basis(e3,lp)*double(fock_space_rep(1+lp,:));
%      c4 = c4 + ex_basis(e4,lp)*double(fock_space_rep(1+lp,:));
%                     end
%       fct_lin_pol = 0;   fct_alpha= 0; fct_CD= 0; fct_CD2= 0; fct_alpha2 =0;
% for j1 = 1:N %I should pre save these but this is quick anyway
%    for j2 = 1:N
%        for j3 = 1:N
%            for j4=1:N
% %[alpha_av,CD_av,CD2_av] = abs_CD_cont_3rd_order(mu,R,[j1,j2,j3,j4],pol(1:2),kk1);   
% % [alpha_L_4,alpha_R_4,alpha_L_5, alpha_R_5] = abs_CD_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol(1:2),kk1);
% % 
% % alpha_term4 = (alpha_L_4+alpha_R_4)/2; alpha_term5 = (alpha_L_5+alpha_R_5)/2;
% % CD_term4 = (alpha_L_4-+alpha_R_4)/2; CD_term5 = (alpha_L_5-alpha_R_5)/2;
% % alpha_av - alpha_term4 -alpha_term5
% % CD_av2 - CD_term4
% % CD_av - CD_term5
% 
%    cont_fc = c1(j1)*c2(j2)*c3(j3)*c4(j4);
%    fct_lin_pol = fct_lin_pol + cont_fc*alpha_lin_pol(j1,j2,j3,j4); 
%    fct_alpha = fct_alpha + cont_fc*alpha_av(j1,j2,j3,j4); 
%    fct_CD = fct_CD + cont_fc*CD_av(j1,j2,j3,j4);
%    
% %    [fct_x1_tmp,fct_y1_tmp,~] = full_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol,kk1);  
% % fct_x1_tmp-2*alpha_av
%    
% % [alpha_av,CD_av] = abs_CD_cont_3rd_order(mu,R,[j1,j2,j3,j4],pol,kk2);        
% % fct_alpha2 = fct_alpha2 + cont_fc*alpha_av;   
% %    fct_CD2 = fct_CD2 + cont_fc*CD_av;
%   %weighting factors for x and y component with each probe ordering
%   %usually kk2 lost in RWA anyway
%            end
%        end
%    end
% end      
%  %fct_x1
% %  fct_y1
% for j = 1:length(om_u_rng)
% 
% %should give a om_r_rng by t_sep_rng
% tmp_ee = rho_tau(:,:,j,e1,e2);
% tmp_gg = rho_tau_p (:,:,j,e1,e2);
% 
% trace_ee = mtimesx(window_op_ee(:,:,e3,e4),tmp_ee);
% trace_ef = mtimesx(window_op_ef(:,:,e3,e4),tmp_ee);  
% 
% trace_gg = mtimesx(window_op_gg(:,:,e3,e4),tmp_gg); 
% 
% %Spp_alpha(:,:,j) = Spp_alpha(:,:,j) + fct_alpha *(trace_ef-(trace_ee+trace_gg));   
% Spp_f(:,:,j) = Spp_f(:,:,j) + fct_alpha *(trace_ef); 
% Spp_e(:,:,j) = Spp_e(:,:,j) + fct_alpha *(trace_ee); 
% Spp_g(:,:,j) = Spp_g(:,:,j) + fct_alpha *(trace_gg); 
% 
% % trace_ee = mtimesx(window_op2_ee(:,:,e3,e4),tmp_ee);
% % trace_ef = mtimesx(window_op2_ef(:,:,e3,e4),tmp_ee);  
% % trace_gg = mtimesx(window_op2_gg(:,:,e3,e4),tmp_gg );  
% % 
% % Spp_CD(:,:,j) = Spp_CD(:,:,j) + fct_CD *(trace_ef-(trace_ee+trace_gg));
% Spp_CDf(:,:,j) = Spp_CDf(:,:,j) + fct_CD *(trace_ef); 
% Spp_CDe(:,:,j) = Spp_CDe(:,:,j) + fct_CD *(trace_ee); 
% Spp_CDg(:,:,j) = Spp_CDg(:,:,j) + fct_CD *(trace_gg); 
% 
% end                %+ fct_y2 * RWA dropped terms
%                 
%             end
%         end
%     end
% end
% t_delay_range_fs = t_sep_rng/convfact*1000;
% %include factor of 2 om_r which just comes from solution to maxwells eqns
% %Spp_alpha = Spp_alpha .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
% %Spp_CD = Spp_CD .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
% 
% Spp_f = Spp_f .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
% Spp_e = Spp_e .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
% Spp_g = Spp_g .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
% 
% Spp_CDf = Spp_CDf .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
% Spp_CDe = Spp_CDe .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
% Spp_CDg = Spp_CDg .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
% 
% Spp_alpha = Spp_alpha + real(Spp_f -(Spp_e+Spp_g)); %save operators
% Spp_CD = Spp_CD + real(Spp_CDf - (Spp_CDe+Spp_CDg));
% 
% pause(1) 
% 
% end
% 
% 
% 
% 
