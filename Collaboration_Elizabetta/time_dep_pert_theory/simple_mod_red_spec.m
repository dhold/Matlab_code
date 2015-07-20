%Used to calculate the third order response function for a dimer using
%Redfield theory and explicitly included, damped Harmonic oscillator modes
%get inverse cm units
Temp = 300; %temp in Kelvin
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
E_r_w = @(om) sech(pi*tau_r.*om)*sqrt(tau_u*pi/2);   
    
end
%%
sys =4; %choose  system
[H_site,mu,R,lambda,gam,om_0,lam_dru, gam_dru,om_vib,displ,pdm,sd_shift] = sys_chooser(sys);
N = length(H_site); sz1 = N+1+N*(N-1)/2;

%set assumptions about how to solve, e.g. if modified Redfield used
use_markov = true; inc_double_ex = true; use_mod_red  =true;
points_to_save = 5000;

g_cnst = zeros(N,1);
for n = 1:N
    g_cnst(n) = line_broad_fn_markov(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n});
end
Gamma_min = 1/min(real(g_cnst));

t1_range = linspace(0,min(10*tau_u,70*Gamma_min),points_to_save);
t3_range = linspace(0,min(10*tau_r,70*Gamma_min),points_to_save); 


t_int_rng = [linspace(0,Gamma_min,10000),...
            linspace(Gamma_min,10*Gamma_min,20000),...
            linspace(10*Gamma_min,50*Gamma_min,20000)];
 t_int_rng = t_int_rng([true,diff(t_int_rng)~=0]);
%calculate lin broadening functions in site basis
g_t_1 = zeros(N,length(t1_range));g_t_3 = zeros(N,length(t3_range));
g_broad = zeros(N,length(t_int_rng)); 
g_deriv = g_broad;  g_sec_der = g_broad;  lam_tot = zeros(N,1);
t = sym('t','real');

for n =1:N

g_broad(n,:)   = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);    
g_deriv(n,:)  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);
g_sec_der(n,:) = line_broad_fn_sec_der(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);
lam_tot(n) = sum(lambda{n}) + sum(lam_dru{n});
g_t_1(n,:)  = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range);
g_t_3(n,:)  = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t3_range);
end

[ex_basis,H0ex] = eig(H_site);
H0 = zeros(length(H_site)+1); H0(2:end,2:end) = H_site;

%% Generate some Gaussian random noise for areas with static disorder

num_realisations = 1 %choose decent number of realisations. 

%If more are needed save old data and add more

%shifts of the site energies
site_shift = randn(num_realisations,N); %randn has variance 1, assumes indep
for j =  1:N
  site_shift(:,j) = site_shift(:,j)*sd_shift(j); %change variance
end
delta_w = sum(site_shift,2)/N; %mean shift to excited state manifold energy

%Could also put some static disorder in the relative dipole moments
%orientations, but I have no idea how and am pretty sure this wouldn't be
%independent of the site energies/interaction energies anyway
%Also if these can fluctuate they can be driven from equilibrium positions
%by interactions with a pump beam or interactions with the bath, leading to
%time dependent site-site interactions.  Leading to complex dynamics.

om_u_rng = 1.0e+04 * [1.1910    1.2210    1.2410    1.2710];

t_sep_rng_fs = 0:2000; %0 fs to  2 ps
t_sep_rng = t_sep_rng_fs/1000*convfact; 

points_to_save = 2^10+1; 
des_freq_r = 1e4*linspace(-2,2,points_to_save);
om_r_rng =des_freq_r(10^7./des_freq_r >=650 & 10^7./des_freq_r <=950);

Spp_alpha = zeros(length(om_r_rng),length(t_sep_rng),length(om_u_rng));
Spp_CD = Spp_alpha; 

%%  Calculate dipole averages for the pump probe geometry
%Assuming no static disorder impacts the dipoles (apart from (random) 
%molecular alignment) we can precompute these values.  The actual values 
%will be averaged over exciton states though, which depend on site energies

 typ_op = 10^7/820; %typical frequency
k_u = 2*pi*abs(typ_op)*kpu; k_r = 2*pi*abs(typ_op)*kpr; %average amplitudes of these para
kk1 = [-k_u;k_u;k_r;-k_r]; kk2 = [k_u;-k_u;k_r;-k_r]; %order of interaction (left to right)
%take probe polarization wavevectors along x

pol_linear = [pol{1};pol{2};[1,0,0];[1,0,0]];
pol_L = [pol{1};pol{2};[1/sqrt(2),+1i/sqrt(2),0];[1/sqrt(2),-1i/sqrt(2),0]];
pol_R = [pol{1};pol{2};[1/sqrt(2),-1i/sqrt(2),0];[1/sqrt(2),+1i/sqrt(2),0]];

 [alpha_av,CD_av,alpha_lin_pol] = ...
   dipole_fourth_order_av(mu,R,pol_L,pol_R,pol_linear,kk1,false);
   
%% Two body averages, required for linear spec and Neq density populations

xx_av = zeros(N,N); yx_av = zeros(N,N);
for j1 = 1:N
    for j2 = 1:N
        xx_av(j1,j2) = dot(mu(j1,:),mu(j2,:))/3;
        yx_av(j1,j2) = 1i*10^7/800*dot(mu(j1,:),cross(mu(j2,:),R(j1,:)-R(j2,:)))/6;
    end
end

   
%% This part is effected by static disorder loops
for real_lp = 1: num_realisations
 %%
H_site_shifted = H_site + diag(site_shift(real_lp,:));
    %Calculate Hamiltonian included double excited and diagonalise
 [H_tot,fock_space_rep,mu_ex] = ...
    generate_ex_vib_ham(H_site_shifted,[],[],[],[],[]) ;   
   [M_prj,H_exciton] = eig(H_tot); 
   %calculate modified tensor Redfield    
 
  %%
  c_nk = M_prj(2:N+1,2:N+1); c_nm_f = zeros(N,N,N*(N-1)/2);
  M_prj_dbl = M_prj(N+2:end,N+2:end); fock_dbl = fock_space_rep(N+2:end,:);

  for f = 1:N*(N-1)/2
      cnt=1;
      for m = 1:N
          for n = 1:m-1             
              c_nm_f(n,m,f) =  M_prj_dbl(f,cnt);
              cnt = cnt+1;
          end
      end 
  end

   E = diag(H_exciton);
  E1 = E(2:N+1); E2 = E(N+2:end); inc_double=true;
  
        [D_t,W_t,Wf_t,D0,Wgg,W1ee,W2ee] = spec_fun_simple(...
          t1_range,t3_range,om_r_rng,om_u_rng,tau_r,tau_u,E1,E2,...
          g_t_3,g_t_1,lam_tot,c_nk,c_nm_f,N,inc_double);
    %% Gen mod_red tensor
    if 1==0
   R_mod_red =  mod_redfield_calc(g_broad,g_deriv,g_sec_der,lam_tot,...
                       M_prj(2:N+1,2:N+1),diag(H_exciton(2:N+1,2:N+1)),t_int_rng);
                   %R_kkkk = -Sum_{k'~=k} R_k'k'kk
              
  R_mod_red = R_mod_red - diag(sum(R_mod_red,1));    
    else %normal Redfield
H_single = H_tot(2:N+1,2:N+1); 
H_double =  blkdiag(H_tot(1,1),H_tot(N+2:end,N+2:end));
[R_red]= redfield_calc_2(H_single,H_double, fock_space_rep,beta,gam_dru,...
            lam_dru,gam,lambda,om_0,true); 
%[R_red]= redfield_calc_2(H_single,H_double, fock_space_rep,beta,gam_dru,...
%            lam_dru,[],[],[],true);         
        R_mod_red = zeros(N,N);
        for k = 1:N
            for kk = 1:N
        R_mod_red(kk,k) = R_red(kk+1,kk+1,k+1,k+1);
            end
        end
       R_mod_red = R_mod_red - diag(sum(R_mod_red,1));
    end
%   for k =1:length(R_mod_red)
%     R_mod_red(k,k) =  R_mod_red(k,k) - sum(R_mod_red(k,:));
%   end     
 %%   Prop density matrix  
t_sep_range_fs = 0:5:2000;  t_sep_range = t_sep_range_fs /1000*convfact;
    D_tev = zeros(length(t_sep_range),N,length(om_u_rng),N);
    D_actual = zeros(length(t_sep_range),N,length(om_u_rng));
      for j = 1:length(om_u_rng)
          for k = 1:N
          init_st = zeros(N,1);  init_st(k) = D0(j,k);
          init_nrm = norm(init_st); init_st = init_st/init_nrm;
          cor_fc = 1e-7;
          den_prop_fn = @(t,v) -R_mod_red*v; %+ cor_fc*(init_nrm-norm(v))*ones(size(v));
          
          options = odeset('AbsTol',1e-15);
         [a,b] = ode45( den_prop_fn,[0,t_sep_range(end)],init_st,options);
          D_tev(:,:,j,k) = init_nrm*interp1(a,b,t_sep_range); 
          D_actual(:,:,j) = D_actual(:,:,j) + xx_av(k,k)*D_tev(:,:,j,k);
          end
          D_actual(:,:,j) = D_actual(:,:,j)/norm(D_actual(1,:,j));
      end
      
      figure
plot(t_sep_range_fs, D_actual(:,:,j))
      figure
plot(t_sep_range_fs, sum(D_actual(:,:,j).^2,2))
  %% Calculate stuff
  
  
end