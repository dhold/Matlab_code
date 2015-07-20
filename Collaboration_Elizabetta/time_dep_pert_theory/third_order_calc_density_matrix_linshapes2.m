%Used to calculate the third order response function for a dimer using
%Redfield theory and explicitly included, damped Harmonic oscillator modes
%get inverse cm units
Temp = 77; %temp in Kelvin
[convfact, beta,speed_unit]= inv_cm_unit_sys(Temp);
B = beta;

kpr = [0,0,1];   %probe unit direction vector, i.e. along z
theta = atan(sqrt(2)); %angle between pump and probe
kpu = [sin(theta),0,cos(theta)];   %pump unit direction vector

%polarization vectors for pump beam
pol_u= [1,0,0;1,0,0]; 
%probe is assumed to be either along x or circularly polarized


type_of_pulse = 'Gaussian';
%either give half width at half max or sd
Hwhm_u_fs = 75; %half width at half maximum, pump
[tau_u,E_u,E_u_w,E_u_inc] = pulsefns(Hwhm_u_fs,type_of_pulse);
Hwhm_r_fs = 75; %half width at half maximum, probe
[tau_r,E_r,E_r_w] = pulsefns(Hwhm_r_fs,type_of_pulse);
%I'm not sure this normalisation is a good idea, if I change the pulse
%width then I change the energy in the pulse... oh well

sys =4; %choose  system
[H_site,mu,R,lambda,gam,om_0,lam_dru, gam_dru,om_vib,displ,pdm,sd_shift] =sys_chooser(sys);
N = length(H_site); sz1 = N+1+N*(N-1)/2;
%set assumptions about how to solve, e.g. if modified Redfield used
use_markov = true; inc_double_ex = true; use_mod_red  =true;
site_dep_bath = false; %sites have same bath / diff baths
%set range of probe frequencies to scan
om_r_rng = linspace(10^7./(700),10^7./(400),1000);
%choose range of pulse seperations to scan
t_sep_range_ps = (0:10:4500)/1000; 
%choose range of pump frequencies to use
om_u_rng = 1.0e+04*[1.9802];%,2.1];
%% Calculate linebroadening function

points_to_save = 150000;
t1_range = linspace(0,15*tau_u,points_to_save);
t3_range = linspace(0,20*tau_r,points_to_save); 
g_cnst = zeros(N,1);
for n = 1:N
    g_cnst(n) = line_broad_fn_markov(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n});
end
Gamma_min = 1/min(real(g_cnst)); %Markovian decay rate

int_calc_point = 200000;
t_early = linspace(0,3*Gamma_min/2,ceil(int_calc_point/2));
t_later = linspace(3*Gamma_min/2,30*Gamma_min,ceil(int_calc_point/2));
t_int_rng = [t_early,t_later(2:end)]; % combine into 

if site_dep_bath
g_t_1 = zeros(N,length(t1_range));g_t_3 = zeros(N,length(t3_range));
g_broad = zeros(N,length(t_int_rng)); 
g_deriv = g_broad;  g_sec_der = g_broad;  lam_tot = zeros(N,1);
for n =1:N
g_broad(n,:) = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);    
g_deriv(n,:) = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng);
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
 tol = 1e-14;
g_broad   = line_broad_fn_full(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);   
g_deriv  = line_broad_fn_deriv(beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);
g_sec_der = gradient(g_deriv,t_int_rng); %due to divergence at zero, calculate from first deriv

lam_tot = sum(lambda{n}) + sum(lam_dru{n});
%enumerate linebroadening functions at times t1 and t3, interp when
%possible to save time
g_t_1  = [interp1(t_int_rng,g_broad,t1_range(t1_range<=max(t_int_rng))),line_broad_fn_full(...
    beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range(t1_range>max(t_int_rng)), tol)];
g_t_3  = [interp1(t_int_rng,g_broad,t3_range(t3_range<=max(t_int_rng))),line_broad_fn_full(...
    beta,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range(t1_range>max(t_int_rng)), tol)];
toc
end

if site_dep_bath
H_site = H_site+diag(lam_tot);    
else
H_site = H_site+lam_tot*eye(N); %renormalise energies due to vib coupling
end
 Den_av = zeros(N,length(t_sep_range_ps));
%% Generate some Gaussian random noise for areas with static disorder

num_realisations = 6; %choose decent number of realisations. 

%If more are needed save old data and add more, check convergence

%shifts of the site energies
site_shift = randn(num_realisations,N); %randn has variance 1, assumes indep
%site_shift = site_shift*0 %uncomment to remove static disorder
for j =  1:N
  site_shift(:,j) = site_shift(:,j)*sd_shift(j); %change variance
end
delta_w = sum(site_shift,2)/N; %mean shift to excited state manifold energy
site_shift_centred = site_shift - repmat(delta_w,1,N); %no centre broadening

%% Precompute all the possible dipole expectation values 
% of the form <mu_4 dot e_out* ... mu_1 dot e_1  exp(i*k dot r_4..)>_{isotropic}
% these are computed as taylor expansions exp(1ix) ~ 1 + ix in terms of 4th
% and 5th order tensor averages

pol_linear = [pol_u;[1,0,0];[1,0,0]];
tmp = [[1,1i,0];[1,-1i,0]]/sqrt(2);
pol_L = [pol_u;tmp]; pol_R = [pol_u;tmp'];
av_lin_pol = zeros(N,N,N,N); 
av_left4 = av_lin_pol;  av_right4 = av_lin_pol; 
av_left5 = av_lin_pol;  av_right5 = av_lin_pol; 
k_set = [-kpr;kpr;-kpu;kpu];
for k1 = 1:N
    for k2 = 1:N
       for k3 = 1:N 
           for k4 = 1:N  
              
              mu_set = [mu(k1,:);mu(k2,:);mu(k3,:);mu(k4,:)].';
              av_lin_pol(k1,k2,k3,k4) = tensor_av(mu_set,pol_linear);
              av_left4(k1,k2,k3,k4) = tensor_av(mu_set,pol_L); 
              av_right4(k1,k2,k3,k4) = tensor_av(mu_set,pol_R); 
              
              for j = 1:4
              av_left5(k1,k2,k3,k4) = av_left5(k1,k2,k3,k4) +... 
                        tensor_av([mu_set;R(k1,:)],[pol_L,k_set(k1)]);     
              av_right5(k1,k2,k3,k4) = av_right5(k1,k2,k3,k4) +... 
                        tensor_av([mu_set;R(k1,:)],[pol_R,k_set(k1)]); 
              end
           end
       end
    end
end
av_circ_pol = (av_left4 + av_right4)/2; %average absorbtion
CD_cont = (av_left4 - av_right4); %Circular dichroism

  
%%

 t_sep_range = t_sep_range_ps*convfact;
D_mat_av = zeros(N,N,length(t_sep_range_ps)); %only for first pump freq
for real_lp = 1: num_realisations
    %% calculate system parameters for this realisation of disorder

    H_site_shifted = H_site + diag(site_shift(real_lp,:));
    [H_full, fock_space_rep] = generate_ex_vib_ham(H_site_shifted);    
   [H_exciton,M_prj] = eig(H_full); %diagonalise
    
    E = diag(H_exciton);  E = E-E(1); %gs energy might not be zero
    E1 = E(2:N+1); E2 = E(N+2:end); %energy by manifold

LL = length(H_exciton);
 manifold_chk = round(sum(double(fock_space_rep),2));
 %Q_j(1+a,1+b,n) = c_nk(a,n).*c_nk(b,n)
  Q_j = zeros(LL,LL,N); %mixes only within the same excitation manifold
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
 end     %this four loop enumerates Q_j
 tmp1 = reshape(Q_j,LL,LL,N,1,1,1);
 QQ = squeeze(mtimesx(permute(tmp1,[4,3,1,2,5,6]),permute(tmp1,[3,4,5,6,1,2])));
%this is the matrix QQ(a,b,c,d) := sum_n <a | Q_n |b> <c|Q_n |d>
%note that |a>_{ex} = sum_n c_nk |k>_{site}
% QQ(1+a,1+b,1+c,1+d) = sum_n c_nk(a,n)*c_nk(b,n)*c_nk(c,n)*c_nk(d,n), a<N     
    

%% calc spec params
use_fft = true;
[D_t,W_t,Wf_t,D0,Wgg,W1ee,W2ee] = spec_fun_simple(...
          t1_range,t3_range,om_r_rng,om_u_rng,tau_r,tau_u,E1,E2,...
          g_t,g_t_door,lam_tot,c_nk,c_nm_f,N,true,use_fft) ;

   
%% use modified Redfield theory


   R_mod_red =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,B,...
                       M_prj(2:N+1,2:N+1)',diag(H_exciton(2:N+1,2:N+1)),t_int_rng,1e-2);

       lg = R_mod_red>0;lg2 = logical(eye(size(R_mod_red)));
       if any(any(lg &~lg2))
           warning('positive elements in Redfield off diag removed')
    R_mod_red(R_mod_red>0) = 0;
    R_mod_red  = R_mod_red  - diag(sum(R_mod_red,1)); %balance condition     
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

