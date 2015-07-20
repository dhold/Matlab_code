%Used to calculate the third order response function for a dimer using
%Redfield theory and explicitly included, damped Harmonic oscillator modes
%get inverse cm units
Temp = 300; %temp in Kelvin
[convfact, beta,speed_unit]= inv_cm_unit_sys(Temp);
kpr = [0,0,1];   epr = [1,0,0];   bpr = [0,1,0];
%probe unit vectors

theta = atan(sqrt(2)); %angle between pump and probe
kpu = [sin(theta),0,cos(theta)];   epu = [1,1,0]/sqrt(2); bpu=  cross(kpu,epu);
%pump unit vectors

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

sys =3; %choose  system
%% PC645 dimer
if sys == 1
%dimer specific parameters

    fle = open('Hamiltonian_save.mat');
        H_site = fle.PC645Hamiltonian;
            E0 = fle.E0_PC645; %these energies are relative 
            %to the lowest exciton state anyway
            N=2;
       H_site = H_site(1:N,1:N);%only these 2
       N=length(H_site);
       H_site = H_site - eye(N)*E0;
      
H0 = zeros(length(H_site)+1); H0(2:end,2:end) = H_site;

mu = fle.PC645_dip_and_pos([4,8],1:3); %4 and 8 are the dimer
take_same_amp = true;
if take_same_amp
%assume allare the same amplitude
for k = 1:N
    mu(k,:) = mu(k,:)/norm(mu(k,:)); %site based mu
end
end
R = fle.PC645_dip_and_pos([4,8],4:6); %4 and 8 are the dimer
%convert units from angstrom to cm^
R = R*1e-8;
     
[ex_basis,H0ex] = eig(H0);

mu_ex = [mu(1,:)*ex_basis(2,2) + mu(2,:)*ex_basis(3,2);...
        mu(1,:)*ex_basis(3,2) + mu(2,:)*ex_basis(2,3)];

%Bath parameters 

lambdaD =100;
omegaD=100;
omegaU = 650;  %omegaU = 1034; 
gammaU = 5.30884;
lambdaU = 44;

om_0 = {omegaU,omegaU};   
lambda ={lambdaU ,lambdaU }; %reorganisation energy
gamma = {gammaU,gammaU}; 
%include these explicitly in Hamiltonian
om_vib = [om_0{:}];  numvib = [4,4]
displ = [[0;0],eye(2)*sqrt(2*lambdaU/omegaU)];

%over damped, include in Redfield only
 lam_dru = {lambdaD,lambdaD}; %reorganisation energy of drude modes
gam_dru = {omegaD,omegaD};   
sz1 =3; %total number of electronic states

elseif sys == 2
%% B820_bacteriochlorophyll_dimer
flenme ='B820_bacteriochlorophyll_dimer_params.mat'; 
%insert name of file containing the system parameters
    fle = open(flenme);
        H_site = fle.H_site; N = length(H_site);
        [ex_basis,H0ex] = eig(H_site);
        
mu = fle.mu;  pdm = fle.pdm; %set to empty if not known
take_same_amp = true;
if take_same_amp
%assume all are the same amplitude
for k = 1:N
    mu(k,:) = mu(k,:)/norm(mu(k,:)); %site based mu
end
end
%mu(2,:) = mu(1,:); %test with them the same
R = fle.R;  R = R*1e-8;   %convert units from angstrom to cm^-1
if isempty(pdm)
    pdm_n_pos = [];
else
   pdm_n_pos = [pdm;R]; 
end
%read out bath parameters
lam_dru  = fle.lam_dru; gam_dru = fle.gam_dru;
%lam_dru = {250,250}; %gam_dru = {40,40}; %stronger Drude
%underdamped
om_0 = fle.om_0; lambda = fle.lambda;  gamma = fle.gamma;
%explicitly included modes

om_vib = [om_0{:}];
numvib = [3,3];
displ = [[0;0],[sqrt(2*lambda{1}/om_0{1}) ; 0],...
            [0 ; sqrt(2*lambda{2}/om_0{2})],...
    [sqrt(2*lambda{1}/om_0{1});sqrt(2*lambda{2}/om_0{2})]];
%displacements of states, ground, excited states and finally double excited
sz1 =4; %total number of electronic states

sd_shift = [200,200]; %standard deviation of site energy fluctuations 
else  %calculate parameters for general set up
 pdm = [];
 R = [0,0,-1/2;0,0,1/2];
 %take these mu to be align at an angle theta to R
 theta1 = 45*pi/180;  phi1 = 20*pi/180;
 theta2 = 45*pi/180;  phi2 = -20*pi/180;
 ry = @(a) [cos(a),0,sin(a);0,1,0;-sin(a),0,cos(a)]; %a runs from 0 to pi
rz = @(b) [cos(b),-sin(b),0;sin(b),cos(b),0;0,0,1];
 
 mu(1,:) = (R(2,:)-R(1,:)) * ry(theta1)*rz(phi1);
 mu(2,:) = (R(2,:)-R(1,:)) * ry(theta2)*rz(phi2);
 R = R*10^(-8);
 
V=200; E0 = 10^7/800;       H_site = [E0, V; V, E0];
 lam_dru  = {35,35}; gam_dru = {50,50};
%lam_dru = {250,250}; %gam_dru = {40,40}; %stronger Drude
%underdamped
om_0 = {400,400}; lambda = {V/2,V/2};  gamma = {1/convfact,1/convfact}; %about ps damping
%explicitly included modes

om_vib = [om_0{:}];
numvib = [3,3];
displ = [[0;0],[sqrt(2*lambda{1}/om_0{1}) ; 0],...
            [0 ; sqrt(2*lambda{2}/om_0{2})],...
    [sqrt(2*lambda{1}/om_0{1});sqrt(2*lambda{2}/om_0{2})]];
%displacements of states, ground, excited states and finally double excited
sz1 =4; %total number of electronic states

sd_shift = 0*[200,200];   %no disorder
end
sz2 = prod(numvib(numvib~=0));  


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

t_sep_rng_fs = 0:3000; %0 fs to  3 ps
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
%pol{3} = [1,0,0]; % pol_left = [1;-1i;0]/sqrt(2);  pol_right = [1;1i;0]/sqrt(2);

pol{1}=  [1,0,0];  pol{2} = pol{1} ; %others along x
%pol{1} = [1/sqrt(2),1/sqrt(2),0]; pol{2} = pol{1} ; %45 degrees to probe apparently
%pol_linear = [pol{1};pol{2};[1,0,0];[1,0,0]];
pol_L = [pol{1};pol{2};[1/sqrt(2),+1i/sqrt(2),0];[1/sqrt(2),-1i/sqrt(2),0]];
pol_R = [pol{1};pol{2};[1/sqrt(2),-1i/sqrt(2),0];[1/sqrt(2),+1i/sqrt(2),0]];

        % alpha_av = zeros(N,N,N,N); CD_av = zeros(N,N,N,N);
         alpha_L_4 = zeros(N,N,N,N); alpha_R_4 = zeros(N,N,N,N);
         alpha_L_5 = zeros(N,N,N,N); alpha_R_5 = zeros(N,N,N,N);   
         alpha_lin_pol = zeros(N,N,N,N);
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
        tensor_av([mu(jj,:);R(jj(j),:)],[pol_L;1i*kk1(j,:)]);   
    alpha_R_5(j1,j2,j3,j4) = alpha_R_5(j1,j2,j3,j4) +...
        tensor_av([mu(jj,:);R(jj(j),:)],[pol_R;1i*kk1(j,:)]);       
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
   
   sz  = sz1*sz2; %total size of matrix
%% This part is effected by static disorder loops
for real_lp = 1: num_realisations
%% 


rho0 = zeros(size(H_site)); rho0(1)=1;  rho0 = reshape(rho0,numel(rho0),1);

H_site_shifted = H_site + diag(site_shift(real_lp,:));
%H0_renorm = H_site - diag(cellfun(@sum, lam_dru) - cellfun(@sum, lambda));

if sys==1
[H_ex_vib,fock_space_rep,mu_ex,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(H_site_shifted,om_vib,numvib,displ,mu) ;    
else
[H_ex_vib,fock_space_rep,mu_ex,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(H_site_shifted,om_vib,numvib,displ,mu,[]) ;
end
%M_prj projects to the exciton basis

H_el = indiv_op{1};  H_vib = indiv_op{2};  
H_el_ex = indiv_op{4};  M_prj = indiv_op{5}; %projector to exciton basis
sz1 = length(indiv_op{1}) ; sz2 = length(indiv_op{2}); %H_vib
mu_ex_sym = indiv_op{6};

%have to generate this every loop as it is in the interaction basis, I
%could just project it each time but that probably isn't much faster
Lindblad_op = gen_lindblad(beta,om_vib,numvib,M_prj,gamma,sz1,sz2);

if real_lp == 1 %first realisation only
    %% Construct interaction operators (unit)
for j=1:N
    temp = zeros(sz1*sz2); 
    temp(1:sz2,(1+j*sz2):((j+1)*sz2)) = diag(ones(sz2,1));

    %populate all the ground -> excitation of site j states with vib included
V_ge{j} = temp  - diag(diag(temp));
V_eg{j} = temp' - diag(diag(temp));

%me_ge_ex{j} =ex_basis'*mu_ge{j}*ex_basis;

end
%now the coherences to the double exciton statesc if included
if sz1 == 1+N+(N-1)*N/2
tmp = 1:size(fock_space_rep,1);
 for j=1:N  %prefactor is mu_j, so states which it can mix to must have an
     temp = H_ex_vib*0; 
     for k =1:N  
         if k~=j
             lg = fock_space_rep(:,k); %states with ex at kth
    tmp2 = tmp(lg); %to these elements
    for kk = 1:length(tmp2)
       temp((1+k*sz2):((k+1)*sz2), 1+(tmp2(kk)-1)*sz2:tmp2(kk)*sz2) = ...
           diag(ones(sz2,1));
    end
         end
     end
    %populate all the ground -> excitation of site j states with vib included
V_ef{j} = temp-diag(diag(temp)); 
V_fe{j} = temp'-diag(diag(temp)); 
%in principle there can also be doubley excited states
%mu_ef{j} = mu_ef{j}-diag(diag(mu_ef{j}));     
mu_hilb{j} = V_ef{j} + V_fe{j} + V_eg{j} + V_ge{j};
 end 
else
    for j = 1:N
V_ef{j} = zeros(size(H_ex_vib));
V_fe{j} = zeros(size(H_ex_vib));
%mu_ef{j} = mu_ef{j}-diag(diag(mu_ef{j}));     
mu_hilb{j} = V_ef{j} + V_fe{j} + V_eg{j} + V_ge{j};    
    end
end
end


%% Generate redfield prop op

use_markov = true;
 %Calculate Markovian parameters from redfield
 %H0_renorm = H_site - diag(cellfun(@sum, lam_dru) + cellfun(@sum, lambda));
 if sys==1
[R_red]= redfield_calc(H_el(2:N+1,2:N+1),beta,gam_dru,...
            lam_dru,{[],[]},{[],[]},{[],[]},use_markov);
 else %includes doubley excited
     %H_el is already renormalised as it were
H_single = H_el(2:N+1,2:N+1); H_double =  blkdiag(H_el(1,1),H_el(2+N:end,2+N:end));    
[R_red]= redfield_calc_2(H_single,H_double, fock_space_rep,beta,gam_dru,...
            lam_dru,{[],[]},{[],[]},{[],[]},use_markov);     
 end
 %generate full op in Liouville space, sparse of courese
[R_red_op_full,ham_rs,R_red_sc]  = gen_redfield_Liouville(R_red,sz1,sz2,false);
%last parameter decides if I should bother rescaling my Hamiltonian with
%imaginary elements from the Redfield tensor

 H_exciton = H_exciton + diag(ham_rs);

%decoherence_op = Lindblad_op - R_red_op_full; 
tmp = sparse(1:length(H_exciton),1:length(H_exciton),ones(length(H_exciton),1));
L_op = -1i*(kron(tmp,sparse(H_exciton))-kron(sparse(H_exciton.'),tmp));    

supop = sparse(L_op + Lindblad_op - R_red_op_full); 

de_fn = @(t,v) supop*v;

%%  Get out specfic reduced operators
%reduced operators acting only on p_gg' elements

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
des_freq_u = linspace(-2e+04,2e+04,points_to_save); 
%des_freq_u = linspace(0.7e+04,1.8e+04,points_to_save); %angular frequency
%pic about 8 different frequencies from this which are of interest for some
%reason, e.g. resonant with certain transitions

pump_freq = [H_exciton(sz2+1,sz2+1)-H_exciton(2,2),...
    H_exciton(sz2+1,sz2+1)-H_exciton(1,1),...
H_exciton(2*sz2+1,2*sz2+1)-H_exciton(2,2),...
    H_exciton(2*sz2+1,2*sz2+1)-H_exciton(1,1)];   

mid_freq = mean( des_freq_u );
om_s = -1i*mid_freq ; %~frequency of oscillation to scale out   


t_step = 2*pi/(des_freq_u(end)-des_freq_u(1)); %set by range
t_end = 2*pi/(des_freq_u(2)-des_freq_u(1)); %set by frequency spacing
t1_range = 0:t_step:t_end ;
t_end_ps = t_end / convfact; %end time in ps, useful to know I think
fft_sh = repmat(exp(om_s*t1_range),length(supop_ge),1);


%length of the ground + 1st excited manifold, range of time_seperations and
%number of transitions
rho_neq = zeros([sz2*N,sz2*N,length(t1_range),N,N]); %this ends up in the gs manifold
rho_neq_p = zeros([sz2,sz2,length(t1_range),N,N]); %end in ex-state manifold

rho_tau = zeros([sz2^2*N^2,length(t_sep_rng),length(pump_freq),N,N]);
rho_tau_p = zeros([sz2^2,length(t_sep_rng),length(pump_freq),N,N]);
%take only linearly indep elements
%rho_tau = zeros([sz2^2*N^2/2+sz2*N,length(t_sep_rng),length(pump_freq),N,N]);
%rho_tau_p = zeros([sz2^2/2+sz2,length(t_sep_rng),length(pump_freq),N,N]);
for e1 = 1:N %first interation
        
    op1 = V_eg{e1}; 
    temp1 = op1*rho_0; temp1 = temp1(sz2+1:sz2*(N+1),1:sz2);
    temp1  =reshape(temp1,(numel(temp1)),1);
    
      op2 = V_ge{e1};%commutator term
    temp2 = rho_0*op2; temp2 = temp2(sz2+1:sz2*(N+1),1:sz2);
    temp2  =reshape(temp2,(numel(temp2)),1);    
        
    Gamma_s = real(R_red_sc(1,1+e1,1,1+e1)); %~coherence decay rate
   % num_scaling = exp(-Gamma_s*t1_range); 
    %keep the om_s scaling for when I perform the Fourier transform
    
    %additional rescaling om_eg om_ge beyond om_s, om_s is simply the shift
    %of the actual frequency range from the Fourier transform, this will
    %need to be shifted
    
    om_shift = H_exciton(1+(1+e1)*sz2,1+(1+e1)*sz2)-H_exciton(1,1); 
    om_eg =  om_shift;
    om_ge = -om_shift;  %shifted wrong way, lost in RWA

    num_scaling_eg = exp((1i*om_eg-Gamma_s)*t1_range); 
    num_scaling_ge = exp((1i*om_ge-Gamma_s)*t1_range); 
    
    %scale out oscillations and decay from electronic DOF
    supop_ge_scaled = supop_ge  + (1i*om_ge+Gamma_s)*sparse(eye(length(supop_ge)));
    supop_eg_scaled = supop_eg  + (1i*om_eg+Gamma_s)*sparse(eye(length(supop_eg)));
     de_fn_ge_scaled = @(t,v) supop_ge_scaled*v; 
     de_fn_eg_scaled  = @(t,v) supop_eg_scaled*v;   
     
    output_DE_fun(length(t1_range),temp1,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);

    ode45(de_fn_eg_scaled,[0,t_end],temp1,options);


    
    [tmp1,tmp2]  =  output_DE_fun(length(t1_range),temp1,'get_data+clean');
 %   figure
%plot(tmp1,phase(tmp2(:,1)))
    %clean gets rid of empty points with no data
    tmp3 = (interp1(tmp1,tmp2,t1_range).').*repmat(num_scaling_eg,length(temp1),1); 
    %flip matrix so time is in dimension two
    
    %Fourier transform into frequency space, range will be around om_char(e1)
    %I really want to do a laplace transform but only if it is analytic
    tmp4 = fftshift(fft(tmp3.*fft_sh,[],2),2)/length(t1_range);
        %reshape this to a matrix with 3rd dimension time
    tmp4 = reshape(tmp4,sz2*N,sz2,length(t1_range));
    
    for e2 = 1:N
        
        op_2 = V_ge{e2}; op_2 = op_2(1:sz2,sz2+1:sz2*(N+1)); 
        %take only the important section
        
            rho_neq (:,:,:,e1,e2) = mtimesx(tmp4,op_2);        
            rho_neq_p (:,:,:,e1,e2) = mtimesx(op_2,tmp4); 
            %prop these in time for selected frequencies to get the density
            %matrix at any time step.  For shorter pulse seperations this
            %will not be true
tic
        for lp =1 : length(pump_freq)
            
     E_fct = reshape(abs(E_u_w(des_freq_u-pump_freq(lp))).^2,1,1,length(des_freq_u));
     %integrate the frequency domain function with respect to this     
     %E_fct2 =reshape(abs(E_u_w(des_freq_u+pump_freq(lp))).^2,1,1,length(des_freq_u));
     %this term is usually dropped in the rotating wave approximation
   
    temp2  = rho_neq(:,:,:,e1,e2)+conj(permute(rho_neq(:,:,:,e1,e2),[2,1,3]));
     %add complex conjugate from other half of the density matrix
     
    temp2  = trapz(des_freq_u,temp2.*repmat(E_fct,size(temp2,1),size(temp2,2),1 ),3);
    temp2 = reshape(temp2, numel(rho_neq (:,:,1,e1,e2)),1);  %reshape into Louiville space ket
    temp2_red = temp2(LI_ele_e,:);
    
   temp3  =rho_neq_p(:,:,:,e1,e2)+conj(permute(rho_neq_p(:,:,:,e1,e2),[2,1,3]));
   temp3  = trapz(des_freq_u,temp3.*repmat(E_fct,size(temp3,1),size(temp3,2),1 ),3);  
   
   temp3  =  reshape(temp3, numel(rho_neq_p (:,:,1,e1,e2)),1);
    temp3_red  = temp3(LI_ele_g,:);
   
   
     %could scale out electronic coherence decay
    output_DE_fun(length(t_sep_rng),temp2,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-18);
    %tic
    ode45(de_fn_ee ,[t_sep_rng(1),t_sep_rng(end)],temp2,options);
    %toc
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_rng),temp2,'get_data+clean');
    tmp2 = interp1(tmp1,tmp2,t_sep_rng).';
     rho_tau(:,:,lp,e1,e2) =  tmp2;        %no scaling here     
          
    
    output_DE_fun(length(t_sep_rng),temp3,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-18);
    %tic
    ode45(de_fn_gg ,[t_sep_rng(1),t_sep_rng(end)],temp3,options);
    %toc
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_rng),temp3,'get_data+clean');             
    tmp2 = interp1(tmp1,tmp2,t_sep_rng).';     %no scaling here     
             
     rho_tau_p(:,:,lp,e1,e2) =  tmp2;         
        end
    toc       
            
    end
end
%%

%test these matricies to check the trace behaves properly
test = reshape(rho_tau(:,:,1,1,1),sz2*N,sz2*N,length(t_sep_rng));
test = test./repmat(reshape(diagsum(test,1,2),1,1,length(diagsum(test,1,2))),size(test,1),size(test,2));
test2a = diagsum(mtimesx(test,test),1,2);
test = reshape(rho_tau_p(:,:,1,1,1),sz2,sz2,length(t_sep_rng));
test = test./repmat(reshape(diagsum(test,1,2),1,1,length(diagsum(test,1,2))),size(test,1),size(test,2));
test3a = diagsum(mtimesx(test,test),1,2);

trace_test_ee = reshape(eye(length(rho_neq (:,1,1,1,1))),1,length(rho_neq (:,1,1,1,1))^2);
trace_test_gg = reshape(eye(length(rho_neq_p (:,1,1,1,1))),1,length(rho_neq_p (:,1,1,1,1))^2);

test1 = mtimesx(trace_test_ee,rho_tau);
test2 = mtimesx(trace_test_gg,rho_tau_p);
test3 = reshape(rho_tau_p(:,:,1,1,1),sz2,sz2,length(t_sep_rng));
test3 = test3 - conj(permute(test3,[2,1,3]));
herm_t1 = max(abs(test3(:)));
test4 = reshape(rho_tau(:,:,1,1,1),sz2*N,sz2*N,length(t_sep_rng));
test4 = test4 - conj(permute(test4,[2,1,3]));
herm_t2 = max(abs(test4(:))); %this one is not Hermitian :/


tr_av1 = sum(sum(sum(sum(abs(diff(test1))))))/numel(test1) ;
tr_av1i = sum(sum(sum(sum(abs(imag(test1))))))/numel(test1) ;
tr_av2 = sum(sum(sum(sum(abs(diff(test2))))))/numel(test2) ;
tr_av2i = sum(sum(sum(sum(abs(imag(test2))))))/numel(test2) ;
if tr_av1>eps(100) || tr_av2>eps(100)
warning('possible issue with trace of density matrix, decaying trace')
end
if tr_av1i >eps(100) ||  tr_av2i >eps(100)
 warning('possible issue with trace of density matrix, imag comp not at eps level')
end
if herm_t1 >eps||  herm_t2 >eps
 warning('possible issue with Herminicity')
 
end
if any(isnan([tr_av1,tr_av2,tr_av1i,tr_av2i,herm_t1 ,herm_t2]))
    error('NaN values in density matrix')
end
should_be_negative_lt_30 = log(abs([tr_av1,tr_av2,tr_av1i,tr_av2i,herm_t1 ,herm_t2]))

clear rho_neq rho_neq_p test test2 test3 test4%these take a lot of memory
%%  Calculate the window function
% Assume the pulse is short enough that no time evolution occurs during the
% interaction.
points_to_save = 2^10+1; 

%des_freq_r = 1e4*linspace(0.9,1.5,points_to_save); %angular frequency
des_freq_r = 1e4*linspace(-2,2,points_to_save);
mid_freq_r = mean( des_freq_r );
om_sr = -1i*mid_freq_r ;
om_f = om_sr; %lack of a better idea for the scaling this transition

%choose range of probe frequencies to test for
%~700 - 900nm probe range expected, 
%om_r_rng =des_freq_r(10^7./des_freq_r >=650 & 10^7./des_freq_r <=950);

t_step = 2*pi/(des_freq_r(end)-des_freq_r(1)); %set by range
t_end = 2*pi/(des_freq_r(2)-des_freq_r(1)); %set by frequency spacing
t3_range = 0:t_step:t_end ;
t3_end_ps = t_end / convfact; %end time in ps, useful to know I think


fft_sc_e = repmat(exp(-om_sr*t3_range),length(supop_ge),1);     %scale middle freq back to centre fft
fft_sc_f = repmat(exp(-om_f*t3_range),length(supop_ef),1);     %scale middle freq back to centre fft


%the window op will have elements in either the ground or excited state
%which contribute, the double excited don't as the neq dm doesn't have
%these
window_op_ee = zeros(length(om_r_rng),length(supop_e),N,N);
window_op_ef = window_op_ee; %could include this in window_ee but makes 
%analysis easier in some ways to seperate these contributions
window_op_gg = zeros(length(om_r_rng),length(supop_g),N,N);

%these are associated with the y component, assuming it is heterodyned with
%the first order signal, which we will also calculate in this loop via
% P_alpha^(1) (omega) =  E_r(omega) <<V|G(omega) V^X |rho_eq>>
% we calculate <<V|G(omega) here so this is simple enough
window_op2_ee = zeros(length(om_r_rng),length(supop_e),N,N);
window_op2_ef = window_op_ee; %could include this in window_ee but makes 
%analysis easier in some ways to seperate these contributions
window_op2_gg = zeros(length(om_r_rng),length(supop_g),N,N);

%calculate two site averages, site basis
xx_av = zeros(N,N); yx_av = zeros(N,N);
for j1 = 1:N
    for j2 = 1:N
        xx_av(j1,j2) = dot(mu(j1,:),mu(j2,:))/3;
        yx_av(j1,j2) = 1i*10^7/800*dot(mu(j1,:),cross(mu(j2,:),R(j1,:)-R(j2,:)))/6;
    end
end
xx_ex_av = zeros(N,N); yx_ex_av = zeros(N,N);%exciton basis
for j1 = 1:N
    for j2 = 1:N

     for e1 = 1:N
         for e2=1:N
     xx_ex_av(j1,j2) = xx_ex_av(j1,j2) + ex_basis(e1,j1)*ex_basis(e2,j2)*xx_av(e1,e2) ;
     yx_ex_av(j1,j2) = yx_ex_av(j1,j2) + ex_basis(e1,j1)*ex_basis(e2,j2)*yx_av(e1,e2) ;
         end
     end
     
    end
end
lin_op = zeros(length(supop_g),length(des_freq_r),2); %x and y components dim 3

for e4 = 1:N
    
    Vge_L = reshape(V_ge{e4},numel(V_ge{e4}),1)'; %conjugated
    Vge_L = Vge_L(tmpge); %reduced to elements that it can be mixed to
    Veg_L = reshape(V_eg{e4},numel(V_eg{e4}),1)'; %conjugated
    Veg_L = Veg_L(tmpeg); %reduced to elements that it can be mixed to
    %prop in time as usual but with backwards acting super op

    om_shift = H_exciton(1+e4*sz2,1+e4*sz2)-H_exciton(1,1); 
    om_eg = om_shift; %shifted wrong way, lost in RWA
    om_ge = -om_shift;  

    Gamma_e = real(R_red_sc(1,1+e4,1,1+e4)); %~coherence decay rate

    num_scaling_eg = exp((-1i*om_eg-Gamma_e)*t3_range);  
    num_scaling_ge = exp((-1i*om_ge-Gamma_e)*t3_range); 

    
    %num_scaling = exp(-Gamma_s*t3_range); 
    
    %scale out oscillations and decay from electronic DOF
    supop_eg_scaled = supop_eg  + (+1i*(om_eg)+Gamma_e)*sparse(eye(length(supop_eg)));    
    supop_ge_scaled = supop_ge  + (+1i*(om_ge)+Gamma_e)*sparse(eye(length(supop_ge)));

 
    de_fn_bk_ge = @(t,v) mtimesx(v,'T',supop_ge_scaled).';       
    de_fn_bk_eg = @(t,v) mtimesx(v,'T',supop_eg_scaled).';  
    %ode45 only takes column so can't do left acting directly
    
    output_DE_fun(length(t3_range),Vge_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    tic
    ode45(de_fn_bk_ge ,[t3_range(1),t3_range(end)],Vge_L,options);
    out1a = toc;
    
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vge_L,'get_data+clean');             
    tmp3 = (interp1(tmp1,tmp2,t3_range).') ;    
%     figure
%     plot(t3_range,phase(tmp3(1,:)))
%      ylabel('phase1')
%     figure
%     plot(t3_range,abs(tmp3(1,:)))
%      ylabel('abs1')
     tmp3 = tmp3  .*repmat(num_scaling_ge,size(tmp2,2),1);
     
    %initial window state, time dimension two

    %fourier transform to freq space, shift zero freq to centre, note this
    %is ofset by om_s due to scaling
    tmp3 = fftshift(fft(tmp3.*fft_sc_e,[],2),2)/length(t3_range);   
%     figure
% plot(des_freq_r,tmp3)
     tmp3 = reshape(tmp3,sz2,sz2*N,length(t3_range));
    
     output_DE_fun(length(t3_range),Vge_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    tic
    ode45(de_fn_bk_eg ,[t3_range(1),t3_range(end)],Veg_L,options);
     out1b =toc;
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vge_L,'get_data+clean');             
    tmp3a = (interp1(tmp1,tmp2,t3_range).').*repmat(num_scaling_eg ,size(tmp2,2),1);      
%     figure
%     plot(t3_range,phase(tmp3a(1,:)))
%     ylabel('phase2')          
    tmp3a = fftshift(fft(tmp3a.*fft_sc_e,[],2),2)/length(t3_range);  
    tmp3a = reshape(tmp3a,sz2*N,sz2,length(t3_range));    

    if N==2 %only one double exciton
        
     Vef_L = reshape(V_ef{e4},numel(V_ef{e4}),1)'; %conjugated
     Vef_L = Vef_L(tmpef); %reduced to elements that it can be mixed to        
     Vfe_L = reshape(V_fe{e4},numel(V_fe{e4}),1)'; %conjugated
     Vfe_L = Vfe_L(tmpfe); %reduced to elements that it can be mixed to  
     
        Gamma_f = real(R_red_sc(e4+1,N+2,e4+1,N+2))/2; %~coherence decay rate
            num_scaling = exp(-Gamma_f*t3_range); 
     om_shift = H_exciton(1+(N+1)*sz2,1+(N+1)*sz2)-H_exciton(1+e4*sz2,1+e4*sz2); 
    om_ef = +om_shift;      
    om_fe = -om_shift;
      
    num_scaling_ef = exp((-1i*om_ef-Gamma_f)*t3_range); 
    num_scaling_fe = exp((-1i*om_fe-Gamma_f)*t3_range); 
 
     %I can scale these elements with a different om
        supop_ef_scaled = supop_ef  + (om_sr - 1i*om_ef+Gamma_f)*sparse(eye(length(supop_fe))); 
        supop_fe_scaled = supop_fe  + (om_sr - 1i*om_fe+Gamma_f)*sparse(eye(length(supop_fe))); 
        de_fn_bk_ef =  @(t,v) mtimesx(v,'T',supop_ef_scaled).';  
        de_fn_bk_fe =  @(t,v) mtimesx(v,'T',supop_fe_scaled).';          
        
            output_DE_fun(length(t3_range),Vef_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    tic
    ode45(de_fn_bk_ef ,[t3_range(1),t3_range(end)],Vef_L,options);
    out2 = toc;
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vef_L,'get_data+clean');             
    tmp4 = (interp1(tmp1,tmp2,t3_range,'pchip').')  ;
    

    tmp4 = tmp4.*repmat(num_scaling_fe,size(tmp2,2),1);
    
    tmp4 = fftshift(fft(tmp4.*fft_sc_f,[],2),2)/length(t3_range);
    tmp4 = reshape(tmp4,sz2*N,sz2*N*(N-1)/2,length(t3_range));    
 
            output_DE_fun(length(t3_range),Vfe_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    tic
    ode45(de_fn_bk_fe ,[t3_range(1),t3_range(end)],Vfe_L,options);
    out2a = toc;
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vfe_L,'get_data+clean');             
    tmp4a = (interp1(tmp1,tmp2,t3_range,'pchip').')  ;    

    tmp4a = tmp4a.*repmat(num_scaling_ef,size(tmp2,2),1);
    
    tmp4a = fftshift(fft(tmp4a.*fft_sc_f,[],2),2)/length(t3_range);
    tmp4a = reshape(tmp4a,sz2*N*(N-1)/2,sz2*N,length(t3_range));     
    
    else
        warning('write the more general one')
    end
    for e3 = 1:N  %apply the other operator, left acting commutator
        
        op_3  = V_ge{e3};       op_3 = op_3(1:sz2,sz2+1:sz2*(N+1)); 
        op_3a = V_eg{e3};      op_3a = op_3a(sz2+1:sz2*(N+1),1:sz2); 
        %pick upper right sections
        op_3f  = V_ef{e3};    op_3f = op_3f(sz2+1:sz2*(N+1),sz2*(N+1)+1:end);   
        op_3fa = V_fe{e3};   op_3fa = op_3fa(sz2*(N+1)+1:end,sz2+1:sz2*(N+1));           
        
        tmp_gg = mtimesx(tmp3,op_3a)-mtimesx(op_3,tmp3a);
        %tmp_gg = 1i*(mtimesx(tmp3,op_3a) - conj(permute(mtimesx(tmp3,op_3a),[2,1,3])));
        tmp_ee = mtimesx(tmp3a,op_3)-mtimesx(op_3a,tmp3);  
        %tmp_ee = 1i*( mtimesx(tmp3a,op_3) - conj(permute( mtimesx(tmp3a,op_3),[2,1,3])));
        tmp_ef = mtimesx(tmp4,op_3fa)-mtimesx(op_3f,tmp4a); 
        %tmp_ef = 1i*(mtimesx(tmp4,op_3fa) - conj(permute(mtimesx(tmp4,op_3fa),[2,1,3])));
        
       
        tmp_gg_sav{e3,e4} = reshape(tmp_gg,length(supop_g),length(des_freq_r));  
        tmp_ee_sav{e3,e4}  = reshape(tmp_ee,length(supop_e),length(des_freq_r));
        tmp_ef_sav{e3,e4}  = reshape(tmp_ef,length(supop_e),length(des_freq_r));
        
        %these are terms with just <<V|G(omega) V^(x) section, can use for
        %linear spec, we only need tmp_gg for this as system is in ground
        %state

        lin_op(:,:,1) = lin_op(:,:,1) + tmp_gg_sav{e3,e4}*xx_ex_av(e4,e3);
        lin_op(:,:,2) = lin_op(:,:,2) + tmp_gg_sav{e3,e4}*yx_ex_av(e4,e3);
        
    end
    
end
 clear fft_sc_e    fft_sc_f    fft_sh   

%% Linear spec stuff
%use this to calculate the y component of the pump probe signal
rho_gg = reshape(rho_fl,sqrt(length(rho_fl)),sqrt(length(rho_fl)));
rho_gg = rho_gg(1:sz2,1:sz2);  rho_gg = reshape(rho_gg,numel(rho_gg),1)';
S1_omega = mtimesx(rho_gg,lin_op);  %first order response function in freq space
S1_omega = permute(S1_omega,[2,3,1]); %put frequency values first

lg=(10^7./des_freq_r >=700 & 10^7./des_freq_r <=980);
%small enough that I can save all
alpha = des_freq_r'.*real(S1_omega(:,1));
CD = des_freq_r'.*imag(S1_omega(:,2));
OR = des_freq_r'.*real(S1_omega(:,2)); 
%use these parameters to calculate the output electric field

lam_nm = 10^7./des_freq_r(lg);
% 
%  figure
%  plot(lam_nm,des_freq_r(lg)'.*imag(S1_omega(lg,1))) %ref index
%figure
%plot(lam_nm,alpha(lg)) %abs
% 
% figure
% plot(lam_nm,CD(lg)) %CD
% figure
% plot(lam_nm,des_freq_r(lg)'.*real(S1_omega(lg,2))) %OR

%% find local maxima and minima in spectra
peaks_abs = localMaximum(alpha(lg));
peaks_abs = peaks_abs(peaks_abs ~= 1 &  peaks_abs ~= length(lam_nm));
peaks = localMaximum(CD(lg));
troughs = localMaximum(-CD(lg));
peaks = peaks(peaks ~= 1 & peaks ~= length(lam_nm));
troughs = troughs(troughs ~= 1 & troughs ~= length(lam_nm));

max_CD = CD(lg); min_CD = max_CD(troughs);  max_CD = max_CD(peaks);
max_abs = alpha(lg);  max_abs = max_abs(peaks_abs);
max_save{real_lp} = {max_CD,min_CD,max_abs};
alpha_CD_or_save{real_lp} = {alpha,CD,OR};
%%
for e4=1:N
    for e3 = 1:N
        for lp =1:length(om_r_rng) %loop over all probe frequencies
            
     E_fct = (E_r_w(des_freq_r-om_r_rng(lp)).*conj(E_r_w(om_r_rng(lp)-des_freq_r)));
     
     %E_fct_y is the y electric field output (divided by the sample length)
     %in the limit of a very short sample (beta L <<1 )  (alpha L <<1 )
 %    E_fct_y = 2*real(E_r_w(des_freq_r-om_r_rng(lp)).*conj((CD/2+1i*OR)'.*E_r_w(om_r_rng(lp)-des_freq_r)));    
     
       %integrate over the probe wavepacket to produce the correct input
        tmp_gg2 = trapz(des_freq_r, tmp_gg_sav{e3,e4}.*repmat(E_fct,length(supop_g),1),2);      
        tmp_ee2  = trapz(des_freq_r, tmp_ee_sav{e3,e4}.*repmat(E_fct,length(supop_e),1),2);             
        tmp_ef2  = trapz(des_freq_r, tmp_ef_sav{e3,e4}.*repmat(E_fct,length(supop_e),1),2);     
        
       window_op_gg(lp,:,e3,e4) = tmp_gg2; %multiply with door_gg
       
       window_op_ee(lp,:,e3,e4) = tmp_ee2; %multiply with door_ee
       window_op_ef(lp,:,e3,e4) = tmp_ef2; %multiply with door_ee
     
        end
    end
end
clear   tmp_gg_sav tmp_ee_sav tmp_ef_sav
%  want imag(int_-inf^inf E_r(t) P^(3)(t)), transform to real by multi i
  window_op_gg = -1i*window_op_gg; %multiply with door_gg  
  window_op_ee = -1i*window_op_ee; %multiply with door_ee
  window_op_ef = -1i*window_op_ef;
%%  Calculate the (Louiville space) inner product of the Window and Doorway
om_u_rng = pump_freq; 


%final signal is x only in this scheme (assuming CD effects small etc)
Spp_g = Spp_alpha*0; Spp_e = Spp_alpha*0; Spp_f = Spp_alpha*0; %sep conts
 Spp_CDf = Spp_CD*0;Spp_CDe =Spp_CD*0; Spp_CDg = Spp_CD*0;


   
for e1=1:N  %loop over all possible interactions
    for e2 = 1:N
        for e3 = 1:N
            for e4 = 1:N
                c1 = zeros(1,N); c2 = c1; c3=c1;c4=c1;
                %work out average based on occupation of each different
                %site
                    for lp = 1:N
                        %calculate effective occupancy of sites
                        %i.e. decomposition of ex dipole op into site 
                        %somewhat general to be adapted for DE states
     c1 = c1 + ex_basis(e1,lp)*double(fock_space_rep(1+lp,:));
     c2 = c2 + ex_basis(e2,lp)*double(fock_space_rep(1+lp,:));
     c3 = c3 + ex_basis(e3,lp)*double(fock_space_rep(1+lp,:));
     c4 = c4 + ex_basis(e4,lp)*double(fock_space_rep(1+lp,:));
                    end
      fct_lin_pol = 0;   fct_alpha= 0; fct_CD= 0; fct_CD2= 0; fct_alpha2 =0;
for j1 = 1:N %I should pre save these but this is quick anyway
   for j2 = 1:N
       for j3 = 1:N
           for j4=1:N
%[alpha_av,CD_av,CD2_av] = abs_CD_cont_3rd_order(mu,R,[j1,j2,j3,j4],pol(1:2),kk1);   
% [alpha_L_4,alpha_R_4,alpha_L_5, alpha_R_5] = abs_CD_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol(1:2),kk1);
% 
% alpha_term4 = (alpha_L_4+alpha_R_4)/2; alpha_term5 = (alpha_L_5+alpha_R_5)/2;
% CD_term4 = (alpha_L_4-+alpha_R_4)/2; CD_term5 = (alpha_L_5-alpha_R_5)/2;
% alpha_av - alpha_term4 -alpha_term5
% CD_av2 - CD_term4
% CD_av - CD_term5

   cont_fc = c1(j1)*c2(j2)*c3(j3)*c4(j4);
   fct_lin_pol = fct_lin_pol + cont_fc*alpha_lin_pol(j1,j2,j3,j4); 
   fct_alpha = fct_alpha + cont_fc*alpha_av(j1,j2,j3,j4); 
   fct_CD = fct_CD + cont_fc*CD_av(j1,j2,j3,j4);
   
%    [fct_x1_tmp,fct_y1_tmp,~] = full_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol,kk1);  
% fct_x1_tmp-2*alpha_av
   
% [alpha_av,CD_av] = abs_CD_cont_3rd_order(mu,R,[j1,j2,j3,j4],pol,kk2);        
% fct_alpha2 = fct_alpha2 + cont_fc*alpha_av;   
%    fct_CD2 = fct_CD2 + cont_fc*CD_av;
  %weighting factors for x and y component with each probe ordering
  %usually kk2 lost in RWA anyway
           end
       end
   end
end      
 %fct_x1
%  fct_y1
for j = 1:length(om_u_rng)

%should give a om_r_rng by t_sep_rng
tmp_ee = rho_tau(:,:,j,e1,e2);
tmp_gg = rho_tau_p (:,:,j,e1,e2);

trace_ee = mtimesx(window_op_ee(:,:,e3,e4),tmp_ee);
trace_ef = mtimesx(window_op_ef(:,:,e3,e4),tmp_ee);  

trace_gg = mtimesx(window_op_gg(:,:,e3,e4),tmp_gg); 

%Spp_alpha(:,:,j) = Spp_alpha(:,:,j) + fct_alpha *(trace_ef-(trace_ee+trace_gg));   
Spp_f(:,:,j) = Spp_f(:,:,j) + fct_alpha *(trace_ef); 
Spp_e(:,:,j) = Spp_e(:,:,j) + fct_alpha *(trace_ee); 
Spp_g(:,:,j) = Spp_g(:,:,j) + fct_alpha *(trace_gg); 

% trace_ee = mtimesx(window_op2_ee(:,:,e3,e4),tmp_ee);
% trace_ef = mtimesx(window_op2_ef(:,:,e3,e4),tmp_ee);  
% trace_gg = mtimesx(window_op2_gg(:,:,e3,e4),tmp_gg );  
% 
% Spp_CD(:,:,j) = Spp_CD(:,:,j) + fct_CD *(trace_ef-(trace_ee+trace_gg));
Spp_CDf(:,:,j) = Spp_CDf(:,:,j) + fct_CD *(trace_ef); 
Spp_CDe(:,:,j) = Spp_CDe(:,:,j) + fct_CD *(trace_ee); 
Spp_CDg(:,:,j) = Spp_CDg(:,:,j) + fct_CD *(trace_gg); 

end                %+ fct_y2 * RWA dropped terms
                
            end
        end
    end
end
t_delay_range_fs = t_sep_rng/convfact*1000;
%include factor of 2 om_r which just comes from solution to maxwells eqns
%Spp_alpha = Spp_alpha .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
%Spp_CD = Spp_CD .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));

Spp_f = Spp_f .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_e = Spp_e .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_g = Spp_g .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));

Spp_CDf = Spp_CDf .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_CDe = Spp_CDe .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_CDg = Spp_CDg .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));

Spp_alpha = Spp_alpha + real(Spp_f - (Spp_e+Spp_g)); %save operators
Spp_CD = Spp_CD + real(Spp_CDf - (Spp_CDe+Spp_CDg));

pause(1) 

end




