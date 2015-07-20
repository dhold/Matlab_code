%Used to calculate the third order response function for a dimer using
%Redfield theory and explicitly included, damped Harmonic oscillator modes
%get inverse cm units
Temp = 300; %temp in Kelvin
[convfact, beta,speed_unit]= inv_cm_unit_sys(Temp);

kpr = [0,0,1];   epr = [1,0,0];   bpr = [0,1,0];
%probe unit vectors

theta = atan(sqrt(2)); %angle between pump and probe
kpu = [sin(theta),0,cos(theta)];   epu = [1,0,0]; bpu=  cross(kpu,epu);
%pump unit vectors

tau_u = 150/1000*convfact; %pumpp pulse SD in inverse cm
tau_r = 150/1000*convfact; %probe pulse SD in inverse cm

%I'm not sure this normalisation is a good idea, if I change the pulse
%width then I change the energy in the pulse... oh well
E_u = @(t) exp(-t.^2/tau_u^2/2) / tau_u / sqrt(2*pi); %init e field env
E_r = @(t) exp(-t.^2/tau_r^2/2) / tau_r / sqrt(2*pi);
E_u_w = @(om) exp(-tau_u^2.*om.^2/2);
E_u_inc = @(om,t) exp(-tau_u^2.*om.^2/2)*(1+1i*erfi(tau_u^2*om-1i*t)/sqrt(2)/tau_u)/2;
E_r_w = @(om) exp(-tau_r^2.*om.^2/2);

sys =2; %choose  system
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

lambdaD =100;omegaD=100;
omegaU = 650;  %omegaU = 1034; 
gammaU = 5.30884;
lambdaU = 44;

om_0 = {omegaU,omegaU};   
lambda ={lambdaU ,lambdaU }; %reorganisation energy
gamma = {gammaU,gammaU}; 
%include these explicitly in Hamiltonian
om_vib = [om_0{:}];  numvib = [2,2]
displ = [[0;0],eye(2)*sqrt(2*lambdaU/omegaU)];

%over damped, include in Redfield only
 lam_dru = {lambdaD,lambdaD}; %reorganisation energy of drude modes
gam_dru = {omegaD,omegaD};   


else
%% B820_bacteriochlorophyll_dimer
flenme ='B820_bacteriochlorophyll_dimer_params.mat'; 
%insert name of file containing the system parameters
    fle = open(flenme);
        H_site = fle.H_site; N = length(H_site);
        [ex_basis,H0ex] = eig(H_site);
delta_w = 0*[-200,200]; %static disorder can be added in this way, shifts diag
      H_site = H_site + diag(delta_w);
        
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
om_vib = [om_0{:}];  numvib = [2,2];
displ = [[0;0],[sqrt(2*lambda{1}/om_0{1}) ; 0],...
            [0 ; sqrt(2*lambda{2}/om_0{2})]...
    [sqrt(2*lambda{1}/om_0{1});sqrt(2*lambda{2}/om_0{2})]];
%displacements of states, ground, excited states and finally double excited

end
%% Add in vibrations
rho0 = zeros(size(H0ex)); rho0(1)=1;  rho0 = reshape(rho0,numel(rho0),1);
%this includes NO decay or decoherence with the ground state, not sure if
%this is really sensible at all to assume that the decoherence times are
%longer than the ground state dephasing times

H0_renorm = H_site - diag(cellfun(@sum, lam_dru) + cellfun(@sum, lambda));
if sys==1
[H_ex_vib,fock_space_rep,mu_ex,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(H0_renorm,om_vib,numvib,displ,mu) ;    
else
[H_ex_vib,fock_space_rep,mu_ex,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(H0_renorm,om_vib,numvib,displ,mu,[]) ;
end
%M_prj projects to the exciton basis
H_ex_vib_been_rs = false; %this makes sure I don't rescale this again

H_el = indiv_op{1};  H_vib = indiv_op{2};  
H_el_ex = indiv_op{4};  M_prj = indiv_op{5}; %projector to exciton basis
sz1 = length(indiv_op{1}) ; sz2 = length(indiv_op{2}); %H_vib
mu_ex_sym = indiv_op{6};
clear indiv_op

%% Construct interaction operators (unit)
for j=1:N
    temp = H_ex_vib*0; 
    temp(1:sz2,(1+j*sz2):((j+1)*sz2)) = diag(ones(sz2,1));

    %populate all the ground -> excitation of site j states with vib included
V_ge{j} = temp  - diag(diag(temp));
V_eg{j} = temp' - diag(diag(temp));

%me_ge_ex{j} =ex_basis'*mu_ge{j}*ex_basis;

end
%now the coherences to the double exciton statesc
if sys~=1
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

 %%   Calculate standard decoherence term on vibration    

%L(a) V = a V a^dag - 1/2* (a^dag a V + V a^dag a)
%Lind = gamma{k}*((1+nav)*L(a) + nav*L(a^dag))
% reshape(A * C, NN_1*NN_2,1) = kron(eye(NN_2),A)*reshape(C, NN_1*NN_2,1)
% reshape(C * B, NN_1*NN_2,1) = kron(A.',eye(NN_1))*reshape(C,NN_1*NN_2,1)
% reshape(A * C * B, NN_1*NN_2,1) = kron(B.',A)*reshape(C, NN_1*NN_2,1)
L_so = @(aa) (kron((aa').',aa) - 1/2 *  kron(eye(size(aa)),(aa')*aa)...
               - 1/2 *  kron(((aa')*aa).',eye(size(aa))))  ;
Lindblad_op =  zeros(length(H_ex_vib)^2);

for k = 1:length(numvib)
nav = exp(-beta*om_vib(k))/(1-exp(-beta*om_vib(k)));

%express the annihilation op in the full hilbert space
aa = kron(sparse(eye(sz1)),sparse(eye(prod(numvib(1:k-1)))));
aa = kron(aa, diag(sqrt(1:numvib(k)-1),1)); 
aa = kron(aa,sparse(eye(prod(numvib(k+1:end))))); 
%project into exciton basis (not ex vib basis!) as Redfield op in ex basis
aa = M_prj'*aa*M_prj; adag = aa';

Lindblad_op_indiv{k} =  gamma{k}*((nav+1)*L_so(aa) + nav*L_so(adag));  
Lindblad_op = Lindblad_op  + Lindblad_op_indiv{k};
end    

% %get init state from this
% rho_eq = zeros(length(Lindblad_op),1); rho_eq(1) = 1;
% %populate the vibrations thermally by allowing system to equilibrate
% de_fn = @(t,v)  (Lindblad_op + Ham_op)*v;
% 
% options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-60);
% output_DE_fun(10,rho_eq ,'notsavingnayway'); 
% ode45(de_fn,[0,10],rho_eq ,options);
% [tmp1,tmp2]  =  output_DE_fun(points_to_save,rho_eq ,'get_data');
% rho_eq  = (tmp2(end,:)).';
% reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_eq

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
 %take out imaginary elements of R_red(a,b,a,b) which are just rescales of
 %the transition frequencies omega_{ab}, due to redundencies need only
 %rescale omega_{0 e} and omega_{0 f} etc etc as otheres can be expressed
 %in terms of these
 R_red_sc  = R_red; %scaled redfield
 for a=1:length(R_red_sc)
     for b=1:length(R_red_sc)  
         if a~=b
        R_red_sc = R_red_sc - imag(R_red_sc(a,b,a,b));
         end
     end
 end
         ham_rs = zeros(1,sz2);
         for b= 2:size(R_red,1)
            ham_rs = [ham_rs,ones(1,sz2)*imag(R_red(1,b,1,b))]; 
         end
        H_ex_vib = H_ex_vib + diag(ham_rs);
        %if I run this cell multiple times it will keep rescaling this, this
        %tests whether I have regenerated it before doing this
        if H_ex_vib_been_rs
            warning('H_ex_vib has been rescaled twice or more')
        else
            H_ex_vib_been_rs = true;
        end
 
        R_red2 = conj(permute(reshape(R_red_sc,sz1,sz1,sz1^2),[2,3,1]));
        %reshape(mtimesx(R_redd2,rho_louiville),sz1,sz1) is application
%Expand Redfield operator to include vibrational manifold

% R_red(4,:,:,:) = 0; R_red(:,4,:,:) = 0; 
% R_red(:,:,4,:) = 0; R_red(:,:,:,4) = 0; 

R_red_op_full = zeros(length(H_ex_vib)^2); 
%pad this out to the full thing
temp = zeros(length(H_ex_vib)); 

cnt1vib = 0;  cnt2vib = 1; cnt1 = 1; cnt2 = 1;
%tic
for lp = 1:sz1^2*sz2^2
    
    %iterate along one vibrational lvl
    cnt1vib = cnt1vib+1;
    if cnt1vib > sz2 
        cnt1vib=1; cnt1 = cnt1+1;
        if cnt1 > sz1
            cnt1 = 1;  cnt2vib = cnt2vib+1;
            if cnt2vib > sz2
                cnt2vib = 1; cnt2 = cnt2+1;
            end
        end
    end
    temp = temp*0;
   for a = 1:sz1 
        for b = 1:sz1 
            temp((a-1)*sz2+cnt1vib,(b-1)*sz2+cnt2vib) = R_red_sc(cnt1,cnt2,a,b);
        end
   end

    R_red_op_full(lp,:) = reshape(temp,1,sz1^2*sz2^2 ).';

end
%toc
%N2*(a-1)+A+ N*N2*(N2*(b-1)+B) 

decoherence_op = Lindblad_op - R_red_op_full; 

L_op = -1i*(kron(eye(length(H_ex_vib)),H_exciton)-...
            kron(H_exciton.',eye(length(H_ex_vib))));    

supop = sparse(L_op + decoherence_op);

de_fn = @(t,v) supop*v;

%%  Get out specfic reduced operators
%reduced operators acting only on p_gg' elements

tmpg = zeros(sz1*sz2); tmpg(1:sz2,1:sz2)=1;
    [ag,bg,sg] = find(tmpg);
	tmpg = logical(reshape(tmpg,sz2^2*sz1^2,1));
    supop_g = supop(tmpg,tmpg);

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

%reduced operators acting only on 1st ex state manifold, 

    tmpe  = zeros(sz1*sz2); tmpe(sz2+1:sz2*(N+1),sz2+1:sz2*(N+1))=1;
    [ae,be,se] = find(tmpe);
	tmpe = logical(reshape(tmpe,sz2^2*sz1^2,1));

    supop_e = supop(tmpe,tmpe);
    
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
    
    %important for second
    de_fn_gg = @(t,v) supop_g*v; de_fn_ee = @(t,v) supop_e*v;

    %only important for third (window)
    de_fn_ef = @(t,v) supop_ef*v; de_fn_fe = @(t,v) supop_fe*v;

    
    
%% Calculate initial condition from the operator, should be sensible thermal dist
points_to_save = 40; t_end_ps = 30; t_end = t_end_ps*convfact;
%t_range = linspace(0,t_end,points_to_save);
rho_fl = zeros(length(supop),1); rho_fl(1) = 1;
%populate the vibrations thermally by allowing system to equilibrate
% 

% output_DE_fun(points_to_save,rho_fl,'notsavingnayway'); 
% ode45(de_fn,[0,t_end],rho_fl,options);
% [tmp1,tmp2]  =  output_DE_fun(points_to_save,rho_fl,'get_data');
% rho_fl = (tmp2(end,:)).';
% rho_0 = reshape(rho_fl,sqrt(length(rho_fl)),sqrt(length(rho_fl)));
% 

%use only reduced, should agree
options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-16);
output_DE_fun(points_to_save,[1;zeros(length(supop_g)-1,1)],'notsavingnayway'); 
ode45(de_fn_gg,[0,t_end],[1;zeros(length(supop_g)-1,1)],options);
[tmp1,tmp2]  =  output_DE_fun(points_to_save,rho_fl,'get_data');

rho_0 =    full(sparse(ag,bg,tmp2(end,:),sz1*sz2,sz1*sz2));
rho_fl = reshape(rho_0,numel(rho_0),1);
%should be the same as rho_fl
%can check convergence if I want, but this should be a sufficient amount of
%time and close enough for spectroscopy!
if abs(1-reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl) > eps(100)
    warning('eq density matrix evaluated to have trace neq 1')
    trace_discrep = 1-reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl
    rho_fl = rho_fl/(reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl);
    rho_0 = rho_0/(reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl);
end
%% Calculate to second order the density matrix before the probe beam hits
%A density matrix associated to each path in Louiville space for
%averaging purposes
points_to_save = 2^10+1; 

des_freq_u = linspace(0.9e+04,1.4e+04,points_to_save); %angular frequency
%pic about 8 different frequencies from this which are of interest for some
%reason, e.g. resonant with certain transitions

pump_freq = [H_ex_vib(sz2+1,sz2+1)-H_ex_vib(2,2),...
    H_ex_vib(sz2+2,sz2+2)-H_ex_vib(1,1),H_ex_vib(sz2+1,sz2+1)-H_ex_vib(1,1),...
H_ex_vib(2*sz2+1,2*sz2+1)-H_ex_vib(2,2),H_ex_vib(2*sz2+2,2*sz2+2)-H_ex_vib(1,1),...
    H_ex_vib(2*sz2+1,2*sz2+1)-H_ex_vib(1,1)];   

 freq_pckr = zeros(size(pump_freq));
 for lp = 1:length(pump_freq)
     [~,freq_pckr(lp)] = min(abs(pump_freq(lp)  - des_freq_u));
 end
pump_freq = des_freq_u (freq_pckr);

mid_freq = mean( des_freq_u );
om_s = -1i*mid_freq ; %~frequency of oscillation to scale out   

t_step = 2*pi/(des_freq_u(end)-des_freq_u(1)); %set by range
t_end = 2*pi/(des_freq_u(2)-des_freq_u(1)); %set by frequency spacing
t1_range = 0:t_step:t_end ;
t_end_ps = t_end / convfact; %end time in ps, useful to know I think


t_sep_rng_fs = 0:3000; %0 fs to  3 ps
t_sep_rng = t_sep_rng_fs/1000*convfact; 

%length of the ground + 1st excited manifold, range of time_seperations and
%number of transitions
rho_neq = zeros([sz2*N,sz2*N,length(t1_range),N,N]); %this ends up in the gs manifold
rho_neq_p = zeros([sz2,sz2,length(t1_range),N,N]); %end in ex-state manifold

rho_tau = zeros([sz2^2*N^2,length(t_sep_rng),length(pump_freq),N,N]);
rho_tau_p = zeros([sz2^2,length(t_sep_rng),length(pump_freq),N,N]);

for e1 = 1:N %first interation
        
    op1 = V_eg{e1}; 
    temp1 = op1*rho_0; temp1 = temp1(sz2+1:sz2*(N+1),1:sz2);
    temp1  =reshape(temp1,(numel(temp1)),1);
    
%         op2 = V_ge{e1};%commutator term
%     temp2 = rho_0*op2;   temp2 = temp2(sz2+1:sz2*(N+1),1:sz2);
%     temp2  =reshape(temp2,(numel(temp2)),1);    
        
    Gamma_s = real(R_red_sc(1,1+e1,1,1+e1)); %~coherence decay rate
    num_scaling = exp(-Gamma_s*t1_range); 
    %keep the om_s scaling for when I perform the Fourier transform
%    num_scaling_eg = exp((om_s-Gamma_s)*t1_range); 
%    num_scaling_ge = exp((-om_s-Gamma_s)*t1_range); 
    
    %scale out oscillations and decay from electronic DOF
%    supop_ge_scaled = supop_ge  + (om_s+Gamma_s)*sparse(eye(length(supop_ge)));
    supop_eg_scaled = supop_eg  + (-om_s+Gamma_s)*sparse(eye(length(supop_eg)));
%     de_fn_ge_scaled = @(t,v) supop_ge_scaled*v; 
     de_fn_eg_scaled  = @(t,v) supop_eg_scaled*v;   
     
    output_DE_fun(length(t1_range),temp1,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
   % tic
    ode45(de_fn_eg_scaled,[0,t_end],temp1,options);
  % toc

    
    [tmp1,tmp2]  =  output_DE_fun(length(t1_range),temp1,'get_data+clean');
 %   figure
%plot(tmp1,phase(tmp2(:,1)))
    %clean gets rid of empty points with no data
    tmp3 = (interp1(tmp1,tmp2,t1_range).').*repmat(num_scaling,length(temp1),1); 
    %flip matrix so time is in dimension two
    
    %Fourier transform into frequency space, range will be around om_char(e1)
    %I really want to do a laplace transform but only if it is analytic
    tmp4 = fftshift(fft(tmp3,[],2),2)/length(t1_range);
        %reshape this to a matrix with 3rd dimension time
    tmp4 = reshape(tmp4,sz2*N,sz2,length(t1_range));
    
%     output_DE_fun(length(t1_range),temp2,'notsavingnayway'); 
%     options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
%     tic
%     ode45(de_fn_eg_scaled,[0,t_end],temp2,options);
%     toc
%     
%     [tmp1,tmp2]  =  output_DE_fun(length(t1_range),temp2,'get_data+clean');
%     %clean gets rid of empty points with no data
%     tmp4 = (interp1(tmp1,tmp2,t1_range)).'; %tmp2 output in the wrong shape for this    
    
    
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
            %reshape into Louiville space ket
    temp2  =reshape(temp2,numel(rho_neq (:,:,1,e1,e2)),1);
   
   temp3  =rho_neq_p(:,:,:,e1,e2)+conj(permute(rho_neq_p(:,:,:,e1,e2),[2,1,3]));
   temp3  = trapz(des_freq_u,temp3.*repmat(E_fct,size(temp3,1),size(temp3,2),1 ),3);  
   temp3  =reshape(temp3, numel(rho_neq_p (:,:,1,e1,e2)),1);
 
     %could scale out electronic coherence decay
    output_DE_fun(length(t_sep_rng),temp2,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-30);
    %tic
    ode45(de_fn_ee ,[t_sep_rng(1),t_sep_rng(end)],temp2,options);
    %toc
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_rng),temp2,'get_data+clean');
    tmp2 = interp1(tmp1,tmp2,t_sep_rng).';
     rho_tau(:,:,lp,e1,e2) =  tmp2;        %no scaling here     
             
    output_DE_fun(length(t_sep_rng),temp3,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-30);
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

%test these matricies to check the trace behaves properly
trace_test_ee = reshape(eye(length(rho_neq (:,1,1,1,1))),1,length(rho_neq (:,1,1,1,1))^2);
trace_test_gg = reshape(eye(length(rho_neq_p (:,1,1,1,1))),1,length(rho_neq_p (:,1,1,1,1))^2);

test1 = mtimesx(trace_test_ee,rho_tau);
test2 = mtimesx(trace_test_gg,rho_tau_p);
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
[tr_av1,tr_av2,tr_av1i,tr_av2i]

%%  Calculate the window function
% Assume the pulse is short enough that no time evolution occurs during the
% interaction.
points_to_save = 2^10+1; 

des_freq_r = 1e4*linspace(0.8,1.6,points_to_save); %angular frequency
%des_freq_r = 1e4*linspace(-1.6,1.6,points_to_save);
mid_freq_r = mean( des_freq_r );
om_sr = -1i*mid_freq_r ;
om_f = om_sr; %lack of a better idea for the scaling

%~700 - 900nm probe range expected, 
om_r_rng =des_freq_r(10^7./des_freq_r >=700 & 10^7./des_freq_r <=980);

t_step = 2*pi/(des_freq_r(end)-des_freq_r(1)); %set by range
t_end = 2*pi/(des_freq_r(2)-des_freq_r(1)); %set by frequency spacing
t3_range = 0:t_step:t_end ;
t3_end_ps = t_end / convfact; %end time in ps, useful to know I think

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
        yx_av(j1,j2) = 1i*mid_freq_r*dot(mu(j1,:),cross(mu(j2,:),R(j1,:)-R(j2,:)))/6;
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
%         for f4 =  1:N*(N-1)/2
%             
%     Vef_L = reshape(V_ef{e4,f4},numel(V_ef{e4,f4}),1)'; %conjugated
%     Vef_L = Vef_L(tmpef); %reduced to elements that it can be mixed to
%             
%             Gamma_f = real(R_red_sc(e4+1,N+1+f4,e4+1,N+1+f4)); %~coherence decay rate
%             num_scaling = exp(-Gamma_f*t3_range); 
%             %in principle I can scale these elements with a different om
%         supop_fe_scaled = supop_fe  + (om_f+Gamma_f)*sparse(eye(length(supop_eg)));           
%         de_fn_bk_ef = @(v,t) v*supop_fe_scaled;  
%         
%     output_DE_fun(length(t3_range),Vef_L,'notsavingnayway'); 
%     options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-30);
%     tic
%     ode45(de_fn_bk_ef ,[t3_range(1),t3_range(end)],Vef_L,options);
%     toc
%     
%     [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vef_L,'get_data+clean');             
%     tmp2 = interp1(tmp1,tmp2,t3_range);     
%     
%     %still not sure If I should take every transition seperately
%     Vef_L(:,:,e4,f4) = tmp2;
%         
%         end
        
   
        Gamma_e = real(R_red_sc(1,1+e4,1,1+e4)); %~coherence decay rate

    num_scaling = exp(-Gamma_s*t3_range); 
    
    %scale out oscillations and decay from electronic DOF
    supop_ge_scaled = supop_ge  + (om_sr+Gamma_s)*sparse(eye(length(supop_ge)));
    supop_eg_scaled = supop_eg  + (-om_sr+Gamma_s)*sparse(eye(length(supop_eg)));
 
    de_fn_bk_ge = @(t,v) mtimesx(v,'T',supop_ge_scaled).';       
    de_fn_bk_eg = @(t,v) mtimesx(v,'T',supop_eg_scaled).';  
    %ode45 only takes column so can't do left acting directly
    
    output_DE_fun(length(t3_range),Vge_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-30);
    tic
    ode45(de_fn_bk_ge ,[t3_range(1),t3_range(end)],Vge_L,options);
    out1a = toc
    
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vge_L,'get_data+clean');             
    tmp3 = (interp1(tmp1,tmp2,t3_range).').*repmat(num_scaling,size(tmp2,2),1);     
    %figure
    %plot(t3_range,phase(tmp3(1,:)))
    %initial window state, time dimension two
    
    %fourier transform to freq space, shift zero freq to centre, note this
    %is ofset by om_s due to scaling
    tmp3 = fftshift(fft(tmp3,[],2),2)/length(t3_range);
    
    tmp3 = reshape(tmp3,sz2,sz2*N,length(t3_range));
    
     output_DE_fun(length(t3_range),Vge_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-30);
    tic
    ode45(de_fn_bk_eg ,[t3_range(1),t3_range(end)],Veg_L,options);
     out2 =toc
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vge_L,'get_data+clean');             
    tmp3a = (interp1(tmp1,tmp2,t3_range).').*repmat(num_scaling.*exp(2*om_sr*t3_range) ,size(tmp2,2),1);      
    tmp3a = reshape(tmp3a,sz2*N,sz2,length(t3_range));    
    
    if N==2 %only one double exciton
        
     Vef_L = reshape(V_ef{e4},numel(V_ef{e4}),1)'; %conjugated
     Vef_L = Vef_L(tmpef); %reduced to elements that it can be mixed to        
        
        Gamma_f = real(R_red_sc(e4+1,N+2,e4+1,N+2)); %~coherence decay rate
            num_scaling = exp(-Gamma_f*t3_range); 
            %in principle I can scale these elements with a different om
        supop_ef_scaled = supop_ef  + (om_f+Gamma_f)*sparse(eye(length(supop_fe))); 
        de_fn_bk_ef =  @(t,v) mtimesx(v,'T',supop_ef_scaled).';  
        
            output_DE_fun(length(t3_range),Vef_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-30);
    tic
    ode45(de_fn_bk_ef ,[t3_range(1),t3_range(end)],Vef_L,options);
    toc
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vef_L,'get_data+clean');             
    tmp4 = (interp1(tmp1,tmp2,t3_range,'pchip').').*repmat(num_scaling,size(tmp2,2),1);      
    %initial window state, fourier transform to freq space, time dim 2
    % figure
   % plot(t3_range,phase(tmp4(1,:)))   
    
    tmp4 = fftshift(fft(tmp4,[],2),2)/length(t3_range);
    tmp4 = reshape(tmp4,sz2*N,sz2*N*(N-1)/2,length(t3_range));    
        
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
        tmp_ee = mtimesx(tmp3a,op_3);%-mtimesx(op_3a,tmp3);                         
        tmp_ef = mtimesx(tmp4,op_3fa);%-mtimesx(op_3f,tmp4a); 
        
       % tmp_gg = reshape(mtimesx(tmp3,op_3,'C')-mtimesx(op_3,tmp3,'C')...
        %                ,length(supop_g),length(des_freq_r));                         
           %complex conjugates are from lower left sections of the
                  %operator, which are conjugates by Hermitivity.
                    
       % tmp_ee = mtimesx(tmp3,'C',op_3)-mtimesx(op_3,'C',tmp3);
       % tmp_ef = mtimesx(tmp4,op_3f,'C')-mtimesx(op_3f,tmp4,'C');
             %must add the V_fe(t) V_ef - V_fe V_ef(t) cont as well
             % this corresponds to e->f->e transitions
        tmp_gg = reshape(tmp_gg,length(supop_g),length(des_freq_r));  
        tmp_ee = reshape(tmp_ee,length(supop_e),length(des_freq_r));
        tmp_ef = reshape(tmp_ef,length(supop_e),length(des_freq_r));
        
        %these are terms with just <<V|G(omega) V^(x) section, can use for
        %linear spec, we only need tmp_gg for this as system is in ground
        %state

        lin_op(:,:,1) = lin_op(:,:,1) + tmp_gg*xx_ex_av(e4,e3);
        lin_op(:,:,2) = lin_op(:,:,2) + tmp_gg*yx_ex_av(e4,e3);
        
        for lp =1:length(om_r_rng) %loop over all probe frequencies
            
     E_fct = E_r_w(des_freq_r-om_r_rng(lp)).*conj(E_r_w(om_r_rng(lp)-des_freq_r));
             
       %integrate over the probe wavepacket to produce the correct input
        tmp_gg2 = trapz(des_freq_r, tmp_gg.*repmat(E_fct,length(supop_g),1),2);      
        tmp_ee2  = trapz(des_freq_r, tmp_ee.*repmat(E_fct,length(supop_e),1),2);             
        tmp_ef2  = trapz(des_freq_r, tmp_ef.*repmat(E_fct,length(supop_e),1),2);     
        
       window_op_gg(lp,:,e3,e4) = tmp_gg2; %multiply with door_gg
       
       window_op_ee(lp,:,e3,e4) = tmp_ee2; %multiply with door_ee
       window_op_ef(lp,:,e3,e4) = tmp_ef2; %multiply with door_ee
        
        end
    end
    
end
%% Linear spec stuff
rho_gg = reshape(rho_fl,sqrt(length(rho_fl)),sqrt(length(rho_fl)));
rho_gg = rho_gg(1:sz2,1:sz2);  rho_gg = reshape(rho_gg,numel(rho_gg),1)';
S1_omega = mtimesx(rho_gg,lin_op);  %first order response function in freq space
S1_omega = permute(S1_omega,[2,3,1]); %put frequency values first

figure
plot(des_freq_r,imag(S1_omega(:,1)))
figure
plot(des_freq_r,real(S1_omega(:,1)))

figure
plot(des_freq_r,imag(S1_omega(:,2)))
figure
plot(des_freq_r,real(S1_omega(:,2)))

%%
window_op_gg = 1i*window_op_gg; %multiply with door_gg  
window_op_ee = 1i*window_op_ee; %multiply with door_ee
window_op_ef = 1i*window_op_ef;
%%  Calculate the (Louiville space) inner product of the Window and Doorway
om_u_rng = pump_freq; 

Spp_x = zeros(length(om_r_rng),length(t_sep_rng),length(om_u_rng));
%final signal is x only in this scheme (assuming CD effects small etc)
Spp_g = Spp_x; Spp_e = Spp_x; Spp_f = Spp_x; %sep conts
Spp_y = Spp_x; %has almost no physical meaning in that this is heterodyned 
%with the rest of the signal... but may be illustrative somehow
k_u = 2*pi*om_s*kpu; k_r = 2*pi*om_s*kpr; %average amplitudes of these para
kk1 = [-k_u;k_u;k_r]; kk2 = [k_u;-k_u;k_r]; %order of interaction (left to right)
%polarization wavevectors take all along x
pol{1} = [1;0;0]; pol{2} = [1;0;0]; pol{3} = [1;0;0];  
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
                    fct_x1 = 0; fct_x2 = 0; fct_y1= 0; fct_y2 =0;
for j1 = 1:N %I should pre save these but this is quick anyway
   for j2 = 1:N
       for j3 = 1:N
           for j4=1:N
   [fct_x1_tmp,fct_y1_tmp,~] = full_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol,kk1);
  [fct_x2_tmp,fct_y2_tmp,~] = full_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol,kk2); 
  
  cont_fc = c1(j1)*c2(j2)*c3(j3)*c4(j4);
  fct_x1 = fct_x1 + cont_fc*fct_x1_tmp; fct_x2 = fct_x2 + cont_fc*fct_x2_tmp;
  fct_y1 = fct_y1 + cont_fc*fct_y1_tmp; fct_y2 = fct_y2 + cont_fc*fct_y2_tmp;
  %weighting factors for x and y component with each probe ordering
  %usually kk2 lost in RWA anyway
           end
       end
   end
end      
for j = 1:length(om_u_rng)

%should give a om_r_rng by t_sep_rng
tmp_ee = rho_tau(:,:,j,e1,e2);
tmp_gg = rho_tau_p (:,:,j,e1,e2);

trace_ee = mtimesx(window_op_ee(:,:,e3,e4),tmp_ee);
trace_ef = mtimesx(window_op_ef(:,:,e3,e4),tmp_ee);  

trace_gg = mtimesx(window_op_gg(:,:,e3,e4),tmp_gg); 

Spp_x(:,:,j) = Spp_x(:,:,j) + fct_x1*(trace_ee+trace_ef+trace_gg);   
Spp_f(:,:,j) = Spp_f(:,:,j) + fct_x1*(trace_ef); 
Spp_e(:,:,j) = Spp_e(:,:,j) + fct_x1*(trace_ee); 
Spp_g(:,:,j) = Spp_g(:,:,j) + fct_x1*(trace_gg); 
% trace_ee = mtimesx(window_op_ee_y(:,:,e3,e4),tmp_ee);
% trace_ef = mtimesx(window_op_ef_y(:,:,e3,e4),tmp_ee);  
% 
% trace_gg = mtimesx(window_op_gg_y(:,:,e3,e4),tmp_gg );  
% 
% Spp_y(:,:,j) = Spp_y(:,:,j) + fct_y1*(trace_ee+trace_ef+trace_gg);


end                %+ fct_y2 * RWA dropped terms
                
            end
        end
    end
end
t_delay_range_fs = t_sep_rng/convfact*1000;
%include factor of 2 om_r which just comes from solution to maxwells eqns
Spp_x = Spp_x .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_y = Spp_y .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));


%% General pseudo colour figures
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);

figure
pcolor(t_delay_range_fs,om_r_rng,real(Spp_x(:,:,4))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
ylabel('probe frequency \omega, cm^{-1}')  
 colormap(CMRmap)
colorbar
figure
pcolor(t_delay_range_fs,om_r_rng,real(Spp_x(:,:,6))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
ylabel('probe frequency \omega, cm^{-1}')  
 colormap(CMRmap)
colorbar
% 
% figure
% pcolor(t_delay_range_fs,om_r_rng,imag(Spp_y(:,:,1))) 
% shading flat
% set(gcf, 'renderer', 'zbuffer');
% figure
% pcolor(t_delay_range_fs,om_r_rng,imag(Spp_y(:,:,3))) 
% shading flat
% set(gcf, 'renderer', 'zbuffer');

%% Cut throughs
for j = 1:2:3
%lg=om_r_rng<1.71e4 & om_r_rng>1.7e4;
%tmp = max(om_r_rng(lg));
%tmp2 = find(lg,1);
%pnts = tmp2-6:3:tmp2+6;
%pnts = 1:20:size(Spp_x,1);
%pnts = lg;

pnts = 1:5:size(Spp_x,1);

figure

plot(t_delay_range_fs,real(Spp_x(1:5:end,:,j)))
xlabel('Time delay (fs)')
ylabel('Signal')  
end
%%
figure
to_plot2=zeros(length(pnts),floor(length(t_delay_range_fs)/2));
j=3;
for k = 1:length(pnts)
    tmp = real(Spp_x(pnts(k),:,j))- mean(real(Spp_x(pnts(k),:,j)));
%         tmp2 = linspace(0,max(t_delay_range),length(t_delay_range)*10);
%     tmp = interp1(t_delay_range,tmp,tmp2,'pchip');
[to_plot1,to_plot2(k,:)] = ezfft(t_sep_rng,tmp);
end
plot(to_plot1,to_plot2)
xlabel('angular frequency, cm^{-1}')
ylabel('Fourier component amplitude')  

%% Plot frequency slices at delay times

t_pnts = [40,100,150,200,250,300,400,500,1000]; pickr = t_pnts*0;
for k =1 : length(t_pnts)
    [~,pickr(k)] = min(abs(t_delay_range_fs-t_pnts(k)));
end
figure
plot(10^7./om_r_rng,real(Spp_x(:,pickr,1)))
xlabel('Wavelength (nm)')
ylabel('Signal,(au)')  


%%
if 1==0
%   figure
%plot(f_plot2,f_plot2.'.*real(S_1_Neq_FT(f_sel2,1,1))) %ref index
figure
plot(f_plot2,f_plot2.'.*imag(S_1_Neq_FT(f_sel2,1,1))) %abs
[~,b]=max(abs(imag(S_1_Neq_FT(f_sel2,1,1))));
ylabel('\Delta \alpha au')
xlabel('frequency, cm^{-1}')  
  figure
plot(f_plot2,f_plot2.'.*real(S_1_Neq_FT(f_sel2,1,2))) %cd
xlabel('\Delta \eta, cm^{-1}')  
ylabel('CD au')
%figure
%plot(f_plot2,f_plot2.'.*imag(S_1_Neq_FT(f_sel2,1,2))) %or
%%
om_pad = repmat(f_plot2.',1,length(t_sep_rng));
Dalpha_J= (8*pi^2.*om_pad).*imag(S_1_Neq_FT(f_sel2,:,1));
Deta_J  = (16*pi^2.*om_pad).*real(S_1_Neq_FT(f_sel2,:,2));

figure
plot(t_sep_rng/convfact,Dalpha_J([b-5,b,b+5],:)) %abs
ylabel('\Delta \alpha au')
xlabel('t_sep, ps')  
  figure
plot(t_sep_rng/convfact,Deta_J([b-5,b,b+5],:) ) %cd
ylabel('\Delta \eta')  
xlabel('t_sep, ps')

%%
figure
[toplot1,toplot2] = ezfft(t_sep_rng-t_sep_rng(1),Dalpha_J(67,:)-mean(Dalpha_J(67,:)));

plot(t_sep_rng/convfact,Dalpha_J(67,:)) %abs
ylabel('\Delta \alpha au')
xlabel('t_sep, ps')  
  figure
plot(t_sep_rng/convfact,Deta_J(67,:) ) %cd
ylabel('\Delta \eta')  
xlabel('t_sep, ps')

figure
plot(toplot1,toplot2)
%%  2D plots
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);
om_pad = repmat(f_plot2.',1,length(t_sep_rng));

  Dalpha_J= (8*pi^2.*om_pad).*imag(S_1_Neq_FT(f_sel2,:,1));
  p_rng = t_sep_rng>0.25*convfact & t_sep_rng<1*convfact;
  figure
  pcolor(f_plot2 ,t_sep_rng(p_rng) ,Dalpha_J(:,p_rng).')
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('Absorption change \Delta \alpha, for different probe wavelengths and pulse seperations');
 colormap(CMRmap)
colorbar

Deta_J  = (16*pi^2.*om_pad).*real(S_1_Neq_FT(f_sel2,:,2));
  figure
  pcolor( f_plot2,t_sep_rng(p_rng) ,Deta_J(:,p_rng).')
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('CD shift \Delta \eta, for different probe wavelengths and pulse seperations');
colormap(CMRmap)
colorbar


end