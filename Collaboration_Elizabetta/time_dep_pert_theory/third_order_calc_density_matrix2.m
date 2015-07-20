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
Hwhm_u_fs = 18; %half width at half maximum, pump
Hwhm_r_fs = 18; %half width at half maximum, probe

%I'm not sure this normalisation is a good idea, if I change the pulse
%width then I change the energy in the pulse... oh well
gaussian_pulses = true;
if gaussian_pulses
tau_u = (Hwhm_u_fs/sqrt(2*log(2)))/1000*convfact; %pumpp pulse SD in inverse cm
tau_r = (Hwhm_r_fs/sqrt(2*log(2)))/1000*convfact; %probe pulse SD in inverse cm    
E_u = @(t) exp(-t.^2/tau_u^2/2) / tau_u / sqrt(2*pi); %init e field env
E_r = @(t) exp(-t.^2/tau_r^2/2) / tau_r / sqrt(2*pi);
E_u_w = @(om) exp(-tau_u^2.*om.^2/2);
E_u_inc = @(om,t) exp(-tau_u^2.*om.^2/2)*(1+1i*erfi(tau_u^2*om-1i*t)/sqrt(2)/tau_u)/2;
E_r_w = @(om) exp(-tau_r^2.*om.^2/2);
else  %sech pulses
tau_u = Hwhm_u_fs/asech(1/sqrt(2))/2/1000*convfact; 
tau_r = Hwhm_r_fs/asech(1/sqrt(2))/2/1000*convfact;
E_u = @(t) sech(t/2/tau_u).^2 / tau_u / 4; %init e field env
E_r = @(t)  sech(t/2/tau_r).^2 / tau_r / 4;
E_u_w = @(om) om.*csch(-pi*tau_u.*om)*sqrt(pi/2);
%E_u_inc = @(om,t) don't know
E_r_w = @(om) om.*csch(-pi*tau_r.*om)*sqrt(pi/2);    
    
end
    

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
om_0={640,640}
om_vib = [om_0{:}];
numvib = [4,4];
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

H0_renorm = H_site - diag(cellfun(@sum, lam_dru) - cellfun(@sum, lambda));
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
%H_vib_op = zeros(length(H_ex_vib));
for k = 1:length(numvib)
nav = exp(-beta*om_vib(k))/(1-exp(-beta*om_vib(k)));

%express the annihilation op in the full hilbert space
aa = kron(sparse(eye(sz1)),sparse(eye(prod(numvib(1:k-1)))));
aa = kron(aa, diag(sqrt(1:numvib(k)-1),1)); 
aa = kron(aa,sparse(eye(prod(numvib(k+1:end))))); 
adag = aa'; %conjugate
%project into exciton basis (not ex vib basis!) as Redfield op in ex basis
aa = M_prj'*aa*M_prj;  adag = M_prj'*adag*M_prj;

Lindblad_op_indiv{k} =  gamma{k}*((nav+1)*L_so(aa) + nav*L_so(adag));  
%H_vib_op  = H_vib_op  + (adag*aa+eye(size(aa))/2)*om_0{1};
Lindblad_op = Lindblad_op  + Lindblad_op_indiv{k};
end    
Lindblad_op = sparse(Lindblad_op);
clear Lindblad_op_indiv
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
 %R_red = real(R_red);
 %take out imaginary elements of R_red(a,b,a,b) which are just rescales of
 %the transition frequencies omega_{ab}, due to redundencies need only
 %rescale omega_{0 e} and omega_{0 f} etc etc as otheres can be expressed
 %in terms of these
 R_red_sc  = R_red; %scaled redfield
ham_rs = zeros(1,sz2);  %will rescale Hamiltonian
if 1 == 0  %change to remove
 for k = 1:length(R_red_sc)
    R_red_sc(k,k,k,k) =  real(R_red(k,k,k,k)); %remove eps level imag parts
 end
 % Now also remove all the eps level mix matches that mean R(a,b,a,b) ~=
 % conj(R(b,a,b,a)) by averaging
 for k = 1:length(R_red_sc)
     for j= (k+1):length(R_red_sc)

    R_red_sc(k,j,k,j) = (R_red(k,j,k,j) + conj(R_red(j,k,j,k)))/2; 
    R_red_sc(j,k,j,k) = conj(R_red_sc(k,j,k,j)); 
         
     end
 end 
 
 
 for a=2:N+1
        R_red_sc(1,a,1,a) = R_red_sc(1,a,1,a) - 1i*imag(R_red(1,a,1,a));    
        R_red_sc(a,1,a,1) = R_red_sc(a,1,a,1) - 1i*imag(R_red(a,1,a,1)); 
        
       for b=N+1+(1:(N-1)*N/2)
        R_red_sc(1,b,1,b) = R_red_sc(1,b,1,b) - 1i*imag(R_red(1,b,1,b));    
        R_red_sc(b,1,b,1) = R_red_sc(b,1,b,1) - 1i*imag(R_red(b,1,b,1));           
%          actually don't rescale these, just ground->whatever
%         R_red_sc(a,b,a,b) = R_red_sc(a,b,a,b) - 1i*imag(R_red(a,b,a,b));
%         R_red_sc(b,a,b,a) = R_red_sc(b,a,b,a) - 1i*imag(R_red(b,a,b,a));        
       end                    
 end

for b= 2:size(R_red,1)
            ham_rs = [ham_rs,ones(1,sz2)*imag(R_red(1,b,1,b))]; 
end            
        
        H_ex_vib = H_ex_vib + diag(ham_rs);
        H_exciton = H_exciton + diag(ham_rs);
        %if I run this cell multiple times it will keep rescaling this, this
        %tests whether I have regenerated it before doing this
        if H_ex_vib_been_rs
            warning('H_ex_vib has been rescaled twice or more')
        else
            H_ex_vib_been_rs = true;
        end
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
 R_red_op_full = sparse( R_red_op_full);
%toc
%N2*(a-1)+A+ N*N2*(N2*(b-1)+B) 

decoherence_op = Lindblad_op - R_red_op_full; 

L_op = -1i*(kron(eye(length(H_exciton)),H_exciton)-...
            kron(H_exciton.',eye(length(H_exciton))));    

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
    %note that if I add spont emmision this is not the case

%reduced operators acting only on 1st ex state manifold

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

options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-16);
output_DE_fun(points_to_save,[1;zeros(length(supop_g)-1,1)],'notsavingnayway'); 
ode45(de_fn_gg,[0,t_end],[1;zeros(length(supop_g)-1,1)],options);
[tmp1,tmp2]  =  output_DE_fun(points_to_save,rho_fl,'get_data');

rho_0 =    full(sparse(ag,bg,tmp2(end,:),sz1*sz2,sz1*sz2));
rho_fl = reshape(rho_0,numel(rho_0),1);

if abs(1-reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl) > eps(100)
    warning('eq density matrix evaluated to have trace neq 1')
    trace_discrep = 1-reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl
    rho_fl = rho_fl/(reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl);
    rho_0 = rho_0/(reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl);
end
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
   
    %temp2  = 1i*(rho_neq(:,:,:,e1,e2)-conj(permute(rho_neq(:,:,:,e1,e2),[2,1,3]))); 
    temp2  = rho_neq(:,:,:,e1,e2)+conj(permute(rho_neq(:,:,:,e1,e2),[2,1,3]));

    
     %add complex conjugate from other half of the density matrix
    temp2  = trapz(des_freq_u,temp2.*repmat(E_fct,size(temp2,1),size(temp2,2),1 ),3);
            
%     tmpee  = zeros(sz1*sz2,sz1*sz2); 
%     tmpee(sz2+1:sz2*(N+1),sz2+1:sz2*(N+1))=temp2;
% 	tmpee = reshape(tmpee,sz2^2*sz1^2,1);      
%     tmpee = tmpee(tmpe,:); %take section corresponding to the first excited manifold density matrix            
  temp2 = reshape(temp2, numel(rho_neq (:,:,1,e1,e2)),1);  %reshape into Louiville space ket
%   test = abs(temp2-tmpee); sum(test(:))
    
    
    %temp3  =1i*(rho_neq_p(:,:,:,e1,e2)-conj(permute(rho_neq_p(:,:,:,e1,e2),[2,1,3])));
   temp3  =rho_neq_p(:,:,:,e1,e2)+conj(permute(rho_neq_p(:,:,:,e1,e2),[2,1,3]));
   temp3  = trapz(des_freq_u,temp3.*repmat(E_fct,size(temp3,1),size(temp3,2),1 ),3);  
   
%   tmpgg  = zeros(sz1*sz2,sz1*sz2); 
%     tmpgg(1:sz2,1:sz2)=temp3;
% 	tmpgg = reshape(tmpgg,sz2^2*sz1^2,1);      
%     tmpgg = tmpgg(tmpg,:); %take section corresponding to the gs manifold density matrix            

   temp3  =reshape(temp3, numel(rho_neq_p (:,:,1,e1,e2)),1);
   % test = abs(temp3-tmpgg); sum(test(:))   
   
   
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
herm_t2 = max(abs(test4(:))); %this one is not Hermitian sometimes, 
%I really don't know why


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
[tr_av1,tr_av2,tr_av1i,tr_av2i,herm_t1 ,herm_t2]

clear rho_neq rho_neq_p %these take a lot of memory
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
om_r_rng =des_freq_r(10^7./des_freq_r >=700 & 10^7./des_freq_r <=980);

t_step = 2*pi/(des_freq_r(end)-des_freq_r(1)); %set by range
t_end = 2*pi/(des_freq_r(2)-des_freq_r(1)); %set by frequency spacing
t3_range = 0:t_step:t_end ;
t3_end_ps = t_end / convfact; %end time in ps, useful to know I think

%tmp=(0:length(des_freq_r)-1);
%phase_shift=exp(-1i*2*pi*(length(des_freq_r)-1)/2*tmp/length(des_freq_r)); 
% calculating phase shift from the fact the FFT doesn't centre the data on
% zero
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
    %initial window state, fourier transform to freq space, time dim 2
%     figure
%     plot(t3_range,phase(tmp4(1,:)))
%      ylabel('phase2')
%     figure
%     plot(t3_range,abs(tmp4(1,:)))
%      ylabel('abs2')

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
        
       % tmp_gg = reshape(mtimesx(tmp3,op_3,'C')-mtimesx(op_3,tmp3,'C')...
        %                ,length(supop_g),length(des_freq_r));                         
           %complex conjugates are from lower left sections of the
                  %operator, which are conjugates by Hermitivity.
                    
       % tmp_ee = mtimesx(tmp3,'C',op_3)-mtimesx(op_3,'C',tmp3);
       % tmp_ef = mtimesx(tmp4,op_3f,'C')-mtimesx(op_3f,tmp4,'C');
             %must add the V_fe(t) V_ef - V_fe V_ef(t) cont as well
             % this corresponds to e->f->e transitions
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
%lg2 = 10^7./des_freq_r >=500 & 10^7./des_freq_r <=1110;

alpha = des_freq_r'.*real(S1_omega(:,1));
CD = des_freq_r'.*imag(S1_omega(:,2));
OR = des_freq_r'.*real(S1_omega(:,2)); 
%use these parameters to calculate the output electric field

lam_nm = 10^7./des_freq_r(lg);
% 
%  figure
%  plot(lam_nm,des_freq_r(lg)'.*imag(S1_omega(lg,1))) %ref index
figure
plot(lam_nm,alpha(lg)) %abs
% 
 figure
 plot(lam_nm,CD(lg)) %CD
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
% 
%        %integrate over the probe wavepacket to produce the correct input
%         tmp_gg2 = trapz(des_freq_r, tmp_gg_sav{e3,e4}.*repmat(E_fct_y,length(supop_g),1),2);      
%         tmp_ee2  = trapz(des_freq_r, tmp_ee_sav{e3,e4}.*repmat(E_fct_y,length(supop_e),1),2);             
%         tmp_ef2  = trapz(des_freq_r, tmp_ef_sav{e3,e4}.*repmat(E_fct_y,length(supop_e),1),2);         
%        
%         window_op2_gg(lp,:,e3,e4) = tmp_gg2;       
%         %y-comp is heterodyned with first order signal and so is different
%         window_op2_ee(lp,:,e3,e4) = tmp_ee2;      
%         window_op2_ef(lp,:,e3,e4) = tmp_ef2;
%        
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

Spp_alpha = zeros(length(om_r_rng),length(t_sep_rng),length(om_u_rng));
%final signal is x only in this scheme (assuming CD effects small etc)
Spp_g = Spp_alpha; Spp_e = Spp_alpha; Spp_f = Spp_alpha; %sep conts
Spp_CD = Spp_alpha;  Spp_CDf = Spp_CD;Spp_CDe =Spp_CD; Spp_CDg = Spp_CD;

typ_op = 10^7/820;
k_u = 2*pi*abs(typ_op)*kpu; k_r = 2*pi*abs(typ_op)*kpr; %average amplitudes of these para
kk1 = [-k_u;k_u;k_r;-k_r]; kk2 = [k_u;-k_u;k_r;-k_r]; %order of interaction (left to right)
%take probe polarization wavevectors along x
%pol{3} = [1,0,0]; % pol_left = [1;-1i;0]/sqrt(2);  pol_right = [1;1i;0]/sqrt(2);
pol{1}=  [1,0,0];  pol{2} = pol{1} ; %others along x
%pol{1} = [1/sqrt(2),1/sqrt(2),0]; pol{2} = pol{1} ; %45 degrees to probe apparently
pol_L = [pol{1};pol{2};1/sqrt(2),+1i/1/sqrt(2),0;1/sqrt(2),-1i/1/sqrt(2),0];
pol_R = [pol{1};pol{2};1/sqrt(2),-1i/1/sqrt(2),0;1/sqrt(2),+1i/1/sqrt(2),0];
        %dipole averaging functions, random orientation
%    xxxx_av = @(mu,jj) (1/15)*(dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),mu(jj(4),:))+...
%                 dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),mu(jj(4),:))+...
%                 dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),mu(jj(4),:))); 
%             %also yyyy etc                                     
%    yyxx_av = @(mu,jj) (1/30)*(4*dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),mu(jj(4),:))-...
%                 dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),mu(jj(4),:))-...
%                 dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),mu(jj(4),:)));    
%             %also xxzz, xyxy etc            
%    xxxyz_av = @(mu,jj,R,k)  (1/30)*...
%             (dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),cross(mu(jj(4),:),R(jj(k),:)))+...
%              dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),cross(mu(jj(4),:),R(jj(k),:)))+...
%              dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),cross(mu(jj(4),:),R(jj(k),:))));
         
        % alpha_av = zeros(N,N,N,N); CD_av = zeros(N,N,N,N);
         alpha_L_4 = zeros(N,N,N,N); alpha_R_4 = zeros(N,N,N,N);
         alpha_L_5 = zeros(N,N,N,N); alpha_R_5 = zeros(N,N,N,N);         
for j1 = 1:N %pre save these
   for j2 = 1:N
       for j3 = 1:N
           for j4=1:N
         jj = [j1,j2,j3,j4];
               
% [alpha_av(j1,j2,j3,j4),CD_av(j1,j2,j3,j4),CD2_av] = ...
%     abs_CD_cont_3rd_order(mu,R,[j1,j2,j3,j4],pol(1:2),kk1);  
%     if CD2_av~=0
%        warning('other contirbution to CD found, consider this term')
%     end
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
                    fct_alpha= 0; fct_CD= 0; fct_CD2= 0; fct_alpha2 =0;
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

Spp_alpha(:,:,j) = Spp_alpha(:,:,j) + fct_alpha *(trace_ef-(trace_ee+trace_gg));   
Spp_f(:,:,j) = Spp_f(:,:,j) + fct_alpha *(trace_ef); 
Spp_e(:,:,j) = Spp_e(:,:,j) + fct_alpha *(trace_ee); 
Spp_g(:,:,j) = Spp_g(:,:,j) + fct_alpha *(trace_gg); 

% trace_ee = mtimesx(window_op2_ee(:,:,e3,e4),tmp_ee);
% trace_ef = mtimesx(window_op2_ef(:,:,e3,e4),tmp_ee);  
% trace_gg = mtimesx(window_op2_gg(:,:,e3,e4),tmp_gg );  
% 
 Spp_CD(:,:,j) = Spp_CD(:,:,j) + fct_CD *(trace_ef-(trace_ee+trace_gg));
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
Spp_alpha = Spp_alpha .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_CD = Spp_CD .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_f = Spp_f .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_e = Spp_e .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_g = Spp_g .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_CDf = Spp_CDf .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_CDe = Spp_CDe .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
Spp_CDg = Spp_CDg .* 2.*repmat(om_r_rng.',1,length(t_sep_rng),length(om_u_rng));
%% General pseudo colour figures
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);


lambda_nm = 10^7./om_r_rng;

figure
pcolor(t_delay_range_fs,lambda_nm,real(Spp_alpha(:,:,1))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
%ylabel('probe frequency \omega, cm^{-1}')  
ylabel('probe wavelength \lambda, nm')  
 colormap(CMRmap)
colorbar

figure
pcolor(t_delay_range_fs,lambda_nm,real(Spp_CD(:,:,1))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
%ylabel('probe frequency \omega, cm^{-1}') 
ylabel('probe wavelength \lambda, nm')  
 colormap(CMRmap)
colorbar
if 1==0
%%

figure
pcolor(t_delay_range_fs,lambda_nm,real(Spp_CDf(:,:,1))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
%ylabel('probe frequency \omega, cm^{-1}')  
ylabel('probe wavelength \lambda, nm')  
 colormap(CMRmap)
colorbar

figure
pcolor(t_delay_range_fs,lambda_nm,real(Spp_CDe(:,:,1))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
%ylabel('probe frequency \omega, cm^{-1}')  
ylabel('probe wavelength \lambda, nm')  
 colormap(CMRmap)
colorbar

figure
pcolor(t_delay_range_fs,lambda_nm,real(Spp_CDg(:,:,1)))
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
%ylabel('probe frequency \omega, cm^{-1}')  
ylabel('probe wavelength \lambda, nm')  
 colormap(CMRmap)
colorbar
%%

figure
pcolor(t_delay_range_fs,lambda_nm,real(Spp_f(:,:,1))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
%ylabel('probe frequency \omega, cm^{-1}')  
ylabel('probe wavelength \lambda, nm')  
 colormap(CMRmap)
colorbar

figure
pcolor(t_delay_range_fs,lambda_nm,real(Spp_e(:,:,1))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
%ylabel('probe frequency \omega, cm^{-1}')  
ylabel('probe wavelength \lambda, nm')  
 colormap(CMRmap)
colorbar

figure
pcolor(t_delay_range_fs,lambda_nm,real(Spp_g(:,:,1)))
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
%ylabel('probe frequency \omega, cm^{-1}')  
ylabel('probe wavelength \lambda, nm')  
 colormap(CMRmap)
colorbar
end
% figure
% pcolor(t_delay_range_fs,om_r_rng,real(Spp_x(:,:,6))) 
% shading flat
% set(gcf, 'renderer', 'zbuffer');
% xlabel('Time delay (fs)')
% ylabel('probe frequency \omega, cm^{-1}')  
%  colormap(CMRmap)
% colorbar
% % 
% figure
% pcolor(t_delay_range_fs,om_r_rng,imag(Spp_y(:,:,1))) 
% shading flat
% set(gcf, 'renderer', 'zbuffer');
% figure
% pcolor(t_delay_range_fs,om_r_rng,imag(Spp_y(:,:,3))) 
% shading flat
% set(gcf, 'renderer', 'zbuffer');
                                    

%% Cut throughs
for j = 1:1%:length(pump_freq)
%lg=om_r_rng<1.71e4 & om_r_rng>1.7e4;
%tmp = max(om_r_rng(lg));
%tmp2 = find(lg,1);
%pnts = tmp2-6:3:tmp2+6;
%pnts = 1:20:size(Spp_x,1);
%pnts = lg;
pnts = [peaks,peaks+8,troughs,troughs+8]; 
%pnts = 1:5:size(Spp_alpha,1);

figure1=figure;

plot1 = plot(t_delay_range_fs,real(Spp_alpha(pnts,:,j)));
%plot1 = plot(t_delay_range_fs,real(Spp_f(pnts,:,j)));
for el = 1:length(pnts)
set(plot1(el),'DisplayName',strcat('\lambda = ',num2str(lam_nm(pnts(el))),'nm'));
end

xlabel('Time delay (fs)');
ylabel('Signal');
legend('show');

figure

plot2 = plot(t_delay_range_fs,real(Spp_CD(pnts,:,j)));
%plot2 = plot(t_delay_range_fs,real(Spp_CDf(pnts,:,j)));
for el = 1:length(pnts)
set(plot2(el),'DisplayName',strcat('\lambda = ',num2str(lam_nm(pnts(el))),'nm'));
end
xlabel('Time delay (fs)')
ylabel('Signal')  
legend('show');

to_plot2=zeros(length(pnts),floor(length(t_delay_range_fs)/2));
to_plot2y = to_plot2;

for k = 1:length(pnts)
    tmp = real(Spp_alpha(pnts(k),:,j))- mean(real(Spp_alpha(pnts(k),:,j)));
    tmp2 = real(Spp_CD(pnts(k),:,j))- mean(real(Spp_CD(pnts(k),:,j)));    
    %tmp = real(Spp_f(pnts(k),:,j))- mean(real(Spp_f(pnts(k),3*end/4:end,j)));
  %  tmp2 = real(Spp_CDf(pnts(k),:,j))- mean(real(Spp_CDf(pnts(k),3*end/4:end,j))); 
[to_plot1,to_plot2(k,:)] = ezfft(t_sep_rng,tmp);
[to_plot1y,to_plot2y(k,:)] = ezfft(t_sep_rng,tmp2);
end
figure
rng = to_plot1<3000;
plot3 = plot(to_plot1(rng),sqrt(to_plot2(:,rng)));
for el = 1:length(pnts)
set(plot3(el),'DisplayName',strcat('\lambda = ',num2str(lam_nm(pnts(el))),'nm'));
end
xlabel('angular frequency, cm^{-1}')
ylabel('Fourier component amplitude')  
legend('show');
figure
plot4 = plot(to_plot1y(rng),sqrt(to_plot2y(:,rng)));
for el = 1:length(pnts)
set(plot4(el),'DisplayName',strcat('\lambda = ',num2str(lam_nm(pnts(el))),'nm'));
end
xlabel('angular frequency, cm^{-1}')
ylabel('Fourier component amplitude')  
legend('show');
end
%% Plot frequency slices at delay times

t_pnts = [150,250,1000]; 

pickr = t_pnts*0;
for k =1 : length(t_pnts)
    [~,pickr(k)] = min(abs(t_delay_range_fs-t_pnts(k)));
end
figure
plot(10^7./om_r_rng-5,real(-Spp_alpha(:,pickr,1)))
xlabel('Wavelength (nm)')
ylabel('Signal,(au)')  
%%
figure
plot(10^7./om_r_rng,real(Spp_CD(:,pickr,1)))
xlabel('Wavelength (nm)')
ylabel('Signal,(au)')  
