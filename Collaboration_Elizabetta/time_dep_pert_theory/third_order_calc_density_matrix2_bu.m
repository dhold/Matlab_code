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

E_u = @(t) exp(-t.^2/tau_u^2/2) / tau_u / sqrt(2*pi); %init e field env
E_r = @(t) exp(-t.^2/tau_r^2/2) / tau_r / sqrt(2*pi);

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
displ = [[0;0],[sqrt(2*lambda{1}*om_0{1}) ; 0],...
            [0 ; sqrt(2*lambda{2}*om_0{2})]...
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
        R_red2 = conj(permute(reshape(R_red,sz1,sz1,sz1^2),[2,3,1]));
        %reshape(mtimesx(R_redd2,rho_louiville),sz1,sz1) is application
%Expand Redfield operator to include vibrational manifold

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
            temp((a-1)*sz2+cnt1vib,(b-1)*sz2+cnt2vib) = R_red(cnt1,cnt2,a,b);
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

%% Calculate initial condition from the operator, should be sensible thermal dist
points_to_save = 40; t_end_ps = 30; t_end = t_end_ps*convfact;
%t_range = linspace(0,t_end,points_to_save);
rho_fl = zeros(length(supop),1); rho_fl(1) = 1;
%populate the vibrations thermally by allowing system to equilibrate

options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-12);
output_DE_fun(points_to_save,rho_fl,'notsavingnayway'); 
ode45(de_fn,[0,t_end],rho_fl,options);
[tmp1,tmp2]  =  output_DE_fun(points_to_save,rho_fl,'get_data');
rho_fl = (tmp2(end,:)).';
rho_0 = reshape(rho_fl,sqrt(length(rho_fl)),sqrt(length(rho_fl)));
%can check convergence if I want, but this should be a sufficient amount of
%time and close enough for spectroscopy!
if abs(1-reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl) > eps(100)
    warning('eq density matrix evaluated to have trace neq 1')
    trace_discrep = 1-reshape(eye(sz1*sz2),1,sz1^2*sz2^2)*rho_fl
end
%% Calculate linear order "density matrix" 
%(not density matrix until Efield integral is performed)

%set number of points to take and how many to save
points_to_save = 2^7; t_end_ps = tau_u/convfact*6; t_end = t_end_ps*convfact;
t1_range = linspace(0,t_end,points_to_save);
t_step = (t1_range(2)-t1_range(1));
rel_freq_range1 = pi/t_step*linspace(-1,1,points_to_save);
f1_spacing = rel_freq_range1 (2)-rel_freq_range1 (1);

points_to_save = 2^7; t_end_ps = tau_u/convfact*6; t_end = t_end_ps*convfact;
t3_range = linspace(0,t_end,points_to_save);
t_step = (t3_range(2)-t3_range(1));
rel_freq_range3 = pi/t_step*linspace(-1,1,points_to_save);
f3_spacing = rel_freq_range3 (2)-rel_freq_range3 (1);

rho_1D = zeros(length(rho_fl),length(t1_range),N);
rho_3W = zeros(length(rho_fl),length(t3_range),N);

rho_1cD = rho_1D;  rho_3cW = rho_3W; %think these might be ccs


om_char = diag(H_el_ex);om_char= om_char(2:N+1);
%characteristic frequencies of exciton transitons
%cc_map = rho_0*0; cc_map = 

lg_mat1 = H_ex_vib*0; lg_mat2 = H_ex_vib*0;
lg_mat2(1:sz2,sz2+(1:sz2*N)) = 1;
lg_mat1(sz2+(1:sz2*N),1:sz2) = 1;
    mat_s1 = sparse(1:length(H_ex_vib)^2,1:length(H_ex_vib)^2,...
                reshape(lg_mat1,length(H_ex_vib)^2,1));
    mat_s2 = sparse(1:length(H_ex_vib)^2,1:length(H_ex_vib)^2,...
                reshape(lg_mat2,length(H_ex_vib)^2,1));            

for e1 = 1:N   
    
    %Calculate important bit for doorway
    op1 = V_eg{e1};
    
    temp = op1*rho_0;  %op1*rho_0 - rho_0*op1; is first commutator
    % other half will be cc of this
    temp = reshape(temp,numel(temp),1);
    
    Gamma_s = real(R_red(1,1+e1,1,1+e1)); %~coherence decay rate
    om_s = -1i*om_char(e1); %~frequency of oscillation

    
    num_scaling = exp((om_s-Gamma_s)*t1_range); 
    %scaling, also partially scale out exponential decay from the
    %electronic dof
    supop_scaled = supop + (-om_s+Gamma_s)*mat_s1;
    de_fn_scaled = @(t,v) supop_scaled*v;
    
    output_DE_fun(length(t1_range),temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    tic
    ode45(de_fn_scaled,[0,t_end],temp,options);
    toc
    
    [tmp1,tmp2]  =  output_DE_fun(length(t1_range),temp,'get_data');
    % remove bad points
    
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    %interpolate to linear order
    tmp2 = (interp1(tmp1,tmp2,t1_range)).'; %tmp2 output in the wrong shape for this
    
    rho_1D(:,:,e1) = tmp2.*repmat(num_scaling,size(rho_1D,1),1); 
    %note this is rho_1 in the exciton basis
    %doorway function is this times V_{ge}{e2}
    
    %calculate window bit, note operator acts left
    temp = V_eg{e1};  %op1*rho_0 - rho_0*op1; is first commutator
    % other half will be cc of this
    temp = reshape(temp,numel(temp),1);
    
    num_scaling = exp((om_s-Gamma_s)*t3_range); 
    %scaling, also partially scale out exponential decay from the
    %electronic dof
    supop_scaled = supop + (-om_s+Gamma_s)*mat_s1 ;
    %de_fn_scaled = @(t,v) (supop_scaled'*conj(v))'; % <<V| L = (L' |V>>)'
    %note that v should be a row vec but is input as column in ode solvers
    %so just conjugate rather than conjugate transpose twice
    de_fn_scaled = @(t,v) conj(mtimesx(supop_scaled,'C',v,'G'));
    
    output_DE_fun(length(t3_range),temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    tic
    ode45(de_fn_scaled,[0,t_end],temp,options);
    toc
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),temp,'get_data');
    % remove bad points
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    %interpolate to linear order

    tmp2 = (interp1(tmp1,tmp2,t3_range)).'; %tmp2 output in the wrong shape for this
    
    rho_3W(:,:,e1) = tmp2.*repmat(num_scaling,size(rho_3W,1),1); 
    %note this is rho_1 in the exciton basis
    %doorway function is this times V_{ge}{e2}    
    
  %**********************************************  
 % CC, should be anyway??
 
   %this calculates the commutator type term, think this should be cc
    op1 = V_ge{e1};
    
    temp = rho_0*op1; 
    temp = reshape(temp,numel(temp),1);
    
    num_scaling = exp((-om_s-Gamma_s)*t1_range); 
    supop_scaled = supop + (om_s+Gamma_s)*mat_s2 ;
    de_fn_scaled = @(t,v) supop_scaled*v;
    
    output_DE_fun(points_to_save,temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn_scaled,[0,t_end],temp,options);
    
    [tmp1,tmp2]  =  output_DE_fun(points_to_save,temp,'get_data');
    % remove bad points
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    %interpolate to linear order
    tmp2 = (interp1(tmp1,tmp2,t1_range)).'; %tmp2 output in the wrong shape for this
    
    rho_1cD(:,:,e1) = tmp2.*repmat(num_scaling,size(rho_1D,1),1); 
    %note this is rho_1 in the exciton basis
    %doorway function is this times V_{ge}{e2}
    
    %calculate window bit, note operator acts left
    temp = V_ge{e1};  
    temp = reshape(temp,numel(temp),1);
    
    num_scaling = exp((-om_s-Gamma_s)*t3_range); 
    supop_scaled = supop + (om_s+Gamma_s)*mat_s2 ;
    de_fn_scaled = @(t,v) conj(mtimesx(supop_scaled,'C',v,'G'));
    
    output_DE_fun(points_to_save,temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn_scaled,[0,t_end],temp,options);
 
    [tmp1,tmp2]  =  output_DE_fun(points_to_save,temp,'get_data');
    % remove bad points
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    %interpolate to linear order
    tmp2 = (interp1(tmp1,tmp2,t3_range)).'; %tmp2 output in the wrong shape for this
    
    rho_3cW(:,:,e1) = tmp2.*repmat(num_scaling,size(rho_1D,1),1); 
    %note this is rho_1 in the exciton basis
    %window function is this times V_{ge}{e2}    
 
    
end
%%
test = 1i*(reshape(rho_3W(:,:,1), sz1*sz2,sz1*sz2,length(t3_range))-...
           reshape(rho_3cW(:,:,1), sz1*sz2,sz1*sz2,length(t3_range))) ;
test2 = 1i*(reshape(rho_1D(:,:,1), sz1*sz2,sz1*sz2,length(t1_range))-...
           reshape(rho_1cD(:,:,1), sz1*sz2,sz1*sz2,length(t1_range))) ;  
test1a =  test - conj(permute(test,[2,1,3])); %should be zero
test2a =  test2 - conj(permute(test2,[2,1,3])); %should be zero
%they generally aren't... this is concerning

 test4 =  rho_1D(reshape(logical(lg_mat1),length(H_ex_vib)^2,1),:,1) -...
             conj(rho_1cD(reshape(logical(lg_mat2),length(H_ex_vib)^2,1),:,1));      
%% Calculate linear spec stuff based on rho_1D

om_char = mean(eig(H_el(2:N+1,2:N+1))); %characteristic frequency
k_pr = 2*pi*om_char; %assumed to go along z%
shift_param = om_char; %shift centre of fft 
ft_scale_shift =exp(1i*shift_param*t1_range); %shift ft along
S_1 = zeros(2,length(t1_range));

[ex_basis,~] = eig(H_el(2:N+1,2:N+1));

tic
for j1 = 1:N
            
        %averages are much easier to calculate in the site basis so
        %take superposition of rho_1D(:,:,j1) corresponding to site basis
        %mu_ex_k = sum_j c_kj mu_j, express as Liouville op
        % mu_j rho = sum_k b_{kj} mu_ex_k 
       temp = ex_basis(j1,1) *  rho_1D(:,:,1) ;
       %temp = ex_basis(j,j1) *  rho_1D(:,:,1) ;
        for j = 2:N
           temp = temp + ex_basis(j1,j) * rho_1D(:,:,j);
        %temp = temp + ex_basis(j,j1) * rho_1D(:,:,j);
        end
       tempc = ex_basis(j1,1) *  rho_1cD(:,:,1) ;
        for j = 2:N
            tempc = tempc + ex_basis(j1,j) * rho_1cD(:,:,j);
        end      
        temp = reshape(temp,sz1*sz2,sz1*sz2,length(t1_range));
        %this is currently still in the exciton basis, project back
        temp = mtimesx(M_prj,'C',temp);
        temp = mtimesx(temp,M_prj);
        
    for j2 = 1:N
        
        xx_av = dot(mu(j1,:),mu(j2,:))/3;
        yx_av = 1i*k_pr*dot(mu(j1,:),cross(mu(j2,:),R(j1,:)-R(j2,:)))/6;
%         op2 = ex_basis(1,j2) * reshape(V_ge{1}+V_eg{1},numel(V_ge{1}),1);
%                for j = 2:N
%         op2 = op2+  ex_basis(j,j2) * reshape(V_ge{j}+V_eg{j},numel(V_ge{j}),1);
%                end
        op2 = (V_ge{j2}+V_eg{j2});
               temp2 = mtimesx(op2,temp);
               %temp2c = mtimesx(temp,'C',op2));
%                
%         %need to perform a trace over the matrix in the first two
%         %dimensions, sadly matlab doesn't let you do this so I use a custom
%         %function
         R1 = reshape(diagsum(temp2,1,2),size(t1_range));
        
               %R1 = mtimesx(op2,'C',temp);  %this performs a trace
               %R1c = mtimesx(op2,'C',tempc); 
        
        temp3 = 1i*(R1 - conj(R1));
        
        if xx_av ~=0
        S_1(1,:) = S_1(1,:) + ft_scale_shift.*temp3.*xx_av;
        end
        if yx_av ~=0
        S_1(2,:) = S_1(2,:) + ft_scale_shift.*temp3.*yx_av; 
        end
        
    end
end
   toc     
     
   S_1_ft = fftshift(fft(S_1,[],2),2); %shift zero frequency to the centre
 %%  plot to check etc
   
   lf = length(rel_freq_range1);
   %f_sel = [((lf - floor(lf/20)) : lf),1:floor(lf/20)]; %range to take from fft
   %f_plot = [-rel_freq_range1(floor(lf/20)+1:-1:2),rel_freq_range1(1:floor(lf/20)+1)] + shift_param;
   f_plot = rel_freq_range1 + shift_param;
   lambda_nm = 10^7./f_plot;
   f_sel = lambda_nm>500 & lambda_nm < 700;
    f_plot = f_plot(f_sel);
   %shift freq range
   %plot the figures of important things like 
   
   figure
plot(f_plot,f_plot.*real(S_1_ft(1,f_sel))) %ref index
figure
plot(f_plot,f_plot.*imag(S_1_ft(1,f_sel))) %abs
   
  figure
plot(f_plot,f_plot.*real(S_1_ft(2,f_sel))) %cd
figure
plot(f_plot,f_plot.*imag(S_1_ft(2,f_sel))) %or


%% Calculate doorway and window function

rho_door = zeros([size(rho_1D),N]); rho_door_p = rho_door;
rho_window = zeros([size(rho_3W),N]); rho_window_p = rho_window;
%for e1 = 1:N
    for e2 = 1:N
        %(:,:,e1)
        %doorway wave packet
rho_door(:,:,:,e2) = reshape(mtimesx(reshape(rho_1D,sz1*sz2,sz1*sz2,...
                    length(t1_range),N),V_ge{e2}),sz1^2*sz2^2,...
                    length(t1_range),N);
rho_door_p(:,:,:,e2) = reshape(mtimesx(V_ge{e2},reshape(rho_1D,sz1*sz2,sz1*sz2,...
                    length(t1_range),N)),sz1^2*sz2^2,length(t1_range),N);                
                
rho_window(:,:,:,e2) = reshape(mtimesx(reshape(rho_3W,sz1*sz2,sz1*sz2,...
                    length(t3_range),N),V_ge{e2}),sz1^2*sz2^2,...
                    length(t3_range),N);
                
rho_window_p(:,:,:,e2) = reshape(mtimesx(V_ge{e2},reshape(rho_3W,sz1*sz2,sz1*sz2,...
                    length(t3_range),N)),sz1^2*sz2^2,length(t3_range),N);                
  %p denotes "prime"              
    end
%end

%% solve for signal in X, assumes normal detection heterodyned with itself


t_prime_max = 8*tau_u; int_points = 4000;
t_max  = 8*tau_r;
De = zeros(size(rho_door)); Dg = De;
We = zeros(size(rho_window)); Wg = We;
%operator is primed due to the fact that it is time reversed
no_evo_approx = true; %ignore evolution during t and t' 


if no_evo_approx %exp(+/- i H_m t/t') ~ 1
        
        t_int_rng = t_prime_max*linspace(-1,1,int_points);
    
     for tlp = 1:length(t1_range); 
             t1 = t1_range(tlp);
    Eu_fc = trapz(t_int_rng,E_u(t_int_rng - t1) .* conj(E_u(t_int_rng)));         
    De(:,tlp,:,:) = Eu_fc*rho_door(:,tlp,:,:); %particle prop
    Dg(:,tlp,:,:) = Eu_fc*rho_door_p(:,tlp,:,:); %hole prop
     end
      t_int_rng = t_max*linspace(-1,1,int_points);
     for tlp = 1:length(t3_range); 
             t3 = t3_range(tlp);
    Eu_fc = trapz(t_int_rng,E_r(t_int_rng) .* conj(E_r(t_int_rng) + t3));         
    We(:,tlp,:,:) = Eu_fc*rho_window(:,tlp,:,:);
    Wg(:,tlp,:,:) = Eu_fc*rho_window_p(:,tlp,:,:);
     end     
     
else %I'm not convinced I have this correct, not sure how to deal with the
    %-ve time evo, DO NOT trust results using this as it is
    
for e1=1:N
    for e2 = 1:N

    %Gamma_s = real(R_red(1,1+e2,1,1+e2)); %don't include decay scaling
    %om_s = -1i*om_char(e1); %don't include freqscaling
    %num_scaling = exp((-om_s)*t_prime_range); 
    %supop_scaled = supop + (om_s)*mat_s2 ;
    de_fn1 = @(t,v) mtimesx(supop,'C',v);  %right acting conj op     
    de_fn2 = @(t,v) conj(mtimesx(supop,'C',v,'G'));  %left acting op   
 tic       
for tlp = 1:length(t1_range);
   
    %prop system to t_prime_max
    
    temp = rho_door(:,tlp,e1,e2);
    %solve ODE
    
     output_DE_fun(int_points,temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn1,[0,t_prime_max],temp,options);
    time_save(tlp,e1,e2) = toc;
    %get data and filter points with no data (too many points requested)
    [tmp1,tmp2]  =  output_DE_fun(points_to_save,temp,'get_data');   
     tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
     tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);   %time dim 2
     %no need to interpolate as I am integrating over this range
     %***** This might actually cause problems if the number of points is
     %low and so seperation is sig frac of pulse width****
     temp = rho_door_p(:,tlp,e1,e2);
     % same for rho_D'
     
      output_DE_fun(int_points,temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn1,[0,t_prime_max],temp,options);
    time_save2(tlp,e1,e2) = toc;
    %get data and filter points with no data (too many points requested)
    [tmp1g,tmp2g]  =  output_DE_fun(points_to_save,temp,'get_data');   
     tmp3g = [true; diff(tmp1g)~=0] & [true; tmp1g(2:end)~=0];
     tmp1g = tmp1g(tmp3g); tmp2g = tmp2g(tmp3g,:);   %time dim 1   
        
    %calculate electric fields from laser pulse at these times
     t1 = t1_range(tlp);
    Eu_fc = E_u(tmp1 - t1) .* conj(E_u(tmp1));
    Eu_fc = repmat(Eu_fc,1,size(tmp2,2));
    Eu_fcg = E_u(tmp1g - t1) .* conj(E_u(tmp1g));
    Eu_fcg = repmat(Eu_fcg,1,size(tmp2g,2));     
    %defined on page 374 of Mukamel, I need the Right-FT of this eventually
    
    %particle prop
    De(:,tlp,e1,e2) = trapz(tmp1,Eu_fc.*tmp2,1); %int over dim 1
    %hole prop
    Dg(:,tlp,e1,e2) = trapz(tmp1g,Eu_fcg.*tmp2g,1) ;
    
end
toc
    tic
for tlp = 1:length(t3_range);
   
    %prop system to t_prime_max
    
    temp = rho_window(:,tlp,e1,e2); %this one is propagated in ground state
    %solve ODE

     output_DE_fun(int_points,temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn2,[0,t_max],temp,options);
    time_save3(tlp,e1,e2) = toc;
    %get data and filter points with no data (too many points requested)
    [tmp1g,tmp2g]  =  output_DE_fun(points_to_save,temp,'get_data');   
     tmp3g = [true; diff(tmp1g)~=0] & [true; tmp1g(2:end)~=0];
     tmp1g = tmp1g(tmp3g); tmp2g = tmp2g(tmp3g,:);   %time dim 2
     %no need to interpolate as I am integrating over this range
     temp = rho_window_p(:,tlp,e1,e2);
     % same for rho_W'

      output_DE_fun(int_points,temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn2,[0,t_max],temp,options);
    time_save4(tlp,e1,e2) = toc;
    %get data and filter points with no data (too many points requested)
    [tmp1,tmp2]  =  output_DE_fun(points_to_save,temp,'get_data');   
     tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
     tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);   %time dim 2    
        
    %calculate electric fields from laser pulse at these times
     t3 = t3_range(tlp);
    Er_fc = conj(E_r(tmp1 + t3) ).* E_r(tmp1);
    Er_fc = repmat(Er_fc,1,size(tmp2,2));
    Er_fcg = conj(E_r(tmp1g + t3) ).* E_r(tmp1g);
    Er_fcg = repmat(Er_fcg,1,size(tmp2g,2));
    
    %defined on page 374 of Mukamel, I need the Right-FT of this eventually
    We(:,tlp,e1,e2) = trapz(tmp1,Er_fc.*tmp2,1) ; %int over dim 1
    Wg(:,tlp,e1,e2) = trapz(tmp1g,Er_fcg.*tmp2g,1) ;
    
end
toc
    end
end
end
%% Ft into frequency space, scaling to central freq range

%fft_scaleD = repmat(exp(1i*omchar*t1_range),size(De,1),1); 
W_sc_fc = om_char;
fft_scaleW = repmat(exp(1i*W_sc_fc*t3_range),size(We,1),1);

%consider only a few discrete frequencies for om_u as I need to time prop
%this and such like, also don't want lots to plot
[~,b] =eig(H_ex_vib); %eigenvalues of the system
b = diag(b);

om_u_rng = [b(sz1+1)-b(1),b(sz1+1)-b(2),b(2*sz1+1)-b(1),b(2*sz1+1)-b(2),1.7065e+04];
%corresponding to ground vib -> 1ex no vib and 1st vib -> 1ex etc
We_om = zeros(size(We)); Wg_om = We_om;
De_om = zeros(size(De,1),length(om_u_rng),N,N);  Dg_om = De_om;

for e1=1:N
    for e2=1:N
       % De(:,:,e1,e2) = fft(De(:,:,e1,e2).*fft_scaleD,[],2);
       % Dg(:,:,e1,e2) = fft(Dg(:,:,e1,e2).*fft_scaleD,[],2);
        We_om(:,:,e1,e2) = fftshift(fft(We(:,:,e1,e2).*fft_scaleW,[],2),2);
        Wg_om(:,:,e1,e2) = fftshift(fft(Wg(:,:,e1,e2).*fft_scaleW,[],2),2);
        
        for j=1:length(om_u_rng)
            om_pu_fc = exp(1i*om_u_rng(j)*t1_range);
            om_pu_fc = reshape(om_pu_fc,1,size(De,2));
            om_pu_fc = repmat(om_pu_fc,size(De,1),1);
            
            tmp = reshape(De(:,:,e1,e2).*om_pu_fc,...
                sqrt(size(De,1)),sqrt(size(De,1)),length(t1_range));
            tmp = 1i*trapz(t1_range,tmp-conj(permute(tmp,[2,1,3])),3);
       De_om(:,j,e1,e2) = reshape((tmp+tmp')/2,size(De,1),1);  %ensure hermit    
            
            tmp2 = reshape(Dg(:,:,e1,e2).*om_pu_fc,...
                sqrt(size(Dg,1)), sqrt(size(Dg,1)),length(t1_range));
            tmp2 = 1i*trapz(t1_range,tmp2-conj(permute(tmp2,[2,1,3])),3);
       Dg_om(:,j,e1,e2) = reshape((tmp2+tmp2')/2,size(De,1),1);      
                
        end
        
    end
end

%% Calculate transient pump probe signal at a range of delays
t_delay_range_fs = 0:5:3000;
t_delay_range = t_delay_range_fs/1000*convfact;
%also do the dipole averaging, which is why I took the tedious decision to
%calculate these things for each state individually
De_om_tau = zeros(size(De_om,1),length(t_delay_range),length(om_u_rng),N,N);
Dg_om_tau = zeros(size(De_om,1),length(t_delay_range),length(om_u_rng),N,N);

de_fn = @(t,v) mtimesx(supop,v); 
for e1=1:N
    for e2=1:N
        for j=1:length(om_u_rng)

    init = De_om(:,j,e1,e2); 
    output_DE_fun(length(t_delay_range),init,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-5);
    ode45(de_fn,[0,max(t_delay_range)],init,options);
    [tmp1,tmp2]  =  output_DE_fun(points_to_save,temp,'get_data');   
     tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
     tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);   %time dim 2
     %any(any(isnan(tmp2)))
     tmp2 = interp1(tmp1,tmp2,t_delay_range);
      %any(any(isnan(tmp2)))
    De_om_tau(:,:,j,e1,e2) = tmp2.';
    
     init = Dg_om(:,j,e1,e2); 
    output_DE_fun(length(t_delay_range),init,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-5);
    ode45(de_fn,[0,max(t_delay_range)],init,options);
    [tmp1,tmp2]  =  output_DE_fun(points_to_save,temp,'get_data');   
     tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
     tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);   %time dim 2
     tmp2 = interp1(tmp1,tmp2,t_delay_range);   
        
    Dg_om_tau(:,:,j,e1,e2) =  tmp2.';
        end  
    end
end


%% Calculate final signal
t3_step = t3_range(2) - t3_range(1);
om_r_rng = linspace(-1/(t3_step),1/(t3_step),length(t3_range))*pi;
om_r_rng = om_r_rng+om_char; %shift due to fft

Spp_x = zeros(length(om_r_rng),length(t_delay_range),length(om_u_rng));
%final signal is x only in this scheme (assuming CD effects small etc)
Spp_y = Spp_x; %has almost no physical meaning in that this is heterodyned 
%with the rest of the signal... but may be illustrative somehow
k_u = 2*pi*om_char*kpu; k_r = 2*pi*om_char*kpr; %average amplitudes of these para
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
    
    temp1 = reshape(We_om(:,:,e3,e4),sz1*sz2,sz1*sz2,length(om_r_rng),1);
    temp1 = 1i*(temp1 - conj(permute(temp1,[2,1,3]))); %subtract HC
    temp2 = reshape(De_om_tau(:,:,j,e1,e2),sz1*sz2,sz1*sz2,1,length(t_delay_range));   
   tmp1 =  mtimesx(repmat(temp1,1,1,1,length(t_delay_range)),...
                   repmat(temp2,1,1,length(om_r_rng),1));
   tmp1 = diagsum(tmp1,1,2); %take trace along den mat dimensions
   
    temp1 = reshape(Wg_om(:,:,e3,e4),sz1*sz2,sz1*sz2,length(om_r_rng),1);
    temp1 = 1i*(temp1 - conj(permute(temp1,[2,1,3]))); %subtract HC
    temp2 = reshape(Dg_om_tau(:,:,j,e1,e2),sz1*sz2,sz1*sz2,1,length(t_delay_range));     
   tmp2 =  mtimesx(repmat(temp1,1,1,1,length(t_delay_range)),...
                   repmat(temp2,1,1,length(om_r_rng),1));
   tmp2 = diagsum(tmp2,1,2); %take trace along den mat dimensions
   
         Spp_x(:,:,j) = Spp_x(:,:,j) + fct_x1*(tmp1+tmp2);
         Spp_y(:,:,j) = Spp_y(:,:,j) + fct_y1*(tmp1+tmp2);
end                %+ fct_y2 * RWA dropped terms
                
            end
        end
    end
end

Spp_x = Spp_x .* 2.*repmat(om_r_rng.',1,length(t_delay_range),length(om_u_rng));
Spp_y = Spp_y .* 2.*repmat(om_r_rng.',1,length(t_delay_range),length(om_u_rng));


%% General pseudo colour figures
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);

figure
pcolor(t_delay_range_fs,om_r_rng,real(Spp_x(:,:,1))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
ylabel('probe frequency \omega, cm^{-1}')  
 colormap(CMRmap)
colorbar
figure
pcolor(t_delay_range_fs,om_r_rng,real(Spp_x(:,:,3))) 
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
lg=om_r_rng<1.71e4 & om_r_rng>1.7e4;
tmp = max(om_r_rng(lg));
tmp2 = find(lg,1);
pnts = tmp2-6:3:tmp2+6;
%pnts = 1:20:size(Spp_x,1);
%pnts = lg;

figure

plot(t_delay_range_fs,real(Spp_x(pnts,:,j)))
xlabel('Time delay (fs)')
ylabel('Signal')  
end
%%
figure
to_plot2=zeros(length(pnts),floor(length(t_delay_range)/2));
j=3;
for k = 1:length(pnts)
    tmp = real(Spp_x(pnts(k),:,j))- mean(real(Spp_x(pnts(k),:,j)));
%         tmp2 = linspace(0,max(t_delay_range),length(t_delay_range)*10);
%     tmp = interp1(t_delay_range,tmp,tmp2,'pchip');
[to_plot1,to_plot2(k,:)] = ezfft(t_delay_range,tmp);
end
plot(to_plot1,to_plot2)
xlabel('Time delay (fs)')
ylabel('frequency, cm^{-1}')  

%%

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
