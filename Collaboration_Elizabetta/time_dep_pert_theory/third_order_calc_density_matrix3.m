%Used to calculate the third order response function for a dimer using
%Redfield theory and explicitly included, damped Harmonic oscillator modes
%get inverse cm units
Temp = 300; %temp in Kelvin
[convfact, beta,speed_unit]= inv_cm_unit_sys(Temp);

kpr = [0,0,1];   epr = [1,0,0];   bpr = [0,1,0];
%probe unit vectors

theta = atan(sqrt(2)); %angle between pump and probe
kpu = [sin(theta),0,cos(theta)];   epu = [1/sqrt(2),1/sqrt(2),0]; bpu=  cross(kpu,epu);
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
om_vib = [om_0{:}];  numvib = [3,3]
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

        %static disorder parameters for each site, gaussian widths
        s1_static = 200;
        s2_static = 200;

        wav = sqrt(s1_static^2+s2_static^2)/sqrt(2); %shift to average frequency
        %this can be dealt with in a convolution of the central frequency
        delta_w = sqrt(s1_static^2+s2_static^2); %shift to w_2-w_1

        
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
%explicitly included modes, could include static disorder here?
om_vib = [om_0{:}];  numvib = [3,3];
displ = [[0;0],[sqrt(2*lambda{1}/om_0{1}) ; 0],...
            [0 ; sqrt(2*lambda{2}/om_0{2})]...
    [sqrt(2*lambda{1}/om_0{1});sqrt(2*lambda{2}/om_0{2})]];
%displacements of states, ground, excited states and finally double excited

end
%% Static disorder terms

centre_shift =  wav*randn*[1,1]; 
diff_shift = delta_w*randn*[-1,1]; %static disorder can be added in this way, shifts diag
H_site = H_site + diag(diff_shift+centre_shift);


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
    
%     %important for first interaction propogation and window
%     de_fn_ge = @(t,v) supop_ge*v; de_fn_eg = @(t,v) supop_eg*v;     
%     %important for second
%     de_fn_gg = @(t,v) supop_g*v; de_fn_ee = @(t,v) supop_e*v;
%     %only important for third (window)
%     de_fn_ef = @(t,v) supop_ef*v; de_fn_fe = @(t,v) supop_fe*v;
    
    
%%  Calculate the window function
% Assume the pulse is short enough that no time evolution occurs during the
% interaction.
points_to_save = 2^10+1; 

des_freq_r = 1e4*linspace(0.9,1.5,points_to_save); %angular frequency

%choose range of probe frequencies to test for
%~700 - 900nm probe range expected, 
om_r_rng =des_freq_r(10^7./des_freq_r >=700 & 10^7./des_freq_r <=980);

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

%calculate the linear operator that describes linear order CD
lin_op = zeros(length(supop_g),length(des_freq_r),2); %x and y components dim 3
[tmp_gg_sav, tmp_ee_sav, tmp_ef_sav] = window_fun(...
          des_freq_r,H_exciton,R_red_sc,supop_ge ,V_ge,V_eg,tmpge,supop_ef,V_ef,V_fe,tmpef);
      
for e4 = 1:N
    for e3 = 1:N      

        lin_op(:,:,1) = lin_op(:,:,1) + tmp_gg_sav{e3,e4}*xx_ex_av(e4,e3);
        lin_op(:,:,2) = lin_op(:,:,2) + tmp_gg_sav{e3,e4}*yx_ex_av(e4,e3);
        
    end    
end
rho_gg = reshape(rho_fl,sqrt(length(rho_fl)),sqrt(length(rho_fl)));
rho_gg = rho_gg(1:sz2,1:sz2);  rho_gg = reshape(rho_gg,numel(rho_gg),1)';
S1_omega = mtimesx(rho_gg,lin_op);  %first order response function in freq space
S1_omega = permute(S1_omega,[2,3,1]); %put frequency values first

alpha = des_freq_r'.*real(S1_omega(:,1)) ;
CD = des_freq_r'.*imag(S1_omega(:,2));
OR = des_freq_r'.*real(S1_omega(:,2)); 
%use these parameters to calculate the output electric field to first order

%% Calculate window function Heterodyned with the input electric field 
% for left and right circularly polarized light this will be the same
% elements besides averaging of the dipole terms
    
for e4=1:N
    for e3 = 1:N
        for lp =1:length(om_r_rng) %loop over all probe frequencies
            
     E_fct = 2*real(E_r_w(des_freq_r-om_r_rng(lp)).*conj(E_r_w(om_r_rng(lp)-des_freq_r)));
      
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
points_to_save = 2^10+1; 

des_freq_u = linspace(0.9e+04,1.4e+04,points_to_save); %angular frequency
%pic about 8 different frequencies from this which are of interest for some
%reason, e.g. resonant with certain transitions

%calculate the operator for the doorway response function after the
%interaction
[rho_neq_g,rho_neq_e] = doorway_fun(...
          des_freq_u,H_exciton,R_red_sc,supop_ge ,V_ge,V_eg,sz2,V_fe);
      
 %for a range of pump carrier frequencies, calculate the full (integrated 
 % over electric field) doorway wavepacket and propagate in time
pump_freq = [H_ex_vib(sz2+1,sz2+1)-H_ex_vib(2,2),...
    H_ex_vib(sz2+2,sz2+2)-H_ex_vib(1,1),H_ex_vib(sz2+1,sz2+1)-H_ex_vib(1,1),...
H_ex_vib(2*sz2+1,2*sz2+1)-H_ex_vib(2,2),H_ex_vib(2*sz2+2,2*sz2+2)-H_ex_vib(1,1),...
    H_ex_vib(2*sz2+1,2*sz2+1)-H_ex_vib(1,1)];   
 
[rho_neq_g_t,rho_neq_e_t,rho_neq_fg_t] = doorway_prop(...
          om_rng,E_u_w,om_u_range,N,sz2,t_sep_rng,...
          supop_g,supop_e,rho_neq_g,rho_neq_e,supop_fg,rho_neq_fg) ;


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
            %reshape into Louiville space ket
    temp2  =reshape(temp2,numel(rho_neq (:,:,1,e1,e2)),1);
   
    %temp3  =1i*(rho_neq_p(:,:,:,e1,e2)-conj(permute(rho_neq_p(:,:,:,e1,e2),[2,1,3])));
   temp3  =rho_neq_p(:,:,:,e1,e2)+conj(permute(rho_neq_p(:,:,:,e1,e2),[2,1,3]));
   temp3  = trapz(des_freq_u,temp3.*repmat(E_fct,size(temp3,1),size(temp3,2),1 ),3);  
   temp3  =reshape(temp3, numel(rho_neq_p (:,:,1,e1,e2)),1);
 
     %could scale out electronic coherence decay
    output_DE_fun(length(t_sep_rng),temp2,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-20);
    %tic
    ode45(de_fn_ee ,[t_sep_rng(1),t_sep_rng(end)],temp2,options);
    %toc
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_rng),temp2,'get_data+clean');
    tmp2 = interp1(tmp1,tmp2,t_sep_rng).';
     rho_tau(:,:,lp,e1,e2) =  tmp2;        %no scaling here     
             
    output_DE_fun(length(t_sep_rng),temp3,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-20);
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
test3 = reshape(rho_tau_p(:,:,1,1,1),sz2,sz2,length(t_sep_rng));
test3 = test3 - conj(permute(test3,[2,1,3]));
herm_t1 = max(abs(test3(:)));
test4 = reshape(rho_tau(:,:,1,1,1),sz2*N,sz2*N,length(t_sep_rng));
test4 = test4 - conj(permute(test4,[2,1,3]));
herm_t2 = max(abs(test4(:)));


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
[tr_av1,tr_av2,tr_av1i,tr_av2i,herm_t1 ,herm_t2]




%% Linear spec stuff
%use this to calculate the y component of the pump probe signal
rho_gg = reshape(rho_fl,sqrt(length(rho_fl)),sqrt(length(rho_fl)));
rho_gg = rho_gg(1:sz2,1:sz2);  rho_gg = reshape(rho_gg,numel(rho_gg),1)';
S1_omega = mtimesx(rho_gg,lin_op);  %first order response function in freq space
S1_omega = permute(S1_omega,[2,3,1]); %put frequency values first

lg=(10^7./des_freq_r >=700 & 10^7./des_freq_r <=980);
lg2 = 10^7./des_freq_r >=500 & 10^7./des_freq_r <=1110;

alpha = des_freq_r'.*real(S1_omega(:,1));
CD= des_freq_r'.*imag(S1_omega(:,2));
OR = des_freq_r'.*real(S1_omega(:,2)); 
%use these parameters to calculate the output electric field

lam_nm = 10^7./des_freq_r(lg);

figure
plot(lam_nm,des_freq_r(lg)'.*imag(S1_omega(lg,1)))
figure
plot(lam_nm,des_freq_r(lg)'.*real(S1_omega(lg,1)))

figure
plot(lam_nm,des_freq_r(lg)'.*imag(S1_omega(lg,2)))
figure
plot(lam_nm,des_freq_r(lg)'.*real(S1_omega(lg,2)))



%%
for e4=1:N
    for e3 = 1:N
        for lp =1:length(om_r_rng) %loop over all probe frequencies
            
     E_fct = 2*real(E_r_w(des_freq_r-om_r_rng(lp)).*conj(E_r_w(om_r_rng(lp)-des_freq_r)));
     
     %E_fct_y is the y electric field output (divided by the sample length)
     %in the limit of a very short sample (beta L <<1 )  (alpha L <<1 )
     E_fct_y = 2*real(E_r_w(des_freq_r-om_r_rng(lp)).*conj((CD/2+1i*OR)'.*E_r_w(om_r_rng(lp)-des_freq_r)));    
     
       %integrate over the probe wavepacket to produce the correct input
        tmp_gg2 = trapz(des_freq_r, tmp_gg_sav{e3,e4}.*repmat(E_fct,length(supop_g),1),2);      
        tmp_ee2  = trapz(des_freq_r, tmp_ee_sav{e3,e4}.*repmat(E_fct,length(supop_e),1),2);             
        tmp_ef2  = trapz(des_freq_r, tmp_ef_sav{e3,e4}.*repmat(E_fct,length(supop_e),1),2);     
        
       window_op_gg(lp,:,e3,e4) = tmp_gg2; %multiply with door_gg
       
       window_op_ee(lp,:,e3,e4) = tmp_ee2; %multiply with door_ee
       window_op_ef(lp,:,e3,e4) = tmp_ef2; %multiply with door_ee

       %integrate over the probe wavepacket to produce the correct input
        tmp_gg2 = trapz(des_freq_r, tmp_gg_sav{e3,e4}.*repmat(E_fct_y,length(supop_g),1),2);      
        tmp_ee2  = trapz(des_freq_r, tmp_ee_sav{e3,e4}.*repmat(E_fct_y,length(supop_e),1),2);             
        tmp_ef2  = trapz(des_freq_r, tmp_ef_sav{e3,e4}.*repmat(E_fct_y,length(supop_e),1),2);         
       
        window_op2_gg(lp,:,e3,e4) = tmp_gg2;       
        %y-comp is heterodyned with first order signal and so is different
        window_op2_ee(lp,:,e3,e4) = tmp_ee2;      
        window_op2_ef(lp,:,e3,e4) = tmp_ef2;
       
        end
    end
end

% %%
% window_op_gg = -1i*window_op_gg; %multiply with door_gg  
% window_op_ee = -1i*window_op_ee; %multiply with door_ee
% window_op_ef = -1i*window_op_ef;
%%  Calculate the (Louiville space) inner product of the Window and Doorway
om_u_rng = pump_freq; 

Spp_x = zeros(length(om_r_rng),length(t_sep_rng),length(om_u_rng));
%final signal is x only in this scheme (assuming CD effects small etc)
Spp_g = Spp_x; Spp_e = Spp_x; Spp_f = Spp_x; %sep conts
Spp_y = Spp_x; %has almost no physical meaning in that this is heterodyned 
%with the rest of the signal... but may be illustrative somehow
k_u = 2*pi*abs(om_s)*kpu; k_r = 2*pi*abs(om_s)*kpr; %average amplitudes of these para
kk1 = [-k_u;k_u;k_r]; kk2 = [k_u;-k_u;k_r]; %order of interaction (left to right)
%take probe polarization wavevectors along x
pol{3} = [1;0;0];  
pol{1}=  [1;0;0];  pol{2} = pol{1} ; %others along x
pol{1} = [1/sqrt(2);1/sqrt(2);0]; pol{2} = pol{1} ; %45 degrees to probe apparently
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
 %fct_x1
%  fct_y1
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

trace_ee = mtimesx(window_op2_ee(:,:,e3,e4),tmp_ee);
trace_ef = mtimesx(window_op2_ef(:,:,e3,e4),tmp_ee);  
trace_gg = mtimesx(window_op2_gg(:,:,e3,e4),tmp_gg );  
% 
 Spp_y(:,:,j) = Spp_y(:,:,j) + fct_y1*(trace_ee+trace_ef+trace_gg);


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
pcolor(t_delay_range_fs,om_r_rng,real(Spp_y(:,:,4))) 
shading flat
set(gcf, 'renderer', 'zbuffer');
xlabel('Time delay (fs)')
ylabel('probe frequency \omega, cm^{-1}')  
 colormap(CMRmap)
colorbar

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
figure

plot(t_delay_range_fs,real(Spp_y(1:5:end,:,j)))
xlabel('Time delay (fs)')
ylabel('Signal')  
end
%%

to_plot2=zeros(length(pnts),floor(length(t_delay_range_fs)/2));
j=3;
for k = 1:length(pnts)
    tmp = real(Spp_x(pnts(k),:,j))- mean(real(Spp_x(pnts(k),:,j)));
    tmp2 = real(Spp_y(pnts(k),:,j))- mean(real(Spp_y(pnts(k),:,j)));    
%         tmp2 = linspace(0,max(t_delay_range),length(t_delay_range)*10);
%     tmp = interp1(t_delay_range,tmp,tmp2,'pchip');
[to_plot1,to_plot2(k,:)] = ezfft(t_sep_rng,tmp);
[to_plot1y,to_plot2y(k,:)] = ezfft(t_sep_rng,tmp);
end
figure
plot(to_plot1,to_plot2)
xlabel('angular frequency, cm^{-1}')
ylabel('Fourier component amplitude')  
figure
plot(to_plot1y,to_plot2y)
xlabel('angular frequency, cm^{-1}')
ylabel('Fourier component amplitude')  
%% Plot frequency slices at delay times

t_pnts = [40,100,150,200,300,500,1000]; pickr = t_pnts*0;
for k =1 : length(t_pnts)
    [~,pickr(k)] = min(abs(t_delay_range_fs-t_pnts(k)));
end
figure
plot(10^7./om_r_rng,-real(Spp_x(:,pickr,1)))
xlabel('Wavelength (nm)')
ylabel('Signal,(au)')  

figure
plot(10^7./om_r_rng,real(Spp_y(:,pickr,1)))
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