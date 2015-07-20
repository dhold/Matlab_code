%% Generally useful constants
clear rho_0 %either do this or set an init cond
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m^-1
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}
% to run the code not as a function just uncomment and use each module
% Energy in Joules,  extra factor of 10 is because thats speed in cm/s
convfact = 2*pi * light_speed *length_unit  *10^(-12); %2 pi * "c" * 100m^-1 * 1ps 
use_reduced_mat = false;

%% Parameters


%relating to pulses
tsep = 400*convfact/1000;
om1 = 1134; om2 = 1066;
tau1 = 100*convfact/1000/sqrt(log(2)); tau2 = 150*convfact/1000/sqrt(log(2));
t01 = 5*tau1;  %assume first pulse hits at 5 sd away
t02 = t01 + tsep;
%relating to time computed and what is saved
 tendps = 2; %end time in pico seconds
tend =convfact*tendps; %end time in units you get by putting everything in wavenumbers
numpoints = [100, linspace(0,tendps,20).*convfact]; 
save_file_name = 'dimer_data_density_resvib.mat';
%save in ~1 fs time steps

%system properties

 E1 = 1042; %Energy of site one
 E2 = 1042; %Energy of site two
mu = ...
[ -0.823072771604003,	-0.5246284653639453,    -0.21752284012024087;...
-0.21108138821716999,	-0.9529594677893477,	-0.2175152875063172];
N = 3; 
R_disp = 10^(-8)*[72.062,-12.081,90.095; 72.913	,9.924,90.122];
R_disp = R_disp - repmat(mean(R_disp),size(R_disp,1),1);
%R_disp = R_disp*10^(-8); % convert from angstrom to inverse cm

%pi/2 rotations about various axes
R_z = [0,-1,0;1,0,0;0,0,1]; R_y = [0,0,1;0,1,0;-1,0,0]; R_x = [1,0,0;0,0,-1;0,1,0];
%use to average over orientations

 V = 92; %Coupling should be able to calculate this from the dipole maagnitude
 
 omegavib =  184; %Frequency of quantised vibrations / level spacing,
% %assumed to be the same in each site, choose (omega_1+omega_2)/2 if not
% the case
 viblvls = 5; %Number of vibrational levels to include, integer
 coupg = 20; %coupling term between electronic and vibrational DOF
 Temp = 300; %temperature in Kelvin
 beta = (2* pi *hbar * light_speed * length_unit)/ ( Temp * boltz_const);

            use_dru = true; %whether to include drude spec
            ham_for_int_pic = 0;
        Kap2 = 2;  %tier truncation
        Kap1 = inf; %freq truncation
        
 if use_dru 
    gamma_dru = {[],100,100}; %drude decay constant
    lambda_dru = {[],6,6}; % weighting of distrubtion
 % %J_j(omega) = 2 lamda_j gamma_j /hbar * (omega/(omega^2+lambda_j^2)
    Kappa = 0; %truncation parameter satisfying Kappa >>omega_0 beta hbar/2pi
 %beyond this e^(-v_(k>kappa) t)*v_k ~delta(t) is made
 %second truncation parameter, there are Kappa + 1 explicitly 
 %treated exponentials, this number is the maximum order the sum of the
 %number of powers can take, e.g. n_1 + n_2 +... + n_Kappa <= Kap2
 saveuptotier = 1;
 else
   saveuptotier = 0;  
  end



%% Hamiltonian terms considered explicitly
%Excitonic Hamilonian
Hel = zeros(3); 
 
Hel(2,2) = E1; Hel(3,3) = E2;%On diagonal terms
Hel(3,2) = V;  Hel(2,3) = V;  % Coupling terms from dipole dipole interaction
syms E_1 E_2 V_12
Helsym = sym(zeros(3));
Helsym(2,2) = E_1; Helsym(3,3) = E_2;%On diagonal terms
Helsym(3,2) = V_12;  Helsym(2,3) = V_12; % Coupling terms from dipole dipole interaction

%Vibrational Hamilonian, relative degrees of freedom only
Hexvib = 0;
if viblvls > 1 %zero/1 means no quantised vibration considered
tmp = ones(viblvls,1);
tmp(1) = 0;
for k=2:viblvls
    tmp(k) = tmp(k-1) + omegavib;
end
Hvib = diag(tmp); clear tmp 

[basis_proj,newE] = eig( Hel);
basis_proj = kron(basis_proj,eye(viblvls));
%projector to delocalised basis

Htot = kron(Hel,eye(size(Hvib)))+kron(eye(size(Hel)),Hvib); 
%total time independent Hamiltonian
%which will be explicitly considered QM

%Hamilonian coupling electronic and vibrational levels
%Hexvib = -g/sqrt(2)  * (n_1 - n_2)*(b^dag_rel + b_rel)

Hexvib = diag(-coupg/sqrt(2) * sqrt(1:viblvls-1),1);
Hexvib = Hexvib+Hexvib.';
Hexvib = kron([0,0,0;0,-1,0;0,0,1],Hexvib); %ground state doesn't drive vibrations

Htot = Htot + Hexvib;  %Full Hamiltonian with electric + vib + coupling
%in chromophore basis
%Htot = basis_proj'*Htot*basis_proj; %project to exciton basis
else
    [basis_proj,~] = eig(Hel); 
    Htot = Hel; Hvib = 0;
   %rho Htot = basis_proj'*Hel*basis_proj;
end

%coefficients used in the HOM sum
 if use_dru
     
 N2 = sum(cellfun(@length,gamma_dru));
 QQ = zeros(N,2); cnt=0;
 cc = zeros(1,sum(N2 + Kappa));
 clear  cc_com vv cc_acom
    for j=1:N
        [cc1,cc2R,cc2I,vv1,vv2,QQ(j,1),tmp] = ...
    coeffients_from_brownian_new([],[],[],Temp,Kappa,lambda_dru{j},gamma_dru{j});
 % QQ is in general length two for each point as I could (potentially)
 % included high frequency decay terms from poles in J(omega) which would
 % have anti commutator parts
        cc(cnt+1:cnt+length(tmp)) = tmp; cnt =cnt+length(tmp);
        cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
         cc_acom{j}= [cc2I,cc1*0];
     end
     
 end
 %% Pulse parameters
%1-2 microjoule for pump pulse, pulse F width HM 100-150 fs
k1 = [0,0,1]; k2 = sqrt([0,1,2]/3);
E_01 = 100*[0,0,1]; E_02=0*sqrt([0,1,2]/3);  %direction is important

%Work out positions needed to extract polarization in direction k_2
 %this is only true to third order
 %dot(k1,rpos(j,:)) = 0;, dot(k2,rpos(j,:)) = pi/2*(j-1);
 
 phi2 = [0,0,0,0]; phi1 = (0:3) * pi;
 rpos = zeros(4,3); %can take second to be all zeros
 
 if k1(1) == 0 && k1(2) == 0
 %k1 is along z, take this component to be zero    
 %take other positions along y, this is enough to solve
 for j =1:3
     rpos(j+1,2) = pi*j/(2*k2(2));  
 end
 else
 error('you have not written anything for this case yet David')
 end


 
tt =sym('t','real'); omega_1 =sym('w1','real'); omega_2 =sym('w2','real');

syms Env_1 Env_2
 E_t_1 = E_01*Env_1*(exp(-1i*omega_1*tt)+exp(+1i*omega_1*tt));
 % it is generally possible to always pick positions such that the k1.r
 % factors are zero, so I won't include this here
 %Env_1 = exp(-(t-t01)^2/tau1^2);
 
 E_t_2 = E_02*Env_2*(exp(-1i*omega_2*tt+1i*dot(k2,rpos(1,:)))...
                    +exp(+1i*omega_2*tt-1i*dot(k2,rpos(1,:))) );  
 
 %Env_2 = exp(-(t-t02)^2/tau2^2);
 %set to zero for non test runs
 %ideally calculate the density matrix from the impact of the first pulse 
 % for a long time, then calculate polarization from induced from the
 % second pulse for each seperation.  Linearity makes this much more
 % efficient.  This should also be useable to 2D spec via rephasing pulses
 B_t_1 = cross(k1,E_t_1);
 B_t_2 = cross(k2,E_t_2);
 V_int = sym(zeros(N));
 for k = 1:N-1
 
    V_int(1,k+1) = -dot(mu(k,:),E_t_1+E_t_2)... %interactions
                   +dot(Hel(k+1,k+1)*cross(R_disp(k,:),mu(k,:))...
                    *1i/2 , B_t_1 + B_t_2);
                %only include mag dipole moment not quadrupole yet
        
 end
 V_int = V_int + V_int'; % just diagonal anyway


 %% Project into interaction picture basis


if ham_for_int_pic ==0 %include electronic
 U_tev = expm(-1i*Helsym*tt); %U_tev_full = kron(U_tev,eye(size(Hvib)));
elseif ham_for_int_pic == 1 %include viblvls as well as electronic
    Helvibshym = kron(Helsym,eye(size(Hvib)))+kron(eye(size(Hel)),Hvib);
      U_tev = expm(-1i*Helvibshym*tt);
 V_int = kron(V_int,eye(size(Hvib)));   %U_tev_full =U_tev;
elseif ham_for_int_pic == 2 %include viblvls as well as electronic and coupling
  Htotsym = kron(Helsym,eye(size(Hvib)))+kron(eye(size(Hel)),Hvib) + Hexvib;   
    U_tev = expm(-1i*Htotsym*tt);
 V_int = kron(V_int,eye(size(Hvib)));  %U_tev_full =U_tev;
else
     U_tev = expm(-1i*Helsym*tt); %U_tev_full = kron(U_tev,eye(size(Hvib)));
end
 V_int_IP = U_tev' * V_int * U_tev;
  
 %apply RWA by removing fast oscillating terms
 
 max_osc_freq = 200; %remove all frequencies which are faster than this
 V_int_IP2 = sym(zeros(size(V_int_IP)));
for k = 1:length(V_int_IP)
    for kk = k+1:length(V_int_IP) %no diag elements and hermitian
        
 coeffmatsym = express_exp_series(V_int_IP(k,kk),tt);
    tmp = sym(0); 

for k3 = 1:size(coeffmatsym ,2)
    if coeffmatsym{1,k3}~=0
        
        freq = coeffmatsym{2,k3};
        freqsub = subs(freq,[tt, omega_1, omega_2, E_1, E_2, V_12]...
                                ,[1 om1 om2 E1 E2 V ]);
    
    if abs(imag(freqsub)) < max_osc_freq 
        %in general there may be other resonances than this obvious one
        %which should be included, like those with vibrational levels

         %tmp = tmp + coeffmatsym{1,k3}*exp(freq);
         %use this to obtain an expression where w1 and w2 can be varied
         %etc but not sufficiently to change the RWA made
         
         tmp = subs(coeffmatsym{1,k3},...
                [Env_1, Env_2,omega_1, omega_2, E_1, E_2, V_12],...
                 [exp(-(tt-t01)^2/tau1^2),exp(-(tt-t02)^2/tau2^2),...
                 om1 om2 E1 E2 V ])*exp(freqsub*tt) + tmp;   
    end
    end
end
        V_int_IP2(k,kk) = tmp;
    end
end
V_int_IP2 = V_int_IP2 + V_int_IP2' ; %add CC
 
%V_int_full = kron(V_int_IP2,ones(size(Hvib)));


            
%% Construct Ltot, the total Louville for self coupling density matricies

% For N X N A and B we have the following
% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho


if ham_for_int_pic ==0 %include electronic
 Htot_m_H0 = Htot - kron(Hel,eye(size(Hvib))); 
elseif ham_for_int_pic == 1 %include viblvls as well as electronic
 Htot_m_H0 = Htot - kron(Hel,eye(size(Hvib))) - kron(eye(size(Hel)),Hvib);  
elseif ham_for_int_pic == 2 %include viblvls as well as electronic and coupling
  Htot_m_H0 = Htot - Htot; %all of it is gone!
else
 Htot_m_H0 = Htot - kron(Hel,eye(size(Hvib))); 
end

if isempty(Kap1)
    Kap1 = inf; %no truncation based on Kap1
end
    
if Kap2>0 && use_dru   %else trivial

numwithn = zeros(1,Kap2+1); %size(cc2,2) is number of poles in J(omega)
tot_poles = sum(cellfun(@length,cc_com));
for kk = 0:Kap2
numwithn(kk+1) = nchoosek(tot_poles-1+kk,kk);
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+size(cc2,2))) 
end

if length(Kap1)>1 %this is an overloading to allow you to just pass nn 
    %straight to the function if it has already been calculated before
   nn = Kap1; totfreq = horzcat(vv{:}).';
else
    
nn = zeros(sum(numwithn),tot_poles); %if Kap1 <inf then less elements may be present
totfreq = horzcat(vv{:}).';
count1=0; tmp3 = 0;
for k =2:Kap2+1

        count2 = count1 + size(tmp3,1); 
        
        tmp3 = zeros(length(count1+1:count2)*size(nn,2),size(nn,2));
       for j = count1+1:count2
        tmp = repmat(nn(j,:),size(nn,2),1); %select element from tier below
        tmp2 = tmp+eye(size(nn,2)); %add one to every position
        tmp3(1+(j-(count1+1))*size(nn,2):(j-count1)*size(nn,2) ,:) = tmp2;       
       end
       %remove any elements of tmp3 which fail to satisfy the high
       %frequency cut off condition
       lg = true(size(tmp3,1),1);
       if isfinite(Kap1)
            for j =1:length(count1+1:count2)*size(nn,2)
                lg(j) = tmp3(j,:)*real(totfreq) < Kap1;
            end
       end
       tmp3 = tmp3(lg,:);
       %now remove duplicate entries
       jj=1;
       while jj < size(tmp3,1)-1
           lg = any(tmp3(jj+1:end,:)  ~= repmat(tmp3(jj,:),size(tmp3,1)-jj,1),2);
           tmp3 = tmp3([true(jj,1);lg],:);
           jj=jj+1;
       end
       count1 = count2;
       nn(count1+1:count1 + size(tmp3,1),:) = tmp3;
end
end

%would be good to include a complicated truncation condition such that
%weakly coupled modes could be taken only to lower orders

%clean any nn values which are all zero besides the first
nn = nn([true;any(nn(2:end,:)~=0,2)],:);
tierofnn = sum(nn,2);
for k  =2:Kap2+1    
    numwithn(k) = sum(tierofnn == (k-1)); %set actual number included
end
else
   nn=0; numwithn=1; %trivial case, DO NOT ENTER Kap1 as nn and also Kap2 =0;
end
%% Set initial condition
if ~exist('rho_0','var') %no initial condition, have to create standard one
    if viblvls > 1
        rho_0 = zeros(size(Hel)); rho_0(1,1) = 1;
        vibpop = 1./(exp((1:viblvls).*beta*omegavib)-1); %BE stats
        vibpop = vibpop/sum(vibpop); %normalise, sum not sum squared
        rho_0 = kron(rho_0,diag(vibpop));
        %1 excitation in exiton basis (excited state) and random vib pops
        rho_0 = basis_proj*rho_0*basis_proj'; %excited state in chromophore basis
    else
        rho_0 = zeros(N);
        rho_0(1) = 1; %ground state
        %rho_0 = basis_proj*rho_0*basis_proj'; %should be same in chromophore basis
    end
end
rho_vec = zeros(size(nn,1)*numel(rho_0),1);
%initial condition
rho_vec(1:numel(rho_0)) = reshape(rho_0,numel(rho_0),1);
clear rhoharch
% Can reshape to zeros(size(rho_heir,1)*numel(rho),1); 
%Column vector of ever possible element of every order of density matrix.
% Ordered as you would get if you called the density matrix by one number
% e.g. if rho is an N X N matrix, rho(n) is the [n-1,mod(N)]+1th row and 
%the floor((n-1)/N)+1 th column.

%% Calculate the coupling from the heirarchy

%general up and down coupling terms
coup_com = zeros(size(nn,1));  %coefficient will be the sqrt of number it couples to (upward)
%coup_p1_save= zeros(size(nn,1),size(nn,1),N);  %takes up too much memory

coup_com_save{N}= sparse(coup_com);

%must also have the term that accounts for the fact that down coupling
%terms are different except for the matsubara terms
%coup_m1 = coup_com ;
coup_m1_acom = coup_com ;  %anti commutator part
coup_acom_save= coup_com_save;
totpoles = sum(cellfun(@length,cc_com));
tierofnn = sum(nn,2);
totfreq = horzcat(vv{:}).';
rng_j =0; %full_rng = 1:totpoles; 

for j = 1:N
    if ~isempty(cc_com{j}) %no bath at this mode for some reason
    rng_j = rng_j(end)+1:length(cc_com{j})+rng_j(end); 
    rng_rest = [1:rng_j(1)-1 , rng_j(end)+1:totpoles];
    cc = cc_com{j};  ccI = cc_acom{j};
for k =1:sum(numwithn(1:end-1))
    
    currenttier = tierofnn(k);
    tierabove = currenttier+1;

    tierpickr = abs(tierofnn-tierabove)<eps(10);

    nnsec = nn(k,rng_j); %section of interest
    nnsec2 = nn(k,rng_rest); %rest
    %which obviously must still match, but must be element of j section
    %changed

    temp0 = repmat(nnsec,numwithn(tierabove+1),1); 
    %make same size as nn in tier above in dimension 1
    temp02 = repmat(nnsec2,numwithn(tierabove+1),1); 
    
    temp = temp0 - nn(tierpickr,rng_j); 
    %take away elements in the tier above

    temp4 = temp02 - nn(tierpickr,rng_rest);
    
    temp2 = sum(temp~=0,2) ==1; %only one element +/- 1 diff
    
    temp3 = sum(temp,2) < -1/2 & sum(temp,2) > -3/2; %ones with a -1 net 
    temp4 = all(temp4==0,2); %match rest of the terms
    
    comb_lg = temp2 & temp3 & temp4; %full condition required
    tierpickr(tierpickr) = comb_lg ; 
    %use to put elements in the right row position in coupling matrix,
    %comb_lg is displaced by the number of tiers below
    
    if any(comb_lg)
        
    [~,temp44] = find(temp(comb_lg,:)); %find row position, i.e. which coefficient
    temp7 = sqrt(sqrt(abs(cc(temp44)).^2+abs(ccI(temp44)).^2));
    coup_com(k,tierpickr)=temp7.*sqrt(1+nnsec(temp44));
    
    %elements with minus are in the positions which are the transpose of this 
    %due to the rescaling the coefficients should be the same up and down
    %for the matsubara terms and the commutator part of the drude terms
    %now second down coupling operator
    
    coup_com(tierpickr,k) = cc(temp44).*sqrt(1+nnsec(temp44))./temp7;%...

    coup_m1_acom(tierpickr,k) = ccI(temp44).*sqrt(1+nnsec(temp44))./temp7;%...

    end

end
    %Note I am scaling things w.r.t. (abs(a)^2+abs(b)^2)^(1/4)
    %for a term with operator a V^O + b V^X

coup_com_save{j} = sparse(coup_com); coup_com = 0*coup_com;
coup_acom_save{j}= sparse(coup_m1_acom); coup_m1_acom=coup_com;
     
    end 
end
%down coupling from real frequency and residue modes is just  coup_p1.'


% calculate the sum_{j=1}^N sum_{k=0}^Kappa n_{jk} v_{ik}
if size(nn,1)==1
  const_factor = 0; 
else
    %reshape(vv.',numel(vv),1)
    %([vv(1,:).';vv(2,:).';vv(3,:).';vv(4,:).'])
    
const_factor = nn*totfreq;%([vv(1,:).';vv(2,:).']);
non_hermit = abs(imag(const_factor))>eps(max(imag(const_factor)));
%const_factor = const_factor -imag(sum(const_factor))/length(const_factor);
end

%% Pass function that calculates stuff all the things required

    temp = kron(const_factor,ones(numel(Htot),1));
	decay_term = - sparse(1:length(temp),1:length(temp),temp);
    clear temp
    szmat = length(rho_vec)/numel(Htot);
    spmat = sparse(1:szmat,1:szmat,ones(szmat,1));

  U_tev_2 = subs(U_tev,[E_1,E_2,V_12],[E1,E2,V]);  

topass{1} = {Htot_m_H0,V_int_IP2, QQ,U_tev_2,viblvls}; %operators
topass{2} = {coup_com_save,coup_acom_save,spmat,decay_term};  %Heirarchy stuff 

HEOM_light_DE(1,topass) ;


options = odeset('OutputFcn',@output_fun_ODE);

 output_fun_ODE(numpoints,rho_vec(1:(sum(numwithn(1:(saveuptotier+1)))...
           *numel(rho_0))),save_file_name);

%clearvars -except tend basis_proj convfact rho_vec nn
clear topass 
%save some memory
 

%% Propogate in time and solve system
nntot = nn; %create file and save this stuff
save(save_file_name,'nntot','convfact') 
clear nntot

[time_units,rho_vec] = ode45(@HEOM_light_DE,[0 tend],rho_vec,options);

time_units = time_units/convfact; %convert to ps from cm^-1


