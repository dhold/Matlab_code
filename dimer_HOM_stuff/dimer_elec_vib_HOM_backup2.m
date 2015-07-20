function [Time_units,rho_vec,nn,basis_proj,coup_p1,coup_m1,coup_m1_acom]=...
    dimer_elec_vib_HOM_backup2(E1,V,omegavib,viblvls,...
            coupg,tendps,gamma_dru,lambda_dru,Kappa,Kap2,Temp,rho_0)
%This code is used to model a dimer system of two chromophores with the
%Heirachial orders method.
% Unit setup
%All terms in inverse cm, conversion to sensible units achieved by
%E = hc / lamba => 6.62606957 *10^{-33} 299792458 *10 * "E in code" gives
%
tic
if nargin == 0  %use a select set nargin
%%
 E1 = 1042; %Energy of initial state (compared to excitation in other chromopore)
 V = 92; %Coupling
 omegavib = 1111; %Frequency of vibrations / level spacing,
% %assumed to be the same in each protein, choose (omega_1+omega_2)/2
 viblvls = 3; %Number of vibrational levels to include, integer
 coupg = 267.1; %coupling term between electronic and vibrational DOF
 tendps = 0.5; %end time in pico seconds
 gamma_dru = [100;100]; %drude decay constant
 lambda_dru = [6;6]; % weighting of distrubtion
 % %J_j(omega) = 2 lamda_j gamma_j /hbar * (omega/(omega^2+lambda_j^2)
 Kappa = 1; %truncation parameter satisfying Kappa >>omega_0 beta hbar/2pi
 %beyond this e^(-v_(k>kappa) t)*v_k ~delta(t) is made
 Kap2 = 3; %second truncation parameter, there are Kappa + 1 explicitly 
 %treated exponentials, this number is the maximum order the sum of the
 %number of powers can take, e.g. n_1 + n_2 +... + n_Kappa <= Kap2
 Temp = 300; %temperature in Kelvin
 include_truncation_correction = true; %for debugging
end
N=2; %number of sites, set to two for dimer code.
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}
% to run the code not as a function just uncomment and use each module
% Energy in Joules,  extra factor of 10 is because thats speed in cm/s
convfact = 2*pi * light_speed *length_unit  *10^(-12); %2 pi * "c" in cm/s * 1ps 
tend =convfact*tendps; %end time in units you get by putting everything in wavenumbers

beta = (2* pi *hbar * light_speed * length_unit)/ ( Temp * boltz_const);

%% Hamiltonian terms considered explicitly
%Excitonic Hamilonian
Hel = zeros(2); 

Hel(1,1) = E1; Hel(2,2) = 0;  %On diagonal terms
Hel(1,2) = V;  Hel(2,1) = V;  % Coupling terms from dipole dipole interaction

%Vibrational Hamilonian, relative degrees of freedom only
if viblvls > 1 %zero/1 means no quantised vibration considered
tmp = ones(viblvls,1);
tmp(1) = 0;
for k=2:viblvls
    tmp(k) = tmp(k-1) + omegavib;
end
Hvib = diag(tmp); clear tmp 
%Note this is H_vib = omega_vib * (b_r^dag b_r) not the COM excitation
notquitezero = kron(eye(size(Hel)),diag(E1*(-eps*(viblvls-1):2*eps:eps*(viblvls-1))));
%weighted in order to just keep the vibrational exited states ahead, this
%is quite a lazy solution I could just reorder afterwards
[basis_proj,~] = eig( kron(Hel,eye(viblvls))+notquitezero ); 
%projector to delocalised basis

Htot = kron(Hel,eye(viblvls))+kron(eye(2),Hvib); %total Hamiltonian
%which will be explicitly considered QM

%Hamilonian coupling electronic and vibrational levels
%Hexvib = -g/sqrt(2)  * (n_1 - n_2)*(b^dag_rel + b_rel)

Hexvib = diag(-coupg/sqrt(2) * sqrt(1:viblvls-1),1);
Hexvib = Hexvib+Hexvib.';
Hexvib = kron([-1,0;0,1],Hexvib);

Htot = Htot + Hexvib;  %Full Hamiltonian with electric + vib + coupling
%in chromophore basis
%Htot = basis_proj'*Htot*basis_proj; %project to exciton basis
else
    [basis_proj,~] = eig(Hel); 
    Htot = Hel;
   %rho Htot = basis_proj'*Hel*basis_proj;
end

%coefficients used in the HOM sum
cc = zeros(N,Kappa+1); vv = cc;
vv(:,1) = gamma_dru;  %Drude decay constant

vv(:,2:end) = (2*pi/beta)*kron(ones(N,1),1:Kappa); %These are the Matsuraba frequencies

cc(:,1) = gamma_dru.*lambda_dru.*(cot(beta.*gamma_dru/2)-1i);

cc(:,2:end) = (4/beta) * repmat(gamma_dru.*lambda_dru,1,Kappa)...
                .*   (         vv(:,2:Kappa+1)./ ...
                (vv(:,2:Kappa+1).^2 -repmat(gamma_dru.^2,1,Kappa))   );

%% Construct density matrix and higher orders 

rho = zeros(size(Htot));  % zeroth order in Heirarchy, ACTUAL density matrix

if Kap2>0 %else trivial

numwithn = zeros(1,Kap2+1);
for kk = 0:Kap2
numwithn(kk+1) = nchoosek((Kappa+1)*N-1+kk,kk); 
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+1)) 
end

nn = zeros(sum(numwithn),N*(Kappa+1)); 
count1=0;
for k =2:Kap2+1

        count2 = count1 + numwithn(k-1); 
        
        tmp3 = zeros(length(count1+1:count2)*N*(Kappa+1),N*(Kappa+1));
       for j = count1+1:count2
        tmp = repmat(nn(j,:),N*(Kappa+1),1); %select element from tier 
        tmp2 = tmp+eye(N*(Kappa+1)); %add one to ever position
        tmp3(1+(j-(count1+1))*N*(Kappa+1):(j-count1)*N*(Kappa+1) ,:) = tmp2  ;     
       end
       %remove duplicates
       jj=1;
       while jj < size(tmp3,1)-1
           lg = any(tmp3(jj+1:end,:)  ~= repmat(tmp3(jj,:),size(tmp3,1)-jj,1),2);
           lg = [true(jj,1);lg]; %#ok<AGROW>
           tmp3 = tmp3(lg,:);
           jj=jj+1;
       end
       count1 = count2;
       nn(count1+1:count1 + numwithn(k),:) = tmp3;

    
end

%% Set initial condition
if ~exist('rho_0','var') %no initial condition, have to create standard one
    if viblvls > 1

        vibpop = 1./(exp((1:viblvls).*beta*omegavib)-1); %BE stats
        vibpop = vibpop/sum(vibpop); %normalise, sum not sum squared
        rho_0 = full(sparse(1:length(rho),1:length(rho),[zeros(1,length(rho)/2),vibpop]));
        %1 excitation in exiton basis (excited state) and random vib pops
        rho_0 = basis_proj*rho_0*basis_proj'; %excited state in chromophore basis
    else
        rho_0 = zeros(N);
        rho_0(end) = 1; %excited state in exciton basis (as energy ordered)
        rho_0 = basis_proj*rho_0*basis_proj'; %excited state in chromophore basis
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

%% Calculate the

%general up coupling
coup_p1 = zeros(size(nn,1));  %coefficient will be the sqrt of number it couples to (upward)

%must also have the term that accounts for the fact that down coupling to 
% the k=0 states are different
coup_m1 = coup_p1 ; 
coup_m1_acom = coup_p1 ;  %anti commutator part
for k =1:sum(numwithn(1:end))
    
    temp0 = kron(ones(size(nn,1),1),nn(k,:));
    temp = temp0 - nn;

    temp2 = sum(temp~=0,2) ==1; %only one element +/- 1 diff
    temp3 = sum(temp,2) < -1/2 & sum(temp,2) > -3/2; %ones with a - net about

    [temp4,temp44] = find(temp(temp2 & temp3,:)); %find row position
    temp5 = 1+floor((temp44-1)/(Kappa+1)); %which chromophore
    temp6 = 1+mod(temp44-1,Kappa+1); %which order
    
    temp7 = sqrt(abs(diag(cc(temp5,temp6))));
    
    coup_p1(k,temp2 & temp3)=temp7.*sqrt(1+diag(temp0(temp4,temp44)));
    
    %now second down coupling operator

    temp3 = sum(temp,2) > 1/2 & sum(temp,2) < 3/2; %ones with a + net about

    [temp4,temp44] = find(temp(temp2 & temp3,:)); %find row position
    temp5 = 1+floor((temp44-1)/(Kappa+1)); %which chromophore
    temp6 = 1+mod(temp44-1,Kappa+1); %which order
    
    temp7 = real(diag(cc(temp5,temp6))) ./ sqrt(abs(diag(cc(temp5,temp6))));
    temp8 = 1i*imag(diag(cc(temp5,temp6))) ./ sqrt(abs(diag(cc(temp5,temp6))));
    
    coup_m1(k,temp2 & temp3)= temp7.*sqrt(diag(temp0(temp4,temp44)));
    coup_m1_acom(k,temp2 & temp3) = temp8.*sqrt(diag(temp0(temp4,temp44)));
    
    %elements with minus are the transpose of this 
    
end

%standard down coupling is just  coup_p1.'
else
   nn=0;  coup_p1=0; coup_m1 =0; coup_m1_acom =0; %trivial case
end



%% Construct Ltot, the total Louville for self coupling density matricies

% For N X N A and B we have the following
% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho
L = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));

Q = L*0;

%************ this BIT is playing up
if include_truncation_correction
for j = 1:N
    
    Qj = zeros(N);
    Qj(j,j) = 1;  Qj = kron(Qj,eye(viblvls)); eyeQ = eye(length(Qj));
      
    Q = Q + (lambda_dru(j).*(2./(beta*gamma_dru(j))-cot(beta*gamma_dru(j)/2))...
        - dot(cc(j,2:end),1./vv(j,2:end)))...
        .*(kron(Qj.',eyeQ) + kron(eyeQ,Qj) - 2*kron(Qj.',Qj)); 
    %reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
    %[Q_i,[Q_i,rho_n]] = (1 X Q_i + Q_i^T X 1 -2 Q_i^T X Q_i ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
end
end

Ltot = (L-Q);

%also calculate the sum_{j=1}^N sum_{k=0}^Kappa n_{jk} v_{ik}

const_factor = nn*([vv(1,:).';vv(2,:).']);



%% Propogate in time and solve system
topass{1} = Ltot; topass{2} = coup_p1; topass{3} = coup_m1;
topass{4} = coup_m1_acom; topass{5} =const_factor;
topass{6} = N; topass{7} = viblvls;
rhs_dif_eq(1,topass); %pass parameters to function

clear Ltot Lp1 Lm1 const_factor topass

%no point keeping two copies really, save memory etc

%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
%[T,Y] = ode45(@rhs_dif_eq,tspan,initcond,options);
%split over several time intervals so as not to fill memory in general
toc

[Time_units,rho_vec] = ode45(@rhs_dif_eq,[0 tend],rho_vec);
Time_units = Time_units/convfact; %convert to ps from cm^-1
end
%%
function drho = rhs_dif_eq(t,rho_vc) %Couples to same numwithn value

persistent L0 Qj Uj coup_com coup_d_acom N CF first_time select_ele choose_ele

        if isempty(t) %pass empty t to get persistent vars back
            drho{1} = L0; drho{2}= Qj; drho{3}=Uj; drho{4}=CF; drho{5}= select_ele;
            drho{2} = choose_ele;
            return
        end
        if iscell(rho_vc)
            
            %clear L0 Qj Uj coup_com coup_d_acom N CF setselect select_ele choose_ele
            
            L0 = rho_vc{1}; %basic Louvillian;
            coup_u = rho_vc{2}; %upward coupling
            coup_d = rho_vc{3}; %downward coupling 
            coup_com = coup_u + coup_d; %total coupling
            %w/ correction for zero order
            coup_d_acom = rho_vc{4}; %correction to downward coupling
            
            CF = rho_vc{5}; %constant factor
            
            N = rho_vc{6};
            vibs = rho_vc{7};
            %matrix coupling stuff
            drho =[];
            
        for j = 1:N
    
    Qjsec = zeros(N);  Qjsec(j,j) = 1;  Qjsec = kron(Qjsec,eye(vibs)); 
    eyeQ = eye(length(Qjsec)); 
    Uj{j} = -1i*sparse(kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ)); %anti commutator
    Qj{j} = -1i*sparse(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    
        end
        first_time = true;

            return           
        end
 
        drho = 0*rho_vc;
        oneunit = length(L0);

        %reason for seperation is possible parallelisation after splitting 
        % rho_vec into "number of cores" chunks and running in parallel
     % tic  

     krng = 1:(size(rho_vc,1)/oneunit);
    trueunit = true(1,oneunit);
    
if first_time      %Construct massive (square) sparse matrix representing 
    %all the operator possibilities

         for k = krng
             c_com = coup_com(k,:); c_anticom = coup_d_acom(k,:);
             select_ele{k} =  (kron(trueunit,c_com~=0))~=0;
             choose_ele{k} =  (kron(trueunit,c_anticom~=0))~=0;
         end
         first_time = false;
end
%tic

for k = krng

        %basic louvillian term with coherence decay terms 
        temp = (L0-eye(oneunit)*CF(k))*rho_vc(1+oneunit*(k-1):oneunit*k);
        %temp = (L0)*rho_vc(1+oneunit*(k-1):oneunit*k); 
        
        %Now the terms of the form 
        % sum_{j=1}^N (-1i[Q_j,sum(c_{vec(n)} p_{vec(n)})] + ...
        %               - { Q_j,sum(cc_{vec(n)} p_{vec(n)}) } )
        
        c_com = coup_com(k,:); c_anticom = coup_d_acom(k,:);
        %tic  %moves this as persistent variable as it wasn't quick
%         select_ele =  (kron(trueunit,c_com~=0))~=0;
%         choose_ele =  (kron(trueunit,c_anticom~=0))~=0;
        %toc
        %tic
        for j =1:N
            if any(c_com~=0 |  c_anticom~=0)

           temp = temp + Qj{j}*(c_com(c_com~=0)*reshape...
               (rho_vc(select_ele{k}),oneunit,sum(c_com~=0)).').';
            %this is just a way of doing the sum efficiently (I think)
           temp = temp + Uj{j}*(c_anticom(c_anticom~=0)*...
               reshape(rho_vc(choose_ele{k}),oneunit,sum(c_anticom~=0))).';
           
           %might be an idea just to save c_com(c_com~=0) etc as array
           % and save performing a logic operation and sum every time
            end
        end
        %toc
        drho(1+oneunit*(k-1):oneunit*k) = temp;
        
end
   % toc
end
