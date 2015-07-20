function [Time_units,rho_vec,nn,basis_proj]=dimer_elec_vib_HOM(E1,V,omegavib,viblvls,...
            coupg,tendps,gamma_dru,lambda_dru,Kappa,Kap2,Temp,rho_0)
%This code is used to model a dimer system of two chromophores with the
%Heirachial orders method.
% Unit setup
%All terms in inverse cm, conversion to sensible units achieved by
%E = hc / lamba => 6.62606957 *10^{-33} 299792458 *10 * "E in code" gives
%
if nargin == 0  %use a select set nargin
 E1 = 1042; %Energy of initial state (compared to excitation in other chromopore)
 V = 92; %Coupling
 omegavib = 1111; %Frequency of vibrations / level spacing,
% %assumed to be the same in each protein, choose (omega_1+omega_2)/2
 viblvls = 3; %Number of vibrational levels to include, integer
 coupg = 267.1; %coupling term between electronic and vibrational DOF
 tendps = 2; %end time in pico seconds
 gamma_dru = [100;100]; %drude decay constant
 lambda_dru = [10;10]; % weighting of distrubtion
 % %J_j(omega) = 2 lamda_j gamma_j /hbar * (omega/(omega^2+lambda_j^2)
 Kappa = 2; %truncation parameter satisfying Kappa >>omega_0 beta hbar/2pi
 %beyond this e^(-v_(k>kappa) t)*v_k ~delta(t) is made
 Kap2 = 1; %second truncation parameter, there are Kappa + 1 explicitly 
 %treated exponentials, this number is the maximum order the sum of the
 %number of powers can take, e.g. n_1 + n_2 +... + n_Kappa <= Kap2
 Temp = 300; %temperature in Kelvin
end
N=2; %number of sites, set to two for dimer code.

% to run the code not as a function just uncomment and use each module
% Energy in Joules,  extra factor of 10 is because thats speed in cm/s
convfact = 2*pi * 299792458*100  *10^(-12); %2 pi * "c" in cm/s * 1ps 
tend =convfact*tendps; %end time in units you get by putting everything in wavenumbers

beta = (6.62606957 * 10^(-34) * 299792458 *100)/ (Temp * 1.3806488 * 10^(-23));

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
   %rho Htot = basis_proj'*Hel*basis_proj;
end
%% Construct density matrix and higher orders 

rho = zeros(size(Htot));  % zeroth order in Heirarchy, ACTUAL density matrix

if Kap2>0 && Kappa>0 %else trivial
numwithn = zeros(1,Kap2+1);
for kk = 0:Kap2
numwithn(kk+1) = nchoosek((Kappa+1)*N-1+kk,kk); 
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+1)) 
end
nn = zeros(sum(numwithn),N*(Kappa+1)); 
intlist = intpartgen2(Kap2,N*(Kappa+1)); count1=1; 

for k =2:Kap2+1
   
    tmp = intlist{k};
    tmp2 = zeros(size(tmp,1),N*(Kappa+1));
    tmp2(1:size(tmp,1),1:size(tmp,2)) = tmp;
    tmp3 = perms2(tmp2,'unique');
    
    count2 = count1 + size(tmp3,1);
    
    nn(count1+1:count2,:) = tmp3;
    
    count1=count2;
    
end

%Calculate the coupling to density matricies with order one higher/lower
%Elements with no changes are simply an identity mixing
coup_p1 = false(count2); 
for k =1:count2
    
    temp = kron(ones(count2,1),nn(k,:)) - nn;
    temp2 = sum(logical(temp),2) ==1; %only one element +/- 1 diff
    temp3 = sum(temp,2) < -1/2 & sum(temp,2) > -3/2; %ones with a - net about
    coup_p1(:,k) = temp2 & temp3;
    %elements with minus are the transpose of this 
end
else
   nn=1;  coup_p1=false; %trivial case
end

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

%% Construct Ltot, the total Louville for self coupling density matricies

%  L = kron(Htot,eye(length(Htot)))-kron(eye(length(Htot)),Htot.');
L = kron(Htot.',eye(length(Htot)))-kron(eye(length(Htot)).',Htot);

%coefficients used in the sum
cc = zeros(N,Kappa+1); vv = cc;
vv(:,1) = gamma_dru;  vv(:,2:end) = 2*pi*kron(1:Kappa,ones(Kappa,1))/beta;
cc(:,1) = gamma_dru.*lambda_dru.*(cot(beta.*gamma_dru/2)-1i);

cc(:,2:end) = (4/beta)*kron(gamma_dru.*lambda_dru,ones(1,Kappa))...
     .*(vv(:,2:Kappa+1)./(vv(:,2:Kappa+1).^2 -kron(gamma_dru,ones(1,Kappa))));

Q = L*0;

%************ this BIT is playing up
%if 1==0 %uncomment to add 
for j = 1:N
    
    Qj = zeros(N);
    Qj(j,j) = 1;  Qj = kron(Qj,eye(viblvls));
%     Q = Q + (lambda_dru(j).*(2./(beta*gamma_dru(j))-1i) ...
%         - dot(1./cc(j,:),vv(j,:))).*(kron(Qj,eye(length(Qj)).') + ...
%         kron(eye(length(Qj)),Qj.') - 2*kron(Qj,Qj.'));       
    Q = Q + (lambda_dru(j).*(2./(beta*gamma_dru(j))-1i) ...
        - dot(cc(j,:),1./vv(j,:))).*(kron(Qj.',eye(length(Qj))) + ...
        kron(eye(length(Qj)).',Qj) - 2*kron(Qj.',Qj)); 
    %[Q_i,[Q_i,rho_n]] = (Q_i X 1^T + 1 X Q_i^T -2 Q_i X Q_i^T ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
end
%end
Ltot = -1i.*(L+Q);

%also calculate the sum_{j=1}^N sum_{k=0}^Kappa n_{jk} v_{ik}

const_factor = nn*([vv(1,:).';vv(2,:).']);
if 1==0
%This part was wrong
%% Construct Ltotplus, the Louville coupling to those with coup_p1


Lp1 = {}; %cell array with sparse matrix at the coupled locations
Lp1{sum(double(any(coup_p1)))} = []; 


for kk = 1:size(nn,1) %loop over all
    if any(coup_p1(:,kk))

            Q = 0*Q;
            for j = 1:N
    
    Qj = zeros(N);
    Qj(j,j) = 1;  Qj = kron(Qj,eye(viblvls));
%     Q = Q + sqrt(dot(nn(kk,1+(j-1)*(Kappa+1):j*(Kappa+1))+1,abs(cc(j,:))))...
%         .*(kron(Qj,eye(length(Qj)).') - kron(eye(length(Qj)),Qj.'));   
    Q = Q + sqrt(dot(nn(kk,1+(j-1)*(Kappa+1):j*(Kappa+1))+1,abs(cc(j,:))))...
        .*(kron(Qj.',eye(length(Qj))) - kron(eye(length(Qj)).',Qj));            

    %[Q_i,rho_n] = (Q_i X 1^T + 1 X Q_i^T) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
            end
            Lp1{kk}= -1i*sparse(Q);
                  
    end
end



%% Construct Ltotminus, the Louville down coupling 


Lm1 = {}; %cell array with sparse matrix at the coupled locations
coup_m1 = coup_p1.';%Down couplings
Lm1{sum(double(any(coup_m1)))}  = []; 

for kk = 1:size(nn,1) %loop over all

    if any(coup_m1(:,kk))

            Q = 0*Q;
            for j = 1:N
    
    Qj = zeros(N);
    Qj(j,j) = 1;  Qj = kron(Qj,eye(viblvls));
%     Q = Q + sqrt(dot(nn(kk,1+(j-1)*(Kappa+1):j*(Kappa+1)),1./abs(cc(j,:))))...
%             .*(sum(cc(j,:))*kron(Qj,eye(length(Qj)).') ...
%             - sum(conj(cc(j,:)))*kron(eye(length(Qj)),Qj.')); 
    Q = Q + sqrt(dot(nn(kk,1+(j-1)*(Kappa+1):j*(Kappa+1)),1./abs(cc(j,:))))...
            .*(sum(cc(j,:))*kron(Qj.',eye(length(Qj))) ...
            - sum(conj(cc(j,:)))*kron(eye(length(Qj)).',Qj));  
    %[Q_i,rho_n] = (Q_i X 1^T + 1 X Q_i^T) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
            end
            Lm1{kk}= -1i*sparse(Q);
         
        
    end
end
else %do it correctly
   
    Lp1 = zeros(N,size(rho_vec,1)); %coefficient for each coupling mat
    Lm1 = zeros(N,size(rho_vec,1)); %coefficient for each coupling mat
    coup_m1 = coup_p1.';%Down couplings
    Qjarraydown ={}; %array with all the matricies which will be used
    Qjarrayup ={}; 
    for j=1:N
       Qj = zeros(N);
       Qj(j,j) = 1;  Qj = kron(Qj,eye(viblvls));
       Qjarraydown{j} = kron(Qj.',eye(length(Qj))) - kron(eye(length(Qj)).',Qj);
       for k = 1:viblvls
          
           Qjarrayup{j,k} = (cc(j,k)*kron(Qj.',eye(length(Qj))) ...
            - conj(cc(j,k))*kron(eye(length(Qj)).',Qj));
           
       end
    end
    for j = 1:size(rho_vec,1)
        for k = 1:size(rho_vec,1)
            if coup_p1(j,k)
                
                tmp = nn(j,:) - nn(k,:); %should be one diff with val +1
                tmp2 = find(tmp);
                jj = floor(tmp2/viblvls); %which chromophore
                Lp1(jj,j) = sqrt(nn(j,tmp2)+1)*sqrt(abs(cc(jj,j)));
                
            elseif coup_m1(j,k) 

                tmp = nn(j,:) - nn(k,:); %should be one diff with val -1
                tmp2 = find(tmp);
                jj = floor(tmp2/viblvls); %which chromophore
                Lm1(jj,j) = sqrt(nn(j,tmp2))/sqrt(abs(cc(jj,j)));
                
            end

        end
    end

    
end
%% Propogate in time and solve system
topass{1} = Ltot; topass{2} = Lp1; topass{3} = Lm1;
topass{4} = const_factor; topass{5} = coup_p1;
rhs_dif_eq(1,topass); %pass parameters to function

clear Ltot Lp1 Lm1 const_factor coup_p1 topass
%no point keeping two copies really, save memory etc

%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
%[T,Y] = ode45(@rhs_dif_eq,tspan,initcond,options);
%split over several time intervals so as not to fill memory in general
[Time_units,rho_vec] = ode45(@rhs_dif_eq,[0 tend],rho_vec);
Time_units = Time_units/convfact; %convert to ps from cm^-1
end
%%
function drho = rhs_dif_eq(t,rho_vc) %Couples to same numwithn value

persistent L0 L1 LM1 CF lg

        if isempty(t) %pass empty t to get persistent vars back
            drho{1} = L0; drho{2}= L1; drho{3}=LM1; drho{4}=CF; drho{5}= lg;
            return
        end
        if iscell(rho_vc)
            L0 = rho_vc{1};
            L1 = rho_vc{2}; 
            LM1 = rho_vc{3}; 
            CF = rho_vc{4};
            lg = rho_vc{5};
            drho =[];
            return           
        end

        drho = 0*rho_vc;
        oneunit = length(L0);
        lg2 = lg';
        %reason for seperation is possible parallelisation after splitting 
        % rho_vec into "number of cores" chunks and running in parallel
     % tic  
for k = 1:(size(rho_vc,1)/oneunit)

        %temp = (L0-sparse(1:oneunit,1:oneunit,CF(k)))*rho_vc(1+oneunit*(k-1):oneunit*k);   
        temp = (L0-eye(oneunit)*CF(k))*rho_vc(1+oneunit*(k-1):oneunit*k);  
        if any(lg(:,k))
        for kk = find(lg(:,k))

            temp = temp + L1{k}*rho_vc(1+oneunit*(kk-1):oneunit*kk);
            
        end
        end
        if any(lg2(:,k))
        for kk = find(lg2(:,k))

            temp = temp + LM1{k}*rho_vc(1+oneunit*(kk-1):oneunit*kk);
            
        end
        end
        drho(1+oneunit*(k-1):oneunit*k) = temp;
        
end
   % toc
end


% function out = rhs_dif_eq_0(rho_h_arch,tosave) %Couples to same numwithn value
% 
% persistent Ltot vv
% 
%         if nargin ==2
%             Ltot = tosave{1};
%             vv = tosave{2}; %vv(1) = 
%             out = 0;
%             return           
%         end
% 
%         out = 0*rho_h_arch;
% 
% for k = 0:size(rho_heir,1)-1
% 
%         out(1+length(rho)*k:length(rho)*(k+1),:) = Ltot*rho;        
%         
% 
% end
% 
% end
%  
% function out = rhs_dif_eq_1(rho_h_arch) %Couples to +-1 numwithn value
% 
% 
% end

%create larger matrix of all orders of density matricies old method
% rho_heir = zeros((Kappa+1)^Kap2,Kappa+1); 
% 
% intlist = intpartgen2(Kap2,Kappa+1); count1=1; 
% numwithn = [1,zeros(1,Kap2)]; %number of density operators with each order
% for k =2:length(intlist)
%    
%     tmp = intlist{k};
%     tmp2 = zeros(size(tmp,1),Kappa+1);
%     tmp2(1:size(tmp,1),1:size(tmp,2)) = tmp;
%     tmp3 = perms2(tmp2,'unique');
%     
%     numwithn(k) = size(tmp3,1);
%     
%     count2 = count1 + size(tmp3,1);
%     
%     rho_heir(count1+1:count2,:) = tmp3;
%     
%     count1=count2;
%     
% end
% 
% rho_heir = rho_heir(1:count2,:); %trim extra crap
% nn = zeros((Kappa+1)^Kap2,2*(Kappa+1)) ; %vector index   
% 
% %all included within truncation
% count3 = 0;
% numwithn2 = [0,cumsum(numwithn)];
% for k = 1:Kap2
% 
%     
%     tmp = kron(ones(numwithn2(Kap2+2-k),1),rho_heir(numwithn2(k)+1:numwithn2(k+1),:));
%     tmp2 = kron(rho_heir(1:numwithn2(Kap2-k+2),:),ones(numwithn(k),1));
%     nn(count3+1:count3+size(tmp,1),:) = [tmp,tmp2];
%     count3 = count3 + size(tmp,1);
%  
% end
% clear rho_heir   
% %sort by total order (=sum(nn))
% 
% tmp = sum(nn,2);
% 
% [tmp,tmp2] = sort(tmp);
% 
% nn = nn(tmp2,:);
% 