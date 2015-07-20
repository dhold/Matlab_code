%% Construnct N chromophore Hamiltonian with explicitly considered vibrational modes

%  Hel = sum_{k=1}^N eta_k simga^+_k simga^-_k 
%        +  sum_{k \ne j} V_{kj} (simga^+_k simga^-_j  + c.c.)
%Hel = [ E1 , V12, V13,..., V1N; V21, E2,V23 ,...;....; VN1,...,EN] ; 
% The excitonic Hamilonian with levels and dipole couplings
% It is assumed only one excitation is ever present and so this is N level
% It is possible to assume WLOG that EN or E1 is zero by rescaling
% By Farrrrrrrrrr the most sensible choice is E1 = 0 but ordering is
% reversed in the literature for much the same reason cm^-1 is the unit...
%

Hel= [17113.9, 319.374, 30.4857, 7.5811;319.374, 17033.3, 20.3238, -43.8736;...
 30.4857, 20.3238, 15807.4, 3.3873;7.5811, -43.8736, 3.3873, 16371.9];
Hvib = 1;
%the labels for the sites are: {"DBVd", "DBVc", "PCBc158", "MBVb"}
%Hel = [1023,45,12 ; 0,800,23 ; 0,0,0]; 
N = length(Hel);
%Hel = Hel-diag(diag(Hel)) + Hel'; %I haven't written in the lower diag eles
%eventually I will probably just save these and load whenever...

%the next point to consider is which viblvls to treat QM, construction of a
%vibrational Hamiltonian will reflect this
% Hvib = sum_{k=1}^n omega_k b_k^{dagger} b_k
% we included only "viblvls" vibrational states

% Noteable modes in PE545 in cm^-1 as always..
% sQM 0.0013 0.0072 0.045 0.0578 0.045
% ωQM 207     244     312    372   438
% sQM 0.0924 0.0761 0.0578 0.0313 0.0578
% ωQM   514   718     813    938    1111
% sQM 0.1013 0.0265 0.0072 0.0113
% ωQM  1450    1520   1790    2090

viblvls = 0; %zero means no quantised vibration considered
%Could make this mode specfic, e.g. a vector
omegavib = [];%1108;%[938,1111,1450];  
%omegavib should contain only high energy modes treated as part of the ham
%g = ω*sqrt(s)
%Note that having multiple modes will slow the code down vastly
%sqm = [0.0313,0.0578,0.1013];  
% excitations expected to be important in assisting the transfer of excitions.
% low energy vibrations are occupied thermally and treated as part of the
% bath.
Ecut = inf; %cut off energy used as a truncation condition with multiple modes
om_0 = [1108]; %#ok<*NBRAK> %om_0 is the values of brownian modes included
lambda =44.32; %reorganisation energy
gamma = [5.3]; %damping of modes
sqm = lambda./omegavib;

lam_dru = 100; %amplitude of drude modes
gam_dru = 100; %width

nmodes = length(omegavib);

if viblvls > 1 && ~isempty(omegavib)
    %Not that best way of enumerating these states but I had it written
    %already so just adapted
    numwithn = zeros(1,viblvls);
for kk = 1:(viblvls-1)*nmodes+2
numwithn(kk) = nchoosek(nmodes-2+kk,kk-1); 
end
    occmat = zeros(sum(numwithn),nmodes); count1=0;
for k =2:(viblvls-1)*nmodes+2

        count2 = count1 + numwithn(k-1); 
        
        tmp3 = zeros(length(count1+1:count2)*nmodes,nmodes);
       for j = count1+1:count2
        tmp = repmat(occmat(j,:),nmodes,1); %select element from tier 
        tmp2 = tmp+eye(nmodes); %add one to ever position
        tmp3(1+(j-(count1+1))*nmodes:(j-count1)*nmodes,:) = tmp2  ;     
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
       occmat(count1+1:count1 + numwithn(k),:) = tmp3;
    
end
occmat = occmat(1:count2+numwithn(k),:);
tmp = max(occmat,2);
occmat = occmat(tmp<=viblvls,:);


%%
tmp = sum(occmat.*repmat(omegavib,size(occmat,1),1),2);
%Hexvib = sum_{j=1}^N simga^+_j simga^-_j *...
%       sum_{k =1}^N (g_{kj} b^{dag}_k  + c.c.)
% g = ω*sqrt(s), the values for modes in 
% Hexvib_cmp compontents of Hamilonian coupling electronic 
% and vibrational levels
coupg = sqrt(sqm).*omegavib;
%also equals omegavib*lambda
coupg = repmat(coupg,N,1); %assume coupling the same from each chromophore
%in general factor of exp(1i*dot(k,r_j)) would be included, but for low E
Hexvib_cmp = zeros(size(occmat,1));
for k=1:size(occmat,1) 
    occsec = repmat(occmat(k,:),size(occmat,1),1);
    Diff = occsec-occmat;
    for j= 1:nmodes %sum over each operator (g_{kj} b^{dag}_k  + c.c.)
        %here it is assumed that that g_{kj} = g_{k}
        lg = Diff(:,j)==1 & all(Diff(:,[1:j-1,j+1:end])==0,2);
        Hexvib_cmp(k,lg) = coupg(1,j) * sqrt(occmat(lg,j));
    end
end
Hexvib_cmp = Hexvib_cmp+Hexvib_cmp'; %add cc

%assuming this is the same for all chromophores or you would have to do
%each section individually

tmp = tmp(:);
% may be an idea to include only states that have a total vibrational
% energy below some threshold
[tmp2,tmp3] = sort(tmp);
Hexvib_cmp = Hexvib_cmp(tmp3,tmp3);
tmp = tmp2(tmp2<Ecut); %Ecut is a cut off energy specified before hand
tmp3 = tmp3(tmp2<Ecut); %make sure to be away of the ordering
Hexvib_cmp = Hexvib_cmp(tmp2<Ecut,tmp2<Ecut); %**check this is correct**

Hvib = diag(tmp); 
Hexvib = kron(eye(size(Hel)),Hexvib_cmp);

[tmp5,tmp4] = eig(Hel); 
basis_proj = kron(tmp5,eye(length(Hvib)));
%Hel_ext = kron(tmp,eye(length(Hvib)));

Htot = Hexvib + kron(Hel,eye(length(Hvib)))+kron(eye(N),Hvib); %total Hamiltonian
%which will be explicitly considered QM

occmat = occmat(tmp2<Ecut,:); occmat=occmat(tmp3,:);

else
    [basis_proj,~] = eig(Hel); 
    Htot = Hel;
   %rho Htot = basis_proj'*Hel*basis_proj;
end

%% Generate initial condition if required
Temp =300; %units Kelvin
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); 
beta = (2 * pi * hbar * light_speed * length_unit)/ ( Temp * boltz_const);
if viblvls > 0 && ~isempty(omegavib)
        vibpop = 1./(exp((diag(Hvib)+Hvib(2,2)).*beta)-1); %BE stats
        vibpop = vibpop/sum(vibpop); %normalise, sum not sum squared
        
        % Use the occupation of the diagonalised (electronically)
        % Hamiltonian, specificy the population c_n of the nth Excition 
        % state, note the first value is taken as the ground state occ
        
        if 1==0 %this method only works for totally statistical mixtures 
            %of the states, a more general one is below
        c_n = [0,0,1];
        temp = kron(c_n,vibpop);
        rho_0 = full(sparse(1:length(Htot),1:length(Htot),temp));
        %1 excitation in exiton basis (excited state) and random vib pops
        % so project back to chromophore basis
        else
        %c_n = [0,0,0;0,1/2,1/2;0,1/2,1/2];  %general superposition
        c_n = [0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,1];
        rho_0 = kron(c_n,diag(vibpop));
            
        end
else
    rho_0 = [0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,1];
end
        rho_0 = basis_proj*rho_0*basis_proj';
      
 %% set parameters

 Kappa = 0;
 
 [cc1,cc2R,cc2I,vv1,vv2,QQ] = ...
    coeffients_from_brownian(lambda,gamma,om_0,Temp,Kappa,lam_dru,gam_dru);
 
cc1 = repmat(cc1,N,1); cc2R = repmat(cc2R,N,1); cc2I = repmat(cc2I,N,1); 
 vv1 = repmat(vv1,N,1); vv2 = repmat(vv2,N,1);


 Kap1=inf; %max frequency truncation
 Kap2 =4; %max level truncation
 tendps =0.7; 
 convfact = 2 * pi * light_speed * length_unit * 10^(-12);
 numpoints =[200, [0.15,0.35,0.55,0.7]*convfact]; 
 saveonlyrho00 = false;
 use_reduced_mat = false;
 vibstates = length(Hvib);
 %%  Run the full calculation
 
 [Time_units,rho_vec,nn,total_prop_op]=multi_site_drude_and_brownian_HOM...
            (Htot,QQ,rho_0,cc1,cc2R,cc2I,vv1,vv2,Kap1,Kap2,vibstates,...
            numpoints,tendps,Temp,saveonlyrho00,use_reduced_mat,'saved_data.mat');
        

 
        