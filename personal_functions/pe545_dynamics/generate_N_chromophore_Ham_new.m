%  this version allows a different set of modes at each site
%% Construct N chromophore Hamiltonian with explicitly considered vibrational modes

% PC645
%Hel= [17113.9, 319.374, 30.4857, 7.5811;319.374, 17033.3, 20.3238, -43.8736;...
% 30.4857, 20.3238, 15807.4, 3.3873;7.5811, -43.8736, 3.3873, 16371.9];
%the labels for the sites are: {"DBVd", "DBVc", "PCBc158", "MBVb"} 
% PE545
fle = open('Hamiltonian_save.mat');
if 1==0 % PC645
    
    Hel = fle.PC645Hamiltonian;
    E0 = fle.E0_PC645; N = length(Hel);
       
lambdaD =100;omegaD=100;
omegaU = 1108.0;gammaU = 5.30884;lambdaU = 44.32;

 omegavib = {omegaU,omegaU,[],[]};
num_vib = cellfun(@length,omegavib);
sQM   =  {[0.0578],[0.0578],[],[]};
om_0 = {[],[],omegaU,omegaU};
%om_0 is the values of damped brownian modes included
lambda ={44.32,44.32,44.32,44.32}; %reorganisation energy
gamma = {5.30884,5.30884,5.30884,5.30884};   
    
 lam_dru = {lambdaD,lambdaD,lambdaD,lambdaD}; %reorganisation energy of drude modes
gam_dru = {omegaD,omegaD,omegaD,omegaD};   
    %coupg = sqrt(sQM{j}).*omegavib{j}; 
   viblvls = {[3],[3],[],[] };  
else % PE545

    Hel = fle.PE545Hamiltonian;
    E0 = fle.E0_PE545; N = length(Hel);
   % Hel(5:6,5:6) = [Hel(5,5),92;92,Hel(5,5)-1092]; %modified
    
% Noteable modes in PE545 in cm^-1 as always..
% sQM 0.0013 0.0072 0.045 0.0578 0.045
% wQM 207     244     312    372   438
% sQM 0.0924 0.0761 0.0578 0.0313 0.0578
% wQM   514   718     813    938    1111
% sQM 0.1013 0.0265 0.0072 0.0113
% wQM  1450    1520   1790    2090
%modes 5 and 6 are the most significant for the highest excitons
lambdaD = 100; omegaD =100;

%omegavib = {[],[],[],[],[438,1111,1450],[438,1111,1450],[],[]};
omegavib = {[],[],[],[],[],[],[],[]};
num_vib = cellfun(@length,omegavib);
sQM   =  {[],[],[],[],[0.045,0.0578,0.1013],[0.045,0.0578,0.1013],[],[]};
%om_0 = {[],[],[],[],1111,1111,[],[]}; 
om_0 = {[],[],[],[],[],[],[],[]}; 
%om_0 is the values of damped brownian modes included
lambda ={[],[],[],[],64,64,[],[]}; %reorganisation energy
gamma = {[],[],[],[],6,6,[],[]};   
    
 lam_dru = {lambdaD,lambdaD,lambdaD,lambdaD,lambdaD,lambdaD,lambdaD,lambdaD};
 %reorganisation energy of drude modes
gam_dru = {omegaD,omegaD,omegaD,omegaD,omegaD,omegaD,omegaD,omegaD};   
    %coupg = sqrt(sQM{j}).*omegavib{j}; 
   viblvls = {[],[],[],[],[1,1,1],[1,1,1],[],[] };  

end
flname = 'saved_data_pe545_no_mode.mat';
Ecut = inf; %cut off energy used as a truncation condition with multiple modes
nmodes = length(omegavib);
%%
if sum(num_vib)>0
    %%
    clear H_vib_sec H_ex_vib
     cnt=0;
    for j = 1: N       
        
        coupg = sqrt(sQM{j}).*omegavib{j}; 
        om_rng = omegavib{j};
        for jj = 1:length(om_rng )
    
                        cnt = cnt+1;
          H_ex_vib_sec{cnt} = sparse(diag(coupg(jj)*sqrt(1:viblvls{j}(jj)-1),1)...
                            +diag(coupg(jj)*sqrt(1:viblvls{j}(jj)-1),-1));
          H_vib_sec{cnt} = sparse(1:viblvls{j}(jj),1:viblvls{j}(jj),...
                            om_rng(jj)*(0:(viblvls{j}(jj)-1)));     
        end
    end
H_vib = sparse(1:prod(cellfun(@prod,viblvls)),1:prod(cellfun(@prod,viblvls)),...
                zeros(prod(cellfun(@prod,viblvls)),1)); 

for j1 = 1:sum(num_vib)
    tmp = sparse(1);
    for j2 = 1:sum(num_vib)
        if j2==j1
        tmp = kron(tmp,H_vib_sec{j1} );
        else
        tmp = kron(tmp,eye(size(H_vib_sec{j1})) ) ;
        end
    end
    H_vib = H_vib + tmp;
end
    Htot = kron(eye(size(Hel)),H_vib) +  kron(Hel,eye(size(H_vib)));
cnt = 0;
  for j = 1: N
    om_rng = omegavib{j}; tmp = sparse(Hel*0); tmp(j,j) = 1; 
        for jj = 1:length(om_rng )
            cnt=cnt+1; tmp2 = tmp;
            for j2 = 1:sum(num_vib)
                if j2==cnt  %construct operator via tensor products etc
                tmp2 = kron(tmp2,H_ex_vib_sec{j2});   
                else
                tmp2 = kron(tmp2,eye(size(H_ex_vib_sec{j2})));
                end
            end
      Htot = Htot  + tmp2;
        end
  end
  
  [tmp,~] = eig(Hel); 
basis_proj = kron(tmp,eye(length(H_vib)));
  
    clear H_ex_vib H_vib_sec
%% Trim highest energy states if required, complicates partial traces
if isfinite(Ecut) %cut states with high on diagonal energy
    tmp  = diag(Htot);
    lg = tmp<Ecut;
    Htot = Htot(lg,lg);

end
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
if ~all(cellfun(@isempty,omegavib))
    %assume each site is in thermal equilibrium
%     cnt = 1;
%     for j = 1: N       
%         om_rng = omegavib{j}; tmp=1; 
%         for j1 = 1:length(om_rng)
%             for j2 = 1:length(om_rng)
%                 if j2 == j1
%              tmp = kron(tmp,H_vib_sec{cnt+j1});
%                 else
%              tmp = kron(tmp,eye(size(H_vib_sec{cnt+j1})));   
%                 end
%             end
%         end
%                 cnt = cnt + length(om_rng);
%                 thermpop{j} = expm(-beta*tmp);
%                 thermpop{j} = thermpop{j}/trace(thermpop{j});
%     end
%         
        thermpop = expm(-beta*H_vib); 
        thermpop = thermpop/trace(thermpop);
        c_n = zeros(size(Hel));
        c_n(end) = 1; 
        %electronic initial state, not this is exciton basis stuff
        rho_0 = kron(c_n,thermpop); %tensor product with vib states
        %to get a reasonable thermal initial state.
        
else
            rho_0 = zeros(size(Hel));
        rho_0(end) = 1; 
        
end
        rho_0 = basis_proj*rho_0*basis_proj'; %project to site basis
      
 %% set parameters

 Kappa = 0;
 
 QQ = zeros(N,2); cnt=0;
 cc = zeros(1,sum(cellfun(@length,lambda))+ sum(cellfun(@length,lam_dru)) + Kappa);
 clear cc_com cc_acom vv
 for j=1:N
 [cc1,cc2R,cc2I,vv1,vv2,QQ(j,1),tmp] = ...
    coeffients_from_brownian_new(lambda{j},gamma{j},om_0{j},Temp,Kappa,lam_dru{j},gam_dru{j});
 % QQ is in general length two for each point as I could (potentially)
 % included high frequency decay terms from poles in J(omega) which would
 % have anti commutator parts
 cc(cnt+1:cnt+length(tmp)) = tmp; cnt =cnt+length(tmp);
 cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
 cc_acom{j}= [cc2I,cc1*0];
 
 end

 Kap1= inf; %max frequency truncation
 Kap2 = 5; %max level truncation
 tendps =0.7; 
 convfact = 2 * pi * light_speed * length_unit * 10^(-12);
 numpoints =[40, linspace(0,tendps,50).*convfact]; 
 saveuptotier =0;
 
 use_reduced_mat = false;
 if exist('H_vib','var')
 vibstates = max(length(H_vib),1);
 else
     vibstates = 1;
 end
 %%  Run the full calculation
 
 [Time_units,rho_vec,nn,total_prop_op]=multi_site_drude_and_brownian_HOM_3...
            (Htot,QQ,rho_0,cc_com,cc_acom,vv,Kap1,Kap2,vibstates,...
            numpoints,tendps,saveuptotier,use_reduced_mat,flname);

 
        