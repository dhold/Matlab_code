%% Set constant terms
%All terms in inverse cm, conversion to sensible units achieved by
%E = hc / lamba => 6.62606957 *10^{-33} 299792458 *10 * "E in code" gives
%Energy in Joules,  extra factor of 10 is because thats speed in cm/s
E1 = 1042; %Energy of initial state
V = 92; %Coupling
omegavib = 1111; %Frequency of vibrations / level spacing,
%assumed to be the same in each protein
viblvls = 20; %Number of vibrational levels to include, integer
coupg = 267.1; %coupling term between electronic and vibrational DOF

tendps = 10; %end time in pico seconds
convfact = 2*pi * 299792458 * 10^2 *10^(-12); %2 pi * "c" in cm/s * 1ps 
tend =convfact*tendps; %end time in units you get by putting everything in wavenumbers


%% Compute Matrix of subsystem (treated QM)
%Density matrix, electronic energy levels
rho = zeros(2);  rho(1,1) = 1;

%Excitonic Hamilonian
Hel = zeros(2); 

Hel(1,1) = E1; Hel(2,2)=0;  %On diagonal terms
Hel(1,2) = V;  Hel(2,1) = V;  % Coupling terms from dipole dipole interaction

%Vibrational Hamilonian, relative degrees of freedom only
tmp = ones(viblvls,1);
tmp(1) = 0;
for k=2:viblvls
    tmp(k) = tmp(k-1) + omegavib;
end
Hvib = diag(tmp); clear tmp 
%Note this is H_vib = omega_vib * (b_r^dag b_r) not the COM excitation

Htot = kron(Hel,eye(viblvls))+kron(eye(2),Hvib); %total Hamiltonian
%Each quadrant is vib Ham multiplied by the electric element

%Hamilonian coupling electronic and vibrational levels
%Hexvib = -g/sqrt(2)  * (n_1 - n_2)*(b^dag_rel + b_rel)

Hexvib = diag(-coupg/sqrt(2) * sqrt(1:viblvls-1),1);
Hexvib = Hexvib+Hexvib.';
Hexvib = kron([-1,0;0,1],Hexvib);

Htot = Htot + Hexvib;  %Full Hamiltonian with electric + vib + coupling
%in chromophore basis

   %eig( kron(Hel,eye(viblvls))) is purely electric Ham in full basis
    [basis_proj,Helnew] = eig( kron(Hel,eye(viblvls)) ); %project to delocalised basis
    

%% Time evolution with no decoherence from bath, no thermal init population

use_deloc_basis = true;  %If false will use chromophore site basis 
%if true will use diagonalised electron Ham basis

%exact diagonalisation is fine
init_state = zeros(length(Htot),1); 

%logic to find first excited state of elec w/ zero vib energy...
     tmplg1 = diag(Helnew)>0; %find positive values (excited)
     tmplg2 = abs(basis_proj(1,:))>eps;  %find values with stuff in zero vib state
     tmplg3 = tmplg1' & tmplg2;
     
 if use_deloc_basis
    
    Htot = basis_proj'*Htot*basis_proj;  %project Ham to deloc basis

    init_state(tmplg3) = 1; %First electronic excited state at zero vib ex
 
 else
     
     init_state = basis_proj(:,tmplg3); %#ok<UNRCH> %value in non projected basis

 end
[eigvc,eigvls] = eig(Htot);


init_state = eigvc'*init_state; %project to eigenbasis
state_at_time_t = @(t) exp(-1i*t*diag(eigvls)).*init_state; %time ev func

state_in_basis = @(t) eigvc*state_at_time_t(t); %project to initial basis

trange = linspace(0,tend,1000);

rhoYY = zeros(length(init_state),1000);
lpcnt = 0;
for tt = trange
    lpcnt = lpcnt+1;
    tmp = state_in_basis(tt);
    rhoYY(:,lpcnt) = tmp;
    %E_should_be_constant(lpcnt) = tmp'*Htot*tmp; appears fine
    
end


if use_deloc_basis
    
    rhoYY2 = sum(abs(rhoYY(~tmplg1,:)).^2);
    % Convention takes Hamiltonians with diagonally increasing
    % eigenvalues but this is not used for this work, hence the swapping
    
else
    rhoYY2 = basis_proj * rhoYY; %#ok<UNRCH>
    rhoYY2 = sum(abs(rhoYY2(~tmplg1,:)).^2)  
    
end



%% Plot graphs if I want to

figure
plot(trange/convfact,rhoYY2)  %Graph of occupation


%% Now do it with density matricies and a thermal population of vibrations

Temperature = 200; %in wavenumbers.. 
bose_ein_dist = 1./(exp(-diag(Hvib)/T)-1);

rho = zeros(length(Htot),length(Htot),length(trange));
%density matrix




%d/dt rho(t) = -i/hbar [H_tot, rho(t)]




