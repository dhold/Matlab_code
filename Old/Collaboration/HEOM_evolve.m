%% Construct Hamiltonian 

V = 100; %Coupling
omegavib = 1111; %(angular) Frequency of relative vibrations
viblvls = 0; %Number of vibrational excited states to include, integer
coupg = 267.1; %coupling term between electronic and vibrational DOF
%J = 2*pi*(1/434-1/489)*10^7; %energy gap units cm^-1, because 
J = omegavib;

Temp =300; %units Kelvin
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); 
Beta = (2 * pi * hbar * light_speed * length_unit)/ ( Temp * boltz_const);

%end time in units you get by putting everything in wavenumbers

Hex = [0 , V; V, 2*J]; N=2; 
Hvib = diag(0:omegavib:viblvls*omegavib);

Htot = kron(eye(size(Hex)),Hvib) + kron(Hex,eye(size(Hvib)));

Hex_vib = diag(sqrt(1:viblvls),1);
Hex_vib = Hex_vib + Hex_vib';

Htot = Htot + kron(diag([coupg/sqrt(2),-coupg/sqrt(2)]),Hex_vib);
[tmp,tmp2] = eig(Hex); 
basis_proj = kron(tmp,eye(length(Hvib)));
%setup the initial condition, assume thermally populated vib lvls

vibpops = 1./(exp(Beta*omegavib*(1:(viblvls+1)))-1);
vibpops = vibpops/sum(vibpops);

%rho_0 = kron([1/2,1/2;1/2,1/2],diag(vibpops));  %assume exciton state populated
%rho_0 = kron([0.5,0.3;0.3,0.5],diag(vibpops));
rho_0 = kron([0,0;0,1],diag(vibpops));
rho_0 = basis_proj*rho_0*basis_proj'; %project back to chorophore basis


%% set parameters
Ecut = inf; %cut off energy used as a truncation condition with multiple modes
om_0 = 1108;  %om_0 is the values of brownian modes included
lambda =444.32; %reorganisation energy
gamma = 5.3; %damping of modes
%sqm = lambda./omegavib;

lam_dru = 40; %reorganisation energy of drude modes
gam_dru = 100; %width
 Kappa = 0;
 
 [cc1,cc2R,cc2I,vv1,vv2,QQ] = ...
    coeffients_from_brownian_new(lambda,gamma,om_0,Temp,Kappa,lam_dru,gam_dru);
 
cc1 = repmat(cc1,N,1); cc2R = repmat(cc2R,N,1); cc2I = repmat(cc2I,N,1); 
 vv1 = repmat(vv1,N,1); vv2 = repmat(vv2,N,1);


 Kap1= inf; %max frequency truncation
 Kap2 = 6; %max level truncation
tendps = 1; %end time in pico seconds

 convfact = 2 * pi * light_speed * length_unit * 10^(-12);
 tend = convfact * tendps; 
 numpoints =[200, [0.2,0.4,0.6,0.8,1]*convfact]; 
 %numpoints =200;
saveuptotier = 0;
 use_reduced_mat = false;
 vibstates = length(Hvib);
 %%  Run the full calculation
 %QQ=0 %uncomment to test without this factor
[Time_units,rho_vec,nn,total_prop_op]=multi_site_drude_and_brownian_HOM_2...
            (Htot,QQ,rho_0,cc1,cc2R,cc2I,vv1,vv2,Kap1,Kap2,vibstates,...
            numpoints,tendps,saveuptotier,use_reduced_mat,'saved_data.mat');