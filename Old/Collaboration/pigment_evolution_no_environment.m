function [time_units,rho_vc,basis_proj] = pigment_evolution_no_environment
%This script will be used to calculate the response of a dimer of two 
% identical pigments with a relative vibrational mode and no environmental
% interactions

V = 100; %Coupling
omegavib = 1111; %(angular) Frequency of relative vibrations
viblvls = 6; %Number of vibrational excited states to include, integer
coupg = 267.1; %coupling term between electronic and vibrational DOF
J = 2*pi*(1/434-1/489)*10^7; %energy gap units cm^-1, because 

tendps = 1; %end time in pico seconds
Temp =300; %units Kelvin
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); 
Beta = (2 * pi * hbar * light_speed * length_unit)/ ( Temp * boltz_const);
tend = 2 * pi * length_unit * light_speed * 10^(-12) * tendps; 
%end time in units you get by putting everything in wavenumbers

Hex = [0 , V; V, 2*J];

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

rho_0 = kron(diag([0,1]),diag(vibpops));  %assume exciton state populated
rho_0 = basis_proj*rho_0*basis_proj'; %project back to chorophore basis

%time to evolve in Louiville space

LL = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));
rho_vec= reshape(rho_0,numel(rho_0),1);
propog(-1,LL);

 [time_units,rho_vc] = ode45(@propog,[0,tend],rho_vec);
 time_units = time_units / (2 * pi * length_unit * light_speed) *(10^12);
    function drho = propog(t,rho_vc)
             persistent LL
             if t == -1
                 LL = rho_vc;
             end
                drho =  LL*rho_vc;
