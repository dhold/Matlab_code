% Nonlinear interaction solve
light_speed = 299792458; %units m s^{-1}, NOT cm 
length_unit = 2*pi/100; % units m, this is 1/[1(cm)^{-1}]
hbar = 1.05457173*10^(-34); %units Js
boltz_const = 1.3806488 * 10^(-23);  
convfact = light_speed / length_unit * 10^(-12);
%converts ps to inverse cm 

 flename = 'test_dimer_plus_light.mat';

reorg_E = {[],20,20}; 
drude_gam = {[],53,53};
Temp_kel = 300;
dipole_vec = [4,1,0;3,sqrt(7),1];
syms t
pulse_width_fs = 30; 
pulse_width_icm = pulse_width_fs*convfact/1000; 
transverse_beam_factor = 1; %factor due to scaling from the beam waist
t0 = 5*pulse_width_icm ;
pump_power_J = 1.5e-6; %1-2 micro joule range

tendps = 2*t0;
numpoints =[200, linspace(0,2*t0,20).*convfact]; 

E_peak = sqrt(pump_power_J/pulse_width_icm)/transverse_beam_factor;
E_t = E_peak*exp(-(t-t0)^2/pulse_width_icm^2/2)/2 ; %Envelope for each rotating component
E_dirn = [1,0,0]; 
%R_disp = [0,0,0;0,0,0]; %not needed here as only dipole considerd

 Kappa = 0; %matsubara frequencies
 Kap2 = 1; %max level truncation
  Kap1= inf; %max frequency truncation
  wave_vec_trunc = 3;
  
om_0 = 10000; om_0j = cellfun(@sum,reorg_E)+[0, om_0 , om_0+300]; 

Ham = diag(om_0j); N = length(Ham);

omega = 10000;

MM = 100;  V = MM*(diag(ones(N-2,1),1) + diag(ones(N-2,1),-1));
%periodic boundary conditions mean that first and last are also coupled
%note that M(1,:) relates to the ground state only
V(1,end) = MM;  V(end,1) = MM;
Ham(2:end,2:end) = Ham(2:end,2:end) + V;

%%
dipole_op = Ham*0;
for k = 2:(1+size(dipole_vec,1))
dipole_op(1,k) = dot(E_dirn,dipole_vec(1,:));
end
dipole_op = dipole_op + dipole_op';
%%

 N = length(Ham);

N2 = sum(cellfun(@length,reorg_E));
 QQ = zeros(N,2); cnt=0;
 cc = zeros(1,sum(N2 + Kappa));
 clear  cc_com vv cc_acom
 for j=1:N
 [cc1,cc2R,cc2I,vv1,vv2,QQ(j,1),tmp] = ...
    coeffients_from_brownian_new([],[],[],Temp_kel,Kappa,reorg_E{j},drude_gam{j});
 % QQ is in general length two for each point as I could (potentially)
 % included high frequency decay terms from poles in J(omega) which would
 % have anti commutator parts
 cc(cnt+1:cnt+length(tmp)) = tmp; cnt =cnt+length(tmp);
 cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
  cc_acom{j}= [cc2I,cc1*0];
 end

use_RWA = false; rho_0 = [1,0,0;0,0,0;0,0,0];  
use_reduced_mat = false; saveuptotier=0;

[Time_units,rho_vec,nn,total_prop_op]=HEOM_plus_light...
            (Ham,QQ,rho_0,cc_com,cc_acom,vv,Kap1,Kap2,1,...
            E_t, omega,dipole_op, wave_vec_trunc , use_RWA ,...
            numpoints,tendps,saveuptotier,use_reduced_mat, flename);