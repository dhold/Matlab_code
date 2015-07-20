%% calculate internal things within dimer
om_0 = sym('w_0'); % transition frequency, inverse cm
omega = sym('w');  %frequency of the light used
mu_ab = sym('mu');  %transition dipole matrix element 

dipole_op = [0,mu_ab ;conj(mu_ab),0];
dipole_op_L = (kron( dipole_op,eye(length( dipole_op)))-kron(eye(length( dipole_op)), dipole_op.'));
H0 = [0,0;0,om_0];

 t = sym('t');
U0 = expm(-1i*t* (kron(H0,eye(length(H0)))-kron(eye(length(H0)),H0)) );
%time evolution operator

E0 = sym('E0');
%E = @(t) E0*[cos(omega*t), sin(omega*t),0];
E =  E0*(exp(-1i*omega*t)+exp(1i*omega*t))/2;
Erwa = E0*exp(-1i*omega*t)/2;

C1 = zeros(size(H0)); 
C1(1,2) = 1; gam = sym('gam');
Linblad = -gam/2*(kron(C1'*C1,eye(length(H0)))...
            +kron(eye(length(H0)),C1'*C1)) + gam*(kron(C1,C1)) ;
%% Solve in interaction picture with RWA
rho_0 = zeros(size(H0));
rho_0(1,1) = 1; rho_0  = reshape(rho_0,numel(rho_0),1);

first_order_corr = -1i*int( U0*dipole_op_L*rho_0*E,t,0,tt);
first_order_corr_RWA = -1i*int( U0*dipole_op_L*rho_0*Erwa,t,0,tt);
first_order_corr_disp = int( U0*(Linblad+ -1i*dipole_op_L*Erwa)*rho_0,t,0,tt); 
trace_at_t = trace(dipole_op*reshape(first_order_corr ,2,2) );

trace_at_t2 =subs(trace_at_t,[omega,om_0,mu_ab,E0],[1000,999,1,1]) ;

[a,b] = eig(U0*(Linblad -1i*dipole_op_L*E)*(U0)');

%% subin symbols with different parameters
%recall params
% om_0 = sym('w_0'); % transition frequency, inverse cm
% omega = sym('w');  %frequency of the light used
% lambda = sym('lam'); 
% mu_ab = sym('mu');  %transition dipole matrix element 
% R = sym('R'); %distance between monomers and dimer centre
% theta = [sym('theta1'),sym('theta2')]; %angle between R and x axis
% phi = [sym('phi1'),sym('phi2')]; 
% thetaR = sym('thetaR'); %angles of R relative to the light
% phiR = sym('phiR');
% E0 = sym('E0'); B0 = sym('B0'); gam = sym('Gam'); tend = sym('tt');
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}
% to run the code not as a function just uncomment and use each module
% Energy in Joules,  extra factor of 10 is because thats speed in cm/s
convfact = 2*pi * light_speed *length_unit  *10^(-12); %2 pi * "c" in cm/s * 1ps 
tend =convfact*2; %end time in units you get by putting everything in wavenumbers

 
tobesubed = {om_0,omega,lambda,mu_ab,R, theta1,phi1,theta2,phi2,  thetaR,phiR,E0,B0,gam};
tobesubin = {10000,9999,2*pi/10000,100,1e-2,  pi/3,pi/4,pi/6,pi/10 ,   0,0,1,1,0.1};
rho_at_tt = vpa(subs(rho_first_order,tobesubed,tobesubin));
trange = linspace(0,tend,1000); rho_at_t = zeros(length(rho_0),100);
for k = 1:100
    
rho_at_t(:,k) = double(subs(rho_at_tt ,tt,trange(k)));

end