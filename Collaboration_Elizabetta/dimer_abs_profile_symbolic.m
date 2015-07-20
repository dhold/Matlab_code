%% calculate internal things within dimer
om_0 = sym('w_0'); % transition frequency, inverse cm
omega = sym('w');  %frequency of the light used
lambda = sym('lam'); 
mu_ab = sym('mu');  %transition dipole matrix element 
R = sym('R'); %distance between monomers and dimer centre
%didn't like these as vectors subing back in
%theta = [sym('theta1'),sym('theta2')]; %angle between R and x axis
%phi = [sym('phi1'),sym('phi2')]; %angle between R and y axis
%mu = mu_ab*[sin(theta).*cos(phi);sin(theta).*sin(phi);cos(theta)];
theta1 = sym('theta1'); theta2 = sym('theta2'); %angle between R and x axis
phi1 = sym('phi1'); phi2 = sym('phi2'); %angle between R and y axis
mu = mu_ab*[sin(theta1).*cos(phi1),sin(theta2).*cos(phi2);...
            sin(theta1).*sin(phi1),sin(theta2).*sin(phi2);cos(theta1),cos(theta2)];

R12 = [0;0;R]; %R12 is taken as the z direction

mm(:,1) = (1i*pi/lambda)*cross(R12,mu(:,1));
mm(:,2) = (1i*pi/lambda)*cross(-R12,mu(:,2));

V = (dot(mu(:,1),mu(:,2)) - 3*mu(3,1)*mu(3,2))/R^3;

%% project coordinate system so that light is travelling down the new z axis


thetaR = sym('thetaR'); %angles of R relative to the light
phiR = sym('phiR');
Rotmat = [cos(-thetaR),0,sin(-thetaR);0,1,0;-sin(-thetaR),0,cos(-thetaR)];
Rotmat = Rotmat*[cos(-phiR),-sin(-phiR),0;sin(-phiR),cos(-phiR),0;0,0,1];
%rotated by thetaR about y then phiR about (previous) z axis.  Can point in
%any direction in this way, in general need to sum over all directions
%Rotmat = [1,0,0;0,cos(phiR),-sin(phiR);0,sin(phiR),cos(phiR)]*Rotmat;


mu_nb =  Rotmat.'*mu; %inverse rotation
R12_nb = Rotmat.'*R12;

% H = - mu dot E - m dot B
%   = - mu dot E - pi*1i*/lambda ( R X mu ) dot B
% E = E0*[cos(omega*t), pm sin(omega*t),0]
% B = B0*[pm sin(omega*t) , -cos(omega*t),0]

H0 = [0,0,0;0,om_0,V;0,conj(V),om_0];
[proj_M,ex_E] = eig(H0); %transform to exciton basis (trivial here)
 t = sym('t');
U0 = expm(-1i* (kron(ex_E*t,eye(length(H0)))-kron(eye(length(H0)),ex_E*t)) );
%time evolution operator

E0 = sym('E0'); B0 = sym('B0');%hbar = c = 1 cm^-1 = 1
%E = @(t) E0*[cos(omega*t), sin(omega*t),0];
E =  E0*exp(-1i*omega*t)*[(1+exp(2i*omega*t))/2, 1i*(1+exp(2i*omega*t))/2,0];
%B = @(t) B0*[sin(omega*t), -cos(omega*t),0];
B = B0*exp(-1i*omega*t)*[1i*(1+exp(2i*omega*t))/2,-(1+exp(2i*omega*t))/2,0];
Erwa = E0*exp(-1i*omega*t)*[1, 1i,0]/2;
Brwa = B0*exp(-1i*omega*t)*[1i,-1,0]/2;

v1 =  -(dot(mu(:,1),E) + dot(mm(:,1),B) );
v2 =  -(dot(mu(:,2),E) + dot(mm(:,2),B) );
Hintmat = proj_M.'*[0,v1,v2;v1,0,0;v2,0,0]*proj_M/sqrt(2);
Lint = -1i*kron(eye(length(H0)),Hintmat) - kron(Hintmat,eye(length(H0)));
Lint_I = U0'*Lint*U0;

v1 =  -(dot(mu(:,1),Erwa) + dot(mm(:,1),Brwa) ); %same elements in RWA
v2 =  -(dot(mu(:,2),Erwa) + dot(mm(:,2),Brwa) );
Hrwmat = proj_M.'*[0,v1,v2;v1,0,0;v2,0,0]*proj_M/sqrt(2);
Lintrwa = -1i*kron(eye(length(H0)),Hrwmat) - kron(Hrwmat,eye(length(H0)));
Lintrwa_I = U0'*Lintrwa *U0;

gam = sym('Gam'); %spontanious decay rate of monomer excited state
%construct linblad operator
C1 = zeros(size(H0)); C2 = C1;
C1(1,2) = 1; C2(1,3)=1; C1 = proj_M'*C1*proj_M; C2 = proj_M'*C2*proj_M;
Linblad = -gam/2*(kron(C1'*C1+C2'*C2,eye(length(H0)))...
            +kron(eye(length(H0)),C1*C1'+C2*C2')) + ...
            gam*(kron(C1',C1')+kron(C2',C2')) ;
Linblad_I  = U0'*Linblad*U0;

%% Solve in interaction picture with RWA
rho_0 = zeros(size(H0));
rho_0(1,1) = 1; rho_0  = reshape(rho_0,numel(rho_0),1);
tt = sym('tt');
tic
rho_first_order = - 1i*int(Lint_I*rho_0,t,0,tt); %maybe not a great idea
toc
v1 =  -(dot(mu(:,1),[1, 1i,0]/sqrt(2)) ); %same elements in RWA
v2 =  -(dot(mu(:,2),[1, 1i,0]/sqrt(2)) );
dipole_op = proj_M.'*[0,v1,v2;v1,0,0;v2,0,0]*proj_M/sqrt(2);

trace_at_t = trace(dipole_op*reshape(rho_first_order,3,3) );
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
tobesubin = {10000,9999,2*pi/10000,100,1e-2,  pi/3,pi/4,pi/6,pi/10 ,0,0,1,1,0.1};
rho_at_tt = vpa(subs(rho_first_order,tobesubed,tobesubin));
trace_at_tt = vpa(subs(trace_at_t,tobesubed,tobesubin));
trange = linspace(0,tend,1000); rho_at_t = zeros(length(rho_0),100);
for k = 1:100
    
rho_at_t(:,k) = double(subs(rho_at_tt ,tt,trange(k)));

end