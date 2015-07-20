%%
om_0 = 1000; % transition frequency, inverse cm
omega = 999;  %frequency of the light used
lambda = om_0; 
mu_ab = 10; %transition dipole matrix element 
R = 5; %distance between monomers and dimer centre
theta = [pi/2,pi/3]; %angle between R and x axis
phi = [0,pi/6]; %angle between R and y axis
mu = mu_ab*[sin(theta).*cos(phi);sin(theta).*sin(phi);cos(theta)];

R12 = [0;0;R]; %R12 is taken as the z direction

mm(:,1) = (1i*pi/lambda)*cross(R12,mu(:,1));
mm(:,2) = (1i*pi/lambda)*cross(-R12,mu(:,2));

V = (dot(mu(:,1),mu(:,2)) - 3*mu(3,1)*mu(3,2))/R^3;



%% project coordinate system so that light is travelling down the new z axis


thetaR = pi/6; %angles of R relative to the light
phiR = pi/9;
Rotmat = [cos(thetaR),0,sin(thetaR);0,1,0;-sin(thetaR),0,cos(thetaR)];
Rotmat = [cos(phiR),-sin(phiR),0;sin(phiR),cos(phiR),0;0,0,1]*Rotmat;
%rotated by thetaR about y then phiR about (previous) z axis.  Can point in
%any direction in this way, in general need to sum over all directions

mu_nb =  Rotmat'*mu; %inverse rotation
R12_nb = Rotmat'*R12;

% H = - mu dot E - m dot B
%   = - mu dot E - pi*1i*/lambda ( R X mu ) dot B
% E = E0*[cos(omega*t), pm sin(omega*t),0]
% B = B0*[pm sin(omega*t) , -cos(omega*t),0]


H0 = [0,0,0;0,om_0,V;0,conj(V),om_0];
[proj_M,ex_E] = eig(H0); %transform to exciton basis (trivial here)
% ex_E = proj_M Hint proj_M'
energ = diag(ex_E);

E0 = 100; B0 = 100; %hbar = c = 1 cm^-1 = 1
%E = @(t) E0*[cos(omega*t), sin(omega*t),0];
E = @(t) E0*exp(-1i*omega*t)*[(1-exp(2i*omega*t))/2, 1i*(1-exp(2i*omega*t))/2,0];
%B = @(t) B0*[sin(omega*t), -cos(omega*t),0];
B = @(t) B0*exp(-1i*omega*t)*[1i*(1-exp(2i*omega*t))/2,-(1-exp(2i*omega*t))/2,0];
%Hint = @(t) -proj_M'*[0,dot(mu(:,1),E(t)) + dot(mm(:,1),B(t)),...
%               dot(mu(:,2),E(t)) + dot(mm(:,2),B(t));...
%             dot(mu(:,1),E(t)) + dot(mm(:,1),B(t)),0, 0 ;...
%               dot(mu(:,2),E(t)) + dot(mm(:,2),B(t)),0,0];
detuning = omega-energ(2:3); %take level for which detuning is smallest
tmp = min(abs(detuning)); 

Lint = @(t) kron(eye(length(H0)),Hint(t)) - kron(Hint(t),eye(length(H0)));

gam = 0.1; %spontanious decay rate of monomer excited state
%construct linblad operator
C1 = zeros(size(H0)); C2 = C1;
C1(1,2) = 1; C2(1,3)=1; C1 = proj_M'*C1*proj_M; C2 = proj_M'*C2*proj_M;
Linblad = -gam/2*(kron(C1'*C1+C2'*C2,eye(length(H0)))...
            +kron(eye(length(H0)),C1*C1'+C2*C2')) + ...
            gam*(kron(C1',C1')+kron(C2',C2')) ;

v1 = sym('v1'); v2 = sym('v2'); %analytically construct louivillian
Hintmat = [0,v2-v1,v1+v2;v1-v2,0,0;v2+v1,0,0]/sqrt(2);
%v1 = -(dot(mu(:,1),E(t)) + dot(mm(:,1),B(t)) ) etc etc
Lvil = kron(eye(length(H0)),Hintmat)-kron(Hintmat.',eye(length(H0))) ;

%% Solve in interaction picture with RWA
t = sym('t'); 
if 1==1
    om_0 = sym('om_0'); V=sym('V');  omega = sym('omega'); gam = sym('gam');
   H0 = [0,0,0;0,om_0,V;0,conj(V),om_0];
    [D,ex_E] = eig(H0); 
    C1(1,2) = 1; C2(1,3)=1; C1 = D'*C1*D; C2 = D'*C2*D;
Linblad = -gam/2*(kron(C1'*C1+C2'*C2,eye(length(H0)))...
            +kron(eye(length(H0)),C1*C1'+C2*C2')) + ...
            gam*(kron(C1',C1')+kron(C2',C2')) ;
end


U0 = expm(-1i* (kron(ex_E*t,eye(length(H0)))-kron(eye(length(H0)),ex_E*t)) );
Erwa = E0*exp(-1i*omega*t)*[1/2, 1i/2,0];
Brwa = B0*exp(-1i*omega*t)*[1i/2,-1/2,0];
v1 =  -(dot(mu(:,1),Erwa) + dot(mm(:,1),Brwa) );
v2 =  -(dot(mu(:,2),Erwa) + dot(mm(:,2),Brwa) );
Hintmat = [0,v2-v1,v1+v2;v1-v2,0,0;v2+v1,0,0]/sqrt(2);
Lint_I = U0'* kron(eye(length(H0)),Hintmat)-kron(Hintmat.',eye(length(H0))) *U0;
Linblad_I = U0'* Linblad *U0; %project to interaction picture

rho_0 = reshape([1,0,0;0,0,0;0,0,0],9,1); %ground state or thermal maybe
