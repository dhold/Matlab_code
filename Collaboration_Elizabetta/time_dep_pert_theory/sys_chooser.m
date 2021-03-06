function [H_site,mu,R,lambda,gam,om_0,lam_dru, gam_dru,om_vib,displ,...
             pdm,sd_shift] = sys_chooser(sys)

         %set things that are not always output as empty
        pdm = []; sd_shift = [];   
        lambda = [] ;gam = [] ; om_0 = [] ;
        lam_dru = [] ; gam_dru = [] ;om_vib = [] ;displ = [] ;
         
%% PC645 dimer
if sys == 1
%dimer specific parameters

    fle = open('Hamiltonian_save.mat');
        H_site = fle.PC645Hamiltonian;
            E0 = fle.E0_PC645; %these energies are relative 
            %to the lowest exciton state anyway
            N=2;
       H_site = H_site(1:N,1:N);%only these 2
       N=length(H_site);
       H_site = H_site - eye(N)*E0;
      


mu = fle.PC645_dip_and_pos([4,8],1:3); %4 and 8 are the dimer
take_same_amp = true;
if take_same_amp
%assume allare the same amplitude
for k = 1:N
    mu(k,:) = mu(k,:)/norm(mu(k,:)); %site based mu
end
end
R = fle.PC645_dip_and_pos([4,8],4:6); %4 and 8 are the dimer
%convert units from angstrom to cm^
R = R*1e-8;

%Bath parameters 

lambdaD =100;
omegaD=100;
omegaU = 650;  %omegaU = 1034; 
gammaU = 5.30884;
lambdaU = 44;

om_0 = {omegaU,omegaU};   
lambda ={lambdaU ,lambdaU }; %reorganisation energy
gam = {gammaU,gammaU}; 
%include these explicitly in Hamiltonian
om_vib = [om_0{:}];  
displ = [[0;0],eye(2)*sqrt(2*lambdaU/omegaU)];

%over damped, include in Redfield only
 lam_dru = {lambdaD,lambdaD}; %reorganisation energy of drude modes
gam_dru = {omegaD,omegaD};   
sd_shift = [0,0];
elseif sys == 2
%% B820_bacteriochlorophyll_dimer
flenme ='B820_bacteriochlorophyll_dimer_params.mat'; 
%insert name of file containing the system parameters
    fle = open(flenme);
        H_site = fle.H_site; N = length(H_site);        
        
mu = fle.mu;  pdm = fle.pdm; %set to empty if not known
take_same_amp = true;
if take_same_amp
%assume all are the same amplitude
for k = 1:N
    mu(k,:) = mu(k,:)/norm(mu(k,:)); %site based mu
end
end
%mu(2,:) = mu(1,:); %test with them the same
R = fle.R;  R = R*1e-8;   %convert units from angstrom to cm^-1

%read out bath parameters
lam_dru  = fle.lam_dru; gam_dru = fle.gam_dru;
%lam_dru = {250,250}; %gam_dru = {40,40}; %stronger Drude
%underdamped
om_0 = fle.om_0; lambda = fle.lambda;  gam = fle.gamma;
%explicitly included modes

om_vib = [om_0{:}];

displ = [[0;0],[sqrt(2*lambda{1}/om_0{1}) ; 0],...
            [0 ; sqrt(2*lambda{2}/om_0{2})],...
    [sqrt(2*lambda{1}/om_0{1});sqrt(2*lambda{2}/om_0{2})]];
%displacements of states, ground, excited states and finally double excited

sd_shift = [200,200]; %standard deviation of site energy fluctuations 
elseif sys==3  %calculate parameters for general set up
 pdm = [];
 R = [0,0,-1/2;0,0,1/2];
 %take these mu to be align at an angle theta to R
 theta1 = 45*pi/180;  phi1 = 20*pi/180;
 theta2 = 45*pi/180;  phi2 = -20*pi/180;
 ry = @(a) [cos(a),0,sin(a);0,1,0;-sin(a),0,cos(a)]; %a runs from 0 to pi
rz = @(b) [cos(b),-sin(b),0;sin(b),cos(b),0;0,0,1];
 
 mu(1,:) = (R(2,:)-R(1,:)) * ry(theta1)*rz(phi1);
 mu(2,:) = (R(2,:)-R(1,:)) * ry(theta2)*rz(phi2);
 R = R*10^(-8);
 
V=200; E0 = 10^7/800;       H_site = [E0, V; V, E0];
 lam_dru  = {35,35}; gam_dru = {50,50};
%lam_dru = {250,250}; %gam_dru = {40,40}; %stronger Drude
%underdamped
om_0 = {400,400}; lambda = {V/2,V/2};  gam = {5.3088,5.3088}; %about ps damping
%explicitly included modes

om_vib = [om_0{:}];
displ = [[0;0],[sqrt(2*lambda{1}/om_0{1}) ; 0],...
            [0 ; sqrt(2*lambda{2}/om_0{2})],...
    [sqrt(2*lambda{1}/om_0{1});sqrt(2*lambda{2}/om_0{2})]];
%displacements of states, ground, excited states and finally double excited


sd_shift = 0*[200,200];   %no disorder
elseif sys==4  %p545, NOT dimer, 8 sites
    

%values taken from paper Novoderezkhin et all Excitation dynamics in
%Phycoerythrin

%checked agrees
dipole_mom_and_pos = ...
 [-0.088, -2.605, -4.019, 12.17,  8.038, 25.032, 37.431;...
  -4.506, -2.670,  0.200, 13.32,  5.368, 24.566, 5.915;...
  5.102,   0.280,  1.574, 13.59,  -8.024, 44.047, 42.084;...
  2.898,   1.310, -3.721, 12.44,  10.976, 47.584, 29.106;...
  0.448,  -0.873,  4.691, 12.18,  17.923, 10.988, 18.324;...
  0.126,   2.422,  3.992, 11.87,  -3.471, 17.026, 27.542;...
 -3.900,   2.600,  1.258, 12.33, -16.170, 30.834, 12.607;...
  1.383,  -4.355, -1.445, 12.18, -11.503, 24.266, 50.469];

mu = dipole_mom_and_pos(:,1:3).*repmat(dipole_mom_and_pos(:,4)...
                    ./sum(dipole_mom_and_pos(:,1:3).^2,2),1,3);
R = dipole_mom_and_pos(:,5:7)*1e-8;   %convert units from angstrom to cm^-1

%checked agrees
H_int = ... %not sure if H(3,2) = -4 or +4
[0,    1,   -37, 37,   23, 92,  -16,  12,;...
 1,    0,   -4, -11,   33,  -39, -46,  3;...
-37,  -4,    0,  -45,  3,   2,  -11,  34;...
 37,  -11, -45,   0,   -7,  -17, -3,  6;...
 23,  33,    3,   -7,   0,   18,  7,  6;...
 92,  -39,   2,   -17,  18,  0,  40,  26, ;...
-16,  -46,  -11,  -3,   7,  40,   0,  7;...
 12,   3,     34,  6,   6,  26,   7,  0];                    

Site_E =[ 18532 ,18008,17973,18040 18711 ,19574 ,19050,18960  ]; 
%Site_E = [18641 18043 18030 19111 18499 19588 17909 18868]; %model E47
H_site = H_int  + diag(Site_E);
%checked agrees
Freq_and_HR= ...
[207 ,  0.0013 ;...
 244 , 0.0072 ;...
 312,  0.0450 ;...
 372,   0.0578 ;...
 438,  0.0450  ;...
 514,  0.0924  ;...
 718,  0.0761 ;...
 813 , 0.0578 ;...
 938 ,  0.0313 ;...
 1111 ,  0.0578 ;...
 1450 ,  0.1013 ;...
 1520 ,  0.0265 ;...
 1790 ,  0.0072 ;...
 2090 ,  0.0113 ];
tmp = prod(Freq_and_HR,2);
om_0={}; lambda={}; gam = {}; lam_dru={}; gam_dru = {};
for k = 1:8
om_0{k} = Freq_and_HR(:,1); 
lambda{k} = tmp;
%paper takes all mode damping to be = 20cm^-1
gam{k} = 20*ones(14,1);

%For modeling of the PE545 complex the spectral density (Eq S3) is constructed as a sum 
%of two overdamped Brownian oscillators (characteristic frequencies ?01=30 cm?1
%, ?02=90 cm?1 and couplings ?01=40 cm?1, ?02=70 cm?1

 lam_dru{k}  = [40,70]; gam_dru{k} = [30,90];
end

om_vib = [];
sd_shift = 400/sqrt(8*log(2))*ones(8,1); %ignore static disorder for now
% [1, 'PEB50/61C', 0.088, 2.605, 4.019, 12.17, 8.038, 25.032, 37.431;...
% 2, 'DBVA', 4.506, 2.670, 0.200, 13.32, 5.368, 24.566, 5.915;...
% 3, 'DBVB', 5.102, 0.280, 1.574, 13.59, 8.024, 44.047, 42.084;...
% 4, 'PEB82C', 2.898, 1.310, 3.721, 12.44, 10.976, 47.584, 29.106;...
% 5, 'PEB158C', 0.448, 0.873, 4.691, 12.18, 17.923, 10.988, 18.324;...
% 6, 'PEB50/61D', 0.126, 2.422, 3.992, 11.87, 3.471, 17.026, 27.542;...
% 7, 'PEB82D', 3.900, 2.600, 1.258, 12.33, 16.170, 30.834, 12.607;...
% 8, 'PEB158D', 1.383, 4.355, 1.445, 12.18, 11.503, 24.266, 50.469];


% PEB50/61C, 0, 1, (3), 37, (42), 37, (33), 23, (29), 92, (71), 16, (15), 12, (10)
% DBVA, 1, (3), 0, 4, (7), 11, (17), 33, (20), 39, (44), 46, (69), 3(3)
% DBVB, 37, (42), 4, (7), 0, 45, (69), 3, (4), 2, (4), 11, (17), 34, (18)
% PEB82C, 37, (33), 11, (17), 45, (69), 0, 7, (7), 17, (17), 3, (4), 6, (7)
% PEB158C, 23, (29), 33, (20), 3, (4), 7, (7), 0, 18, (16), 7, (8), 6, (6)
% PEB50/61D, 92, (71), 39, (44), 2, (4), 17, (17), 18, (16), 0, 40, (35), 26, (32)
% PEB82D, 16, (15), 46, (69), 11, (17), 3, (4), 7, (8), 40, (35), 0, 7, (7)
% PEB158D
% PEB158D, 12, (10), 3, (3), 34, (18), 6, (7), 6, (6), 26, (32), 7, (7), 0
end