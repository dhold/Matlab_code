
M = -250; om_0  = 10^7/2*(1/820+1/795);
H_site = [om_0,M;M,om_0];

delta = 0.7;
alpha = pi*30/180; %angle between dipoles in that paper

%parameters taken for PDB struct for dimers in LH1 ring
NA1 = [5.312 -39.092   1.415];NC1 =  [9.315 -39.508   2.422];
NA2 = [4.935, -40.264,  11.738];NC2 = [ 9.105, -40.557,  11.514];
mu = [NA1-NC1;NA2-NC2]; 
%these should be equal in absolute value, and magnitude is not important as
%everything is directly proportional to this 
mu(1,:) = mu(1,:)/norm(mu(1,:)); mu(2,:) = mu(2,:)/norm(mu(1,:));
alpha_mine = acos(dot(mu(1,:),mu(2,:))/norm(mu(1,:))/norm(mu(2,:)));
%if I find the structure data it might be worth actually considering this
R = [7.332,-39.11,1.946;6.953,-40.373,11.769]; %in Angstroms
mu_normal = cross(mu(1,:),mu(2,:)); 
Delta_R = R(1,:) - R(2,:);
theta_mine = acos(dot(mu_normal,Delta_R)/norm(Delta_R)/norm(mu_normal));

pdm =[]; % I don't even know how to work this out
%Also haven't actually found drude spec
lam_dru  = {100,100}; gam_dru = {100,100};
%underdamped mode(s)
OMEGA = 300; %450 and 750 also listed, not sure 
om_0 ={OMEGA,OMEGA}; lambda = {delta^2/2*OMEGA,delta^2/2*OMEGA};   
gamma =  {10,10}; %damping listed in paper
%explicitly included modes
om_vib = [];  numvib = []; %set this to whatever
displ = [];


save('B820_bacteriochlorophyll_dimer_params.mat','H_site','mu','R','pdm',...
        'lam_dru','gam_dru','om_0','lambda','gamma','om_vib','numvib','displ')
    
    
    
    
    