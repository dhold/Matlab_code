%Calculate the 2D spectra following the HEOM method of Xu and Yan
%In mixed Hberg Schrodinger picture, orientational avearging will have to
%be done by averaging over angles

%% constants and scaling stuff
 speed_unit = 2*pi * 299792458; %units m s^{-1}, NOT cm 
 length_unit = 0.01; % units m
 hbar = 6.62606957*10^(-34) ;  %units Js
 boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}
 %temperature dep
 Temp = 300; %in kelvin
 beta = ( hbar * speed_unit / length_unit)/ ( Temp * boltz_const); %in cm^-1
convfact = speed_unit / length_unit  *10^(-12); %2 pi * "c" * 100m^-1 * 10^-12 s 

%% System properties, independent of orientation

flenme =''; %insert name of file containing the system parameters
    fle = open('Hamiltonian_save.mat');
        H_site = fle.Hamiltonian;
        
mu = fle.mu;  pdm = fle.pdm; %set to empty if not known

R = fle.R;  R = R*1e-8;   %convert units from angstrom to cm^-1
if isempty(pdm)
    pdm_n_pos = [];
else
   pdm_n_pos = [pdm;R]; 
end
%read out bath parameters
lam_dru  = fle.lam_dru; gam_dru = fle.gam_dru;
%underdamped
om_0 = fle.om_0; lambda = fle.lambda;  gamma = fle.gamma;
%explicitly included modes
om_vib = fle.om_vib;  numvib = fle.numvib; %set this to whatever
displ = fle.om_vib;
     
[ex_basis,H0ex] = eig(H0);

[H_ex_vib,fock_space_rep,mu_full] = ...
    generate_ex_vib_ham(H_site(2:end,2:end),om_vib,numvib,displ,mu,pdm_n_pos) ;

%% Beam properties of the 3 input beams
%unit vectors in direction
k1 = [1/sqrt(2);0;1/sqrt(2)]; k2 = [-1/sqrt(2);0;1/sqrt(2)];
k3 = [0;0;1];
%polarization wavevectors
pol{1} = [1;0;0]; pol{2} = [1;0;0]; pol{3} = [1;0;0]; 

%% Next Calculate HEOM coefficients

% set reference frequencies for excited and double excited manifolds
% to simplify numerics, 
om_eg = mean(diag(H_ex_vib(2:N+1,2:N+1)));  
om_fg = mean(diag(H_ex_vib(N+2:end,N+2:end)));
om_fe = om_fg-om_eg;

Htot_rescale = Htot;
Htot_rescale(2:N+1,2:N+1) = Htot_rescale(2:N+1,2:N+1) - diag(om_eg*ones(N,1));
Htot_rescale(N+2:end,N+2:end) = Htot_rescale(N+2:end,N+2:end) - diag(om_fg*ones(N*(N-1)/2,1));

[total_prop_op,nn]=HEOM_op_gen(Htot,QQ,cc_com,cc_acom,vv,Kap1,Kap2,viblvls);
% this operator will decide the propogator in all cases

input_state = Htot_rescale(1:N+1,1:N+1)*0; input_state(1,2:N+1)=1;
[~,~,~,~,savehrarch] = reduced_HEOM(input_state,[],...
                    Htot,QQ,cc_com,cc_acom,vv,Kap1,Kap2,1);
%savehrarch = {nn,coup_com_save,coup_acom_save,const_factor}; pass back to
%function when doing the calculations


%% Average over Euler angles
syms a b y real
%a,b,y are the euler angles (normally alpha beta gamma or phi theta xi)

%Z1Y2Z3 choice
r1 = [cos(y),-sin(y),0;sin(y),cos(y),0;0,0,1];
r2 = [cos(a),0,sin(a);0,1,0;-sin(a),0,cos(a)]; %a runs from 0 to pi
r3 = [cos(b),-sin(b),0;sin(b),cos(b),0;0,0,1];

TT = r3*r2*r1;  %Matrix of rotation, subs for diff values
TTfn = matlabFunction(TT,'vars',{y,a,b});
%must sum over a variety of angles, weighted with sin(theta)/8pi^2 dist

max_it = 100;
tmp = rand(max_it,1)>1/2;
theta_abis = asin(rand(max_it,1));
theta_abis(tmp) = pi/2-theta_abis(tmp); %50% will be over pi/4
theta_abis(~tmp) = pi/2+theta_abis(~tmp); 
theta_weight = sin(theta_abis)/2; %probability of that value being taken
y_abs = 2*pi*rand(max_it,1); y_weight =1/2/pi;
b_abs = 2*pi*rand(max_it,1); b_weight =1/2/pi;
%use these points for the Euler angles

%% Calculate R^(3)(t3,t2,t1) for each orientation
%rescale k unit vectors with a typical c/omega factor, this is typically ok
%provided beams aren't scanned that much this value may as well be omega_ge

k1 = k1*2*pi/om_eg; k2 = k2*2*pi/om_eg; k3 = k3*2*pi/om_eg; 
rho_eq = zeros(size(Hex_vib)); 
if all(cellfn(isempty(numvib)))
rho_eq(1)=1;
else
    warning('write something that thermally populates vibrations in GS')
    rho_eq(1)=1;
end

for lp = 1:max_it%should converge far before reaching max it
    TT = TTfn(y_abs,theta_abis,b_abs);
    point_weight = y_weight*theta_weight*b_weight;
    
    mu_rot = TT*mu;  R_rot = TT*R;
    %calculate operator which will be mu, this is only the case if ALL
    %beams are polarized along x, which we assume here
    
    %Calculate first the k_s = k_3 + k_2 - k_1 signal
    k_rp = k3 + k2 - k1;
    k_nr = k3 + k1 - k2;
    % AKA the rephasing signal
    % due to the spatial dependence of each of the chromophores in my
    % formalism each operator has a slightly different phase dependent term
    % if one is not interested in CD this formalism is not required
       mu_w_spat{1} = mu_rot.*exp(1i*R_rot'*(-k1)); 
       mu_w_spat{2} = mu_rot.*exp(1i*R_rot'*k2); 
       mu_w_spat{3} = mu_rot.*exp(1i*R_rot'*k3); 
       mu_w_spat{4} = mu_rot.*exp(1i*R_rot'*(-k_rp)); 

    mu_ge{1} = mu_w_spat{1}*pol{1};   %mu_ef{1} = zeros(N,N*(N-1)/2); %not needed
    de_lg = fock_space_rep(sum(double(fock_space_rep))==2,:); %double exciton manifold
    for j = 2:3
        mu_ge{j} = mu_w_spat{j}*pol{j}; temp = zeros(N,N*(N-1)/2);
        for k =1:N
            lg1 = de_lg(:,k) == 1; %find elements with 1 at k, should be line  
            rng = [1:k-1,k+1:N];
            for kk = rng
                lg2 = de_lg(:,kk) == 1;
            temp(k,lg2 & lg3) =  mu_w_spat{j}(:,kk)*pol{j};
            end
        end
        mu_ef{j} = temp;
    end
    
    %final one is vector in nature
    mu_ge_x = mu_w_spat{4}(1,:); mu_ge_y = mu_w_spat{4}(2,:);  mu_ge_z = mu_w_spat{4}(3,:);
    temp_x = zeros(N,N*(N-1)/2);temp_y = zeros(N,N*(N-1)/2);temp_z = zeros(N,N*(N-1)/2);
        for k =1:N
            lg1 = de_lg(:,k) == 1; %find elements with 1 at k, should be line  
            rng = [1:k-1,k+1:N];
            for kk = rng
                lg2 = de_lg(:,kk) == 1;
            temp_x(k,lg2 & lg3) =  mu_w_spat{j}(1,kk); %there should be no overlap
            temp_y(k,lg2 & lg3) =  mu_w_spat{j}(2,kk);
            temp_z(k,lg2 & lg3) =  mu_w_spat{j}(3,kk);
            end
         mu_ef_x = temp_x; mu_ef_y = temp_y;  mu_ef_z = temp_z;            
        end    
    %first calculate rho_eg(t_1)
    rho_eg = mu_ge{1}*rho_eq;
    %propogate with the HEOM
    [rho_eg_t1_tmp,t_range1] = reduced_HEOM(rho_eg,tend,Htot,QQ,cc_com,cc_acom,vv,savehrarch,Kap2,numpoints);
    %reshape this to an actual density matrix, plus auxilliary orders later    
    %this actually contains all the info required for the linear spectrum
    %if this is of interest, 
    
    %Next calculate the time dependence of the operator last used
    
    mu_ge_x
    
    %Now calculate rho_gg(t_2;t_1) rho_ee(t_2;t_1) NOT for RP rho_fg(t_2;t_1)
    %will have to loop over each t_1 sadly...
    
    for lp2 = 1:numpoints
        
        rho_gg_t2_t1 = 
    
    end
end
