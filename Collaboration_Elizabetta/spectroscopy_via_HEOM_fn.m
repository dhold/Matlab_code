function [V_mat,nn,av_freq,freq_scale] = spectroscopy_via_HEOM_fn(flename, Ham,...
                reorg_E, drude_gam, Temp_kel,dipole_vec,R_disp,Kappa,Kap1,Kap2,tendps,numpoints)
light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23);  
    convfact = 2 * pi * light_speed * length_unit * 10^(-12);
    
if nargin == 0
                
 flename = {'test_specxR.mat','test_specyR.mat','test_speczR.mat',...
     'test_specxL.mat','test_specyL.mat','test_speczL.mat'};

reorg_E = {[],20,20}; 
drude_gam = {[],53,53};
Temp_kel = 300;
dipole_vec = [1,0,0;1,0,0];
R_disp = [0,0,0;0,0,0]; %not important here, but in general is

 Kappa = 0; 
 Kap2 = 10; %max level truncation
  Kap1= inf; %max frequency truncation
  
om_0 = 10000; om_0j = cellfun(@sum,reorg_E)+[0, om_0 , om_0+300]; 

Ham = diag(om_0j); N = length(Ham);

MM = 100;  V = MM*(diag(ones(N-2,1),1) + diag(ones(N-2,1),-1));
%periodic boundary conditions mean that first and last are also coupled
%not that M(1,:) relates to the ground state only
V(1,end) = MM;  V(end,1) = MM;
Ham(2:end,2:end) = Ham(2:end,2:end) + V;
end

%%

 N = length(Ham);

%%%%
tmp = diag(Ham);
av_freq = mean(tmp(2:end));
freq_scale = [0,ones(1,N-1)*av_freq];
Ham = Ham - diag(freq_scale);
%%%%%

rho_0 = zeros(size(Ham)); rho_0(1,1) = 1;
%pi/2 rotations about various axes
R_z = [0,-1,0;1,0,0;0,0,1]; R_y = [0,0,1;0,1,0;-1,0,0]; R_x = [1,0,0;0,0,-1;0,1,0];
% directions of light with cp light
eps_zp = [1;1i;0]/sqrt(2); eps_xp = R_y *eps_zp;  eps_yp = R_x' *eps_zp;
eps_zm = [1;-1i;0]/sqrt(2); eps_xm = R_y *eps_zm;  eps_ym = R_x' *eps_zm;

beta_zp = cross([0;0;1],eps_zp); beta_yp = cross([0;1;0],eps_yp); beta_xp = cross([1;0;0],eps_xp);
beta_zm = cross([0;0;1],eps_zm); beta_ym = cross([0;1;0],eps_ym); beta_xm = cross([1;0;0],eps_xm); 

V_mat = zeros(N,N,6);

for k =1:N-1
    
    V_mat(1,k+1,1) = -dot(eps_xp,dipole_vec(k,:)) + 1i*dot(beta_xp,...
        (Ham(k+1,k+1)+freq_scale(k+1))/2*cross(R_disp(k,:),dipole_vec(k,:)));
    V_mat(1,k+1,2) = -dot(eps_yp,dipole_vec(k,:)) + 1i*dot(beta_yp,...
        (Ham(k+1,k+1)+freq_scale(k+1))/2*cross(R_disp(k,:),dipole_vec(k,:)));
    V_mat(1,k+1,3) = -dot(eps_zp,dipole_vec(k,:)) + 1i*dot(beta_zp,...
        (Ham(k+1,k+1)+freq_scale(k+1))/2*cross(R_disp(k,:),dipole_vec(k,:)));
    
    V_mat(1,k+1,4) = -dot(eps_xm,dipole_vec(k,:)) + 1i*dot(beta_xm,...
        (Ham(k+1,k+1)+freq_scale(k+1))/2*cross(R_disp(k,:),dipole_vec(k,:)));
    V_mat(1,k+1,5) = -dot(eps_ym,dipole_vec(k,:)) + 1i*dot(beta_ym,...
        (Ham(k+1,k+1)+freq_scale(k+1))/2*cross(R_disp(k,:),dipole_vec(k,:)));
    V_mat(1,k+1,6) = -dot(eps_zm,dipole_vec(k,:)) + 1i*dot(beta_zm,...
        (Ham(k+1,k+1)+freq_scale(k+1))/2*cross(R_disp(k,:),dipole_vec(k,:)));
end
V_mat = V_mat+conj(permute(V_mat,[2,1,3]));

%rho_0x+ = V_mat(:,:,1) * rho_0;  etc
%%
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


 if nargin<12
 numpoints =[200, linspace(0,tendps,20).*convfact]; 
 end
 saveuptotier =0;
 
 use_reduced_mat = false;   vibstates = 1; 

%%
for k =1:6
    tic
    rho_init = V_mat(:,:,k)*rho_0;
 [Time_units,rho_vec,nn]=multi_site_drude_and_brownian_HOM_3...
            (Ham,QQ,rho_init,cc_com,cc_acom,vv,Kap1,Kap2,vibstates,...
   numpoints,tendps,saveuptotier,use_reduced_mat,flename{k});
if k == 1
    load(flename{k},'nntot')
Kap1 = nntot; %function is overloaded to allow the passing of precomputed nn values
end
toc
end

