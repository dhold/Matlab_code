%Calculate the 2D spectra following the HEOM method of Xu and Yan
%In mixed Hberg Schrodinger picture, orientational avearging will have to
%be done by averaging over angles

%% constants and scaling stuff
 speed_unit = 2 * pi * 299792458; %units m s^{-1}, NOT cm 
 length_unit = 0.01; % units m
 hbar = 6.62606957*10^(-34) ;  %units Js
 boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}
 %temperature dep
 Temp = 300; %in kelvin
 beta = ( hbar * speed_unit / length_unit)/ ( Temp * boltz_const); %in cm^-1
convfact = speed_unit / length_unit  *10^(-12); %2 pi * "c" * 100m^-1 * 10^-12 s 
% HEOM parameters for convergence and such like
Kappa = 0; %# matsubara frequencies
Kap2 = 2; %Tier truncation
Kap1 = inf; %high freq cut off
% other parameters relating to how long evolution is considered for
t_end_ps = 4;
t_end = t_end_ps*convfact;
%% System properties, independent of orientation

flenme ='B820_bacteriochlorophyll_dimer_params.mat'; 
%insert name of file containing the system parameters
    fle = open(flenme);
        H_site = fle.H_site;

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
     
[H_ex_vib,fock_space_rep,mu_full,H_exciton,indiv_op] = generate_ex_vib_ham...
    (H_site,om_vib,numvib,displ,sqrt(sum(mu.^2,2)),pdm_n_pos) ;

[ex_basis,H0ex] = eig(H_ex_vib);
H_vib = indiv_op{2};

%% Beam properties of the 3 input beams
%unit vectors in direction
k1 = [1/sqrt(2);0;1/sqrt(2)]; k2 = [-1/sqrt(2);0;1/sqrt(2)];
k3 = [0;0;1];
%polarization wavevectors
pol{1} = [1;0;0]; pol{2} = [1;0;0]; pol{3} = [1;0;0];  

%% Next Calculate time dependent dipole operators via HEOM

N = length(H_site); %number of sites

if isempty(H_vib);
    H_vib = 0;
end

%Calculate HEOM coefficients
 QQ = zeros(length(H_site),2); cnt=0;
 cc = zeros(1,sum(cellfun(@length,lambda))+ sum(cellfun(@length,lam_dru)) + Kappa);
 clear cc_com cc_acom vv nn
 for j=1:length(H_site)
     
 [cc1,cc2R,cc2I,vv1,vv2,QQ(j,1),tmp] = coeffients_from_brownian_new(...
    lambda{j},gamma{j},om_0{j},Temp,Kappa,lam_dru{j},gam_dru{j});
 % QQ is in general length two for each point as I could (potentially)
 % included high frequency decay terms from poles in J(omega) which would
 % have anti commutator parts
 cc(cnt+1:cnt+length(tmp)) = tmp; cnt =cnt+length(tmp);
 cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
 cc_acom{j}= [cc2I,cc1*0];
 
 end

for j=1:N
    temp = H_ex_vib*0; 
    temp(1:length(H_vib),(1+j*length(H_vib)):((j+1)*length(H_vib))) =...
            diag(ones(length(H_vib),1));

    %populate all the ground -> excitation of site j states with vib included
mu_ge{j} = temp+temp';
end
%now the coherences to the double exciton states
    lg01 =  sum(double(fock_space_rep),2)==1;
    lg02 =  sum(double(fock_space_rep),2)==2; tmp = 1:size(fock_space_rep,1);
%    rnj = 1:size(fock_space_rep,1); rnj=rnj(lg); 
 for j=1:N  %prefactor is mu_j, so states which it can mix to must have an
     %excitation AT jth site
     lg2 = fock_space_rep(:,j) & lg02;
     temp = H_ex_vib*0; 
     for k =1:N  
         if k~=j
             lg = fock_space_rep(:,k); %states with ex at kth
            %lg1 =  lg & lg01; %just kth %this will be at 1+k in matrix
            lg22 = lg & lg2; %states with ex at kth and jth
    tmp2 = tmp(lg); %to these elements
    for kk = 1:length(tmp2)
       temp((1+k*length(H_vib)):((k+1)*length(H_vib)),...
           1+(tmp2(kk)-1)*length(H_vib):tmp2(kk)*length(H_vib)) = ...
           diag(ones(length(H_vib),1));
    end
         end
     end
    %populate all the ground -> excitation of site j states with vib included
mu_ef{j} = temp+temp'; mu_ef{j} = mu_ef{j}-diag(diag(mu_ef{j}));     
 end   
%project the single site coherences to single exciton coherences
%  for j = 1:N
%  mu_g_1ex{j}  = ex_basis'*mu_ge{j}*ex_basis;
%  end
%% evolve in time with HEOM
to_pass1 = Kap1; to_pass2 = Kap1; 
% [total_prop_op,nn]=HEOM_op_gen(H_ex_vib,fock_space_rep,...
%                         QQ,cc_com,cc_acom,vv,Kap1,Kap2,viblvls);
numpoints = 2^11; %might as well make it a power of 2 for the fft
des_t_range = zeros(numpoints,1); 
des_t_range(:) = linspace(0,t_end_ps,numpoints); %value in ps
des_t_range = des_t_range*convfact; %convert to inverse cm from ps
rel_freq_range = linspace(0,2*pi/(des_t_range(2)-des_t_range(1)),numpoints);

for j=1:N
    H_ex_vib_sc = H_ex_vib;  % scale elements by transition frequency to simplify numerics
    H_ex_vib_sc(2:N+1,2:N+1) = H_ex_vib_sc(2:N+1,2:N+1) - ...
                            diag(H_ex_vib(j+1,j+1)*ones(N,1));
[out1,trange1,out2,trange2,savehrarch] = reduced_HEOM_operator(mu_ge{j},...
    t_end, H_ex_vib_sc,fock_space_rep,QQ,cc_com,cc_acom,vv,to_pass1,Kap2,numpoints);
if j==1
    plot(trange1,real(out1))
    hold on
    plot(trange2,real(out2))
    figure
    plot(trange1,imag(out1))
    hold on
    plot(trange2,imag(out2))
end
%savehrarch = {nn,coup_com_save,coup_acom_save,const_factor,total_prop_op1,total_prop_op2};
if ~exist('nn','var')
nn = savehrarch{1};
end
if ~iscell(to_pass1)
to_pass1 = savehrarch; %no point calculating this every time!  Pass back
end
%interpolate to give desired linearly spaced time range and reintroduce the
%phase factor that was removed
out1 = interp1(trange1,out1,des_t_range)...
    .*kron(exp(-1i*H_ex_vib(j+1,j+1)*des_t_range),ones(1,size(nn,1)*N));
out2 = interp1(trange2,out2,des_t_range)...
    .*kron(exp(+1i*H_ex_vib(j+1,j+1)*des_t_range),ones(1,size(nn,1)*N));
% covert to full size in ALL ITS GLORY!
temp = zeros(length(des_t_range),size(nn,1)*length(H_ex_vib)^2);

temp1 = false(length(H_ex_vib));   temp1(2:N+1,1) = true; 
temp1 = reshape(temp1,length(H_ex_vib)^2,1);
temp1 = repmat(temp1,size(nn,1),1); 
temp2 = false(length(H_ex_vib));   temp2(1,2:N+1) = true; 
temp2 = reshape(temp2,length(H_ex_vib)^2,1);
temp2 = repmat(temp2,size(nn,1),1); 

temp(:,temp1) = out1; temp(:,temp2) = out2;
mu_ge_t{j} = sparse(temp.');
% 
% pad1 = kron(1:length(des_t_range),ones(1,length(pdmat)));
% spmap1 = kron(spmap1,ones(1,length(des_t_range)));
% spmap2 = kron(spmap2,ones(1,length(des_t_range)));
% mu_ge_t{j} = sparse(spmap1,pad1,out1,numel(H_ex_vib)*size(nn,1),length(des_t_range))...
%     + sparse(spmap2,pad1,out2,numel(H_ex_vib)*size(nn,1),length(des_t_range));
%stack of H_ex_vib sized matricies all one top of one another, in Louiville
%space format.  

%actually just stack each H_ex_vib sized matrix on top of each other in
%full format for now...


% now the same for the single - double coherences
if ~iscell(to_pass2)
to_pass2 = savehrarch(1:4); %these will be the same but total prop won't be
end
    H_ex_vib_sc = H_ex_vib;  % scale elements by transition frequency to simplify numerics
    H_ex_vib_sc(N+2:end,N+2:end) = H_ex_vib_sc(N+2:end,N+2:end) -diag(...
                                    H_ex_vib(j+1,j+1)*ones(N*(N-1)/2,1));
[out1,trange1,out2,trange2,savehrarch] = reduced_HEOM_operator2(mu_ef{j} ,...
    t_end, H_ex_vib_sc,fock_space_rep,QQ,cc_com,cc_acom,vv,to_pass2,Kap2,numpoints);
%interpolate to give desired linearly spaced time range, note that in this
%case the out1 and out2 and in fact 2:N+1 by N(N-1)/2 matricies
if length(to_pass2)==4
to_pass2 = savehrarch; %pass propogators back
end
%diag_en = diag(H_ex_vib(2:N+1,2:N+1)); diag_en = kron(diag_en,size(nn,1))
%[EE,tt] = meshgrid(diag_en,des_trange);
out1 = interp1(trange1,out1,des_t_range)...
    .*kron(exp(1i*H_ex_vib(j+1,j+1)*des_t_range),ones(1,size(nn,1)*N*(N-1)));
out2 = interp1(trange2,out2,des_t_range)...
    .*kron(exp(1i*H_ex_vib(j+1,j+1)*des_t_range),ones(1,size(nn,1)*N*(N-1)));
% covert to full size in ALL ITS GLORY!
% 
% pdmat = ((1:size(nn,1))-1)*length(H_ex_vib)^2;
% pdmat = kron(pdmat,ones(1,N)).';
% temp = zeros(length(H_ex_vib));   temp(N+2:end,2:N+1) = 1; 
% temp = find(temp);  spmap1 = repmat(temp,size(nn,1),1) + pdmat;
% temp = zeros(length(H_ex_vib));   temp(2:N+1,N+2:end) = 1; 
% temp = find(temp);  spmap2 = repmat(temp,size(nn,1),1) + pdmat;
% 
% pad1 = kron(1:length(des_t_range),ones(1,length(pdmat)));
% spmap1 = kron(spmap1,ones(1,length(des_t_range)));
% spmap2 = kron(spmap2,ones(1,length(des_t_range)));
% mu_ef_t{j}= sparse(spmap1,pad1,out1,numel(H_ex_vib)*size(nn,1),length(des_t_range))...
%     + sparse(spmap2,pad1,out2,numel(H_ex_vib)*size(nn,1),length(des_t_range));
temp = zeros(length(des_t_range),size(nn,1)*length(H_ex_vib)^2);

temp1 = false(length(H_ex_vib));   temp1(N+2:end,2:N+1) = true; 
temp1 = reshape(temp1,length(H_ex_vib)^2,1);
temp1 = repmat(temp1,size(nn,1),1); 
temp2 = false(length(H_ex_vib));   temp2(2:N+1,N+2:end) = true; 
temp2 = reshape(temp2,length(H_ex_vib)^2,1);
temp2 = repmat(temp2,size(nn,1),1); 

temp(:,temp1) = out1; temp(:,temp2) = out2;
mu_ef_t{j} = sparse(temp.'); %time in dimension 2


%suppose I could just add these things together
end
clear temp
%% Calculate first order response function
%rho_eq = H_ex_vib*0; rho_eq(1,1) = 1;
f_leng = length(H_ex_vib)^2; h_leng = length(H_ex_vib);
rho_eq = zeros(size(H_ex_vib)); rho_eq(1,1) = 1;
JJ_t = zeros(length(des_t_range),2); 
JJ_t_tmp = zeros(length(des_t_range),1);
tmp = diag(H_ex_vib); k_pr = 2*pi*mean(tmp(2:N+1)); %2*pi*av freq
for j1 = 1:N
    init_s = (mu_ge{j1}*rho_eq);    
    for j2 = 1:N    
        xx_av = dot(mu(j1,:),mu(j2,:))/3;
        yx_av = 1i*k_pr*dot(mu(j1,:),cross(mu(j2,:),R(j1,:)-R(j2,:)))/6;
        for tlp = 1:length(des_t_range)
          
JJ_t_tmp(tlp) = trace(full(reshape(mu_ge_t{j2}...
            (1:f_leng,tlp),h_leng,h_leng))*init_s); 
%will need this same thing for third as well
        end
  
     JJ_t(:,1) = JJ_t(:,1) +  JJ_t_tmp.*xx_av;
     JJ_t(:,2) =  JJ_t(:,2) + JJ_t_tmp.*yx_av; 
    end
    
end

S_1 = 1i*(JJ_t - conj(JJ_t)); %S_1_FT = fft(S_1,[],1);  
S_1_FT = -1i*fft([S_1*0;S_1],[],1);  
plot(real(S_1_FT))
figure
plot(imag(S_1_FT))

%% Calculate third order response function
%this will be a load slower
rho_eq = H_ex_vib*0; rho_eq(1,1) = 1;
R1_t_tmp = zeros(length(t1_rng),length(t2_rng),length(t3_rng));
R2_t_tmp = R1_t_tmp;   R3_t_tmp = R1_t_tmp;  R4_t_tmp = R1_t_tmp;



k_I = k_3+k_2-k_1; %rephasing
k_II = k_3-k_2+k_1; %non-rephasing
k_III = -k_3+k_2+k_1; %coherence

%preallocate all required response function, comment out any that are not
%required, e.g. for the rephasing signal in x only comment out all terms 
% but S_3I_x here and in the big loop section
S_3I_x = R2_t_tmp; S_3I_y = R2_t_tmp;S_3I_z = R2_t_tmp;
S_3II_x = R2_t_tmp; S_3II_y = R2_t_tmp;S_3II_z = R2_t_tmp;
S_3III_x = R2_t_tmp; S_3III_y = R2_t_tmp;S_3III_z = R2_t_tmp;

for j1 = 1:N
    %init_s = (mu_ge{j1}*rho_eq);    
    for j2 = 1:N 
        for j3 = 1:N
            for j4 = 1:N %this is the final interaction
         %calculate dipole moment averages      
        xxxx_av = (1/15)*(dot(mu(j1),mu(j2))*dot(mu(j3),mu(j4))+...
                dot(mu(j1),mu(j3))*dot(mu(j2),mu(j4))+...
                dot(mu(j2),mu(j3))*dot(mu(j1),mu(j4))); %also yyyy etc
        
        yyxx_av = (1/30)*(4*dot(mu(j1),mu(j2))*dot(mu(j3),mu(j4))-...
                dot(mu(j1),mu(j3))*dot(mu(j2),mu(j4))-...
                dot(mu(j2),mu(j3))*dot(mu(j1),mu(j4)));    %also xxzz etc
        xyyx_av = (1/30)*(-dot(mu(j1),mu(j2))*dot(mu(j3),mu(j4))+4*...
                dot(mu(j1),mu(j3))*dot(mu(j2),mu(j4))-...
                dot(mu(j2),mu(j3))*dot(mu(j1),mu(j4)));               
        yxyx_av = (1/30)*(-dot(mu(j1),mu(j2))*dot(mu(j3),mu(j4))-...
                dot(mu(j1),mu(j3))*dot(mu(j2),mu(j4))+4*...
                dot(mu(j2),mu(j3))*dot(mu(j1),mu(j4)));             
        %only non zero 4 comp averages with last index x
         
        %now xxxyz averages, also yyyzx and zzzxy, this is very general,
        %possibly overkill but takes little computational time so is
        %included, works 9or should at least) for any input beam 
        %configuration and polarizations
        % diff one for each
        xxxyz_av(1) = (1i/30)*...
            (dot(mu(j1),mu(j2))*dot(mu(j3),cross(mu(j4),R(j1)))+...
                dot(mu(j1),mu(j3))*dot(mu(j2),cross(mu(j4),R(j1)))+...
                dot(mu(j2),mu(j3))*dot(mu(j1),cross(mu(j4),R(j1))));
        xxxyz_av(2) = (1i/30)*...
            (dot(mu(j1),mu(j2))*dot(mu(j3),cross(mu(j4),R(j2)))+...
                dot(mu(j1),mu(j3))*dot(mu(j2),cross(mu(j4),R(j2)))+...
                dot(mu(j2),mu(j3))*dot(mu(j1),cross(mu(j4),R(j2))));        
        xxxyz_av(3) = (1i/30)*...
            (dot(mu(j1),mu(j2))*dot(mu(j3),cross(mu(j4),R(j3)))+...
                dot(mu(j1),mu(j3))*dot(mu(j2),cross(mu(j4),R(j3)))+...
                dot(mu(j2),mu(j3))*dot(mu(j1),cross(mu(j4),R(j3))));        
        xxxyz_av(4) = (1i/30)*...
            (dot(mu(j1),mu(j2))*dot(mu(j3),cross(mu(j4),R(j4)))+...
                dot(mu(j1),mu(j3))*dot(mu(j2),cross(mu(j4),R(j4)))+...
                dot(mu(j2),mu(j3))*dot(mu(j1),cross(mu(j4),R(j4))));
     %this bit would be better to call as a C function or similar             
        for tlp1 = 1:t1_rng  %maybe best to have one 3D time range
            for tlp2 = 1:t2_rng
                for tlp3  = 1:t3_rng
                    
V1 = mu_ge_t{j1}(:,:,1) + mu_ef_t{j1}(:,:,1);
V2 = mu_ge_t{j2}(:,:,tlp1) + mu_ef_t{j2}(:,:,tlp1);
V3 = mu_ge_t{j3}(:,:,tlp1+tlp2) + mu_ef_t{j3}(:,:,tlp1+tlp2);   
V4 = mu_ge_t{j4}(:,:,tlp1+tlp2+tlp3) + mu_ef_t{j4}(:,:,tlp1+tlp2+tlp3);
                    
R1_t_tmp(tlp1,tlp2,tlp3) = trace(V2*V3*V4*V1*rho_eq);
R2_t_tmp(tlp1,tlp2,tlp3) = trace(V1*V3*V4*V2*rho_eq); 
R3_t_tmp(tlp1,tlp2,tlp3) = trace(V1*V2*V4*V3*rho_eq); 
R4_t_tmp(tlp1,tlp2,tlp3) = trace(V4*V3*V2*V1*rho_eq); 
%will need this same thing for third as well
                end
            end
        end
        tmp  = 2*imag(R1_t_tmp+R2_t_tmp+R3_t_tmp+R4_t_tmp);
        %if any additional polarizations aren't actually needed just
        %comment out these and don't preallocate them
        [out_x,out_y,out_z] = full_cont_3rd_order(xxxx_av,yyxx_av,...
                xyyx_av,yxyx_av ,xxxyz_av,pol,k_1,k_2,k_3,k_I);
            %this function works out all the contributions to the first
            %order response function in each direction for each combination
            %of wavevectors and polarisations etc etc
     S_3I_x = S_3I_x + tmp*out_x;
     S_3I_y = S_3I_y + tmp*out_y;
     S_3I_z = S_3I_z + tmp*out_z;
     %also e^(i*(w_3*
                %<--->_zzzxy e3_z e2_z e1_z es_x k_y
        [out_x,out_y,out_z] = full_cont_3rd_order(xxxx_av,yyxx_av,...
                xyyx_av,yxyx_av ,xxxyz_av,pol,k_1,k_2,k_3,k_II);
     S_3II_x = S_3II_x + tmp*out_x;
     S_3II_y = S_3II_y + tmp*out_y;
     S_3II_z = S_3II_z + tmp*out_z;
     
        [out_x,out_y,out_z] = full_cont_3rd_order(xxxx_av,yyxx_av,...
                xyyx_av,yxyx_av ,xxxyz_av,pol,k_1,k_2,k_3,k_III);
     S_3III_x = S_3III_x + tmp*out_x;
     S_3III_y = S_3III_y + tmp*out_y;
     S_3III_z = S_3III_z + tmp*out_z;     
  
            end
        end
    end
end

% do fourier transforms to get out the appropriate components

%rescale by some central frequency to recentre the fourier transforms