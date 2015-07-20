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
Kappa = 1; %# matsubara frequencies
Kap2 = 5; %Tier truncation
tier_save = 0;
Kap1 = inf;%1.5e+03; %high freq cut off (inv cm)
% other parameters relating to how long evolution is considered for
t_end_ps = 4;
t_end = t_end_ps*convfact;

numpoints = 2^12; %might as well make it a power of 2 for the fft

des_t_range = zeros(numpoints,1); 
des_t_range(:) = linspace(0,t_end_ps,numpoints); %value in ps
des_t_range = des_t_range*convfact; %convert to inverse cm from ps
rel_freq_range = linspace(0,2*pi/(des_t_range(2)-des_t_range(1)),numpoints);

%% System properties, independent of orientation

flenme ='B820_bacteriochlorophyll_dimer_params.mat'; 
%insert name of file containing the system parameters
    fle = open(flenme);
        H_site = fle.H_site;
delta_w = 0*[-200,200]; %static disorder can be added in this way, shifts diag
      H_site = H_site + diag(delta_w);
        
mu = fle.mu;  pdm = fle.pdm; %set to empty if not known
%mu(2,:) = mu(1,:); %test with them the same
R = fle.R;  R = R*1e-8;   %convert units from angstrom to cm^-1
if isempty(pdm)
    pdm_n_pos = [];
else
   pdm_n_pos = [pdm;R]; 
end
%read out bath parameters
lam_dru  = fle.lam_dru; gam_dru = fle.gam_dru;
%lam_dru = {250,250}; %gam_dru = {40,40}; %stronger Drude
%underdamped
om_0 = fle.om_0; lambda = fle.lambda;  gamma = fle.gamma;
%explicitly included modes
om_vib = fle.om_vib;  numvib = fle.numvib; %set this to whatever
displ = fle.om_vib;
  %lambda = {300,300};   
%subtract reorganisation energy
%H_site = H_site + diag(cellfun(@sum, lam_dru) + cellfun(@sum, lambda));

[H_ex_vib,fock_space_rep,mu_full,H_exciton,indiv_op] = generate_ex_vib_ham...
    (H_site,om_vib,numvib,displ,sqrt(sum(mu.^2,2)),pdm_n_pos) ;

[ex_basis,H0ex] = eig(H_ex_vib);
H_vib = indiv_op{2};

%% Beam properties of the 3 input beams
%unit vectors in direction
k1 = [1/sqrt(2);0;1/sqrt(2)]; k2 = k1; %pp %k2 = [-1/sqrt(2);0;1/sqrt(2)];
k3 = [0;0;1];
%polarization wavevectors
pol{1} = [1;0;0]; pol{2} = [1;0;0]; pol{3} = [1;0;0];  

%% Next Calculate time dependent dipole operators via HEOM

N = length(H_site); %number of sites

if isempty(H_vib);
    H_vib = 0;
end
%lam_dru = {10,10};  %test without much bath
%lambda = {eps(100),eps(100)};  
%lambda = {500,500};  %bit silly coupling to vib to see strong effect
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
mu_ge{j} = mu_ge{j}-diag(diag(mu_ge{j}));

%me_ge_ex{j} =ex_basis'*mu_ge{j}*ex_basis;

end
%now the coherences to the double exciton statesc
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
mu_ef{j} = temp+temp'; 
mu_ef{j} = mu_ef{j}-diag(diag(mu_ef{j}));     

 end   
%project the single site coherences to single exciton coherences
%  for j = 1:N
%  mu_g_1ex{j}  = ex_basis'*mu_ge{j}*ex_basis;
%  end
%% evolve in time with HEOM
to_pass1 = Kap1; to_pass2 = Kap1; 
% [total_prop_op,nn]=HEOM_op_gen(H_ex_vib,fock_space_rep,...
%                         QQ,cc_com,cc_acom,vv,Kap1,Kap2,viblvls);

for j=1:N
    
    H_ex_vib_sc = H_ex_vib;  % scale elements by transition frequency to simplify numerics
    H_ex_vib_sc(2:N+1,2:N+1) = H_ex_vib_sc(2:N+1,2:N+1) - ...
                            diag(H_ex_vib(j+1,j+1)*ones(N,1));
  %compute in exciton basis and then sum these terms to get the ones in
  %site basis by linearity, this makes numerics easier.
             ex_mu = ex_basis'*mu_ge{j}*ex_basis;      
     tic
[out1,trange1,out2,trange2,savehrarch] = reduced_HEOM_operator2(ex_mu,...
    t_end, H_ex_vib_sc,fock_space_rep,QQ,cc_com,cc_acom,vv,to_pass1,Kap2,numpoints,tier_save);
    toc
 
% t_sym=[flipud(out1(2:end,1));out1(:,1)];
% tmp=(0:length(t_sym)-1);
% phase_shift=exp(-1i*2*pi*(length(t_sym)-1)/2*tmp/length(t_sym));
   
%out1= out1*ex_basis(2:N+1,2:N+1); %put back in exciton basis
%out2= out2*ex_basis(2:N+1,2:N+1); %same for this one due to the way I have orientated them
%this allows for easier fitting

%savehrarch = {nn,coup_com_save,coup_acom_save,const_factor,total_prop_op1,total_prop_op2};
%size(savehrarch{1})
if ~exist('nn','var')
nn = savehrarch{1};
nn_red = nn(sum(nn,2)<=tier_save+10*eps);
end
% if ~iscell(to_pass1)
% to_pass1 = savehrarch; %no point calculating this every time!  Pass back
% end %this is currently breaking the code and I can't work out why...

%interpolate to give desired linearly spaced time range and reintroduce the
%phase factor that was removed
 mu_ge_t{j} = interp1(trange1,out1,des_t_range)...
   .*repmat(exp(1i*H_ex_vib(j+1,j+1)*des_t_range),1,size(nn_red,1)*N);
%out2 = interp1(trange2,out2,des_t_range);%...
     %.*repmat(exp(-1i*H_ex_vib(j+1,j+1)*des_t_range),1,size(nn_red,1)*N);
%  out1_ft = fftshift(fft(out1)); %fit splines to these?
  % asymtotically 1/(w-w0)^2, closer in normal spline
 % out2_ft = fftshift(fft(out2));

% now the same for the single - double coherences
%if 1==0
% if ~iscell(to_pass2)
% to_pass2 = savehrarch(1:4); %these will be the same but total prop won't be
% end
ex_mu = ex_basis'*mu_ef{j}*ex_basis;
    H_ex_vib_sc = H_ex_vib;  % scale elements by transition frequency to simplify numerics
    H_ex_vib_sc(N+2:end,N+2:end) = H_ex_vib_sc(N+2:end,N+2:end) -diag(...
                                    H_ex_vib(j+1,j+1)*ones(N*(N-1)/2,1));
    H_ex_vib_sc(2:end,2:end) =H_ex_vib_sc(2:end,2:end) -diag(...
                                    H_ex_vib(j+1,j+1)*ones(N*(N+1)/2,1));                             
tic
  [out1,trange1,out2,trange2,savehrarch] = reduced_HEOM_operator2(ex_mu ,...
    t_end, H_ex_vib_sc,fock_space_rep,QQ,cc_com,cc_acom,vv,to_pass2,Kap2,numpoints,tier_save);

toc
%interpolate to give desired linearly spaced time range, note that in this
%case the out1 and out2 and in fact 2:N+1 by N(N-1)/2 matricies
% if length(to_pass2)==4
% to_pass2 = savehrarch; %pass propogators back
% end
%out1= out1*ex_basis(2:N+1,2:N+1);
mu_ef_t{j} = interp1(trange1,out1,des_t_range)...
     .*repmat(exp(1i*H_ex_vib(j+1,j+1)*des_t_range),1,size(nn_red,1)*N*(N-1));
%out2 = interp1(trange2,out2,des_t_range)...
%     .*repmat(exp(-1i*H_ex_vib(j+1,j+1)*des_t_range),1,size(nn_red,1)*N*(N-1));


%end
%suppose I could just add these things together
end
clear temp out1 out2 %these actually use quite a bit of memory
%% Calculate first order response function
%rho_eq = H_ex_vib*0; rho_eq(1,1) = 1;

rho_eq = zeros(1+N+N*(N-1)/2); rho_eq(1,1) = 1;
A1 = sym('Fge', [1 N]); 
A2 = sym('Gge', [1 N]);
V1 = sym(zeros(N*(N+1)/2+1)); 
V1(1,2:N+1) = A1; V1=V1+V1';
V2 = sym(zeros(N*(N+1)/2+1)); 
V2(1,2:N+1) = A2; V2=V2+V2';
tmp= 2*imag(trace(V2*V1*rho_eq));
             
first_order_fn = matlabFunction(tmp,'vars',{A1,A2});


f_leng = length(H_ex_vib)^2; h_leng = length(H_ex_vib);
S_1  = zeros(length(des_t_range),2); 

mu = fle.mu;  
mu(1,:) = mu(2,:);
tmp = diag(H_ex_vib); scale_freq = mean(tmp(2:N+1));
k_pr = 2*pi*scale_freq; %2*pi*av freq
%scale_freq =0*scale_freq ;
%ft_scale_shift = exp(-1i*scale_freq*des_t_range);
ft_scale_shift =exp(-1i*0*des_t_range);
tic
for j1 = 1:N

    A1 = ex_basis(2,j1+1)*mu_ge_t{1}(1,:);  
            %time dependent operators are expressed in the exciton basis,
        %convert them back to site basis by taking superpositions
            for ell = 2:N
                A1 = A1 + ex_basis(1+ell,j1+1)*mu_ge_t{ell}(1,:);  
            end
    for j2 = 1:N    

        A2 = ex_basis(2,j2+1)*mu_ge_t{1};
            for ell = 2:N
                A2 = A2 + ex_basis(1+ell,j2+1)*mu_ge_t{ell};  
            end
        
        xx_av = dot(mu(j1,:),mu(j2,:))/3;
        yx_av = 1i*k_pr*dot(mu(j1,:),cross(mu(j2,:),R(j1,:)-R(j2,:)))/6;
          
        tmp = first_order_fn(A1,A2);
  
     S_1(:,1) = S_1(:,1) + ft_scale_shift.*tmp.*xx_av;
     S_1(:,2) = S_1(:,2) + ft_scale_shift.*tmp.*yx_av; 
    end
    
end

S_1_FT = fft(S_1,[],1);  
lambda_rng = 10^7./(rel_freq_range(1:end/2)-scale_freq*0);
red_rng = lambda_rng>700 & lambda_rng<950; %range of interest
% brodnin = 100;
% gauss_broad = exp(-des_t_range.^2*brodnin^2/2)/sqrt(sqrt(2*pi)*brodnin);
% S_1_FT_broad = fft(S_1.*repmat(gauss_broad,1,2) ,[],1);  
toc
%% Plot things
if 1==1

%S_1_FT =S_1_FT(1:end/2,:);
%note S_1 is real so S_1_FT(N-k,:) = conj(S_1_FT(N-k,:))
%S_1_FT = -1i*fft([S_1*0;S_1],[],1);  
figure %ref index change and absorb
plot(lambda_rng(red_rng),[real(S_1_FT(red_rng,1)),-imag(S_1_FT(red_rng,1))])
%figure %OR and CD
%plot(lambda_rng(red_rng),[real(S_1_FT(red_rng,2)),imag(S_1_FT(red_rng,2))])

end

%% Calculate third order response function

%ranges for time will be somewhat different 
t1_rng = (1:5:length(des_t_range)/3);
t2_rng = t1_rng; t3_rng = t1_rng;
%R1_t_tmp = zeros(length(t1_rng),length(t2_rng),length(t3_rng));
%R2_t_tmp = R1_t_tmp;   R3_t_tmp = R1_t_tmp;  R4_t_tmp = R1_t_tmp;
Stmp= zeros(length(t1_rng),length(t2_rng),length(t3_rng));

k_1 = scale_freq*k1; k_2 = scale_freq*k2; k_3 = scale_freq*k3; 

k_I = k_3+k_2-k_1; %rephasing
k_II = k_3-k_2+k_1; %non-rephasing
if sum(abs(k_II-k_I)) < 10*eps(norm(k_II)+norm(k_I)) %pp setup
    pp_setup = true;
else
    pp_setup = false;
end
k_III = -k_3+k_2+k_1; %coherence

%preallocate all required response function, comment out any that are not
%required, e.g. for the rephasing signal in x only comment out all terms 
% but S_3I_x here and in the big loop section
%if pp_setup %only 1 signal of interest
%  S_3_x   = zeros(length(t1_rng),length(t2_rng),length(t3_rng)); 
%  S_3_y   = zeros(length(t1_rng),length(t2_rng),length(t3_rng)); 
%else
S_3I_x = zeros(length(t1_rng),length(t2_rng),length(t3_rng)); 
S_3I_y = S_3I_x;%S_3I_z = S_3I_x;
S_3II_x = S_3I_x; S_3II_y = S_3I_x;%S_3II_z = S_3I_x;
%S_3III_x = S_3I_x; S_3III_y = S_3I_x;S_3III_z = S_3I_x;
%end
%analytically calculate expression for a given contribution
rho_eq = zeros(1+N+N*(N-1)/2); rho_eq(1,1) = 1;
A1 = sym('Fge', [1 N]); B1= sym('Fef', [N,N*(N-1)/2 ]) ;
A2 = sym('Gge', [1 N]); B2= sym('Gef', [N,N*(N-1)/2 ]) ;
A3 = sym('Hge', [1 N]); B3= sym('Hef', [N,N*(N-1)/2 ]) ;
A4 = sym('Ige', [1 N]); B4= sym('Ief', [N,N*(N-1)/2 ]) ;
V1 = sym(zeros(N*(N+1)/2+1)); 
V1(1,2:N+1) = A1; V1(2:N+1,N+2:end) = B1; V1=V1+V1';
V2 = sym(zeros(N*(N+1)/2+1)); 
V2(1,2:N+1) = A2; V2(2:N+1,N+2:end) = B2; V2=V2+V2';
V3 = sym(zeros(N*(N+1)/2+1)); 
V3(1,2:N+1) = A3; V3(2:N+1,N+2:end) = B3;  V3=V3+V3';
V4 = sym(zeros(N*(N+1)/2+1));
V4(1,2:N+1) = A4; V4(2:N+1,N+2:end) = B4;  V4=V4+V4';

tmp= 2*imag(trace(V2*V3*V4*V1*rho_eq + V1*V3*V4*V2*rho_eq+...
                  V1*V2*V4*V3*rho_eq + V4*V3*V2*V1*rho_eq));
              %reshape into lioville space shapes which I store them as
B1 = reshape(B1,1,N^2*(N-1)/2); B2 = reshape(B2,1,N^2*(N-1)/2);  
B3 = reshape(B3,1,N^2*(N-1)/2);  B4 = reshape(B4,1,N^2*(N-1)/2);  
third_order_fn = matlabFunction(tmp,'vars',{A1,A2,A3,A4,B1,B2,B3,B4});
%sub in values of mu_ge_t and mu_ef_t corresponding to the appropriate time
%position.
%% Actual calculation
%this will be a load slower
for j1 = 1:N
        A1 = ex_basis(2,j1+1)*mu_ge_t{1}(1,:);  
        B1 = ex_basis(2,j1+1)*mu_ef_t{1}(1,:);
            %time dependent operators are expressed in the exciton basis,
        %convert them back to site basis by taking superpositions
         for ell = 2:N
                A1 = A1 + ex_basis(1+ell,j1+1)*mu_ge_t{ell}(1,:);  
        B1 = B1+ ex_basis(1+ell,j1+1)*mu_ef_t{ell}(1,:);        
         end
   %only at t=0
    for j2 = 1:N 
        A2 = ex_basis(2,j2+1)*mu_ge_t{1};  
        B2 = ex_basis(2,j2+1)*mu_ef_t{1};
         for ell = 2:N
                A2 = A2 + ex_basis(1+ell,j2+1)*mu_ge_t{ell};  
        B2 = B2+ ex_basis(1+ell,j2+1)*mu_ef_t{ell};        
         end
        for j3 = 1:N   
        A3 = ex_basis(2,j3+1)*mu_ge_t{1};  
        B3 = ex_basis(2,j3+1)*mu_ef_t{1};
         for ell = 2:N
                A3 = A3 + ex_basis(1+ell,j3+1)*mu_ge_t{ell};  
        B3 = B3+ ex_basis(1+ell,j3+1)*mu_ef_t{ell};
         end
            for j4 = 1:N %this is the final interaction
        A4 = ex_basis(2,j4+1)*mu_ge_t{1};  
        B4 =ex_basis(2,j4+1)*mu_ef_t{1};
         for ell = 2:N
                A4 = A4 + ex_basis(1+ell,j4+1)*mu_ge_t{ell};  
        B4 = B4+ ex_basis(1+ell,j4+1)*mu_ef_t{ell}; 
         end
         %calculate dipole moment averages      
   xxxx_av = (1/15)*(dot(mu(j1,:),mu(j2,:))*dot(mu(j3,:),mu(j4,:))+...
                dot(mu(j1,:),mu(j3,:))*dot(mu(j2,:),mu(j4,:))+...
                dot(mu(j2,:),mu(j3,:))*dot(mu(j1,:),mu(j4,:))); %also yyyy etc
        %only non zero 4 comp averages with last index x
                 
        yyxx_av = (1/30)*(4*dot(mu(j1,:),mu(j2,:))*dot(mu(j3,:),mu(j4,:))-...
                dot(mu(j1,:),mu(j3,:))*dot(mu(j2,:),mu(j4,:))-...
                dot(mu(j2,:),mu(j3,:))*dot(mu(j1,:),mu(j4,:)));    %also xxzz etc
        xyyx_av = (1/30)*(-dot(mu(j1,:),mu(j2,:))*dot(mu(j3,:),mu(j4,:))+4*...
                dot(mu(j1,:),mu(j3,:))*dot(mu(j2,:),mu(j4,:))-...
                dot(mu(j2,:),mu(j3,:))*dot(mu(j1,:),mu(j4,:)));               
        yxyx_av = (1/30)*(-dot(mu(j1,:),mu(j2,:))*dot(mu(j3,:),mu(j4,:))-...
                dot(mu(j1,:),mu(j3,:))*dot(mu(j2,:),mu(j4,:))+4*...
                dot(mu(j2,:),mu(j3,:))*dot(mu(j1,:),mu(j4,:)));           
            
            %now xxxyz averages, also yyyzx and zzzxy, this is very general,
        %possibly overkill but takes little computational time so is
        %included, works 9or should at least) for any input beam 
        %configuration and polarizations
        % diff one for each
        
        xxxyz_av(1) = (1i/30)*...
            (dot(mu(j1,:),mu(j2,:))*dot(mu(j3,:),cross(mu(j4,:),R(j1,:)))+...
             dot(mu(j1,:),mu(j3,:))*dot(mu(j2,:),cross(mu(j4,:),R(j1,:)))+...
             dot(mu(j2,:),mu(j3,:))*dot(mu(j1,:),cross(mu(j4,:),R(j1,:))));
        xxxyz_av(2) = (1i/30)*...
            (dot(mu(j1,:),mu(j2,:))*dot(mu(j3,:),cross(mu(j4,:),R(j2,:)))+...
             dot(mu(j1,:),mu(j3,:))*dot(mu(j2,:),cross(mu(j4,:),R(j2,:)))+...
             dot(mu(j2,:),mu(j3,:))*dot(mu(j1,:),cross(mu(j4,:),R(j2,:))));        
        xxxyz_av(3) = (1i/30)*...
            (dot(mu(j1,:),mu(j2,:))*dot(mu(j3,:),cross(mu(j4,:),R(j3,:)))+...
             dot(mu(j1,:),mu(j3,:))*dot(mu(j2,:),cross(mu(j4,:),R(j3,:)))+...
             dot(mu(j2,:),mu(j3,:))*dot(mu(j1,:),cross(mu(j4,:),R(j3,:))));        
        xxxyz_av(4) = (1i/30)*...
            (dot(mu(j1,:),mu(j2,:))*dot(mu(j3,:),cross(mu(j4,:),R(j4,:)))+...
             dot(mu(j1,:),mu(j3,:))*dot(mu(j2,:),cross(mu(j4,:),R(j4,:)))+...
             dot(mu(j2,:),mu(j3,:))*dot(mu(j1,:),cross(mu(j4,:),R(j4,:))));
     %this bit would be better to call as a C function or similar                   
      %Matlab will not like this one little bit
      
      
      tic            
            for tlp2 = 1:length(t2_rng)    
                
                for tlp3  = 1:length(t3_rng)
             
tmp = third_order_fn(A1,A2(t1_rng,:),A3(t1_rng+tlp2,:),A4(t1_rng+tlp2+tlp3,:),...
                     B1,B2(t1_rng,:),B3(t1_rng+tlp2,:),B4(t1_rng+tlp2+tlp3,:));
Stmp(:,tlp2,tlp3) = tmp ;
                end
            end
        toc

        %if any additional polarizations aren't actually needed just
        %comment out these and don't preallocate them
%        if pp_setup
%            
%         [out_x,out_y,out_z] = full_cont_3rd_order(xxxx_av,yyxx_av,...
%                 xyyx_av,yxyx_av ,xxxyz_av,pol,k_1,k_2,k_3,k_I);
%      S_3_x = S_3_x + Stmp*out_x;
%      S_3_y = S_3_y + Stmp*out_y;  
%         [out_x,out_y,out_z] = full_cont_3rd_order(xxxx_av,yyxx_av,...
%                 xyyx_av,yxyx_av ,xxxyz_av,pol,k_1,k_2,k_3,k_II); 
%      S_3_x = S_3_x + Stmp*out_x;
%      S_3_y = S_3_y + Stmp*out_y;             
%             
%         else
        [out_x,out_y,out_z] = full_cont_3rd_order(xxxx_av,yyxx_av,...
                xyyx_av,yxyx_av ,xxxyz_av,pol,-k_1,k_2,k_3,k_I);
            %this function works out all the contributions to the first
            %order response function in each direction for each combination
            %of wavevectors and polarisations etc etc
     S_3I_x = S_3I_x + Stmp*out_x;
     S_3I_y = S_3I_y + Stmp*out_y;
     %S_3I_z = S_3I_z + Stmp*out_z;
     %also e^(i*(w_3*
                %<--->_zzzxy e3_z e2_z e1_z es_x k_y
        [out_x,out_y,out_z] = full_cont_3rd_order(xxxx_av,yyxx_av,...
                xyyx_av,yxyx_av ,xxxyz_av,pol,k_1,-k_2,k_3,k_II);
     S_3II_x = S_3II_x + Stmp*out_x;
     S_3II_y = S_3II_y + Stmp*out_y;
     %S_3II_z = S_3II_z + Stmp*out_z;
     
%         [out_x,out_y,out_z] = full_cont_3rd_order(xxxx_av,yyxx_av,...
%                 xyyx_av,yxyx_av ,xxxyz_av,pol,k_1,k_2,k_3,k_III);
%      S_3III_x = S_3III_x + Stmp*out_x;
%      S_3III_y = S_3III_y + Stmp*out_y;
%      S_3III_z = S_3III_z + Stmp*out_z;     
    %    end
            end
        end
    end
end
%% Analysis
        if pp_setup
            
        end
% do fourier transforms to get out the appropriate components

%rescale by some central frequency to recentre the fourier transforms