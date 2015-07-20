%1 cm?1 = 0.000123986 eV 
%h = 2pi*c = 0.01m = 1 codification
 speed_unit = 2*pi * 299792458; %units m s^{-1}, NOT cm 
 length_unit = 0.01; % units m
 hbar = 6.62606957*10^(-34) ;  %units Js
 boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}
 Temp = 300; %in kelvin
 beta = ( hbar * speed_unit / length_unit)/ ( Temp * boltz_const); %in cm^-1
% to run the code not as a function just uncomment and use each module
% Energy in Joules,  extra factor of 10 is because thats speed in cm/s
convfact = speed_unit / length_unit  *10^(-12); %2 pi * "c" * 100m^-1 * 10^-12 s 
%probe unit vectors

kpr = [0,0,1];   epr = [1,0,0];   bpr = [0,1,0];
%pump unit vectors

theta = atan(sqrt(2)); %angle between pump and probe
%note arctan(sqrt(2)) is not the magic angle from the cho 
%paper but the nmr magic angle so it will not remove quadrupole terms
kpu = [sin(theta),0,cos(theta)];  epu = [1,0,0]; %same as  probe
%might be easier to consider it just polarised along x
bpu = cross(kpu,epu);

%PC645 dimer

%dimer specific parameters

    fle = open('Hamiltonian_save.mat');
        H_site = fle.PC645Hamiltonian;
            E0 = fle.E0_PC645; %these energies are relative 
            %to the lowest exciton state anyway
       H_site = H_site(1:2,1:2);%only these 2
       temp = min(diag(H_site));
       H_site = H_site - eye(2)*temp;
       E0 = E0-temp; %set lowest exciton included to zero energy
      % H_site = diag(diag(H_site)) + fle.PC645_el_coup([4,8],[4,8]);
       %don't take newer version of electronic coupling 
       %that is much bigger for some reason.  Uber excitons!
H0 = zeros(length(H_site)+1); H0(2:end,2:end) = H_site;
H0(1) = +E0; %E0 is negative here

mu = fle.PC645_dip_and_pos([4,8],1:3); %4 and 8 are the dimer
R = fle.PC645_dip_and_pos([4,8],4:6); %4 and 8 are the dimer
%convert units from angstrom to cm^
R = R*1e-8;
R12 = R(1,:) - R(2,:); %relative position
%transform so R12 is pointing down the z axis
uu = [R12(2);-R12(1);0]/sqrt(R12(1)^2+R12(2)^2); %axis off rotation
ucross = [0,-uu(3),uu(2);uu(3),0,-uu(1);-uu(2),uu(1),0];
ukron = kron(uu.',uu);
%cos(theta) = (R3/sqrt(R1^2+R2^2+R3^2));
%sin(theta) = (1 - R3^2/(R1^2 + R2^2 + R3^2))^(1/2);
Trot0 = (R12(3)/(R12(1)^2+R12(2)^2+R12(3)^2)^(1/2))*eye(3) + ...
         (1 - R12(3)^2/(R12(1)^2+R12(2)^2+R12(3)^2))^(1/2)*ucross + ...
         (1- R12(3)/(R12(1)^2+R12(2)^2+R12(3)^2)^(1/2))*ukron;
     
mu = (Trot0*(mu.') ).';   R = (Trot0*(R.') ).';    
R12 = (Trot0*(R12.') ).';    
     

[ex_basis,H0ex] = eig(H0);

mu_ex = [mu(1,:)*ex_basis(2,2) + mu(2,:)*ex_basis(3,2);...
        mu(1,:)*ex_basis(3,2) + mu(2,:)*ex_basis(2,3)];

%% Bath parameters 

Kappa = 0 ; %HEOM param
Kap1 = inf; Kap2 = 1; %More truncation parameters
lambdaD =100;omegaD=100;
omegaU = 1108.0;gammaU = 5.30884;lambdaU = 44.432;

%om_0 is the values of underdamped brownian modes included
om_0 = {[],omegaU,omegaU};   
lambda ={[],lambdaU ,lambdaU }; %reorganisation energy
gamma = {[],gammaU,gammaU}; 

%over damped
 lam_dru = {[],lambdaD,lambdaD}; %reorganisation energy of drude modes
gam_dru = {[],omegaD,omegaD};   

 QQ = zeros(length(H0),2); cnt=0;
 cc = zeros(1,sum(cellfun(@length,lambda))+ sum(cellfun(@length,lam_dru)) + Kappa);
 clear cc_com cc_acom vv
 for j=1:length(H0)
     
 [cc1,cc2R,cc2I,vv1,vv2,QQ(j,1),tmp] = coeffients_from_brownian_new(...
    lambda{j},gamma{j},om_0{j},Temp,Kappa,lam_dru{j},gam_dru{j});
 % QQ is in general length two for each point as I could (potentially)
 % included high frequency decay terms from poles in J(omega) which would
 % have anti commutator parts
 cc(cnt+1:cnt+length(tmp)) = tmp; cnt =cnt+length(tmp);
 cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
 cc_acom{j}= [cc2I,cc1*0];
 
 end
 syms t t1 t2 t3 real
use_markov = true;
 %Calculate Markovian parameters from redfield
[R_red,R_red_op ]= redfield_calc(H_site+eye(2)*(-E0),beta,gam_dru,lam_dru,gamma,lambda,om_0,use_markov);

supop = -1i*(kron(eye(length(H0ex)),H0ex)-kron(H0ex.',eye(length(H0ex))))-R_red_op;
%expm(supop*pi)-V*expm(J*pi)*V^(-1) = 0 or should
 PP = feval(symengine,'linalg::jordanForm',sym(supop),'All');
 Jre =double(real(PP(1))); %not convinced by outputting as double in some ways
 Vre =double(real(PP(2))); Jim =double(imag(PP(1))); Vim =double(imag(PP(2)));
% clear PP

red_propop = (Vre+1i*Vim)*expm((Jre+1i*Jim)*t)*(Vre+1i*Vim)^(-1);
red_propop1 =subs(red_propop,t,t1);
red_propop2 =subs(red_propop,t,t2);
red_propop3 =subs(red_propop,t,t3);


rho0 = zeros(size(H0ex)); rho0(1)=1;  rho0 = reshape(rho0,numel(rho0),1);
%this includes NO decay or decoherence with the ground state, not sure if
%this is really sensible at all to assume that the decoherence times are
%longer than the ground state dephasing times

%% Electric field and such like

syms om_u om_r t_sep om1 om2 om3 real  %om1..om3 are ft variables of t1..3


% u is p(u)mp r is p(r)obe
Efield_u_m = epu*exp(-1i*(om_u*t));    Efield_u_p = epu*exp(1i*(om_u*t));
Efield_r_m = epr*(exp(-1i*(om_r*t)));   Efield_r_p = epr*(exp(1i*(om_r*t)));

tau_u = 150*convfact/1000; tau_r = 150*convfact/1000; %convert from fs to 1/cm

E_0_u = 0.1; %peak electric field Volts / Angstrum
%there is no standard way to convert e*Angstrom into cm^-1, 
% hence dipole moments are not scaled from these units it is simply
%important that one keeps 1 eV = 8062.4 cm^-1, standard units for E field
%are V/m, hence we take the units of E to be 
% E_scaled = (E in Volts / Angstrum) * 8064.2, 
%that way the results stay in cm^(-1)
E_0_u = E_0_u* 8064.2 ;
%normalise envelopes to unity intensity 
Env_r = exp(-(t-t3-t_sep)^2/2/tau_r^2)/sqrt(sqrt(pi)*tau_r/2);  
Env_r_ft = (pi^(-1/4)*exp(-(om3^2*tau_r^2))/2)*sqrt(abs(tau_r));

Env_r_ft_bc = (pi^(-1/4)*exp(-((om3 + om_r)*(t*2*1i - t_sep*2*1i +...
    om3*tau_r^2 + om_r*tau_r^2))/2)*abs(tau_r))/tau_r^(1/2);

Env_r_ft_fwd = (pi^(-1/4)*exp(-((om3 - om_r)*(t*2*1i - t_sep*2*1i +...
    om3*tau_r^2 - om_r*tau_r^2))/2)*abs(tau_r))/tau_r^(1/2);

%freq envelope is a Gaussian centred at om_r with width 1/tau_r
%considered to have unit magnitude as it's linear in this
Env_u2 = E_0_u*exp(-(t-t3-t2)^2/2/tau_u^2)/sqrt(sqrt(pi)*tau_u); 
Env_u1 = E_0_u*exp(-(t-t3-t2-t1)^2/2/tau_u^2)/sqrt(sqrt(pi)*tau_u);



%% Calculate first response functions
% k_s direction = +/- k_1 where k_1 is k_pu or k_pr
%in this case just probe

% recall mu_ex = [mu(1,:)*ex_basis(2,2) + mu(2,:)*ex_basis(3,2);...
 %                   mu(1,:)*ex_basis(3,2) + mu(2,:)*ex_basis(2,3)];
coeffs = ex_basis(2:end,2:end);
mu_unit{1} = zeros(3); mu_unit{1}(2,1) = 1;  mu_unit{1}(1,2) = 1;
mu_unit{2} = zeros(3); mu_unit{2}(3,1) = 1;  mu_unit{2}(1,3) = 1;
for e1 = 1:2
mu_u_L_ex{e1} = -kron(mu_unit{e1},eye(3)); %left acting commutator part
mu_u_R_ex{e1} = kron(eye(3),mu_unit{e1});  %#ok<*SAGROW>%right acting, normal  
end


Rsp = norm(R12);
%|kr| and |ku| are just 2 pi/ lambda, or k = omega/c -> 2*pi*omega in cm^-1 
kr = 2*pi*om_r*kpr; %only beam to hit

syms mu_s_1_r mu_s_2_r mu_s_1_u mu_s_2_u  

%these sym variables represent the dipole operator on each site but also
%the exp(i R_k dot{ k_pr or k_pu})
sblst2 = {mu_s_2_r*conj(mu_s_1_r),mu_s_1_r*conj(mu_s_2_r),...
    mu_s_1_r*conj(mu_s_1_r),mu_s_2_r*conj(mu_s_2_r)};

ab1 = dimer_rot_av_fn_simp([1,2],-kr,Rsp,{mu(1,:),mu(2,:)}); %
ab2 = dimer_rot_av_fn_simp([1,2],+kr,Rsp,{mu(2,:),mu(1,:)}); 
ab3 = dimer_rot_av_fn_simp([1,2],kr,0 ,{mu(1,:),mu(1,:)}); 
ab4 = dimer_rot_av_fn_simp([1,2],kr,0 ,{mu(2,:),mu(2,:)}); 
J1 = sym(zeros(3,1)); %S1_nrwa = sym(zeros(3,1));
numpoints = 3000;  
prefct = sym(zeros(3,2));  %symbolic prefactor related to dipole averages
tintrng = linspace(0,1,numpoints); %needs to be long for CD signal as it isn't heterodyned
tint_ext = [-tintrng(end-1:2);tintrng]; %extended to -ve freq

Fs = 1/(tintrng(2)-tintrng(1));
 NFFT = 2^nextpow2(2*numpoints-1);  %next FFT
freq_rng = Fs/2*linspace(0,1,NFFT/2+1); %frequency range (NOT angular)

for e1 = 1:2

    %first calculate the rotational average of the vector components of the
    %dipole matrix elements corresponding to this Fman diagram
    temp2 = expand((mu_s_1_r*coeffs(e1,1)+mu_s_2_r*coeffs(e1,2))...
            *conj(mu_s_1_r*coeffs(e1,1)+mu_s_2_r*coeffs(e1,2)));
    prefct(:,e1) = subs(temp2,sblst2,{ab1,ab2,ab3,ab4});   
    %this prefactor depends on the wavevector (and seperations between
    %sites) so leave it in symbolic form
    
    %The red_propop works in the EXCITON basis
    temp = prefct(:,e1)*trace(reshape(mu_u_R_ex{e1}*red_propop1*mu_u_R_ex{e1}*rho0,3,3)); 
    J1 = J1 + temp;
    %kept in RWA

end
om_int_rng = abs(E0) + linspace(H0ex(2,2)-30/tau_r,H0ex(3,3)+30/tau_r,numpoints);
S1 = 1i*(J1-conj(J1)); %this factor also must be incluced
J1ft = fourier(J1*heaviside(t1),t1,om3);

%S1ft = fourier(S1*heaviside(t1),t1,om3);
S1ft = -(J1ft+conj(subs(J1ft,om3,-om3)));
%same thing in fourier space


temp = subs(S1ft,om_r,om3); %single limit for pulse (wide pulse limit)
%Jones matrix coefficients are as follows
alpha_J  = -om3*(8*pi^2)*real(temp(1)); %P(omega) = i S(omega)
eta_J  = om3*(16*pi^2)*real(temp(2)); 
delta_J  = -om3*(8*pi^2)*imag(temp(2)); 



%% precompute dimer rot averages note R12 points along z
Rsp = norm(R12);

syms K1 K2 K3 real %magnitudes
k1 = K1*kpu; k2 = K2*kpu; k3 = K3*kpr;

av = sym(zeros(3,2,2,2,2)); %first interaction is conjugated
%av_nr = sym(zeros(3,2,2,2,2)); %second interaction is conjugated
av_FO = sym(zeros(3,2,2,2,2));  
%first order terms in k1,k2,k3


for ell = 1:2
    for e3 = 1:2
        for e2 = 1:2
            for e1 = 1:2
deltaRmat = -1i*[R(ell,:)- R(e1,:);R(ell,:)- R(e2,:);R(ell,:)- R(e3,:)].';          
%mumat = [mu(ell,:);mu(e3,:);mu(e2,:);mu(e1,:)].';  
mumat = [mu(e1,:);mu(e2,:);mu(e3,:);mu(ell,:)].';  
for lp = 1:3 %dimension
av(lp,ell,e3,e2,e1) = tensor_average(mumat,1,lp);  
    for lp2=1:3
    
av_FO(lp,ell,e3,e2,e1) = av_FO(lp,ell,e3,e2,e1)+...
                    tensor_average([mumat,deltaRmat(:,1)],lp,lp2)*k1(lp2)...
                     + tensor_average([mumat,deltaRmat(:,2)],lp,lp2)*k2(lp2)...
                     + tensor_average([mumat,deltaRmat(:,3)],lp,lp2)*k3(lp2); 
                              
                 
    end                
end
            end
        end
    end
end


for e1= 1:2
    temp = zeros(3);
    temp(1,e1+1) = 1; 
mu_op{e1} = ex_basis'*(temp+temp')*ex_basis; %project from site to exciton
%mu_op_R{e1} = (red_propop')*kron(mu_op{e1},eye(3))*red_propop; %time dep
mu_op_L{e1} = (red_propop')*kron(eye(3),mu_op{e1})*red_propop;
%mu_op_LR{e1} = mu_op_L{e1}-mu_op_R{e1};

op0{e1} = subs(mu_op_L{e1},t,0);
op1{e1} = subs(mu_op_L{e1},t,t1);
op2{e1} = subs(mu_op_L{e1},t,t1+t2);
op3{e1} = subs(mu_op_L{e1},t,t1+t2+t3); 

% Op0{e1} = subs(mu_op_LR{e1},t,0);
% Op1{e1} = subs(mu_op_LR{e1},t,t1);
% Op2{e1} = subs(mu_op_LR{e1},t,t1+t2);
% Op3{e1} = subs(mu_op_LR{e1},t,t1+t2+t3); 

end
%IMPORTANT*** L is the operator action FROM the left, R is FROM the right****

R1 = sym([0;0;0]); %both beams assumed to be polarized along x only
R2 = sym([0;0;0]);   R3 = sym([0;0;0]);    R4 = sym([0;0;0]);
R1fo = sym([0;0;0]); %both beams assumed to be polarized along x only
R2fo = sym([0;0;0]);   R3fo = sym([0;0;0]);    R4fo = sym([0;0;0]);
%S3_test = sym([0;0;0]);  S3_FO_test = sym([0;0;0]);  

for ell = 1:2
    for e3 = 1:2
        for e2 = 1:2
            for e1 = 1:2
  
  %Rtmp =  trace(reshape(Op1{ell}*Op2{e3}*Op3{e2}*Op0{e1}*rho0,3,3)); 
  
  Rtmp1 = trace(reshape(op1{ell}*op2{e3}*op3{e2}*op0{e1}*rho0,3,3));
  Rtmp2 = trace(reshape(op0{ell}*op2{e3}*op3{e2}*op1{e1}*rho0,3,3));              
  Rtmp3 = trace(reshape(op0{ell}*op1{e3}*op3{e2}*op2{e1}*rho0,3,3));
  Rtmp4 = trace(reshape(op3{ell}*op2{e3}*op1{e2}*op0{e1}*rho0,3,3)); 
  
  tmp = av(:,ell,e3,e2,e1); 
  tmp2 = av_FO(:,ell,e3,e2,e1);
 
  R1 = R1 + tmp*Rtmp1;     R1fo = R1fo + tmp2*Rtmp1;
  R2 = R2 + tmp*Rtmp2;     R2fo = R2fo + tmp2*Rtmp2;
  R3 = R3 + tmp*Rtmp3;     R3fo = R3fo + tmp2*Rtmp3;
  R4 = R4 + tmp*Rtmp4;     R4fo = R4fo + tmp2*Rtmp4; 
  
   % S3_FO_test =  S3_FO_test + tmp*Rtmp;   S3_test =  S3_test + tmp2*Rtmp;
  
            end
        end
    end
end

clear mu_op mu_op_R mu_op_L op0 op1 op2 op3

S3 = R1-conj(R1) + R2-conj(R2) + R3-conj(R3) + R4-conj(R4);
S3_FO = R1fo-conj(R1fo) + R2fo-conj(R2fo) + R3fo-conj(R3fo) + R4fo-conj(R4fo);

%%

   R1_ft = fourier(subs(R1(1),om_r,om3)*heaviside(t1),t1,om1);
   R1_ft = fourier(R1_ft*heaviside(t2),t2,om1+om2);
   R1_ft = -1i*fourier(R1_ft*heaviside(t3),t3,om1 + om2 + om3);

   R1_ft_fo = fourier(subs(R1fo(2),om_r,om3)*heaviside(t1),t1,om1);
   R1_ft_fo = fourier(R1_ft_fo*heaviside(t2),t2,om1+om2);
   R1_ft_fo = -1i*fourier(R1_ft_fo*heaviside(t3),t3,om1 + om2 + om3); 
   
   R2_ft = fourier(subs(R2(1),om_r,om3)*heaviside(t1),t1,om1);
   R2_ft = fourier(R2_ft*heaviside(t2),t2,om1+om2);
   R2_ft = -1i*fourier(R2_ft*heaviside(t3),t3,om1 + om2 + om3);

   R2_ft_fo = fourier(subs(R2fo(2),om_r,om3)*heaviside(t1),t1,om1);
   R2_ft_fo = fourier(R2_ft_fo*heaviside(t2),t2,om1+om2);
   R2_ft_fo = -1i*fourier(R2_ft_fo*heaviside(t3),t3,om1 + om2 + om3); 
   
   R3_ft = fourier(subs(R3(1),om_r,om3)*heaviside(t1),t1,om1);
   R3_ft = fourier(R3_ft*heaviside(t2),t2,om1+om2);
   R3_ft = -1i*fourier(R3_ft*heaviside(t3),t3,om1 + om2 + om3);

   R3_ft_fo = fourier(subs(R3fo(2),om_r,om3)*heaviside(t1),t1,om1);
   R3_ft_fo = fourier(R3_ft_fo*heaviside(t2),t2,om1+om2);
   R3_ft_fo = -1i*fourier(R3_ft_fo*heaviside(t3),t3,om1 + om2 + om3); 
   
   R4_ft = fourier(subs(R4(1),om_r,om3)*heaviside(t1),t1,om1);
   R4_ft = fourier(R4_ft*heaviside(t2),t2,om1+om2);
   R4_ft = -1i*fourier(R4_ft*heaviside(t3),t3,om1 + om2 + om3);

   R4_ft_fo = fourier(subs(R4fo(2),om_r,om3)*heaviside(t1),t1,om1);
   R4_ft_fo = fourier(R4_ft_fo*heaviside(t2),t2,om1+om2);
   R4_ft_fo = -1i*fourier(R4_ft_fo*heaviside(t3),t3,om1 + om2 + om3); 
   
   %%
   
   S3tilde = -(R1_ft + conj(subs(R1_ft,{om1,om2,om3},{-om1,-om2,-om3}))+...
       R2_ft + conj(subs(R2_ft,{om1,om2,om3},{-om1,-om2,-om3}))+...
       R3_ft + conj(subs(R3_ft,{om1,om2,om3},{-om1,-om2,-om3}))+...
       R4_ft + conj(subs(R4_ft,{om1,om2,om3},{-om1,-om2,-om3})));
  
    S3tilde_fo = -(R1_ft_fo + conj(subs(R1_ft_fo,{om1,om2,om3},{-om1,-om2,-om3}))+...
       R2_ft_fo + conj(subs(R2_ft_fo,{om1,om2,om3},{-om1,-om2,-om3}))+...
       R3_ft_fo + conj(subs(R3_ft_fo,{om1,om2,om3},{-om1,-om2,-om3}))+...
       R4_ft_fo + conj(subs(R4_ft_fo,{om1,om2,om3},{-om1,-om2,-om3})));
   
   %% Calculate third order polarization in x first
syms om real
%after taking the fourier transform of P this introduces a factor
%delta(om-om1-om2-om3)
   S3tilde2 = subs(S3tilde,{om3},{om-om1-om2});
   S3tilde_fo_plus = subs(S3tilde_fo,{om3,K1,K2,K3},{om-om1-om2,...
                    2*pi*(om-om1-om2),2*pi*(om2),-2*pi*(om1)});
   S3tilde_fo_minus = subs(S3tilde_fo,{om3,K1,K2,K3},{om-om1-om2,...
                    2*pi*(om-om1-om2),-2*pi*(om2),2*pi*(om1)});
                
  om_u_range = 1.6714e+04; %resonances near 1.6714e+04,1.7338e+04
   %range of pump wavelengths
   numpoints = 300;
   om_plot_rng =abs(E0) + linspace(H0ex(2,2)-10/tau_u,H0ex(3,3)+10/tau_u,numpoints);
   tsep_range = linspace(0.05,3,250)/convfact;
   
   tdelay_fc = exp(-1i*t_sep*(om2+om1));
  
   quad_order= 37;   tau = tau_u;
   
  [omm, weig] = GaussHermite(quad_order,eps(20));
  omm = omm*sqrt(2)/tau;
  
    prefact = 1/sqrt(pi)/tau; 
  %%  
%   temp_fn_x=matlabFunction(S3tilde2*tdelay_fc,'vars',{om,om1,om2,om_u,t_sep});
%   temp_fn_y_plus=matlabFunction(S3tilde_fo_plus*tdelay_fc,'vars',{om,om1,om2,om_u,t_sep});
%   temp_fn_y_minus=matlabFunction(S3tilde_fo_minus*tdelay_fc,'vars',{om,om1,om2,om_u,t_sep});
%   
%     P3x = zeros(length(om_plot_rng),length(tsep_range),length(om_u_range));
%     P3y = P3x;
%     
%  for lp = 1:length(tsep_range)
%           tsp = tsep_range(lp); 
%       for lp1 = 1:length(om_u_range)
%           omu = om_u_range(lp1);
%           temp1 = om_plot_rng*0; temp2 = temp1;
%       for lp2 = 1:quad_order
%           for lp3 = 1:quad_order
%               %solve with quadrature, this will suck with big seperations
%               weights = weig(lp2)*weig(lp3);
% %       temp1 = temp1 + weights*(temp_fn_x(om_plot_rng,omm(lp2)-omu,omm(lp3)+omu, omu,tsp)...
% %                      + temp_fn_x(om_plot_rng,omm(lp2)+omu,omm(lp3)-omu, omu,tsp));
% %        temp2 = temp2 + weights*(temp_fn_y_plus(om_plot_rng,omm(lp2)-omu,omm(lp3)+omu, omu,tsp)...
% %                      + temp_fn_y_minus(om_plot_rng,omm(lp2)+omu,omm(lp3)-omu, omu,tsp));                               
%       temp1 = temp1 + weights*(temp_fn_x(om_plot_rng,omm(lp2)+omu,omm(lp3)-omu, omu,tsp)...
%                      + temp_fn_x(om_plot_rng,omm(lp2)-omu,omm(lp3)+omu, omu,tsp));
%        temp2 = temp2 + weights*(temp_fn_y_plus(om_plot_rng,omm(lp2)+omu,omm(lp3)-omu, omu,tsp)...
%                      + temp_fn_y_minus(om_plot_rng,omm(lp2)-omu,omm(lp3)+omu, omu,tsp));                               
%                                   
%           end
%       end
%       P3x(:,lp,lp1) = (temp1.*om_plot_rng).';   
%       P3y(:,lp,lp1) = (temp2.*om_plot_rng).';   
%       end
%  end
 
 %%  Evaluate using relative variables
 syms om_rel om_c real
 %relative variables

%   sym_to_fn('temp_fn_x.m',vpa(subs(S3tilde2, {om1,om2} ,...
%       {(om_c+om_rel)/sqrt(2),(om_c-om_rel)/sqrt(2)}) ,16),[om,om_c,om_rel])
%   sym_to_fn('temp_fn_y_plus.m',vpa(subs(S3tilde_fo_plus, {om1,om2} ,...
%       {(om_c+om_rel)/sqrt(2),(om_c-om_rel)/sqrt(2)}) ,16),[om,om_c,om_rel])
%   sym_to_fn('temp_fn_y_minus.m',vpa(subs(S3tilde_fo_minus, {om1,om2} ,...
%       {(om_c+om_rel)/sqrt(2),(om_c-om_rel)/sqrt(2)}) ,16),[om,om_c,om_rel])
  
   temp_fn_x=matlabFunction(subs(S3tilde2, {om1,om2} ,...
      {(om_c+om_rel)/sqrt(2),(om_c-om_rel)/sqrt(2)}),'vars',{om,om_c,om_rel});
  temp_fn_y_plus=matlabFunction(subs(S3tilde_fo_plus, {om1,om2} ,...
      {(om_c+om_rel)/sqrt(2),(om_c-om_rel)/sqrt(2)}),'vars',{om,om_c,om_rel});
  temp_fn_y_minus=matlabFunction(subs(S3tilde_fo_minus, {om1,om2} ,...
      {(om_c+om_rel)/sqrt(2),(om_c-om_rel)/sqrt(2)}),'vars',{om,om_c,om_rel});
 
  
 %%
  %pretty unavoidable, will have to tabulate the entire function...
  %quadrature points for om_r and a range to give the appropriate t_sep
  %range 
  tstep = 0.005; tmax = 3; %change these values according to all the bath parameters
  tsep_range_des = 0:tstep:tmax; %ps
  t_dimless = tsep_range_des*convfact*2/tau;

  om_c_range = (0:2*pi/(tmax*convfact):2*pi/(tstep*convfact))/sqrt(2);
  %om_c_tilde = linspace(-om_c_range(end)/2,om_c_range(end)/2,length(om_c_range));
  %om_c_dimless = om_c_tilde*2/tau;
  %  om_c_range = om_c_tilde  - om_c_tilde(1);
  nextfft = nextpow2(length(om_c_range));
  NN = 2^nextfft; 
  
   quad_order= 37;   tau = tau_u;
   
  [omm, weig] = GaussHermite(quad_order,eps(20));
  omm = omm*sqrt(2)/tau;  weig = weig*tau/2;
 
  P3x2 = zeros(length(tsep_range_des),length(om_plot_rng),length(om_u_range));
  P3y2 = P3x2;

  % I pm = int dw_c int dw_rel exp(i*sqrt(2)*w_c*t_sep/tau - (w_c^2 +
  % w_rel^2)/2) * F(om3 = om-sqrt(2)*om_c/tau,(om_c-om_rel)/sqrt(2)/tau...
  %  +/- sqrt(2) om_u,(om_c +om_rel)/sqrt(2)/tau  -/+ sqrt(2) om_u,)
  
  [omm1,omm2] = meshgrid(om_plot_rng,om_c_range);
  %omm1 has constant value along the vertical, omm2 in dim 2
  %therefore omm2 varies as om_c_range along dim 1, which is the one fft is
  %over.
  exp_fct = exp(-omm2.^2*tau_u^2/2)*sqrt(tau/2/sqrt(pi));
  tic
 for lp = 1:length(om_u_range)
          omu = om_u_range(lp);
          omplus = omm + sqrt(2)*omu;  omminus= omm - sqrt(2)*omu;
      for lp1 = 1:quad_order
          %fft defaults to dimension one, fft is over om_c
    P3x2(:,:,lp) = P3x2(:,:,lp) +weig(lp1)*(fft((temp_fn_x(omm1,omm2,omplus(lp1))...
                     +temp_fn_x(omm1,omm2,omminus(lp1))).*exp_fct)...
                     + ifft((temp_fn_x(omm1,-omm2,omplus(lp1))...
                     +temp_fn_x(omm1,-omm2,omminus(lp1))).*exp_fct) );
%     P3y2(:,:,lp) = P3y2(:,:,lp) +weig(lp1)*(fft(temp_fn_y_minus(omm1,omm2,omplus(lp1))...
%                     +temp_fn_y_plus(omm1,omm2,omminus(lp1)).*exp_fct)+...
%                     ifft(temp_fn_y_minus(-omm1,omm2,omplus(lp1))...
%                     +temp_fn_y_plus(-omm1,omm2,omminus(lp1)).*exp_fct));      
     P3y2(:,:,lp) = P3y2(:,:,lp) +weig(lp1)*(fft((temp_fn_y_plus(omm1,omm2,omplus(lp1))...
                    +temp_fn_y_minus(omm1,omm2,omminus(lp1))).*exp_fct)+...
                    ifft((temp_fn_y_plus(omm1,-omm2,omplus(lp1))...
                    +temp_fn_y_minus(omm1,-omm2,omminus(lp1))).*exp_fct));          
      end
 end 
 toc
P3x2 = P3x2*1i; P3y2 = P3y2*1i;
 %%
 % 9 Level Color Scale Colormap with Mapping to Grayscale for Publications.
%
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);

 lambda_nm = 10^7./om_plot_rng;
 om_pad = kron(ones(size(P3x2,1),1),om_plot_rng);
 %should also include 1/n(omega) factor but that is probably just that of
 %water
  Dalpha_J= (8*pi^2.*om_pad).*imag(P3x2);
  figure
  pcolor( lambda_nm ,tsep_range_des(20:floor(end/2)) ,Dalpha_J(20:floor(end/2),:))
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('Absorption change \Delta \alpha, for different probe wavelengths and pulse seperations');
 colormap(CMRmap)
colorbar

Deta_J  = (16*pi^2.*om_pad).*real(P3y2); 
  figure
  pcolor( lambda_nm ,tsep_range_des(20:floor(end/2))  ,Deta_J(20:floor(end/2),:))
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('CD shift \Delta \eta, for different probe wavelengths and pulse seperations');
colormap(CMRmap)
colorbar

 Ddelta_J  = (8*pi^2.*om_pad).*imag(P3y2);
  figure
  pcolor( lambda_nm ,tsep_range_des(20:floor(end/2))  ,Ddelta_J (20:floor(end/2),:))
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('OR shift \Delta \delta, for different probe wavelengths and pulse seperations');
colormap(CMRmap)
colorbar
%%
temp = abs(Deta_J(1,:));
[a,b] = max(temp);
figure
plot(tsep_range_des(1:floor(length(tsep_range_des)/2)),Deta_J(1 ...
:floor(length(tsep_range_des)/2),1:10:end))

%%
figure
hold on 
plot( lambda_nm,[Deta_J(floor(length(tsep_range_des)/16),:);...
    Deta_J(floor(length(tsep_range_des)/8),:);...
    Deta_J(floor(length(tsep_range_des)/4),:);...
    Deta_J(floor(length(tsep_range_des)/2),:)])

%%
figure
plot(tsep_range_des(1:size(Dalpha_J,1)),Dalpha_J(:,[1:15:end]))
%figure
%plot( lambda_nm ,Dalpha_J(1,:))
%figure
%plot( lambda_nm ,Deta_J(1,:))
%figure
%plot( lambda_nm ,Ddelta_J(1,:))
  %%
    Dalpha_J  = -(8*pi^2)*real(P3x); %P(omega) = i S(omega)
Deta_J  = (16*pi^2)*real(P3y); 
Ddelta_J  = -(8*pi^2)*imag(P3y); 

lambda_nm = 10^7./om_plot_rng;

figure
%tsep_range
plot(lambda_nm,Dalpha_J(:,20),'LineWidth',2)
xlabel('Wavelength (nm)');
ylabel('Absorption (a.u.)');

figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'YColor',[0 0 1]);
box(axes1,'on');
hold(axes1,'all');

plot(lambda_nm,Deta_J(:,20),'Parent',axes1,'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Circular dichorism shift','Color',[0 0 1]);
axes2 = axes('Parent',figure1,'YAxisLocation','right','YColor',[0 0.5 0],...
    'ColorOrder',[0 0.5 0;1 0 0;0 0.75 0.75;0.75 0 0.75;0.75 0.75 0;0.25 0.25 0.25;0 0 1],...
    'Color','none');
hold(axes2,'all');
plot(lambda_nm,Ddelta_J(:,20),'Parent',axes2,'LineWidth',2);
ylabel('Optical rotation shift (units of wavevector)','VerticalAlignment','cap',...
    'Color',[0 0.5 0]);

%%
figure
pcolor(tsep_range,lambda_nm,Deta_J);
set(gcf, 'renderer', 'zbuffer');
shading flat

figure
pcolor(tsep_range,lambda_nm,Ddelta_J);
set(gcf, 'renderer', 'zbuffer');
shading flat

figure
pcolor(tsep_range,lambda_nm,Dalpha_J);
set(gcf, 'renderer', 'zbuffer');
shading flat