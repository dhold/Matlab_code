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



lambdaD =100;omegaD=100;
omegaU = 643.8131;gammaU = 5.30884; 

lambdaUrng = [0,0.5,1,2,3,4,5:5:65];
for lppp =1:length(lambdaUrng) %testing a specific feature

lambdaU = lambdaUrng(lppp);
mu = rand(1,3); 
R = rand(1,3);
%[inf,50,20,10,5,2,1]
%om_0 is the values of underdamped brownian modes included
om_0 = {omegaU};   
lambda ={lambdaU}; %reorganisation energy
gamma = {gammaU}; 
H_site = 10000; E0 = 0;
H0ex = zeros(2); H0ex(2,2) = H_site;
%over damped
 lam_dru = {lambdaD}; %reorganisation energy of drude modes
gam_dru = {omegaD};       
    
%% Generate HEOM stuff and redfield prop op
Kappa = 0 ; %HEOM param
Kap1 = inf; Kap2 = 4; %More truncation parameters
 QQ = zeros(length(H_site),2); cnt=0;
 cc = zeros(1,sum(cellfun(@length,lambda))+ sum(cellfun(@length,lam_dru)) + Kappa);
 clear cc_com cc_acom vv
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
 syms t om t1 t2 t3 real
use_markov = true;
 %Calculate Markovian parameters from redfield
[R_red,R_red_op ]= redfield_calc(H_site,beta,gam_dru,...
            lam_dru,gamma,lambda,om_0,use_markov);

supop = -1i*(kron(eye(length(H0ex)),H0ex)-kron(H0ex.',eye(length(H0ex))))-R_red_op;
%expm(supop*pi)-V*expm(J*pi)*V^(-1) = 0 or should
 PP = feval(symengine,'linalg::jordanForm',sym(supop),'All');
 tic
 TT1 = vpa(PP(1),12); TT2 = vpa(PP(2),12);
 Jre =double(real(TT1)); %not convinced by outputting as double in some ways
 Vre =double(real(TT2)); 
 Jim =double(imag(TT1)); 
 Vim =double(imag(TT2));
 toc
% clear PP

red_propop = (Vre+1i*Vim)*expm((Jre+1i*Jim)*t)*(Vre+1i*Vim)^(-1);

%Check jordan form is diagonal
if sum(sum(Jre-diag(diag(Jre))))~=0
    1
%     tmp = eye(size(Jim))*om;
% red_propop_ft = -1i*(Vre+1i*Vim)*diag(1./diag(-Vre+1i*(tmp-Vim)))*(Vre+1i*Vim)^(-1);
% %same thing in fourier space
end
red_propop1 =subs(red_propop,t,t1);

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
mu_unit = zeros(2); mu_unit(2,1) = 1;  mu_unit(1,2) = 1;

mu_u_L_ex = -kron(mu_unit,eye(2)); %left acting commutator part
mu_u_R_ex = kron(eye(2),mu_unit);  %#ok<*SAGROW>%right acting, normal  

J1 = sym(0);
numpoints = 3000;  
tintrng = linspace(0,1,numpoints); %needs to be long for CD signal as it isn't heterodyned
tint_ext = [-tintrng(end-1:2);tintrng]; %extended to -ve freq

Fs = 1/(tintrng(2)-tintrng(1));
 NFFT = 2^nextpow2(2*numpoints-1);  %next FFT
freq_rng = Fs/2*linspace(0,1,NFFT/2+1); %frequency range (NOT angular)

    temp = 1/3*trace(reshape(mu_u_R_ex*red_propop1*mu_u_R_ex*rho0,2,2)); 
    J1 = J1 + temp;

om_int_rng = abs(E0) + linspace(H0ex(2,2)-1000,H0ex(2,2)+1000,numpoints);
S1 = 1i*(J1-conj(J1)); %this factor also must be incluced
J1ft = -1i*fourier(J1*heaviside(t1),t1,om3);

%S1ft = fourier(S1*heaviside(t1),t1,om3);
S1ft = -(J1ft+conj(subs(J1ft,om3,-om3)));
%same thing in fourier space


temp = subs(S1ft,om_r,om3); %single limit for pulse (wide pulse limit)
%Jones matrix coefficients are as follows
nav =  om3*(8*pi^2)*real(temp(1)); 
alpha_J  = om3*(8*pi^2)*imag(temp(1)); %P(omega) = i S(omega)
%eta_J  = om3*(16*pi^2)*real(temp(2)); 
%delta_J  = om3*(8*pi^2)*imag(temp(2)); 
%****************
alpha_J_tp = matlabFunction(alpha_J,'vars',{om3});
%eta_J_tp = matlabFunction(eta_J,'vars',{om3});
% figure
% plot(om_int_rng,-alpha_J_tp(om_int_rng))
%figure
%plot(om_int_rng,eta_J_tp)
to_plot(lppp,:) = -alpha_J_tp(om_int_rng);
%to_plot_2(lppp,:) = -eta_J_tp(om_int_rng);
%********************
end