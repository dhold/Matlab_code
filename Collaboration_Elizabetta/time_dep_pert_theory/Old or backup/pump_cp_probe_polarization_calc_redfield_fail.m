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
kpu = [0,sin(theta),cos(theta)];  epu = [1,0,0]; %same as  probe
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
omegaU = 1108.0;gammaU = 5.30884;lambdaU = 44.32;

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
 J =PP(1); %not convinced by outputting as double in some ways
 V =PP(2); clear PP

red_propop = double(V)*expm(double(J)*t)*double(V)^(-1);
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
    om3*tau_r^2 - om_r*tau_r^2))/2)*abs(tau_r))/tau_r^(1/4);

%freq envelope is a Gaussian centred at om_r with width 1/tau_r
%considered to have unit magnitude as it's linear in this
Env_u2 = E_0_u*exp(-(t-t3-t2)^2/2/tau_u^2)/sqrt(sqrt(pi)*tau_u/2); 
Env_u2_ft = fourier(fourier(Env_u2,t3,om3),t2,om2);
%Matlab doesn't seem to be dirac delta functions for FTs of e^(i om *k)
Env_u1 = E_0_u*exp(-(t-t3-t2-t1)^2/2/tau_u^2)/sqrt(sqrt(pi)*tau_u/2);
Env_u1_ft = fourier(fourier(fourier(Env_u1,t3,om3),t2,om2),t1,om1);
%syms E1 E2 V reals %these are the dipole matrix elements of each site


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
%Jones matrix coefficients are as follows
temp = subs(S1ft,om_r,om3);
alpha_J  = -om3^3/(2*pi^2)*real(temp(1)); 
eta_J  = om3^3/(pi^2)*real(temp(2)); 
delta_J  = -om3^3/(2*pi^2)*imag(temp(2)); 
if 1==1
alpha_toplot = double(subs(alpha_J,om3,om_int_rng));
temp2 = max(abs(alpha_toplot));
alpha_toplot = alpha_toplot/temp2; %scale
eta_toplot   = double(subs(eta_J ,om3,om_int_rng))/temp2;
delta_toplot = double(subs(delta_J  ,om3,om_int_rng))/temp2;

stick_spec = 0*om_int_rng;
[~,temp] = min(abs(om_int_rng-(H0ex(2,2)-H0ex(1,1))));
[~,temp2] = min(abs(om_int_rng-(H0ex(3,3)-H0ex(1,1))));
stick_spec([temp,temp2])=1;

figure
plot(10^7./om_int_rng,alpha_toplot,'LineWidth',2)
xlabel('Wavelength (nm)');
ylabel('Absorption (a.u.)');

figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'YColor',[0 0 1]);
xlim(axes1,[540.963590489601 639.999481382944]);
ylim(axes1,[-6.21533291325352e-05 6.07653854232304e-05]);
box(axes1,'on');
hold(axes1,'all');
plot(10^7./om_int_rng,eta_toplot,'Parent',axes1,'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Circular dichorism (relative to absolute absorption)','Color',[0 0 1]);
axes2 = axes('Parent',figure1,'YAxisLocation','right','YColor',[0 0.5 0],...
    'ColorOrder',[0 0.5 0;1 0 0;0 0.75 0.75;0.75 0 0.75;0.75 0.75 0;0.25 0.25 0.25;0 0 1],...
    'Color','none');
 xlim(axes2,[540.963590489601 639.999481382944]);
 ylim(axes2,[-6.21533291325352e-05 6.07653854232304e-05]);
hold(axes2,'all');
plot(10^7./om_int_rng,delta_toplot,'Parent',axes2,'LineWidth',2);
ylabel('Optical rotation (units of wavevector)','VerticalAlignment','cap',...
    'Color',[0 0.5 0]);

end
P1_ft = subs(Env_r_ft_fwd ,{t_sep,t},{0,0}).*S1ft;

P1inte = 1i*(S1.*subs(Env_r,{t_sept,t3},{0,t1})); 
%term that must be integrated wrt t1 from 0 to infinity to get the FO pol

%P1 = ifourier(int_arg,om3,t);

%% first order polarization 
clear temp temp2

sym_to_fn('temp_fn_11.m',vpa(P1inte(1),12),[t1,t,om_r])
sym_to_fn('temp_fn_22.m',vpa(P1inte(2),12),[t1,t,om_r])
sym_to_fn('temp_fn_33.m',vpa(subs(Env_r,{t_sep,t3},{0,0})*exp(-1i*om_r*t),12),[t,om_r]) 
%%
om_r_range = abs(E0) + linspace(H0ex(2,2)-200,H0ex(3,3)+200,100);
tintrng2 = linspace(-tau_r*7,tau_r*7,1000);

P1x_fwd = zeros(length(om_r_range),length(tintrng2));  %exp(i k.r)  term
P1y_fwd = zeros(length(om_r_range),length(tintrng));

exp_fc =  exp(-(tintrng2).^2/tau_r^2/2)/sqrt(sqrt(pi)*tau_r/2);
        %can select any frequency from this by selectively shifting it
        %along, but easier just to not to   
 
for lp = 1:length(om_r_range)
    
    omr = om_r_range(lp);
    intvar0 = temp_fn_33(tintrng2,omr);
    tic
    for lp2 = 1:length(tintrng2)
    P1x_fwd(lp,lp2) = trapz(tintrng,temp_fn_11(tintrng,tintrng2(lp2),omr));  
    end
    toc
    for lp2 = 1:length(tintrng)
    P1y_fwd(lp,lp2) = trapz(tintrng,temp_fn_22(tintrng,tintrng(lp2),omr)); 
    end    
       toc   
    %integral over the electric field can be performed as a convolution
    %anyway so do it with ffts
    % int_0^inf dx g(x)f(y-x) = ifft(fft(g(x),x,t)*fft(f(x),x,t),t,y)
    % in our case "x" is t1 and "y" is t, "t"is some dummy omega

Sig_X_1(lp) = (4*pi*omr)*trapz(tintrng2 , imag( intvar0.*P1x_fwd(lp,:)));
Sig_Y_1(lp) = trapz(tintrng,abs(P1y_fwd(lp,:)).^2);
%Sig_X_2(lp) = (4*pi*omr)*trapz(tint_ext, imag(exp(-1i*(omr)*tintrng).*exp_fc.*P1x_fwd(lp,:)));

CD_1(lp) = trapz(tint_ext,real(exp(+1i*(omr)*tintrng).*exp_fc.*P1y_fwd(lp,:)));
end
%% X signal
figure
lambda_rng_nm = 10^7./om_r_range;
plot(lambda_rng_nm,Sig_X_1) %note axis in inverse CM

%% Y signal
figure
plot(lambda_rng_nm,Sig_Y_1)

 %P1 = P1*exp(-1i*om_r*t);
%% CD signal

%split output into left and right CP light along z
% Eout = [Ex,Ey] = A[1,i] + B[1,-i]-> A = (Ex+Ey/i)/2 , B = (Ex-Ey/i)/2
% 2*|A|^2 = (|Ex|^2 + |Ey|^2 -  i(Ey*conj(Ex) - conj(Ey)*Ex))/2
% 2*(|A|^2-|B|^2) = 2*imag(Ey*Ex)  this being the CD signal
% assuming attenuation is weak, then Ex_in ~ Ex_out (not that important) so
% the signal is propto the real polarization in Y dirn
figure
plot(lambda_rng_nm,2*CD_1)

%% Calculate third order polarization

% 4 different F-man diagrams contribute
% R / L for acting to the left (usual) or right (commutator term)
% 2Stim emmission
% mu*_e1_R mu_e2_L mu*_e2_L mu_e1_R
% mu*_e1_R mu_e2_L mu_e1_R mu*_e2_L
%
% 2 GS bleaching
% mu*_e2_R mu_e2_R mu_e1_L mu^*_e1_L
% mu*_e1_R mu_e2_R mu*_e1_R mu_e1_R
%
%Calculate operator terms and then coefficients
%note the factor of exp(- i om_r (t-t3)) exp(+/- i om_u t1) with the +/-
%depending on whether the first (-) or second (+) term is conjugated
% precompute dimer rot averages note R12 points along z
Rsp = norm(R12);
%|kr| and |ku| are just 2 pi/ lambda, or k = omega/c -> 2*pi*omega in cm^-1 
kr = 2*pi*om_r*kpr; %second beam to hit
ku = 2*pi*om_u*kpu; %first beam to hit

a1b1 = dimer_rot_av_fn_simp([3,2],-ku-kr,Rsp,{mu(1,:),mu(2,:),mu(2,:),mu(1,:)}); %
a1b2 = dimer_rot_av_fn_simp([3,2],-ku+kr,Rsp,{mu(2,:),mu(1,:),mu(2,:),mu(1,:)}); 
a1b3 = dimer_rot_av_fn_simp([3,2],-ku,Rsp ,{mu(1,:),mu(1,:),mu(2,:),mu(1,:)}); 
a1b4 = dimer_rot_av_fn_simp([3,2],-ku,Rsp ,{mu(2,:),mu(2,:),mu(2,:),mu(1,:)}); 

a2b1 = dimer_rot_av_fn_simp([3,2],+ku-kr,Rsp,{mu(1,:),mu(2,:),mu(1,:),mu(2,:)}); %
a2b2 = dimer_rot_av_fn_simp([3,2],+ku+kr,Rsp,{mu(2,:),mu(1,:),mu(1,:),mu(2,:)}); 
a2b3 = dimer_rot_av_fn_simp([3,2],+ku,Rsp ,{mu(1,:),mu(1,:),mu(1,:),mu(2,:)}); 
a2b4 = dimer_rot_av_fn_simp([3,2],+ku,Rsp ,{mu(2,:),mu(2,:),mu(1,:),mu(2,:)}); 

a3b1 = dimer_rot_av_fn_simp([3,2],-kr,Rsp,{mu(1,:),mu(2,:),mu(1,:),mu(1,:)}); %
a3b2 = dimer_rot_av_fn_simp([3,2],+kr,Rsp,{mu(2,:),mu(1,:),mu(1,:),mu(1,:)}); 
a3b3 = dimer_rot_av_fn_simp([3,2],[0,0,0],Rsp ,{mu(1,:),mu(1,:),mu(1,:),mu(1,:)}); 
a3b4 = dimer_rot_av_fn_simp([3,2],[0,0,0],Rsp ,{mu(2,:),mu(2,:),mu(1,:),mu(1,:)}); 

a4b1 = dimer_rot_av_fn_simp([3,2],-kr,Rsp,{mu(1,:),mu(2,:),mu(2,:),mu(2,:)}); %
a4b2 = dimer_rot_av_fn_simp([3,2],+kr,Rsp,{mu(2,:),mu(1,:),mu(2,:),mu(2,:)}); 
a4b3 = dimer_rot_av_fn_simp([3,2],[0,0,0],Rsp ,{mu(1,:),mu(1,:),mu(2,:),mu(2,:)}); 
a4b4 = dimer_rot_av_fn_simp([3,2],[0,0,0],Rsp ,{mu(2,:),mu(2,:),mu(2,:),mu(2,:)}); 

rho0 = zeros(3^2,1);  rho0(1)=1;
syms a1 a2 a3 a4 b1 b2 b3 b4

subst3 = {a1*b1,a1*b2,a1*b3,a1*b4,a2*b1,a2*b2,a2*b3,a2*b4,...
            a3*b1,a3*b2,a3*b3,a3*b4,a4*b1,a4*b2,a4*b3,a4*b4};
subst4 = {a1b1,a1b2,a1b3,a1b4,a2b1,a2b2,a2b3,a2b4,...
            a3b1,a3b2,a3b3,a3b4,a4b1,a4b2,a4b3,a4b4};
%%  Calculate contributions from the different Fman diagrams
for e1 = 1:2
    for e2  = 1:2
%mu*_e1_R mu_e2_L mu*_e2_L mu_e1_R
%R1 = mu_unit{e1}*(red_propop3*(red_propop2*((red_propop1*mu_unit{e1}*rho0)*mu_unit{e2})*mu_unit{e2}));
R1{e1,e2} = exp(+1i*om_u*t1)*red_propop1*mu_u_R_ex{e1}*rho0;
R1{e1,e2} = red_propop2*mu_u_L_ex{e2}*R1{e1,e2};
R1{e1,e2} = exp(1i*om_r*t3)*trace(reshape(mu_u_R_ex{e1}*(red_propop3*mu_u_L_ex{e2}*R1{e1,e2}),3,3));

%calculate coefficient from averaging of terms
temp = mu_s_1_u*coeffs(e1,1)+mu_s_2_u*coeffs(e1,2);
temp = expand(temp*conj(mu_s_1_u*coeffs(e2,1)+mu_s_2_u*coeffs(e2,2)));
temp2 = expand((mu_s_1_r*coeffs(e2,1)+mu_s_2_r*coeffs(e2,2))...
            *conj(mu_s_1_r*coeffs(e1,1)+mu_s_2_r*coeffs(e1,2)));
 %this last interaction isn't "really" conjugated but has an 
%exp(i(r-T*R_j)dot k_pr) term which (factoring out the r) looks like a c.c.


% seperate these terms and sub in
% mu_s_2_u*conj(mu_s_1_u)   - (R1-R2)               
% mu_s_1_u*conj(mu_s_2_u)   + (R1-R2)                
% mu_s_1_u*conj(mu_s_1_u)   0 (R1-R2)
% mu_s_2_u*conj(mu_s_2_u)   0 (R1-R2)
sblst1 = {mu_s_2_u*conj(mu_s_1_u),mu_s_1_u*conj(mu_s_2_u),...
    mu_s_1_u*conj(mu_s_1_u),mu_s_2_u*conj(mu_s_2_u)};

temp = subs(temp,sblst1,{a1,a2,a3,a4});
temp2 = subs(temp2,sblst2,{b1,b2,b3,b4});        
coeff1{e1,e2} = subs(expand(temp2*temp),subst3,subst4);

%mu*_e1_R mu_e2_L mu_e1_R mu*_e2_L
%R2{e1,e2} =mu_unit{e1}*(red_propop3*(red_propop2*(mu_unit{e1}*(red_propop1*(rho0*mu_unit{e2})))*mu_unit{e2}));
R2{e1,e2} = exp(-1i*om_u*t1)*red_propop1*(mu_u_L_ex{e2}*rho0);
R2{e1,e2} = red_propop2*(mu_u_R_ex{e1}*R2{e1,e2});
R2{e1,e2} = exp(1i*om_r*t3)*trace(reshape(mu_u_R_ex{e1}*(red_propop3*mu_u_L_ex{e2}*R2{e1,e2}),3,3));

temp = conj(mu_s_1_u*coeffs(e2,1)+mu_s_2_u*coeffs(e2,2));
temp = expand(temp*(mu_s_1_u*coeffs(e1,1)+mu_s_2_u*coeffs(e1,2)));
temp2 = expand((mu_s_1_r*coeffs(e2,1)+mu_s_2_r*coeffs(e2,2))...
            *conj(mu_s_1_r*coeffs(e1,1)+mu_s_2_r*coeffs(e1,2)));
temp = subs(temp,sblst1,{a1,a2,a3,a4});
temp2 = subs(temp2,sblst2,{b1,b2,b3,b4});        
coeff2{e1,e2} = subs(expand(temp2*temp),subst3,subst4);        
   
% seems these terms don't actually contribute to CD signal
% mu*_e2_R mu_e2_R mu_e1_L mu_e1^*_L
R3{e1,e2} = exp(-1i*om_u*t1)*red_propop1*(mu_u_L_ex{e1}*rho0);
R3{e1,e2} = red_propop2*(mu_u_L_ex{e1}*R3{e1,e2});
R3{e1,e2} = exp(1i*om_r*t3)*trace(reshape(mu_u_R_ex{e2}*(red_propop3*(mu_u_R_ex{e2}*R3{e1,e2})),3,3));

temp =conj(mu_s_1_u*coeffs(e1,1)+mu_s_2_u*coeffs(e1,2));
temp = expand(temp*(mu_s_1_u*coeffs(e1,1)+mu_s_2_u*coeffs(e1,2)));
temp2 = expand((mu_s_1_r*coeffs(e2,1)+mu_s_2_r*coeffs(e2,2))...
                *conj(mu_s_1_r*coeffs(e2,1)+mu_s_2_r*coeffs(e2,2)));
temp = subs(temp,sblst1,{a1,a2,a3,a4});
temp2 = subs(temp2,sblst2,{b1,b2,b3,b4});        
coeff3{e1,e2} = subs(expand(temp2*temp),subst3,subst4);                  
            
            
% mu*_e2_R mu_e2_R mu*_e1_R mu_e1_R
R4{e1,e2} = exp(+1i*om_u*t1)*red_propop1*(mu_u_R_ex{e1}*rho0);
R4{e1,e2} = red_propop2*(mu_u_R_ex{e1}*R4{e1,e2});
R4{e1,e2} = exp(1i*om_r*t3)*trace(reshape(mu_u_R_ex{e2}*(red_propop3*(mu_u_R_ex{e2}*R4{e1,e2})),3,3));

temp =mu_s_1_u*coeffs(e1,1)+mu_s_2_u*coeffs(e1,2);
temp = expand(temp*conj((mu_s_1_u*coeffs(e1,1)+mu_s_2_u*coeffs(e1,2))));
temp2 = expand((mu_s_1_r*coeffs(e2,1)+mu_s_2_r*coeffs(e2,2))...
                *conj(mu_s_1_r*coeffs(e2,1)+mu_s_2_r*coeffs(e2,2)));
temp = subs(temp,sblst1,{a1,a2,a3,a4});
temp2 = subs(temp2,sblst2,{b1,b2,b3,b4});        
coeff4{e1,e2} = subs(expand(temp2*temp),subst3,subst4);                     
            
    end
end
%% collate all these terms together into one uber function

S3 = sym(zeros(3,1));
for e1 = 1:2
    for e2  = 1:2
        S3 =  S3  +  R1{e1,e2}*coeff1{e1,e2} + R2{e1,e2}*coeff2{e1,e2} ...
                +R3{e1,e2}*coeff3{e1,e2} +R4{e1,e2}*coeff4{e1,e2} ;
    end
end
%note exp(-1i*om_r*t) factor not included in this expression
%tic %no idea why this took so long...
%P3inte = simplify(P3inte,'IgnoreAnalyticConstraints',true,'Steps',2);
%toc
 %% sub in values and them numerically integrate them
%probably best to use a seperate script for this

%loop over range of om_r and om_u
%P_Y_3_sub = subs(P_Y_3_env,{om_r,om_u} ,{1050,1050});
%P_Z_3_sub = subs(P_Y_3_env,{om_r,om_u} ,{1050,1050});

%only likely significant after first pulse hits, start 5 sigma from pulse
%centrec
t_sep_rng = 1:0.5:100; %range in fs
t_sep_rng = t_sep_rng*convfact/1000; %range in inverse cm

%first evaluate in the simplist possible way, with impulsive limit
tic
sym_to_fn('temp_fn_1.m',vpa(S3(1)*Env_r*Env_u2*Env_u1,12),[t,t1,t2,t3,t_sep,om_r,om_u])
sym_to_fn('temp_fn_2.m',vpa(S3(2)*Env_r*Env_u2*Env_u1,12),[t,t1,t2,t3,t_sep,om_r,om_u])
toc


%temp_fn(0,t_sep,t-t_sep) Theta(t-t_sep) is the result of the integral
% but with Gauss Hermite evaluation of the first functions we have
%temp_fn(-tau(tt1+tt2),t_sep-tau*tt2,t-t_sep) Theta(t-t_sep) 
%Theta(t_sep-tau(tt1+tt2)) Theta(t_sep-tau*tt2)
[tt, wei] = GaussHermite(5,eps(20));
tt = tau_u*tt; wei = tau_u*wei;

om_r_range = linspace(0.6e3,2e3,20); om_u_range = linspace(0.6e3,2e3,20); 
Sig_Y_3_impulsive = zeros(length(om_r_range),length(om_u_range),length(t_sep_rng));
Sig_X_3_impulsive = Sig_Y_3_impulsive;


for tseplp = 1:length(t_sep_rng)
    
tsep = t_sep_rng(tseplp);
trng_x = linspace(sqrt(tsep),sqrt(tsep+40*tau_u),20).^2;
hevifct = double(trng_x-tsep >0) + double(trng_x-tsep ==0)/2;  

for om_r_lp = 1:length(om_r_range)
    for om_u_lp = 1:length(om_u_range)
        om_pr = om_r_range(om_r_lp);
        om_pu = om_u_range(om_u_lp);
     %tic   
    Int_Y_3_inpulsive = trng_x*0; Int_X_3_inpulsive = trng_x*0;
                    
for tlp1 = 1:length(tt) %GH quad points in time
    for tlp2 = 1:length(tt)
        tt1 = tt(tlp1); tt2 = tt(tlp2); wei1 = wei(tlp1)*wei(tlp2); 
        
    hevifct1 = double(tsep - (tt1+tt2) >0) + double(tsep - (tt1+tt2)==0)/2+...
        double(tsep - tau_u*tt2 >0) + double(tsep - tau_u*tt2==0)/2;
%factor from Heviside theta functions
    
    if hevifct1 ~=0
    Int_Y_3_inpulsive = Int_Y_3_inpulsive + wei1*hevifct1*hevifct.*...
        temp_fn_2(trng_x,-tt1-tt2,tsep-tt2,trng_x-tsep,tsep,om_pr,om_pu);
    Int_X_3_inpulsive = Int_X_3_inpulsive + wei1*hevifct1*hevifct.*...
        temp_fn_1(trng_x,-tt1-tt2,tsep-tt2,trng_x-tsep,tsep,om_pr,om_pu);
    end
    end
end
%toc
Sig_X_3_impulsive(om_r_lp,om_u_lp,tseplp) =tau_r*trapz(trng_x,...
    exp(-(trng_x-tsep).^2/tau_r^2).*Int_X_3_inpulsive);
Sig_Y_3_impulsive(om_r_lp,om_u_lp,tseplp) = tau_r*trapz(trng_x,abs(Int_Y_3_inpulsive).^2);
%factor of exp(-i * om_r *t) cancels either way
%this signal is not heterodyne detected as it is orthogonal to the probe
%beam polarization giving a weaker signal
    end
end
end