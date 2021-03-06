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

kpr = [1,0,0];   epr = [0,0,1];   bpr = [0,1,0];
% eL = [0,1,+1i]/sqrt(2); eR = [0,1,-1i]/sqrt(2);
% bL = cross(kpr,eL);  bR = cross(kpr,eR);  

%pump unit vectors
theta = atan(sqrt(2)); %angle between pump and probe
%note arctan(sqrt(2)) is not the magic angle from the cho 
%paper but the nmr magic angle so it will not remove quadrupole terms
kpu = [cos(theta),0,sin(theta)];
epu = [-sin(theta),0,cos(theta)];
bpu = cross(kpu,epu);
usetestdimer = true; %false;%
%dimer specific parameters
%mu is elec dipole
%mu = [mu11,mu12,mu13;mu21,mu22,mu23]; %elec dipole of both sites
%mm = [mm11,mm12,mm13;mm21,mm22,mm23]; %intrinsic mag dipole (can be o)
if usetestdimer
    fle = open('Hamiltonian_save.mat');
        H_site = fle.PC645Hamiltonian;
            E0 = fle.E0_PC645; 
       H_site = H_site(1:2,1:2);%only these 2
       H_site = diag(diag(H_site)) + fle.PC645_el_coup([4,8],[4,8]);
       %take newer version of electronic coupling that is much bigger for
       %some reason.  Uber excitons!
  
else
  
om_mono = 1000;  V = 50;
H_site = [om_mono,V;V,om_mono];
end
H0 = zeros(length(H_site)+1); H0(2:end,2:end) = H_site;

if usetestdimer
    
mu = fle.PC645_dip_and_pos([4,8],1:3); %4 and 8 are the dimer
R = fle.PC645_dip_and_pos([4,8],4:6); %4 and 8 are the dimer
%convert units from angstrom to cm^-1
R = R*10^(-8); %1 A = 0.1 nm = 0.1*10^-3*10^-3*10^(-1) cm

mm1 = [0,0,0;0,0,0]; %assume no intrinsic magnetic moment
else
mod_mu = 10;
mu = mod_mu*[1,0,0; 1/sqrt(3),-1/sqrt(3),1/sqrt(3)];    
%R = [R11,R12,R13;R21,R22,R23]; %vector from centre of molecule site dipole centres
R = 10^(-7)*[2,3,-5;-3,4,3]; %units must be cm
mm1 = [0,0,0;0,0,0]; %assume no intrinsic magnetic moment
end


mm2 = -1i * pi * om_mono/(4*pi^2) * cross(R,mu);

mm = mm1 + mm2;

mu_op = zeros(length(H0),length(H0),3); mm_op = mu_op;

for k =1:3
    mu_op(1,2:end,k) = mu(:,k);
    mm_op(1,2:end,k) = mm(:,k);
end

mu_op = mu_op +conj(permute(mu_op,[2,1,3])); 
mm_op = mm_op +conj(permute(mm_op,[2,1,3]));


[ex_basis,H0ex] = eig(H0);
for k = 1:3
mu_op_ex(:,:,k) = ex_basis'*(mu_op(:,:,k))*ex_basis;
mm_op_ex(:,:,k)  = ex_basis'*(mm_op(:,:,k) )*ex_basis;
end

 
%% Bath parameters 
if usetestdimer
    
lambdaD =100;omegaD=100;
omegaU = 1108.0;gammaU = 5.30884;lambdaU = 44.32;
 omegavib = {[],omegaU,omegaU};
num_vib = cellfun(@length,omegavib);
sQM   =  {0.0578,0.0578};
om_0 = {[],omegaU,omegaU};
%om_0 is the values of damped brownian modes included
lambda ={44.32,44.32,44.32,44.32}; %reorganisation energy
gamma = {5.30884,5.30884,5.30884,5.30884};   
    
 lam_dru = {lambdaD,lambdaD,lambdaD,lambdaD}; %reorganisation energy of drude modes
gam_dru = {omegaD,omegaD,omegaD,omegaD};   
    %coupg = sqrt(sQM{j}).*omegavib{j}; 
  % viblvls = {[3],[3],[],[] };   
else
    lam_dru = {[],60,60}; %reorganisation energy drude
    gam_dru = {[],100,100};  % Damping/width, drude
    
    lambda ={[],44.32,44.32}; %reorganisation energy
    gamma = {[],5.30884,5.30884};   %damping
    om_0 = {[],1674,1674}; %resonant frequency
%HEOM parameters
 Kappa = 0;
end

 QQ = zeros(3,2); cnt=0;
 cc = zeros(1,sum(cellfun(@length,lambda))+ sum(cellfun(@length,lam_dru)) + Kappa);
 clear cc_com cc_acom vv
 for j=1:3
 [cc1,cc2R,cc2I,vv1,vv2,QQ(j,1),tmp] = ...
    coeffients_from_brownian_new(lambda{j},gamma{j},om_0{j},Temp,Kappa,lam_dru{j},gam_dru{j});
 % QQ is in general length two for each point as I could (potentially)
 % included high frequency decay terms from poles in J(omega) which would
 % have anti commutator parts
 cc(cnt+1:cnt+length(tmp)) = tmp; cnt =cnt+length(tmp);
 cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
 cc_acom{j}= [cc2I,cc1*0];
 
 end
 
  use_markov = true;
 %Calculate Markovian parameters from redfield
 [R_red,R_red_op ]= redfield_calc(H_site,beta,gam_dru,lam_dru,gamma,lambda,om_0,use_markov);

pickr = false(length(H_site)+1); pickr(2:end,2:end) = true;  %picks exciton
pickr = reshape( pickr,numel(pickr),1);
R_red_op_full = zeros((length(H_site)+1)^2);
R_red_op_full(pickr,pickr) = R_red_op;
% drho(pickr) /dt = R_red_op*rho(pickr)
supop = -1i*(kron(eye(length(H0ex)),H0ex)-kron(H0ex.',eye(length(H0ex))))-R_red_op_full;
[V,J]=jordan(supop);%expm(supop*pi)-V*expm(J*pi)*V^(-1) = 0 or should
red_propop = V*expm(J*t)*V^(-1); %prop operator within markovian redfield
%this includes NO decay or decoherence with the ground state, not sure if
%this is really sensible at all to assume that the decoherence times are
%longer than the ground state dephasing times


 %********see page 176 of mukamel for this derivation
gam = zeros(2,1);
 
for e1 = 1:2
   
    tmp = gam_dru{e1+1}; tmp2 = lam_dru{e1+1};
    for k =1:length(tmp)
       gam(e1) =gam(e1) + tmp2(k)*(2/beta/tmp(k)-1i);
    end
    tmp = gamma{e1+1}; tmp2 = lambda{e1+1};     tmp3 = om_0{e1+1};
    for k = 1:length(tmp)
        
        if tmp(k) == 0
            eta = sqrt(tmp3(k)^2-tmp(k)^2/4);
       gam(e1) = gam(e1) + tmp2(k)*(2*imag(cot(beta*(tmp(k)/2+1i*eta)/2)...
                        /(tmp(k)/2+1i*eta))-1i);
       vn = 2*pi*(1:20)/beta; %take 20 matsubara, probably more than enough
       gam(e1) = gam(e1) + sum(vn.^2./((tmp3(k)^2+vn.^2).^2-tmp(k)^2*vn.^2));
       
        end
    end
    
end
Gam = zeros(2,2);
for e1 = 1:2
    for e2 = 1:2
        Gam(e1,e2) = (gam(e1)+gam(e2))/2;
    end
end


%% Electric field and such like

syms om_u om_r t reals 
syms t1 t2 t3 t_sep reals

% u is p(u)mp r is p(r)obe
%Efield_u = epu*(exp(1i*(om_u*t-dot(kpu,rpos))) + exp(-1i*(om_u*t-dot(kpu,rpos))));
%Efield_r = eLR*(exp(1i*(om_r*t-dot(kpr,rpos))) + exp(-1i*(om_r*t-dot(kpr,rpos))));
Efield_u_m = epu*exp(-1i*(om_u*t));    Efield_u_p = epu*exp(1i*(om_u*t));
Efield_r_m = epr*(exp(-1i*(om_r*t)));   Efield_r_p = epr*(exp(1i*(om_r*t)));

Bfield_u_m = cross(kpu,Efield_u_m); Bfield_u_p = cross(kpu,Efield_u_p);
Bfield_r_m = cross(kpr,Efield_r_m ); Bfield_r_p = cross(kpr,Efield_r_p );
%consider envelopes seperately

tau_u = 15*convfact/1000; tau_r = 15*convfact/1000; %convert from fs

E_0_u = 0.1; %peak electric field Volts / Angstrum
%there is no standard way to convert e*Angstrom into cm^-1, 
% hence dipole moments are not scaled from these units it is simply
%important that one keeps 1 eV = 8062.4 cm^-1, standard units for E field
%are V/m, hence we take the units of E to be 
% E_scaled = (E in Volts / Angstrum) * 8064.2, 
%that way the results stay in cm^(-1)
E_0_u = E_0_u* 8064.2 ;

Env_r = exp(-(t-t3-t_sep)^2/tau_r^2);  
%considered to have unit magnitude as it's linear in this
Env_u2 = E_0_u*exp(-(t-t3-t2)^2/tau_u^2); 
Env_u1 = E_0_u*exp(-(t-t3-t2-t1)^2/tau_u^2);
%syms E1 E2 V reals %these are the dipole matrix elements of each site

for k=1:3
    ex_basis(:,k) = ex_basis(:,k) / norm(ex_basis(:,k));
end
%project magnetic and electric dipole vectors into the new basis

mu_ex = [mu(1,:)*ex_basis(2,2) + mu(2,:)*ex_basis(3,2);...
        mu(1,:)*ex_basis(3,2) + mu(2,:)*ex_basis(2,3)];

mm_ex = [mm(1,:)*ex_basis(2,2) + mm(2,:)*ex_basis(3,2);...
        mm(1,:)*ex_basis(3,2) + mm(2,:)*ex_basis(2,3)];


%% Analytic expression for first order polarization w/ markovian approx
if 1==0
%Note that int(exp[-(t-t1-ts)^2/T^2]*exp[i*om*t1-Gam*t1],{t1,0,inf})
%   = 1/2 * exp( 1/4 * (Gam - i*om)*(-4*t+(Gam-i*om)*T^2+4*ts))*sqrt(pi)*...
%      *erfc((-2*t+(Gam-i*om)*T^2+2*ts)/2T)
% times this by exp(-i*om_pr*t) and integrate from -inf to inf give
% 1i*pi*T*exp(-om_r^2*T^2-om_r*ts) / (Gam - i(om-om_pr))
% Env_r = E_0_u*exp(-(t-t3-t_sep)^2/tau_u^2);

P1sym = sym(zeros(3,1));
syms G1 G2 reals om1 om2

temp1 = sqrt(pi)/2 * exp( 1/4 * (G1 - 1i*(om1-om_r))*(-4*t+(G1-1i*(om1-om_r))*tau_u^2+4*t_sep))...
      *erfc((-2*t+(G1-1i*(om1-om_r))*tau_u^2+2*t_sep)/(2*tau_u));
temp2 = sqrt(pi)/2 * exp( 1/4 * (G2 - 1i*(om2-om_r))*(-4*t+(G2-1i*(om2-om_r))*tau_u^2+4*t_sep))...
      *erfc((-2*t+(G2-1i*(om2-om_r))*tau_u^2+2*t_sep)/(2*tau_u));  
P1sym(3) = temp1*dot(mu_ex(1,:),mu_ex(1,:))/3 + temp2*dot(mu_ex(2,:),mu_ex(2,:))/3 ;
P1sym(2) = temp1*dot(mu_ex(1,:),mm_ex(1,:))/3 + temp2*dot(mu_ex(2,:),mm_ex(2,:))/3 ;

P1sym = P1sym * exp(-1i*om_r*t);
% in cw case intensity for left |E_0 + p_x - i p_y|^2 - |E_0|^2
% or right |E_0 + p_x + i p_y|^2 - |E_0|^2

%foCD = int(imag(P1sym(2)*subs(Env_r,t3,0)),t,[-inf,inf]);

%% Calculate first order polarization with HEOM
rho0 = zeros(size(H0)); rho0(1,1)=1;
for k = 1:3 %need to average over dimensions
init_cond10{k} = (mu_op(:,:,k)*rho0-rho0*mu_op(:,:,k))/3;
init_cond11{k} = (mm_op(:,:,k)*rho0-rho0*mm_op(:,:,k))/3;
end

t_end_ps = 2;
t_end = convfact*t_end_ps; numpoints =1000;
final_desired_trange = linspace(0,t_end,numpoints);
FO_pol_HEOM = zeros(3,numpoints);
%solve using HEOM in order to propogate each initial condition in time
Kap1 = [];
freq_scale =1000; %numerically easier if you take out a factor out of the
%Hamiltonian to reduce the oscillatory components, this is like adding a
%factor of exp(-1i*freq_scale*t) into it
H0scaled = H0 - eye(size(H0))*freq_scale;
for k = 1:3
    
[out1,t_range1,out2,t_range2,ftharg] = reduced_HEOM(init_cond10{k},t_end,H0scaled,QQ(2:end,:)...
          ,{cc_com{2},cc_com{3}},{cc_acom{2},cc_acom{3}},vv,Kap1,3,numpoints);
                            %reduced HEOM only have to propogate a system
                            %which is proportial to the number of excited
                            %states times the size of the heirarchy where
                            %as normal would be N^2 times this
                            if isempty(Kap1)
                                Kap1=ftharg;
                            end
%out1 is first column, out2 first row     
%trim zeros
lg1 = t_range1~=0; lg1(1) = true; lg2 = t_range2~=0; lg2(1) = true; 
out1 = out1(lg1,:); t_range1 = t_range1(lg1); 
out2 = out2(lg2,:); t_range2 = t_range2(lg2); 

tmp1 = interp1(t_range1,out1,final_desired_trange);
tmp1(isnan(tmp1)) = 0;
FO_pol_HEOM(3,:) = FO_pol_HEOM(3,:) + mu_op(1,2:end,k)*(tmp1.');
tmp2 = interp1(t_range2,out2,final_desired_trange);
tmp2(isnan(tmp2)) = 0;
FO_pol_HEOM(3,:) = FO_pol_HEOM(3,:) + (tmp2*mu_op(2:end,1,k)).';

%figure 
%plot(final_desired_trange,abs([tmp1,tmp2]))

[out1,t_range1,out2,t_range2] = reduced_HEOM(init_cond11{k},t_end,...
 H0scaled,QQ,{cc_com{2},cc_com{3}},{cc_acom{2},cc_acom{3}},vv,Kap1,3,numpoints);
                            
lg1 = t_range1~=0; lg1(1) = true; lg2 = t_range2~=0; lg2(1) = true; 
out1 = out1(lg1,:); t_range1 = t_range1(lg1); 
out2 = out2(lg2,:); t_range2 = t_range2(lg2);                      
                            
tmp1 = interp1(t_range1,out1,final_desired_trange);

tmp1(isnan(tmp1)) = 0;
FO_pol_HEOM(2,:) = FO_pol_HEOM(2,:) + mu_op(1,2:end,k)*(tmp1.');
tmp2 = interp1(t_range2,out2,final_desired_trange);
tmp2(isnan(tmp2)) = 0;
FO_pol_HEOM(2,:) = FO_pol_HEOM(2,:) + (tmp2*mu_op(2:end,1,k)).';                          
                            
%figure 
%plot(final_desired_trange,abs([tmp1,tmp2]))
end

%reintroduce this scaling factor 
FO_pol_HEOM = FO_pol_HEOM.*kron(exp(+1i*freq_scale*final_desired_trange),ones(3,1)); 

%This is actually not yet the polarization, we must ALSO integrate each of
%these signals with E_z(t-t1) and B_y(t-t1), these quantities are the same 
%in our chosen unit system / basis so all that needs to be considered is
%the envelope function
%%
if 1==0
numpoints2 = 200;
sig_t_rng = linspace(-5*tau_u,5*tau_u,numpoints2)+t_sep;
sig_t_rng = sig_t_rng(sig_t_rng>0);
[tt2,tt1] = meshgrid(final_desired_trange,sig_t_rng);
%gaussian quadrature points might actually be easier but cba

tmp = exp(-(tt1-tt2-t_sep).^2/tau_u^2);
%Currently no assumptions at all have been made about the probe beam, we
%must now pic a frequency
om_probe = 1050;
tmp2 = exp(-1i*final_desired_trange*om_probe);

FO_pol_HE = zeros(3, length(sig_t_rng));
%FO_pol_HE_NRW =FO_pol_HE;
for k = 1:length(sig_t_rng) %loop over this range
   
    FO_pol_HE(2,k) = trapz(final_desired_trange,tmp2.*tmp(k,:).*FO_pol_HEOM(2,:));
    FO_pol_HE(3,k) = trapz(final_desired_trange,tmp2.*tmp(k,:).*FO_pol_HEOM(3,:));
    %FO_pol_HE_NRW(2,k) = trapz(final_desired_trange,conj(tmp2).*tmp(k,:).*FO_pol_HEOM(2,:));
    %FO_pol_HE_NRW(3,k) = trapz(final_desired_trange,conj(tmp2).*tmp(k,:).*FO_pol_HEOM(3,:));
    %These really are VERY small, nice to know
end

dt = final_desired_trange(2) - final_desired_trange(1);
NFFT  = 2^nextpow2(numpoints);
Y = fft(FO_pol_HEOM(3,:),NFFT)/numpoints;
plot(linspace(0,1,NFFT/2+1)*2/dt,Y(1:NFFT/2+1))
end
%to get end CD signal, integrate Y component times 
% E_env(t-t_1)* exp(-i*om*(t-t1)) * imag(P_Y(t1))
end
%% Calculate decoherence / phase ev terms

om_ge = [H0ex(2,2),H0ex(3,3)]; %expected value of g/ 1st exciton transitions
om_e1e2 = [0,H0ex(2,2)-H0ex(3,3); H0ex(3,3)-H0ex(2,2),0]; %expected value 
%of differences between exciton states /hbar


%******** GET VALUES FOR THESE  *************%%
%gam = [1,1]; %inverse lifetimes of states  a,b
%Gam = [1,1;1,1]; %markovian decoherence rates#
%<delta om_a(t1) delta_om_b(t2)> ~ Gam[a,b] delta(t2-t1)

for e1 = 1:2
    for e2  = 1:2
        if e1~=e2
            Eta(e1,e2) = (Gam(e2,e2)+Gam(e1,e1)-2*Gam(e2,e1));
        else %e1==e2
            Eta(e1,e2) = gam(e1);            
        end
    end
end


%Markovian approximation derived from Redfield theory
for e1 = 1:2
    for e2  = 1:2
        
rr{1,e1,e2} = exp(1i*(-om_ge(e1)*t3 + om_e1e2(e2,e1)*t2 -om_ge(e1)*t1)) *...
            exp(-Gam(e1,e1)*t3 - Eta(e2,e1)*t2 - Gam(e1,e1)*t1);
rr{2,e1,e2} = exp(1i*(-om_ge(e1)*t3 + om_e1e2(e2,e1)*t2 +om_ge(e2)*t1)) *...
            exp(-Gam(e1,e1)*t3 - Eta(e2,e1)*t2 - Gam(e2,e2)*t1);
rr{3,e1,e2} = exp(1i*(-om_ge(e1)*t3 +om_ge(e2)*t1)) *...
            exp(-Gam(e1,e1)*t3 - gam(e2)*t2 - Gam(e2,e2)*t1);        
rr{4,e1,e2} = exp(1i*(-om_ge(e2)*t3 -om_ge(e2)*t1)) *...
            exp(-Gam(e2,e2)*t3 - gam(e1)*t2 - Gam(e1,e1)*t1);
    end
end

%Non Markovian
% for e1 = 1:2
%     for e2  = 1:2
%         
% r2{1,e1,e2} = exp(1i*(-om_ge(e1)*t3 + om_e1e2(e2,e1)*t2 -om_ge(e1)*t1)) *...
%             exp(-Gam(e1,e1)*t3 - Eta(e2,e1)*t2 - Gam(e1,e1)*t1);
% r2{2,e1,e2} = exp(1i*(-om_ge(e1)*t3 + om_e1e2(e2,e1)*t2 +om_ge(e2)*t1)) *...
%             exp(-Gam(e1,e1)*t3 - Eta(e2,e1)*t2 - Gam(e2,e2)*t1);
% r2{3,e1,e2} = exp(1i*(-om_ge(e1)*t3 +om_ge(e2)*t1)) *...
%             exp(-Gam(e1,e1)*t3 - gam(e2)*t2 - Gam(e2,e2)*t1);        
% r2{4,e1,e2} = exp(1i*(-om_ge(e2)*t3 -om_ge(e2)*t1)) *...
%             exp(-Gam(e2,e2)*t3 - gam(e1)*t2 - Gam(e1,e1)*t1);
%     end
% end


%% Calculate first order polarization 
% k_s direction = +/- k_1 where k_1 is k_pu or k_pr
%in this case just probe
int_10  = 0 ; int_11 = 0;
for e1 = 1:2
int_10 = int_10 + exp(1i*(om_r-om_ge(e1))*t1-Gam(e1,e1)*t1)*dot(mu_ex(e1,:),mu_ex(e1,:))/3;
%  
int_11 = int_11 + exp(1i*(om_r-om_ge(e1))*t1-Gam(e1,e1)*t1)*dot(mu_ex(e1,:),mm_ex(e1,:))/3;
end
P1 = sym(zeros(3,1));


temp = express_exp_series(int_10,t1);
temp2 = express_exp_series(int_11,t1);
for lpvar = 1:size(temp,2)
    
P1(3) = P1(3)+temp{1,lpvar}*exp(-1i*om_r*t)*int(exp(temp{2,lpvar}),t1,[0,inf]);
P1(2) = P1(2)+temp{1,lpvar}*exp(-1i*om_r*t)*int(exp(temp2{2,lpvar}),t1,[0,inf]);

end

 P1 = P1*exp(-1i*om_r*t);
% %% Calculate first order polarization
% P_Y_1 = sym(0);
% for e1 = 1:2
% 
%         P_Y_1 = P_Y_1 + dot(mm(e1,:),mu(e1,:))*exp(-Gam(e1,e1)*t1)/3;
%         
% end
% 
% P_Y_1 = simplify(exp(1i*om_pr*(t1-t))*subs(Env_r,t3,t1)*P_Y_1);


%% Calculate third order polarization
%take only interactions with 2 pump and then 1 probe
% All these terms have a prefactor of i^3 int_0^inf dt1 dt2 dt3 exp(i
% om_pr [t3-t]) env_pr(t-t3) env_pump(t-t2-t3)  env_pump(t-t2-t3-t1)
P_Y_3 = sym(0);
for e1 = 1:2
    for e2  = 1:2
     
        P_Y_3 = P_Y_3 + (1/30)*((dot(mu(e1,:),mu(e1,:))*dot(mm(e2,:),mu(e2,:))...
             - 3*dot(mu(e2,:),mu(e1,:))*dot(mm(e2,:),mu(e1,:)))*(cos(theta)-1)+ ...
             2*cos(theta)*(2*dot(mu(e2,:),mu(e2,:))*dot(mu(e1,:),mm(e1,:))-...
             dot(mu(e1,:),mu(e2,:))*dot(mm(e1,:),mu(e2,:))))*rr{1,e1,e2}*exp(1i*om_u*t1);
         
        P_Y_3 = P_Y_3 + (1/30)*((dot(mu(e1,:),mu(e1,:))*dot(mm(e2,:),mu(e2,:))...
             - 3*dot(mu(e2,:),mu(e1,:))*dot(mm(e1,:),mu(e2,:)))*(cos(theta)-1)+ ...
             2*cos(theta)*(2*dot(mu(e2,:),mu(e2,:))*dot(mu(e1,:),mm(e1,:))-...
             dot(mu(e1,:),mu(e2,:))*dot(mm(e1,:),mu(e2,:))))*rr{2,e1,e2}*exp(-1i*om_u*t1);     

        P_Y_3 = P_Y_3 + (1/15)*(2*dot(mu(e1,:),mu(e1,:))*dot(mm(e2,:),mu(e2,:))...
             - dot(mu(e2,:),mu(e1,:))*dot(mm(e2,:),mu(e1,:)))*rr{3,e1,e2}*exp(1i*om_u*t1);
         
        P_Y_3 = P_Y_3 + (1/15)*(2*dot(mu(e1,:),mu(e1,:))*dot(mm(e2,:),mu(e2,:))...
             - dot(mu(e2,:),mu(e1,:))*dot(mm(e2,:),mu(e1,:)))*rr{4,e1,e2}*exp(-1i*om_u*t1);              
         
    end
end

P_Y_3_exp = express_exp_series(1i*om_r*t3*P_Y_3,t1);

P_Z_3 = sym(0);
for e1 = 1:2
    for e2  = 1:2
     
        tmp1 = dot(mu(e1,:),mu(e1,:))*dot(mu(e2,:),mu(e2,:));
        tmp2 = dot(mu(e1,:),mu(e2,:))*dot(mu(e1,:),mu(e2,:));
        
        P_Z_3 = P_Z_3 + (tmp1*(3*cos(theta)^2-1) + tmp2*(3+cos(theta)^2))...
                    *rr{1,e1,e2}*exp(1i*om_u*t1)/30;
         
        P_Z_3 = P_Z_3 + (tmp1*(3*cos(theta)^2-1) + tmp2*(3+cos(theta)^2))...
                    *rr{2,e1,e2}*exp(-1i*om_u*t1)/30;    

        P_Z_3 = P_Z_3 + (tmp1*(2-cos(theta)^2) + tmp2*(3*cos(theta)^2-1))...
                    *rr{3,e1,e2}*exp(1i*om_u*t1)/15;
         
        P_Z_3 = P_Z_3 + (tmp1*(2-cos(theta)^2) + tmp2*(3*cos(theta)^2-1))...
                    *rr{4,e1,e2}*exp(-1i*om_u*t1)/15;         
         
    end
end

P_Z_3_exp = express_exp_series(1i*om_r*t3*P_Z_3,t1);
%% Remove rapidly oscillating terms to improve convergence
if 1==0 %not yet finished
om_r_range = linspace(800,1200);
om_u_range = om_r_range;
max_freq = 5/tau_u; %max frequency of oscillation to bother with
lp = false;
for k = 1:length(a)
lp(k) = P_Y_3_exp{1,k}==0;
end
for kk = 1:length(lp)
if ~lp(kk)
    tmp = subs(P_Y_3_exp{2,kk},{t1,om_r,om_u},{1,om_r_range(end),om_u_range(end)});
    tmp2 = subs(P_Y_3_exp{2,kk},{t1,om_r,om_u},{1,om_r_range(1),om_u_range(1)});
    tmp3 = sign(imag(tmp))~=sign(imag(tmp2));
    
    
    if abs(imag(tmp))<max_freq || abs(imag(tmp2)) <maxfreq || tmp3
    cnt = cnt+1;
    P_Y_3_exp_new{1,cnt} = P_Y_3_exp{1,kk};
    P_Y_3_exp_new{2,cnt} = P_Y_3_exp{2,kk};
    end
end
end

%now remove rapidly oscillating terms with t2

%now remove rapidly oscillating terms with t3
end

%% sub in values and them numerically integrate them
%probably best to use a seperate script for this



%loop over range of om_r and om_u
%P_Y_3_sub = subs(P_Y_3_env,{om_r,om_u} ,{1050,1050});
%P_Z_3_sub = subs(P_Y_3_env,{om_r,om_u} ,{1050,1050});

%only likely significant after first pulse hits, start 5 sigma from pulse
%centre
t_sep_rng = 10:1:100; %range in fs
t_sep_rng = t_sep_rng*convfact/1000; %range in inverse cm

%first evaluate in the simplist possible way, with impulsive limit
P_Z_3_env = simplify(exp(1i*om_r*(t3-t))*P_Z_3);
P_Y_3_env = simplify(exp(1i*om_r*(t3-t))*P_Y_3);

sym_to_fn('temp_fn_1.m',E_0_u^2*P_Z_3_env,[t,t1,t2,t3,t_sep,om_r,om_u])
sym_to_fn('temp_fn_2.m',E_0_u^2*P_Y_3_env,[t,t1,t2,t3,t_sep,om_r,om_u])
%temp_fn(0,t_sep,t-t_sep) Theta(t-t_sep) is the result of the integral
% but with Gauss Hermite evaluation of the first functions we have
%temp_fn(-tau(tt1+tt2),t_sep-tau*tt2,t-t_sep) Theta(t-t_sep) 
%Theta(t_sep-tau(tt1+tt2)) Theta(t_sep-tau*tt2)
[tt, wei] = GaussHermite(9,eps(20));
tt = tau_u*tt; wei = tau_u*wei;

om_r_range = linspace(950,1150,10); om_u_range = linspace(950,1150,10); 
Sig_Y_3_inpulsive = zeros(length(om_r_range),length(om_u_range),length(t_sep_rng));
Sig_Z_3_inpulsive = Sig_Y_3_inpulsive;


for tseplp = 1:length(t_sep_rng)
    
tsep = t_sep_rng(tseplp);
trng = linspace(sqrt(tsep),sqrt(tsep+40*tau_u),20).^2;
hevifct = double(trng-tsep >0) + double(trng-tsep ==0)/2;  

for om_r_lp = 1:length(om_r_range)
    for om_u_lp = 1:length(om_u_range)
        om_pr = om_r_range(om_r_lp);
        om_pu = om_u_range(om_u_lp);
     %tic   
    Int_Y_3_inpulsive = trng*0; Int_Z_3_inpulsive = trng*0;
                    
for tlp1 = 1:length(tt) %GH quad points in time
    for tlp2 = 1:length(tt)
        tt1 = tt(tlp1); tt2 = tt(tlp2); wei1 = wei(tlp1)*wei(tlp2); 
        
    hevifct1 = double(tsep - (tt1+tt2) >0) + double(tsep - (tt1+tt2)==0)/2+...
        double(tsep - tau_u*tt2 >0) + double(tsep - tau_u*tt2==0)/2;
%factor from Heviside theta functions
    
    if hevifct1 ~=0
    Int_Y_3_inpulsive = Int_Y_3_inpulsive + wei1*hevifct1*hevifct.*...
        temp_fn_2(trng,-tt1-tt2,tsep-tt2,trng-tsep,tsep,om_pr,om_pu);
    Int_Z_3_inpulsive = Int_Z_3_inpulsive + wei1*hevifct1*hevifct.*...
        temp_fn_1(trng,-tt1-tt2,tsep-tt2,trng-tsep,tsep,om_pr,om_pu);
    end
    end
end
%toc
Sig_Y_3_inpulsive(om_r_lp,om_u_lp,tseplp) =tau_r*trapz(trng,...
    exp(-(trng-tsep).^2/tau_r^2).*Int_Y_3_inpulsive);
Sig_Z_3_inpulsive(om_r_lp,om_u_lp,tseplp) = tau_r*trapz(trng,abs(Int_Z_3_inpulsive).^2);

%this signal is not heterodyne detected as it is orthogonal to the probe
%beam polarization giving a weaker signal
    end
end
end
%% Now use monte carlo methods to do integral properly


%include envelopes explicitly
P_Z_3_env = simplify(P_Z_3_env*Env_u1*Env_u2*Env_r);
P_Y_3_env = simplify(P_Y_3_env*Env_u1*Env_u2*Env_r);

sym_to_fn('temp_fn_1.m',P_Z_3_env,[t,t1,t2,t3,om_r,om_u])
sym_to_fn('temp_fn_2.m',P_Y_3_env,[t,t1,t2,t3,om_r,om_u])



Int_Y_3 = trng*0; Int_Z_3 = trng*0;
IInt_Y_3 = trng*0; IInt_Z_3 = trng*0;



int_rel_tol = 1e-4;




for tlp = 1:length(trng)
    tt = trng(tlp);
    tmp1Y = 0;  tmp2Y = 0; tmp1Z = 0;  tmp2Z = 0;
    %montecarlo sampling, non independant randvar, 
    %not 100% sure about the validity of this, but it seems like it should
    %converge to the correct answer
    smp_vls = randn(3,400); %randvar with mean of 1, will be tau_u and tau_r
    prb = prod(exp(-smp_vls.^2)/sqrt(2*pi));
    smp_vls(1,:) = smp_vls(1,:)*tau_r + tt-t_sep;
    lg1 = smp_vls(1,:)>=0;  smp_vls = smp_vls(:,lg1); prb=prb(lg1); %reject < 0 
    smp_vls(2,:) = smp_vls(2,:)*tau_u + tt-smp_vls(1,:);
    lg2 = smp_vls(2,:)>=0;  smp_vls = smp_vls(:,lg2); prb=prb(lg2); 
    smp_vls(3,:) = smp_vls(3,:)*tau_u + tt-smp_vls(1,:)-smp_vls(2,:);
    lg3 = smp_vls(2,:)>=0;  smp_vls = smp_vls(:,lg3); prb=prb(lg3); 

    tot_num = length(prb);
    
    for lp = 1:tot_num
        
        tmp1Y = tmp1Y + temp_fn_1...
            (smp_vls(1,lp),smp_vls(2,lp),smp_vls(3,lp),1050,1050)/prb(lp);
        tmp1Z = tmp1Z + temp_fn_2...
            (smp_vls(1,lp),smp_vls(2,lp),smp_vls(3,lp),1050,1050)/prb(lp);
    end
    %montecarlo sampling, independant randvar
    smp_vls2 = randn(3,400);
    prb2 = prod(exp(-smp_vls.^2)/sqrt(2*pi));
    smp_vls2(1,:) = smp_vls2(1,:)*tau_r + tt-t_sep;
    smp_vls2(2,:) = smp_vls2(2,:)*tau_u + t_sep;
    smp_vls2(3,:) = smp_vls2(3,:)*tau_u; %half will be rejected
    lg = all(smp_vls>0);
    smp_vls2 = smp_vls2(:,all(smp_vls>0)); %goes in dim 1 by default
    prb2 = prb2(lg);
    
    tot_num2 = length(prb2);
    for lp = 1:tot_num2
        
        tmp2Y = tmp2Y + temp_fn_1...
            (smp_vls2(1,lp),smp_vls2(2,lp),smp_vls2(3,lp),1050,1050)/prb2(lp);
        tmp2Z = tmp2Z + temp_fn_2...
            (smp_vls2(1,lp),smp_vls2(2,lp),smp_vls2(3,lp),1050,1050)/prb2(lp);
    end

end
