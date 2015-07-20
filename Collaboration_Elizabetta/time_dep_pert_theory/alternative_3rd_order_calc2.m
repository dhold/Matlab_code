%Used to calculate the third order response function for a dimer using
%Redfield theory and explicitly included, damped Harmonic oscillator modes
%get inverse cm units
Temp = 300; %temp in Kelvin
[convfact, beta,speed_unit]= inv_cm_unit_sys(Temp);

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
            N=2;
       H_site = H_site(1:N,1:N);%only these 2
       N=length(H_site);
       H_site = H_site - eye(N)*E0;
      
H0 = zeros(length(H_site)+1); H0(2:end,2:end) = H_site;

mu = fle.PC645_dip_and_pos([4,8],1:3); %4 and 8 are the dimer
take_same_amp = true;
if take_same_amp
%assume allare the same amplitude
for k = 1:N
    mu(k,:) = mu(k,:)/norm(mu(k,:));
end
end
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

%Bath parameters 

lambdaD =100;omegaD=100;
omegaU = 650;  %omegaU = 1034; 
gammaU = 5.30884; 
%lambdaUrng = [0,0.5,1,2,3,4,5:5:65];
lambdaU = 44;%lambdaUrng(lppp);

%[inf,50,20,10,5,2,1]
%om_0 is the values of underdamped brownian modes included
om_0 = {omegaU,omegaU};   
lambda ={lambdaU ,lambdaU }; %reorganisation energy
gamma = {gammaU,gammaU}; 
%include these explicitly in redfield model

%over damped
 lam_dru = {lambdaD,lambdaD}; %reorganisation energy of drude modes
gam_dru = {omegaD,omegaD};       
    


%% Add in vibrations
rho0 = zeros(size(H0ex)); rho0(1)=1;  rho0 = reshape(rho0,numel(rho0),1);
%this includes NO decay or decoherence with the ground state, not sure if
%this is really sensible at all to assume that the decoherence times are
%longer than the ground state dephasing times
displ = [[0;0],eye(2)*sqrt(2*lambdaU/omegaU)];

H0_renorm = H_site - diag(cellfun(@sum, lam_dru) + cellfun(@sum, lambda));

om_vib = [om_0{:}];  numvib = [2,2];
[H_ex_vib,fock_space_rep,~,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(H0_renorm,om_vib,numvib,displ,ones(1,2)) ;
%indiv_op{5} projects to the exciton basis
sz1 = length(indiv_op{1}) ; sz2 = length(indiv_op{2}); %H_vib


 %%   Calculate standard decoherence term on vibration    
 Lindblad_op = sparse(zeros(length(H_ex_vib)^2));
 %The Lindblad operator into the exciton basis
 %Lindblad equation
%d rho/ dt = mu a rho a^d + v a^d rho a -(mu+v)/2 (N rho +rho N + rho) +
%           (mu-v)/2 rho
% make the RHS operator, note mu > v >= 0 and for us mu = gamma*(n+1) and v
% = gamma*n where n is the mean number of excitations in bath osc
%  Ham_op  = Lindblad_op;
for k = 1:length(numvib)
nav = exp(-beta*om_vib(k))/(1-exp(-beta*om_vib(k)));
muu = gamma{k}*(1+nav); nu = gamma{k}*(nav);
%express the annihilation op in the full hilbert space
aa = kron(sparse(eye(sz1)),sparse(eye(prod(numvib(1:k-1)))));
aa = kron(aa, diag(sqrt(1:numvib(k)-1),1)); 
aa = kron(aa,sparse(eye(prod(numvib(k+1:end))))); 
%project into exciton basis (not ex vib basis!) as Redfield op in ex basis
aa = indiv_op{5}'*aa*indiv_op{5};
Id = eye(size(aa));
% adag = aa'; adaga = diag(0:numvib(k)-1);

Lindblad_op = Lindblad_op  +  muu*kron((aa').',aa) +...
                 nu*(kron(aa.',(aa'))-kron(Id,Id)) -...
                 (muu+nu)/2 *(kron(aa'*aa,Id) + kron(Id,aa'*aa));
%         adaga = om_vib(k)*(aa'*aa+Id/2);
% Ham_op =Ham_op  -1i*(kron(Id,adaga)-kron(adaga,Id));    
end    

%% Generate redfield prop op

use_markov = true;
 %Calculate Markovian parameters from redfield
[R_red]= redfield_calc(H0_renorm,beta,gam_dru,...
            lam_dru,{[],[]},{[],[]},{[],[]},use_markov);
%Expand Redfield operator to include vibrational manifold

R_red_op_full = zeros(length(H_ex_vib)^2); 
%pad this out to the full thing
temp = zeros(length(H_ex_vib)); 

cnt1vib = 0;  cnt2vib = 1; cnt1 = 1; cnt2 = 1;
%tic
for lp = 1:sz1^2*sz2^2
    
    %iterate along one vibrational lvl
    cnt1vib = cnt1vib+1;
    if cnt1vib > sz2 
        cnt1vib=1; cnt1 = cnt1+1;
        if cnt1 > sz1
            cnt1 = 1;  cnt2vib = cnt2vib+1;
            if cnt2vib > sz2
                cnt2vib = 1; cnt2 = cnt2+1;
            end
        end
    end
    temp = temp*0;
   for a = 1:sz1 
        for b = 1:sz1 
            temp((a-1)*sz2+cnt1vib,(b-1)*sz2+cnt2vib) = R_red(cnt1,cnt2,a,b);
        end
   end

    R_red_op_full(lp,:) = reshape(temp,1,sz1^2*sz2^2 ).';

end
%toc
%N2*(a-1)+A+ N*N2*(N2*(b-1)+B) 

decoherence_op = Lindblad_op - R_red_op_full; 

L_op = -1i*(kron(eye(length(H_ex_vib)),H_exciton)-...
            kron(H_exciton.',eye(length(H_ex_vib))));    
% [P,D] = eig(full(supop));
% %prop_op = P*expm(D*t)*(P^(-1)); %this would take too much mem/time
% D = diag(D); Pinv = P^(-1); %~= P'
% prop_op = @(t)  P*diag(exp(D*t))*Pinv;

supop = sparse(L_op + decoherence_op);
  
rho_fl = zeros(length(supop),1); rho_fl(1) = 1;
%populate the vibrations thermally by allowing system to equilibrate

de_fn = @(t,v) supop*v;
points_to_save = 40; t_end_ps = 30; t_end = t_end_ps*convfact;
%t_range = linspace(0,t_end,points_to_save);

options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-12);
output_DE_fun(points_to_save,rho_fl,'notsavingnayway'); 
ode45(de_fn,[0,t_end],rho_fl,options);
[tmp1,tmp2]  =  output_DE_fun(points_to_save,rho_fl,'get_data');
rho_fl = tmp2(end,:);
%can check convergence if I want

%% Calculate time dependence of mu
% Simplified numerical setup only considers the parts that matter and
% scales appropriately the frequency

coeffs = ex_basis(2:end,2:end);
for e1 = 1:length(H_site)
    tmp = zeros(length(H0ex)); tmp(1,e1+1)=1;
    mu_unit{e1} = tmp;
    %same thing in full Liouville space
    mu_hilb{e1} = kron(tmp,eye(sz2));
    mu_lio{e1} = reshape(kron(tmp,eye(sz2)),sz2^2*sz1^2,1);
    
    mu_left{e1} = sparse(eye(length(mu_hilb{e1})),mu_hilb{e1});
    mu_right{e1} = sparse(mu_hilb{e1},eye(length(mu_hilb{e1})));
end
tracer_mat = reshape(logical(eye(sz2*sz1)),sz2^2*sz1^2,1);
%calculate the time dependence of these operators
%take only the elements which relate to g-e electronic coherences in upper 
%right of the matrix
tmp  = zeros(sz1*sz2); tmp(1:sz2,sz2+1:sz2*(sz1))=1;
[aa,bb,ss] = find(tmp);
tmp = logical(reshape(tmp,sz2^2*sz1^2,1));
%deco_op_red = decoherence_op(tmp,tmp);

supop_red = supop(tmp,tmp);
om_char =eig(H0_renorm);  %use to scale for easier numerics
%set number of points to take and how many to save
points_to_save = 2^14; t_end_ps = 2; t_end = t_end_ps*convfact;
t_range = linspace(0,t_end,points_to_save);
t_step = (t_range(2)-t_range(1));
rel_freq_range = linspace(0,2*pi/t_step ,points_to_save);
f_spacing = rel_freq_range (2)-rel_freq_range (1);


%round shifts to values that are included in rel_freq_range
for k = 1:N
    [~,tmp2] = min(abs(rel_freq_range-om_char(k)));
    om_char(k) = rel_freq_range(tmp2);
    shift_pos(k) = tmp2;
end

t_range_ext = [t_range(end:-1:2),t_range].';
om_rng_ext = (-1/t_step : 1/t_end :1/t_step)*2*pi/sqrt(2);
% tmpp=(0:length(om_rng_ext)-1).';
% phase_shift=exp(-1i*2*pi*tmpp*(length(om_rng_ext)-1)/2/length(om_rng_ext)); 
% % calculating phase shift from the fact the FFT doesn't centre the data on
% % zero, this must also be padded out to be the correct size
% phase_shift=repmat(phase_shift,1,sum(tmp));

Decay_Rate = zeros(N,sz2);
for e1 = 1:length(H_site)
    tic
    %rescale this operator to account for the phase fluc

supop_red_scale  = sparse(supop_red - om_char(e1)*1i*eye(length(supop_red)));
de_fn = @(t,v) (v.'*supop_red_scale).'; %this I don't think accurately describes
%the time evolution... 
de_fn2 = @(t,v) (v*supop_red_scale-supop_red_scale*v);

    mu_red = mu_lio{e1}(tmp);
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-70);
output_DE_fun(points_to_save,mu_red,'notsavingnayway'); 
ode45(de_fn,[0,t_end],mu_red,options);
[tmp1,tmp2]  =  output_DE_fun(points_to_save,mu_red,'get_data');
   toc
 tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
 tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
 
 %interpolate to constant t range and include the overarching exponential
tmp2 = interp1(tmp1,tmp2,t_range);
% fit exponential decay rate to weakly oscillating lines as init guess



on_diag_cohers = find(mu_red); 

for lpp = 1:length(on_diag_cohers)
Decay_Rate(e1,lpp) = exp_env_fit(t_range,tmp2(:,on_diag_cohers(lpp)),0);
end
%will be approx width of peaks in Fourier spectrum 
% mu_scaled_t{e1} = tmp2.*repmat(...
%     exp(min(Decay_Rate(e1,:))*t_range.'),1,length(supop_red)) ;  

% mu_full_omega{e1} = fft([zeros(size(tmp2,1)-1,size(tmp2,2));tmp2].*repmat(...
%                exp(1i*om_char(e1)*t_range_ext),1,length(supop_red)),[],1)./phase_shift ;  
%mu_full_omega{e1} = fft([zeros(size(tmp2,1)-1,size(tmp2,2));tmp2])./phase_shift ; 
%test = tmp2.*repmat(exp(1i*om_char(e1)*t_range.'),1,length(supop_red)) ; 
%test = [test(end:-1:2,:)*0;test];
%test2 = fft(test./phase_shift);
            
 mu_full_t{e1} = tmp2.*repmat(...
                exp(1i*om_char(e1)*(t_range.')),1,length(supop_red)) ;  
end
clear tmp1 tmp2 tmp3


%to get mu as a matrix use the fact that it is hermitian

%%  Calculate first order time domain response function


f_leng = length(H_ex_vib)^2; h_leng = length(H_ex_vib);
S_1  = zeros(2,length(t_range)); 

 om_char = mean(eig(H0_renorm));
k_pr = 2*pi*om_char; %assumed to go along z%
ft_scale_shift =exp(1i*om_char*t_range);
tmp = zeros(size(t_range));

for j1 = 1:N

    A1 = ex_basis(2,j1+1)*mu_full_t{1}(1,:);  
    B1 = ex_basis(2,j1+1)*mu_full_omega{1}(1,:);  
            %time dependent operators are expressed in the exciton basis,
        %convert them back to site basis by taking superpositions
            for ell = 2:N
                A1 = A1 + ex_basis(1+ell,j1+1)*mu_full_t{ell}(1,:);  
                B1 = ex_basis(2,j1+1)*mu_full_omega{ell}(1,:);  
            end
            O1 = sparse(aa,bb,A1, h_leng , h_leng );
            init_state= (O1+O1')*reshape(rho_fl,h_leng,h_leng);
    for j2 = 1:N    

        A2 = ex_basis(2,j2+1)*mu_full_t{1};
            for ell = 2:N
                A2 = A2 + ex_basis(1+ell,j2+1)*mu_full_t{ell};  
            end
        
        xx_av = dot(mu(j1,:),mu(j2,:))/3;
        yx_av = 1i*k_pr*dot(mu(j1,:),cross(mu(j2,:),R(j1,:)-R(j2,:)))/6;
          
        for tlp = 1:length(t_range)
            O2 = sparse(aa,bb,A2(tlp,:), h_leng , h_leng );
            tmp(tlp) = 2*imag(trace((O2+O2')*init_state));
        end
        
     S_1(1,:) = S_1(1,:) + ft_scale_shift.*tmp.*xx_av;
     S_1(2,:) = S_1(2,:) + ft_scale_shift.*tmp.*yx_av; 
    end
    
end
%%  FT and plot
S_1_FT = fft(S_1,[],2);  
lambda_rng = 10^7./(rel_freq_range(1:end/2)-om_char);
red_rng = lambda_rng>500 & lambda_rng<750;
figure
scale_fct = max(imag(S_1_FT(1,red_rng )));
plotyy(lambda_rng (red_rng),imag(S_1_FT(1,red_rng ))/scale_fct,...
            lambda_rng (red_rng),real(S_1_FT(2,red_rng ))/scale_fct)
xlabel('wavelength, nm')

%test if CD signal appears to be conservative, should end with ~0
figure
test = cumsum(real(S_1_FT(2,1:end/2)));
plot(lambda_rng (red_rng),test(red_rng)-test(find(red_rng,1)-1))

%% This just is not an efficient way to evaluate time evolution
tic
 rho_eq = reshape(rho_fl,h_leng,h_leng); %rho_eq(2,2) = 1; %simple form

A1 = sym('Fge', [sz2 N*sz2]); %B1= sym('Fef', [N,N*(N-1)/2 ]) ;
A2 = sym('Gge', [sz2 N*sz2]); %B2= sym('Gef', [N,N*(N-1)/2 ]) ;
A3 = sym('Hge', [sz2 N*sz2]); %B3= sym('Hef', [N,N*(N-1)/2 ]) ;
A4 = sym('Ige', [sz2 N*sz2]); %B4= sym('Ief', [N,N*(N-1)/2 ]) ;
V1 = sym(zeros(size(H_ex_vib))); %V1(2:N+1,N+2:end) = B1;
V1(1:sz2,sz2+1:(N+1)*sz2) = A1;  V1=V1+V1';
V2 = sym(zeros(size(H_ex_vib))); %V2(2:N+1,N+2:end) = B2; 
V2(1:sz2,sz2+1:(N+1)*sz2) = A2;  V2=V2+V2';
V3 = sym(zeros(size(H_ex_vib))); %V3(2:N+1,N+2:end) = B3; 
V3(1:sz2,sz2+1:(N+1)*sz2) = A3;  V3=V3+V3';
V4 = sym(zeros(size(H_ex_vib)));%V4(2:N+1,N+2:end) = B4; 
V4(1:sz2,sz2+1:(N+1)*sz2) = A4;  V4=V4+V4';

tmp= 2*imag(trace(V2*V3*V4*V1*rho_eq + V1*V3*V4*V2*rho_eq+...
                  V1*V2*V4*V3*rho_eq + V4*V3*V2*V1*rho_eq));
              %reshape into lioville space shapes which I store them as
A1 = reshape(A1,N*sz2^2,1); A2 = reshape(A2,N*sz2^2,1);
A3 = reshape(A3,N*sz2^2,1); A4 = reshape(A4,N*sz2^2,1);
              
third_order_fn = matlabFunction(tmp,'vars',{A1,A2,A3,A4});
toc
clear A1 A2 A3 A4 B1 B2 B3 B4 

   %% Calculate third order polarization components S_{xxxx} & S_{yxxx}
   %by this I mean S_{xxxx;k_r k_r k_u -k_u} + S_{xxxx;k_r k_r -k_u k_u}
% etc as these things are wavevector dependent, just use the characteristic
% wavevector for k_r and k_u as this doesn't make a big difference
 for k = 1:3
    pol{k} = [1,0,0]; % all beams along x, pol{1} = pol{2} as its PP 
 end
 om_char = mean(eig(H0_renorm));
 k_r = 2*pi*om_char*kpr; 
 k_u = 2*pi*om_char*kpu;

 %input pulse widths in fs, sd NOT FWHM
  tau_u_fs = 150;  tau_r_fs = 150;
 tau_u = 150*convfact/1000; tau_r = 150*convfact/1000; %convert from fs to 1/cm

 
 t1_rng = t_range(t_range<4*tau_u); %assume this is sufficient
 numpoints1 = 100; %not exact number, will in general be less
 t1_lp_rng = 1:ceil(length(t1_rng)/numpoints1):length(t1_rng);
 %using this in this way will map this to the actual points of t_range
 %where I have evaluated the original function
 t1_rng = t1_rng(t1_lp_rng); 
 t1_vc = 1:length(t1_rng);
 
 t_sep_max_ps = 1;  t_sep_max =  t_sep_max_ps * convfact;
 
 t2_rng = t_range(t_range<t_sep_max+tau_u); numpoints2 =200;
 t2_lp_rng =  1:floor(length(t2_rng)/numpoints2):length(t2_rng);  
 t2_rng = t2_rng(t2_lp_rng);
 t3_rng = t2_rng;  t3_lp_rng =  t2_lp_rng;
  % Could parallelise this so the loop runs over multiple nodes
 Stmp = zeros(length(t1_rng),length(t2_rng),length(t3_rng));
 
 Sx1 = Stmp;   Sy1 = Stmp;   Sx2 = Stmp;   Sy2 = Stmp;  

for j1 = 1:N
        A1 = ex_basis(2,j1+1)*mu_full_t{1}(1,:);  
            %time dependent operators are expressed in the exciton basis,
        %convert them back to site basis by taking superpositions
         for ell = 2:N
                A1 = A1 + ex_basis(1+ell,j1+1)*mu_full_t{ell}(1,:);      
         end
   %only at t=0
   A1 = A1.';
    for j2 = 1:N 
        A2 = ex_basis(2,j2+1)*mu_full_t{1};  
         for ell = 2:N
                A2 = A2 + ex_basis(1+ell,j2+1)*mu_full_t{ell};   
         end
         A2 = A2.'; 
        for j3 = 1:N   
        A3 = ex_basis(2,j3+1)*mu_full_t{1};  
         for ell = 2:N
                A3 = A3 + ex_basis(1+ell,j3+1)*mu_full_t{ell};  
         end
         A3 = A3.';
            for j4 = 1:N %this is the final interaction
        A4 = ex_basis(2,j4+1)*mu_full_t{1};  
         for ell = 2:N
                A4 = A4 + ex_basis(1+ell,j4+1)*mu_full_t{ell};  
         end
                A4 = A4.'; %flip these so they are correct shape
      tic            
            for tlp2 = 1:length(t2_lp_rng)    
                t2 =  t2_lp_rng(tlp2)-1; %this is the position along t_range
                for tlp3  = 1:length(t3_lp_rng)
                    t3 =  t2_lp_rng(tlp2)-1;
                    %take a reduced range of t1 to support this, if it is
                    %empty then don't both doing anything
                    inc_lg = (t1_lp_rng+t2+t3<size(A4,2));
                    t1_lp_red = t1_vc(inc_lg);
  %tmp = third_order_fn(A1,A2(tlp1,:),A3(tlp1+tlp2,:),A4(t3_rng+tlp1+tlp2,:));     
  if ~isempty(t1_lp_red) %
tmp = third_order_fn(A1,A2(:,t1_lp_red),A3(:,t1_lp_red+t2),A4(:,t1_lp_red+t2+t3));
Stmp(t1_lp_red,tlp2,tlp3) = tmp ;
  end
                end
            end
        toc
        
          kk1 = [-k_u;k_u;k_r]; kk2 = [k_u;-k_u;k_r];
  [fct_x1,fct_y1,~] = full_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol,kk1);
  [fct_x2,fct_y2,~] = full_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol,kk2);    
             
  Sx1 = Sx1 + fct_x1*Stmp;   Sy1 = Sy1 + fct_y1*Stmp;    
  Sx2 = Sx2 + fct_x2*Stmp;   Sy2 = Sy2 + fct_y2*Stmp;    
  
            end
        end
    end
end

% %% Electric field and such like
% 
% 
% E_0_u = 0.1; %peak electric field Volts / Angstrum
% %there is no standard way to convert e*Angstrom into cm^-1, 
% % hence dipole moments are not scaled from these units it is simply
% %important that one keeps 1 eV = 8062.4 cm^-1, standard units for E field
% %are V/m, hence we take the units of E to be 
% % E_scaled = (E in Volts / Angstrum) * 8064.2, 
% %that way the results stay in cm^(-1)
% E_0_u = E_0_u* 8064.2 ;
% %normalise envelopes to unity intensity 
% Env_r = exp(-(t-t3-t_sep)^2/2/tau_r^2)/sqrt(sqrt(pi)*tau_r/2);  
% Env_r_ft = (pi^(-1/4)*exp(-(om3^2*tau_r^2))/2)*sqrt(abs(tau_r));
% 
% Env_r_ft_bc = (pi^(-1/4)*exp(-((om3 + om_r)*(t*2*1i - t_sep*2*1i +...
%     om3*tau_r^2 + om_r*tau_r^2))/2)*abs(tau_r))/tau_r^(1/2);
% 
% Env_r_ft_fwd = (pi^(-1/4)*exp(-((om3 - om_r)*(t*2*1i - t_sep*2*1i +...
%     om3*tau_r^2 - om_r*tau_r^2))/2)*abs(tau_r))/tau_r^(1/2);
% 
% %freq envelope is a Gaussian centred at om_r with width 1/tau_r
% %considered to have unit magnitude as it's linear in this
% Env_u2 = E_0_u*exp(-(t-t3-t2)^2/2/tau_u^2)/sqrt(sqrt(pi)*tau_u); 
% Env_u1 = E_0_u*exp(-(t-t3-t2-t1)^2/2/tau_u^2)/sqrt(sqrt(pi)*tau_u);


 %%  Evaluate signal in short probe limit using relative variables
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
  tstep = 0.01; tmax = 4; %change these values according to all the bath parameters
  tsep_range_des = 0:tstep:tmax; %ps
  t_dimless = tsep_range_des*convfact;

  om_c_range = (-1/(tstep*convfact):1/(tmax*convfact):1/(tstep*convfact))*2*pi/sqrt(2);
  %om_c_tilde = linspace(-om_c_range(end)/2,om_c_range(end)/2,length(om_c_range));
  %om_c_dimless = om_c_tilde*2/tau;
  %  om_c_range = om_c_tilde  - om_c_tilde(1);
 % nextfft = nextpow2(length(om_c_range));  NN = 2^nextfft; 
   tmp=(0:length(om_c_range)-1).';
phase_shift=exp(-1i*2*pi*(length(om_c_range)-1)/2*tmp/length(om_c_range)); 
% calculating phase shift from the fact the FFT doesn't centre the data on
% zero, this must also be padded out to be the correct size
phase_shift=repmat(phase_shift,1,length(om_plot_rng));

  
   quad_order= 1;   tau = tau_u;
   om_plot_rng =abs(E0) + linspace(H0ex(2,2)-10/tau_u,H0ex(3,3)+10/tau_u,numpoints);
  % om_plot_rng = om_plot_rng(60:2:80); %take smaller range to test  
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
  tmp1x = zeros(length(om_c_range),length(om_plot_rng));
  %tmp2x = tmp1x;  tmp1y = tmp1x; tmp2y = tmp1x;
  tic

 for lp = 1:length(om_u_range)-1
          omu = om_u_range(lp);
          omplus = omm + sqrt(2)*omu;  omminus= omm - sqrt(2)*omu;
          %quadrature points for each of the exponentials
          tmp1x = 0*tmp1x; %tmp2x = tmp1x;  
          tmp1y = tmp1x; %tmp2y = tmp1x;
      for lp1 = 1:quad_order
          %fft defaults to dimension one, fft is over om_c
    %P3x2(:,:,lp) = P3x2(:,:,lp) +weig(lp1)*(fft((temp_fn_x(omm1,omm2,omplus(lp1))...
    %                 +temp_fn_x(omm1,omm2,omminus(lp1))).*exp_fct)...
    %                 + ifft((temp_fn_x(omm1,-omm2,omplus(lp1))...
    %                 +temp_fn_x(omm1,-omm2,omminus(lp1))).*exp_fct) );
%     P3y2(:,:,lp) = P3y2(:,:,lp) +weig(lp1)*(fft(temp_fn_y_minus(omm1,omm2,omplus(lp1))...
%                     +temp_fn_y_plus(omm1,omm2,omminus(lp1)).*exp_fct)+...
%                     ifft(temp_fn_y_minus(-omm1,omm2,omplus(lp1))...
%                     +temp_fn_y_plus(-omm1,omm2,omminus(lp1)).*exp_fct));      
   %  P3y2(:,:,lp) = P3y2(:,:,lp) +weig(lp1)*(fft((temp_fn_y_plus(omm1,omm2,omplus(lp1))...
   %                 +temp_fn_y_minus(omm1,omm2,omminus(lp1))).*exp_fct)+...
   %                 ifft((temp_fn_y_plus(omm1,-omm2,omplus(lp1))...
   %                 +temp_fn_y_minus(omm1,-omm2,omminus(lp1))).*exp_fct));          
     tmp1x = tmp1x +weig(lp1)*(temp_fn_x(omm1,omm2,omplus(lp1))...
                    +temp_fn_x(omm1,omm2,omminus(lp1)));
     %tmp2x = tmp2x   +weig(lp1)*( temp_fn_x(omm1,-omm2,omplus(lp1))...
     %                +temp_fn_x(omm1,-omm2,omminus(lp1)));
     tmp1y = tmp1y +weig(lp1)*(temp_fn_y_plus(omm1,omm2,omplus(lp1))...
                    +temp_fn_y_minus(omm1,omm2,omminus(lp1)));
     %tmp2y = tmp2y   +weig(lp1)*( temp_fn_y_plus(omm1,-omm2,omplus(lp1))...
     %                +temp_fn_y_minus(omm1,-omm2,omminus(lp1)));     
                 %these functions are what is obtained after integrating 
      end
      tmp1x = fftshift(fft(tmp1x.*exp_fct)./phase_shift,1)/length(om_c_range);
      tmp2x = fftshift(fft(tmp1y.*exp_fct)./phase_shift,1)/length(om_c_range);
      P3x2(:,:,lp) = tmp1x(floor(end/2)+1:end,:); 
      P3y2(:,:,lp) = tmp2x(floor(end/2)+1:end,:); 
 end 
 toc
%P3x2 = 1i*P3x2; P3y2 = 1i*P3y2;
 %%  Plot quantities of interest
 % 9 Level Color Scale Colormap with Mapping to Grayscale for Publications.
%
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);

 lambda_nm = 10^7./om_plot_rng;
 om_pad = kron(ones(size(P3x2,1),1),om_plot_rng);
 %should also include 1/n(omega) factor but that is probably just that of
 %water
lp = 1;
 
  Dalpha_J= (8*pi^2.*om_pad).*imag(P3x2(:,:,lp));
  p_rng = tsep_range_des>0.25 & tsep_range_des<0.75;
  figure
  pcolor( lambda_nm ,tsep_range_des(p_rng) ,Dalpha_J(p_rng,:))
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('Absorption change \Delta \alpha, for different probe wavelengths and pulse seperations');
 colormap(CMRmap)
colorbar

Deta_J  = (16*pi^2.*om_pad).*real(P3y2(:,:,lp)); 
  figure
  pcolor( lambda_nm ,tsep_range_des(p_rng)  ,Deta_J(p_rng,:))
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('CD shift \Delta \eta, for different probe wavelengths and pulse seperations');
colormap(CMRmap)
colorbar
%%
 lp=1
  Dn= (8*pi^2.*om_pad).*real(P3x2(:,:,lp));
  figure
  pcolor( lambda_nm ,tsep_range_des(1:floor(end/5)) ,Dn(1:floor(end/5),:))
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('Change in refractive index (AU), for different probe wavelengths and pulse seperations');
 colormap(CMRmap)
colorbar

 Ddelta_J  = (8*pi^2.*om_pad).*imag(P3y2(:,:,lp));
  figure
  pcolor( lambda_nm ,tsep_range_des(1:floor(end/5))  ,Ddelta_J (1:floor(end/5),:))
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('OR shift \Delta \delta, for different probe wavelengths and pulse seperations');
colormap(CMRmap)
colorbar
%% Cut through of the CD signal
temp = abs(Deta_J(1,:)); temp2 = diff(temp); temp3 = sign(temp2);
temp4 = diff(temp3);  temp5 = temp4~=0; 
%pad these local maxima out with a point either side
temp5([temp5(2:end),false]) = true; temp5([false,temp5(1:end-1)]) = true; 

[a,b] = max(temp);
gap1 = H0ex(2,2)-H0ex(1,1);   [a1,b1] = min(abs(gap1-om_plot_rng));
%find the max of temp near here

gap2 = H0ex(3,3)-H0ex(1,1);   [a2,b2] = min(abs(gap2-om_plot_rng));
%gap3 = om_0{1}; dom_rng = om_plot_rng(2)-om_plot_rng(1); bb = floor(om_0{1}/dom_rng);

figure
plot(tsep_range_des(1:floor(length(tsep_range_des)/2)),Deta_J(1 ...
:floor(length(tsep_range_des)/2),temp5))
% time cut on resonance
 %%
figure
plot(lambda_nm,Deta_J([1,11,21],:))

%% Freq Cut through line of the absorption signal


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