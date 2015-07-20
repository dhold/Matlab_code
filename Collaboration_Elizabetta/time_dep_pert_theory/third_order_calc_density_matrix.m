%Used to calculate the third order response function for a dimer using
%Redfield theory and explicitly included, damped Harmonic oscillator modes
%get inverse cm units
Temp = 300; %temp in Kelvin
[convfact, beta,speed_unit]= inv_cm_unit_sys(Temp);

kpr = [0,0,1];   epr = [1,0,0];   bpr = [0,1,0];
%probe unit vectors

theta = atan(sqrt(2)); %angle between pump and probe
kpu = [sin(theta),0,cos(theta)];   epu = [1,0,0]; bpu=  cross(kpu,epu);
%pump unit vectors

tau_u = 150/1000*convfact; %pulse SD in inverse cm

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
rho_fl = (tmp2(end,:)).';
rho_0 = reshape(rho_fl,sqrt(length(rho_fl)),sqrt(length(rho_fl)));
%can check convergence if I want, but this should be a sufficient amount of
%time and close enough for spectroscopy!

%% Calculate linear order density matrix

%set number of points to take and how many to save
points_to_save = 2^14; t_end_ps = tau_u/convfact*10; t_end = t_end_ps*convfact;
t_range = linspace(0,t_end,points_to_save);
t_step = (t_range(2)-t_range(1));
rel_freq_range = linspace(0,2*pi/t_step ,points_to_save);
f_spacing = rel_freq_range (2)-rel_freq_range (1);

rho_1 = zeros(length(rho_fl),length(t_range),N);


for e1 = 1:N
    temp = zeros(length(H0ex)); temp(1,e1+1)=1;
    mu_unit{e1} = temp;
    %same thing in full Hilbert space
    mu_hilb{e1} = kron(temp,eye(sz2));mu_hilb{e1} = mu_hilb{e1}+mu_hilb{e1}';
    %project this to the exciton basis which is the basis the operator is
    %in
    op1 = indiv_op{5}'*mu_hilb{e1}*indiv_op{5};
    
    temp = op1*rho_0 - rho_0*op1; %first commutator
    %reduced system to consider
    temp = reshape(temp,numel(temp),1);
    
    output_DE_fun(points_to_save,temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-70);
    [tmpp1,tmpp2]=ode45(de_fn,[0,t_end],temp,options);
    
    
    [tmp1,tmp2]  =  output_DE_fun(points_to_save,temp,'get_data');
    
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    tmp2 = interp1(tmp1,tmp2,t_range);
    
    rho_1(:,:,e1) = tmp2.'; %tmp2 output in the wrong shape for this
    %note this is rho_1 in the exciton basis
end

%% Calculate linear spec terms alternative way
if 1==0
om_char =eig(H0_renorm);  %use to scale for easier numerics
%set number of points to take and how many to save
points_to_save = 2^14; t_end_ps = 2; t_end = t_end_ps*convfact;
t_range = linspace(0,t_end,points_to_save);
t_step = (t_range(2)-t_range(1));
rel_freq_range = linspace(0,2*pi/t_step ,points_to_save);
f_spacing = rel_freq_range (2)-rel_freq_range (1);

rho_1 = zeros(numel(rho_fl),length(t_range),N);

    tmpa  = zeros(sz1*sz2); tmpa(1:sz2,sz2+1:sz2*(sz1))=1;
    [aa,ba,sa] = find(tmpa);
	tmpa = logical(reshape(tmpa,sz2^2*sz1^2,1));
%deco_op_red = decoherence_op(tmp,tmp);

    supop_a = supop(tmpa,tmpa);
    
    tmpb  = zeros(sz1*sz2); tmpb(1:sz2,sz2+1:sz2*(sz1))=1;
    [ab,bb,sb] = find(tmpb);
	tmpb = logical(reshape(tmpb,sz2^2*sz1^2,1));
%deco_op_red = decoherence_op(tmp,tmp);

    supop_b = supop(tmpb,tmpb);
    
    rho_1a = zeros(numel(tmpa),length(t_range),N);
    rho_1b = zeros(numel(tmpb),length(t_range),N);
    
for e1 = 1:N
    

    temp = zeros(length(H0ex)); temp(1,e1+1)=1;
    mu_unit{e1} = temp; %#ok<*SAGROW>
    %same thing in full Hilbert space
    mu_hilb{e1} = kron(temp,eye(sz2));
    mu_hilb{e1} = mu_hilb{e1}+mu_hilb{e1}';

    temp1 = mu_hilb{e1}*rho_0; %first commutat
    temp2 = - rho_0*mu_hilb{e1};
    %reduced system to consider

    
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-70);
    
    supop_red_scale  = sparse(supop_red - om_char(e1)*1i*eye(length(supop_red)));
    de_fn2 = @(t,v) supop_red_scale*v; %scaled for smaller phase oscillations
    
    temp1_red = temp1(tmp);
    output_DE_fun(points_to_save,temp1_red,'notsavingnayway'); 
    
    ode45(de_fn2,[0,t_end],reshape(temp1_red,numel(temp1_red),1),options);
    
    [tmp1,tmp2]  =  output_DE_fun(points_to_save,temp1_red,'get_data');
    
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    tmp2a = interp1(tmp1,tmp2,t_range);
    
    rho_1a(:,:,e1) = (tmp2a.') .*repmat(exp(om_char(e1)*1i.*t_range),numel(tmp),1) ; 
    %tmp2 output in the wrong shape for this, have to include oscillation
    %back in 

    %same for other half
    
    supop_red_scale  = sparse(supop_red - om_char(e1)*1i*eye(length(supop_red)));
    de_fn2 = @(t,v) supop_red_scale*v;
    
    temp2_red = temp2(tmp);
    output_DE_fun(points_to_save,temp2_red,'notsavingnayway'); 
    
    ode45(de_fn2,[0,t_end],reshape(temp2_red,length(temp2_red),1),options);
    
    [tmp1,tmp2]  =  output_DE_fun(points_to_save,temp2_red,'get_data');
    
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    tmp2b = interp1(tmp1,tmp2,t_range);    
    
    rho_1b(:,:,e1) = tmp2b.' .*repmat(-exp(om_char(e1)*1i.*t_range),numel(tmp),1);
%rho_1a and rho_1b make up the non zero sections of the FO density matrix
end
end
%% Calculate linear spec stuff based on this

om_char = mean(eig(H0_renorm)); %characteristic frequency
k_pr = 2*pi*om_char; %assumed to go along z%
shift_param = om_char;
ft_scale_shift =exp(-1i*shift_param*t_range); %shift ft along
S_1 = zeros(2,length(t_range));
tic
for j1 = 1:N
    for j2 = 1:N
        xx_av = dot(mu(j1,:),mu(j2,:))/3;
        yx_av = 1i*k_pr*dot(mu(j1,:),cross(mu(j2,:),R(j1,:)-R(j2,:)))/6;
        
        %project reshaped rho_1 back to site basis       
        temp = mtimesx(reshape(rho_1(:,:,j1),sz1*sz2,sz1*sz2,...
                    length(t_range)),indiv_op{5});
        temp = mtimesx(indiv_op{5},'T',temp);
        temp2 = mtimesx(mu_hilb{j2},temp); %system after application of V_j2
        %need to perform a trace over the matrix in the first two
        %dimensions, sadly matlab doesn't let you do this so I use a custom
        %function
        temp3 = reshape(2*imag(diagsum(temp2,1,2)),size(t_range));
        
        if xx_av ~=0
        S_1(1,:) = S_1(1,:) + ft_scale_shift.*temp3.*xx_av;
        end
        if yx_av ~=0
        S_1(2,:) = S_1(2,:) + ft_scale_shift.*temp3.*yx_av; 
        end
        
    end
end
   toc     
     
   S_1_ft = fft(S_1,[],2);
 %%  plot to check etc
   
   lf = length(rel_freq_range);
   f_sel = [((lf - floor(lf/20)) : lf),1:floor(lf/20)]; %range to take from fft
   f_plot = [-rel_freq_range(floor(lf/20)+1:-1:2),rel_freq_range(1:floor(lf/20)+1)] + shift_param; %shift freq range
   %plot the figures of important things like 
   
   figure
plot(f_plot,f_plot.*real(S_1_ft(1,f_sel))) %ref index
figure
plot(f_plot,f_plot.*imag(S_1_ft(1,f_sel))) %abs
   
  figure
plot(f_plot,f_plot.*real(S_1_ft(2,f_sel))) %cd
figure
plot(f_plot,f_plot.*imag(S_1_ft(2,f_sel))) %or

%% Calculate two time state of density matrix for each combo of dipole elements
%this effectively describes the 

leng1 = 2^4;  %this length must statisfy  
%t_end * numpoints /t1_max*leng1 = integer and t1_max <= t_end_ps
%so these points are included in the range rho_1 is calculated for!

leng2 = 200;  %doesn't need to be power of 2 as no FTing
%leng2 will relate to number of pump-probe seperation times I want
t1_max = 10*tau_u; %this is the time between first and second 
%interactions, and should be of the order of the pulse width, this is quite
%short for the most part, might be better just to calculate it again.

t2_max = 1*convfact; %this is the time between the second and third 
%interaction, choose the length of delays we include

t1_rng = t_range(t_range<=t1_max);
t1_map = 1:length(t1_rng)/leng1:length(t1_rng); leng1 = length(t1_map);
t1_rng = t_range(t1_map); t1_step = (t1_rng(2)-t1_rng(1));

t2_rng = linspace(0,t2_max,leng2); t2_step = (t2_rng(2)-t2_rng(1));

w1_rng = linspace(0,2*pi/t1_step ,leng1);
w1_step = w1_rng(2)-w1_rng(1);

rho_2 = zeros(numel(rho_fl),leng2,leng1,N,N);  %not actually the second 
%order density matrix, would be when integrated over t1 and t2 from 0 to
%infinity weighted by the electric field from the pump pulse
%rho_2(:,t2,t1,j2,j1)
coeffs = ex_basis(2:end,2:end);


for j1 = 1:N
    for j2 = 1:N
        op2 = indiv_op{5}'*mu_hilb{j2}*indiv_op{5};
        %calculate initial condition to be propogates
        rho_2_init = mtimesx(op2,reshape(rho_1(:,t1_map,j1),...
            sz1*sz2,sz1*sz2,leng1));
        %next take away commutator term
        rho_2_init = rho_2_init-mtimesx(reshape(rho_1(:,t1_map,j1),...
            sz1*sz2,sz1*sz2,leng1),op2);        
        
        %this is the initial condition for all the time evolutions
        for tlp = 1:leng1
         %loop over t1 values and calculate evolution over t2   
        
        temp = reshape(rho_2_init(:,:,tlp),sz1^2*sz2^2,1);
    output_DE_fun(leng2 ,temp,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-70);
    ode45(de_fn,[0,t2_max],temp,options);
  
    [tmp1,tmp2]  =  output_DE_fun(leng2 ,temp,'get_data');
    
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    tmp2 = interp1(tmp1,tmp2,t2_rng);
    
    rho_2(:,:,tlp,j2,j1) = tmp2.';
        end
    end
end

%% Calculate actual second order density matrix at time t for given pump

tt_range = (0:0.002:2)*convfact;
tmp = 1:numel(rho_fl);
[tt1,tt2,~] = meshgrid(t1_rng,t2_rng,tmp);

%exp_fct = exp(-((tt-tt2-tt1).^2+(tt-tt2).^2)/tau_u^2/2)/tau_u/sqrt(2*pi);
om_u = om_char; %probe frequency
k_u = 2*pi*om_u*kpu; k_r = 2*pi*om_char*kpr; %average amplitudes of these para

rho_sec_order_p = zeros(numel(rho_fl),length(tt_range),N,N);
rho_sec_order_m = zeros(numel(rho_fl),length(tt_range),N,N);
for j1 = 1:N
    for j2 = 1:N
for t_lp = 1:length(tt_range)
    
    tt = tt_range(t_lp);
    exp_fct = exp(-((tt-tt2-tt1).^2+(tt-tt2).^2)/tau_u^2/2)/tau_u/sqrt(2*pi);
    %exp_fct(:,n,any) is exp(-((tt-t2_rng-t1_rng(n)).^2+(tt-t2_rng).^2)/tau_u^2/2)/tau_u/sqrt(2*pi)
    %exp_fct(n,:,any) is exp(-((tt-t2_rng(n)-t1_rng).^2+(tt-t2_rng(n)).^2)/tau_u^2/2)/tau_u/sqrt(2*pi)
    rho_sec_order_p(:,t_lp,j1,j2) = trapz(t2_rng,trapz(t1_rng,...
        permute(rho_2(:,:,:,j2,j1),[2,3,1]).*exp_fct.*exp(1i*om_u*tt1),2),1);
    rho_sec_order_m(:,t_lp,j1,j2) = trapz(t2_rng,trapz(t1_rng,...
        permute(rho_2(:,:,:,j2,j1),[2,3,1]).*exp_fct.*exp(-1i*om_u*tt1),2),1);    
end
    end
end
%% Calculate the signal from this in the limit of a short probe (delta fn)

%t_sep_rng_fs = 40:10:1000;
%t_sep_rng = t_sep_rng_fs /1000 *convfact;
pol = {epu,epu,epr};
%choose sep range that matches points in tt_range;
t_sep_map = tt_range>0.04*convfact & tt_range<1*convfact;
t_sep_rng = tt_range(t_sep_map);
len_sep = 100;
t_spc = (1:ceil(length(t_sep_rng)/len_sep):length(t_sep_rng));
t_sep_rng = t_sep_rng(t_spc);
len_sep = length(t_sep_rng);
t_sep_map = find(t_sep_map,1,'first') + t_spc;

 kk1 = [-k_u;k_u;k_r]; kk2 = [k_u;-k_u;k_r];
S_1_Neq = zeros(length(tt_range),length(t_sep_rng),2);
%temp = S_1_Neq;

for j1 = 1:N
    for j2 = 1:N
        for j3 = 1:N    
    op3 = indiv_op{5}'*mu_hilb{j3}*indiv_op{5}; %->exciton basis
    tic
        for t_lp = 1:len_sep 
%loop over seperations range 
    %t_sep = t_sep_rng(t_lp);
    t_map = t_sep_map(t_lp);
%propogate the matrix in time after the action of the 
temp2 = op3*reshape(rho_sec_order_p(:,t_map,j1,j2),sz1*sz2,sz1*sz2)-...
         reshape(rho_sec_order_p(:,t_map,j1,j2),sz1*sz2,sz1*sz2)*op3;
%reshape the section to a matrix shape

t_red = tt_range(t_map:end)-tt_range(t_map); %only reduced t_range from after pulse hits

output_DE_fun(length(t_red),rho_sec_order_p(:,t_map,j1,j2),'notsavingnayway');    
%tic
ode45(de_fn,[0,max(t_red)],temp2,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_red),temp2,'get_data');
    
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    tmp2p = interp1(tmp1,tmp2,t_red).';    
%tim1 =  toc   
output_DE_fun(length(t_red),rho_sec_order_m(:,t_map,j1,j2),'notsavingnayway');    

temp2 = op3*reshape(rho_sec_order_m(:,t_map,j1,j2),sz1*sz2,sz1*sz2)-...
         reshape(rho_sec_order_m(:,t_map,j1,j2),sz1*sz2,sz1*sz2)*op3;
%reshape the section to a matrix shape

ode45(de_fn,[0,max(t_red)],temp2,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_red),temp2,'get_data');
    
    tmp3 = [true; diff(tmp1)~=0] & [true; tmp1(2:end)~=0];
    tmp1 = tmp1(tmp3); tmp2 = tmp2(tmp3,:);
    tmp2m = interp1(tmp1,tmp2,t_red).';      
%tic
        for j4=1:N
           op4 = indiv_op{5}'*mu_hilb{j4}*indiv_op{5};  
           
  [fct_x1,fct_y1,~] = full_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol,kk1);
  [fct_x2,fct_y2,~] = full_cont_3rd_order2(mu,R,[j1,j2,j3,j4],pol,kk2); 
          %S_1_Neq(1:t_map-1,t_lp,:) = 0; %no signal before pulse hits
          
          tmp3p = mtimesx(op4,reshape(tmp2p,sz1*sz2,sz1*sz2,length(t_red)));
          tmp3p = diagsum(tmp3p,1,2); %trace
          tmp3m = mtimesx(op4,reshape(tmp2m,sz1*sz2,sz1*sz2,length(t_red)));
          tmp3m = diagsum(tmp3m,1,2); %trace          
          
          S_1_Neq(t_map:end,t_lp,1) = S_1_Neq(t_map:end,t_lp,1)+...
              fct_x1*tmp3p + fct_x2*tmp3m;
          S_1_Neq(t_map:end,t_lp,2) = S_1_Neq(t_map:end,t_lp,2)+...
              fct_y1*tmp3p + fct_y2*tmp3m;
          
        end
        %tim2 = toc
        end
        toc
        end
    end
end
S_1_Neq = -1i*S_1_Neq;
%% Take FFT of S_1_Neq as usual to get the response in freq space

shift_param = om_char*1.0;
ft_scale_shift =exp(-1i*shift_param*tt_range); %shift ft along
ft_scale_shift = reshape(ft_scale_shift,length(tt_range),1);
ft_scale_shift = repmat(ft_scale_shift,1,length(t_sep_rng),2);

S_1_Neq_FT = fftshift(fft(S_1_Neq.*ft_scale_shift),1); 
%fft automatically takes dimension 1 fftshift does not so need to say so

om_rng = shift_param + 2*pi/max(tt_range)*(-(length(tt_range)-1)/2:(length(tt_range)-1)/2);
%due to the fact that S_1_Neq is zero when ever t<t_sep the ft will start
%with a phase of e^i om t_sep, which will cancel when the factor of E_r(om)
%which we divide this by.  
%om2_rng = 2*pi/max(tt_range)*(0:length(tt_range)-1);
[tt_t,om_om] = meshgrid(t_sep_rng,om_rng); %tt_t varies horizontally (dim 2)
other_fct = repmat(exp(1i*tt_t.*om_om),1,1,2);

S_1_Neq_FT = S_1_Neq_FT.*other_fct;

% lf = length(om_rng );
%   f_sel2 = [((lf - floor(lf/2)) : lf),1:floor(lf/2)]; %range to take from fft
%   f_plot2 = [-om_rng(floor(lf/2)+1:-1:2),om_rng(1:floor(lf/2)+1)] + shift_param; %shift freq range
  f_sel2 = om_rng> om_char - 3000 & om_rng < om_char + 3000;
  f_plot2 = om_rng(f_sel2);
%% Plot things!
%frequency cut through (inverse cm) for the minimum seperation time considered
%note that this path only includes those with the pump interacting first

%   figure
%plot(f_plot2,f_plot2.'.*real(S_1_Neq_FT(f_sel2,1,1))) %ref index
figure
plot(f_plot2,f_plot2.'.*imag(S_1_Neq_FT(f_sel2,1,1))) %abs
[~,b]=max(abs(imag(S_1_Neq_FT(f_sel2,1,1))));
ylabel('\Delta \alpha au')
xlabel('frequency, cm^{-1}')  
  figure
plot(f_plot2,f_plot2.'.*real(S_1_Neq_FT(f_sel2,1,2))) %cd
xlabel('\Delta \eta, cm^{-1}')  
ylabel('CD au')
%figure
%plot(f_plot2,f_plot2.'.*imag(S_1_Neq_FT(f_sel2,1,2))) %or
%%
om_pad = repmat(f_plot2.',1,length(t_sep_rng));
Dalpha_J= (8*pi^2.*om_pad).*imag(S_1_Neq_FT(f_sel2,:,1));
Deta_J  = (16*pi^2.*om_pad).*real(S_1_Neq_FT(f_sel2,:,2));

figure
plot(t_sep_rng/convfact,Dalpha_J([b-5,b,b+5],:)) %abs
ylabel('\Delta \alpha au')
xlabel('t_sep, ps')  
  figure
plot(t_sep_rng/convfact,Deta_J([b-5,b,b+5],:) ) %cd
ylabel('\Delta \eta')  
xlabel('t_sep, ps')

%%
figure
[toplot1,toplot2] = ezfft(t_sep_rng-t_sep_rng(1),Dalpha_J(67,:)-mean(Dalpha_J(67,:)));

plot(t_sep_rng/convfact,Dalpha_J(67,:)) %abs
ylabel('\Delta \alpha au')
xlabel('t_sep, ps')  
  figure
plot(t_sep_rng/convfact,Deta_J(67,:) ) %cd
ylabel('\Delta \eta')  
xlabel('t_sep, ps')

figure
plot(toplot1,toplot2)
%%  2D plots
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);
om_pad = repmat(f_plot2.',1,length(t_sep_rng));

  Dalpha_J= (8*pi^2.*om_pad).*imag(S_1_Neq_FT(f_sel2,:,1));
  p_rng = t_sep_rng>0.25*convfact & t_sep_rng<1*convfact;
  figure
  pcolor(f_plot2 ,t_sep_rng(p_rng) ,Dalpha_J(:,p_rng).')
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('Absorption change \Delta \alpha, for different probe wavelengths and pulse seperations');
 colormap(CMRmap)
colorbar

Deta_J  = (16*pi^2.*om_pad).*real(S_1_Neq_FT(f_sel2,:,2));
  figure
  pcolor( f_plot2,t_sep_rng(p_rng) ,Deta_J(:,p_rng).')
set(gcf, 'renderer', 'zbuffer');
shading flat
xlabel('wavelength, nm')
ylabel('time seperation, ps')
title('CD shift \Delta \eta, for different probe wavelengths and pulse seperations');
colormap(CMRmap)
colorbar

error('break code')
%% Plot oscillations in
%% Old code
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