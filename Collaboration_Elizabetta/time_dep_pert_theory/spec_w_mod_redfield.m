Temp = 300; %temp in Kelvin
[convfact, B,speed_unit]= inv_cm_unit_sys(Temp);

%% B820_bacteriochlorophyll_dimer
flenme ='B820_bacteriochlorophyll_dimer_params.mat'; 
%insert name of file containing the system parameters
    fle = open(flenme);
        H_site = fle.H_site; N = length(H_site);
        [ex_basis,H0ex] = eig(H_site);
        
mu = fle.mu;  pdm = fle.pdm; %set to empty if not known
take_same_amp = true;
if take_same_amp
%assume all are the same amplitude
for k = 1:N
    mu(k,:) = mu(k,:)/norm(mu(k,:)); %site based mu
end
end
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
om_0 = fle.om_0; lambda = fle.lambda;  gam_rng = fle.gamma;
%explicitly included modes

om_vib = [om_0{:}];

displ = [[0;0],[sqrt(2*lambda{1}/om_0{1}) ; 0],...
            [0 ; sqrt(2*lambda{2}/om_0{2})],...
    [sqrt(2*lambda{1}/om_0{1});sqrt(2*lambda{2}/om_0{2})]];
%displacements of states, ground, excited states and finally double excited
sz =4; %total number of electronic states

sd_shift = [200,200]; %standard deviation of site energy fluctuations 
H_el = zeros(4);  H_el(2:3,2:3) = H_site;  H_el(4,4) = H_el(2,2)+H_el(3,3);
H_single = H_el(2:N+1,2:N+1); H_double =  blkdiag(H_el(1,1),H_el(2+N:end,2+N:end));    
%also can include interaction between ground and double excited
[ex_basis_2,H_doub_ex] = eig(H_double);
%%  Calculate dipole averages for the pump probe geometry
%Assuming no static disorder impacts the dipoles (apart from (random) 
%molecular alignment) we can precompute these values.  The actual values 
%will be averaged over exciton states though, which depend on site energies

 typ_op = 10^7/820; %typical frequency
 theta = atan(1/sqrt(2));
 kpr = [0,0,1]; kpu = [sin(theta),0,cos(theta)];
 
k_u = 2*pi*abs(typ_op)*kpu; k_r = 2*pi*abs(typ_op)*kpr; %average amplitudes of these para
kk1 = [-k_u;k_u;k_r;-k_r]; kk2 = [k_u;-k_u;k_r;-k_r]; %order of interaction (left to right)
%take probe polarization wavevectors along x
%pol{3} = [1,0,0]; % pol_left = [1;-1i;0]/sqrt(2);  pol_right = [1;1i;0]/sqrt(2);

pol{1}=  [1,0,0];  pol{2} = pol{1} ; %others along x
%pol{1} = [1/sqrt(2),1/sqrt(2),0]; pol{2} = pol{1} ; %45 degrees to probe apparently
pol_linear = [pol{1};pol{2};[1,0,0];[1,0,0]];
pol_L = [pol{1};pol{2};[1/sqrt(2),+1i/sqrt(2),0];[1/sqrt(2),-1i/sqrt(2),0]];
pol_R = [pol{1};pol{2};[1/sqrt(2),-1i/sqrt(2),0];[1/sqrt(2),+1i/sqrt(2),0]];

        % alpha_av = zeros(N,N,N,N); CD_av = zeros(N,N,N,N);
         alpha_L_4 = zeros(N,N,N,N); alpha_R_4 = zeros(N,N,N,N);
         alpha_L_5 = zeros(N,N,N,N); alpha_R_5 = zeros(N,N,N,N);   
         alpha_lin_pol = zeros(N,N,N,N);
for j1 = 1:N %pre save these
   for j2 = 1:N
       for j3 = 1:N
           for j4=1:N
               
         jj = [j1,j2,j3,j4];             

    alpha_lin_pol(j1,j2,j3,j4) = tensor_av(mu(jj,:),pol_linear); 
    alpha_L_4(j1,j2,j3,j4) = tensor_av(mu(jj,:),pol_L);
    alpha_R_4(j1,j2,j3,j4) = tensor_av(mu(jj,:),pol_R);
    
    for j= 1:4
    alpha_L_5(j1,j2,j3,j4) = alpha_L_5(j1,j2,j3,j4) +...
        tensor_av([mu(jj,:);R(jj(j),:)],[pol_L;1i*kk1(j,:)]);   
    alpha_R_5(j1,j2,j3,j4) = alpha_R_5(j1,j2,j3,j4) +...
        tensor_av([mu(jj,:);R(jj(j),:)],[pol_R;1i*kk1(j,:)]);       
    end
           end
       end
   end
end
    alpha_av = (alpha_L_4 + alpha_R_4)/2; 
    lg = alpha_L_5 + alpha_R_5; 
    if any(lg(:)~=0)
        warning('5th order contibution to absorbtion')
      alpha_av =  alpha_av+ lg /2;
    end


    CD_av = (alpha_L_5 - alpha_R_5); 
   lg2 = alpha_L_4 - alpha_R_4;  
   if any(lg2(:)~=0)
       warning('4th order contibution to CD')
    CD_av =  CD_av+ lg2;
   end
   
%% Use function which gives doorway and window functions
tau_u_fs = 150; tau_r_fs = 150;
om_u_rng = 1.0e+04 * [1.1910    1.2210    1.2410    1.2710];
om_r_rng = 1.0e+04 *linspace(0.8,1.4);
t_sep_rng_fs = 0:5:3000; t_sep_rng = t_sep_rng_fs/1000*convfact;

params{1} = B; %beta  = 1/k_B T 
params{2}=gam_rng; params{3}=om_0 ; params{4} = lambda; %underdamped
params{5}=gam_dru; params{6}=lam_dru; %overdamped on each site

tau = [tau_u_fs,tau_r_fs]/1000*convfact;
tic

[D_kk,D_t,W_gg,W_kk,W_qq,R_red_markov,fock_rep]...
    = doorway_fun_simple(H_single,H_double,tau,om_u_rng,om_r_rng,t_sep_rng,params) ;

%[D_full,D_t_full,W_gg_full,W_kk_full,W_qq_full,D_markov...
% ,W_gg_markov,W_kk_markov,W_qq_markov,R_red,D_markov_test,fock_rep] = ...
%        doorway_fun_simple(...
%          H_single,H_double,tau,om_u_rng,om_r_rng,t_sep_rng,params) ;
toc
 partic_ex = zeros(size(fock_rep)); %participation of sites in exciton states
 partic_ex(2:N+1,:) = ex_basis;
tmp = double(fock_rep);

for lp =N+2:size(fock_rep,1) %state loop
    for lp2 = N+2:size(fock_rep,1)
    partic_ex(lp,:) =partic_ex(lp,:) + ex_basis_2(lp-N-1,lp2-N-1)*tmp(lp2,:);
    end
end

%% Calculate sums of dipole averages for each set of exciton transitions
%this is FAR more general that is required here but I wanted to do the most
%general case
fct_alpha_1 = zeros(N,N,N,N); %gg->ge~ge'-> gg or ee then -> eg~e'g -> e'e' or gg
%~ represents transfer between coherences by coupling to vibrations
fct_alpha_2 = zeros(N,N,N^2*(N-1)/2,N^2*(N-1)/2); 
%gg-> ge -> ee then -> ef~e'f' -> e'e'
fct_CD_1 = zeros(N,N,N,N); %gg->ge-> gg or ee then -> eg -> ee or gg
fct_CD_2 = zeros(N,N,N^2*(N-1)/2,N^2*(N-1)/2); 

for lp1=1:N %always gg->ge
    for lp2 = 1:N %ge -> gg or ee 
        for lp3 = 1:N
            for lp3a = 1:length(H_double) %can be ge or ef
                
            for lp4 =1:N
                for lp4a  = 1:length(H_double)
                
                    if lp3a==1 && lp4a == 1 %kg kg
     for e1 = 1:N
         for e2=1:N
             for e3=1:N
                 for e4 = 1:N
fctor = partic_ex(lp1,e1)*partic_ex(lp2,e2)*partic_ex(lp3,e3)*partic_ex(lp4,e4);
    fct_alpha_1(lp1,lp2,lp3,lp4) = fct_alpha_1(lp1,lp2,lp3,lp4)+...
                                    fctor*alpha_av(e1,e2,e3,e4);

    fct_CD_1(lp1,lp2,lp3,lp4) = fct_CD_1(lp1,lp2,lp3,lp4)+...
                                    fctor*CD_av(e1,e2,e3,e4);                                
                 end
             end
         end
     end
                    elseif lp3a~=1  && lp4a~=1 %kf kf
     for e1 = 1:N
         for e2=1:N
             for e3=1:N
                 for e3a=e3+1:N
                    for e4 = 1:N
                        for e4a = e4+1:N
fctor1 = partic_ex(lp1,e1)*partic_ex(lp2,e2)*partic_ex(lp3,e3)...
        *partic_ex(lp4,e4);
fctor2 = partic_ex(lp1,e1)*partic_ex(lp2,e2)*partic_ex(lp3,e3)*partic_ex(lp4,e4);

    fct_alpha_1(lp1,lp2,lp3,lp4) = fct_alpha_1(lp1,lp2,lp3,lp4)+...
                                    fctor*alpha_av(e1,e2,e3,e4);

    fct_CD_1(lp1,lp2,lp3,lp4) = fct_CD_1(lp1,lp2,lp3,lp4)+...
                                    fctor*CD_av(e1,e2,e3,e4);    
                        end
                    end
                 end
             end
         end
     end                        
                    
                    end
                end
            end
            end
        end                         
    end   
end

 %% Calculate the appropriate quantities     
  
 PB_sig = zeros(length(om_r_rng),length(om_u_rng));
 SE_sig = zeros(length(t_sep_rng),length(om_r_rng),length(om_u_rng));
 ESA_sig = zeros(length(t_sep_rng),length(om_r_rng),length(om_u_rng));

  PB_CD_sig = PB_sig;
 SE_CD_sig = SE_sig ;
 EsA_CD_sig =  ESA_sig ;
 
 for e1=1:N  %loop over all possible interactions
    for e2 = 1:N %only exciton-exciton transfer possible here
                c1 = zeros(1,N); c2 = c1; c3=c1;c4=c1;
                %work out average based on occupation of each different
                %site
                    for lp = 1:N
                        %calculate effective occupancy of sites
                        %i.e. decomposition of ex dipole op into site 
                        %somewhat general to be adapted for DE states
     c1 = c1 + ex_basis(e1,lp)*double(fock_rep(1+lp,:));
     c2 = c2 + ex_basis(e1,lp)*double(fock_rep(1+lp,:));
     c3 = c3 + ex_basis(e2,lp)*double(fock_rep(1+lp,:));
     c4 = c4 + ex_basis(e2,lp)*double(fock_rep(1+lp,:));
                    end
      fct_lin_pol = 0;   fct_alpha= 0; fct_CD= 0; fct_CD2= 0; fct_alpha2 =0;
for j1 = 1:N %I should pre save these but this is quick anyway
   for j2 = 1:N
       for j3 = 1:N
           for j4=1:N

   cont_fc = c1(j1)*c2(j2)*c3(j3)*c4(j4);
   fct_lin_pol = fct_lin_pol + cont_fc*alpha_lin_pol(j1,j2,j3,j4); 
   fct_alpha = fct_alpha + cont_fc*alpha_av(j1,j2,j3,j4); 
   fct_CD = fct_CD + cont_fc*CD_av(j1,j2,j3,j4);
   
           end
       end
   end
end      
for j = 1:length(om_u_rng)

PB_sig(:,j) = PB_sig(:,j) -om_r_rng.*W_gg(e2,:)*D_kk(e1,j)*fct_alpha;
SE_sig(:,:,j) = SE_sig(:,:,j) -repmat(om_r_rng.*W_kk(e2,:),length(t_sep_rng)).*D_t(:,e1,j)*fct_alpha;

PB_CD_sig(:,j) =PB_sig;
SE_CD_sig(:,:,j) = SE_sig ;


for q = 1:N2
ESA_sig(:,:,j) = ESA_sig(:,:,j)-repmat(om_r_rng.*W_qq(e2,:),length(t_sep_rng)).*D_t(:,e1,j)*fct_alpha;
EsA_CD_sig(:,:,j) =  EsA_CD_sig(:,:,j);
end
      
end
    end
 end