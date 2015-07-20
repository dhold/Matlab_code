function [time_units,P_t,rho_vc]= Polarization_via_HEOM...
            (Htot,QQ,mu,R,pu_angle,FWHM,E0,cc_com,cc_acom,rho_0,vv,Kap1,Kap2,...
            numpoints,trangeps,saveuptotier,save_file_name)
%  This code uses the HEOM with a time dependent Liouvillian in order to
%  calculate the density matrix of the system (molecular angle dependent)
%  at 4 different points in space (which can be summed to get the 
%  non spatially phase dependent component of the density matrix) 
%
% Unit setup
%All terms in inverse cm, conversion to sensible units achieved by
%E = hc / lamba => 6.62606957 *10^{-33} 299792458 *10 * "E in code" gives
%
%if nargin == 0  %use a select set nargin
%%
 % LL is the total Louivillian for the system of N chromophores and 
 % whatever vibrations are treated quantum mechanically.
 % QQ is the extra convergence parameter made up of all the terms of the 
 % form int_0^t d tau exp(-v (t-tau)) * ( c_1 V^X(tau) + c_2 V^O(tau) )
 %  which are evaluated via exp(-v (t-tau)) ~ delta(t-tau)/v, as v->inf
 % and are thus just markovian terms.
 % cc1 = residues at all the poles of 1/(1-e^{-beta omega}) 
 % cc2R = residues at all the poles of J(omega), contribution to V^x
 % cc2I = residues at all the poles of J(omega), contribution to V^o
 % i.e. the anti commutator part
 % each of these is paired with a frequency
 % vv1 = position of poles of 1/(1-e^{-beta omega})  times -1i
 % vv2 = position of poles of J(omega)  times -1i
 % N is the number of sites / Chromophores
 N = length(cc_com); 

    % gsinc = length(Htot) - N; %one if ground state is included

 if saveuptotier > Kap2
     warning('ordered to save up to a tier above the truncation threshold, will reduce')
     saveuptotier = Kap2
 end
 %Kappa = size(cc1,2); now delt with before the function
% Kappa e.g.= 1; %truncation parameter satisfying Kappa >>omega_0 beta hbar/2pi
 %beyond this e^(-v_(k>kappa) t)*v_k ~delta(t) is made
 % Kap1 e.g. = 4*omega_0 (frequency of significant transition), this
 % truncation parameter is used instead of / as well as Kap2, this
 % truncates the heirarchy based on cc dot vv via the approximation
 % exp(-cc.vv*t) ~ delta(t)/(cc.vv)
% Kap2 e.g.= 4; %second truncation parameter, there are Kappa + 1 explicitly 
 %treated exponentials, this number is the maximum order the sum of the
 %number of powers can take, e.g. n_1 + n_2 +... + n_Kappa <= Kap2
 % numpoints is the number of time points to save of the solved equation
 % for all the auxilary density matricies
 % tendps e.g.= 2 is the end time in pico seconds
% rho_0 is some initial condition
% saveonlyrho00 is a flag which makes the solver save only rho00
if ~exist('saveuptotier','var')
    saveuptotier =0;
end
if ~exist('use_reduced_mat','var')
    use_reduced_mat =false;
end
if ~exist('save_file_name','var')
   save_file_name  = 'saved_data.mat'; 
end
 %use_reduced_mat =1; %include only linearly independent elements of rho
 %i.e. for 2 X 2 rho: rho_vec_red rho_11, rho_21, p_22
 %tested to work ONLY for drude
 
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
% to run the code not as a function just uncomment and use each module
% Energy in Joules,  extra factor of 10 is because thats speed in cm/s
convfact = 2*pi * light_speed *length_unit  *10^(-12); %2 pi * "c" in cm/s * 1ps 
tend =  convfact*trangeps; %end time in units you get by putting everything in wavenumbers
totpoles = sum(cellfun(@length,cc_com));
tic
%% symbolic variables
%this code was written to be used for calculate the state of a density
%matrix after a pump pulse needs parameters
% tstart, polvec, mu,R, L0, 
% L0 = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));

kk = sym('k','real');  %pump mean wavevector
kvec = kk*[sin(pu_angle),0,cos(pu_angle)];
t=sym('t','real'); 
om = sym('om','real'); %pump resonant freq
phi = sym('phi','real');  %phase, occurs from k dot r =phi terms taken

syms a b c %Euler angles
ry = [cos(a),0,sin(a);0,1,0;-sin(a),0,cos(a)];
rz = [cos(b),-sin(b),0;sin(b),cos(b),0;0,0,1];
ry2 = [cos(c),0,sin(c);0,1,0;-sin(c),0,cos(c)];
rot_op = ry2*rz*ry; 

mu = Tmat*mu;  R = Tmat*R; %transform with specified rotation

%to select the correct phase component of the output polarization

sigma = FWHM/(2*sqrt(2*ln(2)));
E_env = E0*exp(-t^2/2/sigma^2)/sigma/sqrt(2*pi); %normalised envelope

LE_fwd = sym(zeros(N+1)); LE_bak = LE_fwd;
for k = 1 : N
LE_fwd(1,k+1) = dot(mu(:,k),polvec)*exp( 1i*(dot(R(:,k),kvec))); 
LE_bak(1,k+1) = dot(mu(:,k),polvec)*exp(-1i*(dot(R(:,k),kvec)));
end
LE_fwd = LE_fwd+LE_fwd.'; LE_bak=LE_bak+LE_bak.';
%note these are not actually hermitian but the sum is
LE_fwd = kron(eye(length(LE_fwd)),LE_fwd)-kron(LE_fwd.',eye(length(LE_fwd)));
LE_bak = kron(eye(length(LE_fwd)),LE_bak)-kron(LE_bak.',eye(length(LE_fwd)));
%commutator super operators in Liouville space


%% Construct the total Louville for self coupling density matricies

% For N X N A and B we have the following
% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho

 L0 = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));
% QQ as input should be the constants in the correction 
if size(QQ) == size(L0) %assume matrix given that is Q
    Q = QQ; clear QQ
else
if size(QQ,1) == 1 %assume same for every site
    QQ = repmat(QQ,N*viblvls,1);
end
 Q = L0*0;
% 
% %if include_truncation_correction
for j = 1:N
    
    Qj = zeros(N);
    Qj(j,j) = 1;  Qj = kron(Qj,eye(viblvls));eyeQ = eye(length(Qj));
    
      %note Qj^2 = Qj so I don't need to worry about this

    Q = Q + QQ(j,1).*(kron(Qj.',eyeQ) + kron(eyeQ,Qj) - 2*kron(Qj.',Qj)); 
    Q = Q  + QQ(j,2) .*(kron(Qj.',eyeQ) - kron(eyeQ,Qj)); 
    
    %reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
    %[Q_i,[Q_i,rho_n]] = (1 X Q_i*Q_i + Q_i^T*Q_i^T X 1 - 2 Q_i^T X Q_i ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
    
    % .*(kron(Vj.',eyeQ) + kron(eyeQ,Vj)); 
end
% %end
end

%in principle Q can also include high frequency decaying components from
% poles in J(omega) 

%%
TT = diag(exp(-diag(L0)*t)); %diagonal elements of system Liouvillian
%this is the projector to the interaction picture essentially

barescalefct = double(TT^(-1)*diff(TT,t)); %should be just a number!
L0_test = L0 - barescalefct; %should have no diagonal elements
if sum(diag(L0_test))>eps(10*N)
    warning('diagonal elements found in interaction pic L0')
    L0_test
end
L0 = L0 - diag(diag(L0));
%note this is not the actual interaction picture liouvillian for
%computational reasons with the HEOM

temp = TT*(exp( 1i*(phi -om*t))*ones(size(L0)))*TT^(-1);
temp2 = TT*(exp(-1i*(phi -om*t))*ones(size(L0)))*TT^(-1);
LE_tmp1 = sym(zeros(size(temp))); LE_tmp2 = LE_tmp1;
%remove rapidly oscillating terms, i.e. oscillating fast than RWA_cutoff
for lp1 = 1:N+1  %must be a less fucking stupid way to do this....
    for lp2 = 1:N+1
        tmp = express_exp_series(temp(lp1,lp2));
        tmpp1 = sym(0);
        for lp3 = 1:size(tmp,2)
            tmpp = subs(diff(tmp{2,lp3},t),t,0); %picks linear comp
            if abs(imag(tmpp)) < RWA_cutoff
                tmpp1 = tmp{1,lp3}*exp(tmp{2,lp3});
            end
        end
        LE_tmp1(lp1,lp2) = tmpp1;        
        tmp2 = express_exp_series(temp2(lp1,lp2));
        tmpp1 = sym(0);
        for lp3 = 1:size(tmp,2)
            tmpp = subs(diff(tmp2{2,lp3},t),t,0); %picks linear comp
            if abs(imag(tmpp)) < RWA_cutoff
                tmpp1 = tmp2{1,lp3}*exp(tmp2{2,lp3});
            end
        end
        LE_tmp2(lp1,lp2) = tmpp1;
    end
end
LE_fwd = E_env*LE_fwd.*LE_tmp1;%this actually IS in the interaction picture as 
%LE is symbolic and time dependent 
LE_bak = E_env*LE_bak.*LE_tmp2;
%% Construct labels for auxilliary density matrix

if isempty(Kap1)
    Kap1 = inf; %no truncation based on Kap1
end
    
if Kap2>0  %else trivial

numwithn = zeros(1,Kap2+1); %size(cc2,2) is number of poles in J(omega)
tot_poles = sum(cellfun(@length,cc_com));
for kk = 0:Kap2
numwithn(kk+1) = nchoosek(tot_poles-1+kk,kk);
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+size(cc2,2))) 
end

if length(Kap1)>1 %this is an overloading to allow you to just pass nn 
    %straight to the function if it has already been calculated
   nn = Kap1; totfreq = horzcat(vv{:}).';
else
    
nn = zeros(sum(numwithn),tot_poles); %if Kap1 <inf then less elements may be present
totfreq = horzcat(vv{:}).';
count1=0; tmp3 = 0;
for k =2:Kap2+1

        count2 = count1 + size(tmp3,1); 
        
        tmp3 = zeros(length(count1+1:count2)*size(nn,2),size(nn,2));
       for j = count1+1:count2
        tmp = repmat(nn(j,:),size(nn,2),1); %select element from tier below
        tmp2 = tmp+eye(size(nn,2)); %add one to every position
        tmp3(1+(j-(count1+1))*size(nn,2):(j-count1)*size(nn,2) ,:) = tmp2;       
       end
       %remove any elements of tmp3 which fail to satisfy the high
       %frequency cut off condition
       lg = true(size(tmp3,1),1);
       if isfinite(Kap1)
            for j =1:length(count1+1:count2)*size(nn,2)
                lg(j) = tmp3(j,:)*real(totfreq) < Kap1;
            end
       end
       tmp3 = tmp3(lg,:);
       %now remove duplicate entries
       jj=1;
       while jj < size(tmp3,1)-1
           lg = any(tmp3(jj+1:end,:)  ~= repmat(tmp3(jj,:),size(tmp3,1)-jj,1),2);
           tmp3 = tmp3([true(jj,1);lg],:);
           jj=jj+1;
       end
       count1 = count2;
       nn(count1+1:count1 + size(tmp3,1),:) = tmp3;
end
end

%would be good to include a complicated truncation condition such that
%weakly coupled modes could be taken only to lower orders

%clean any nn values which are all zero besides the first
nn = nn([true;any(nn(2:end,:)~=0,2)],:);
tierofnn = sum(nn,2);
for k  =2:Kap2+1    
    numwithn(k) = sum(tierofnn == (k-1)); %set actual number included
end
else
   nn=0; numwithn=1; %trivial case, DON'T ENTER Kap1 as nn and also Kap2 =0;
end

toc
tic
%% Set initial condition as lioville space
rho_vec = zeros(size(nn,1)*numel(rho_0),1);
%initial condition
rho_vec(1:numel(rho_0)) = reshape(rho_0,numel(rho_0),1);

% Can reshape to zeros(size(rho_heir,1)*numel(rho),1); 
%Column vector of ever possible element of every order of density matrix.
% Ordered as you would get if you called the density matrix by one number
% e.g. if rho is an N X N matrix, rho(n) is the [n-1,mod(N)]+1th row and 
%the floor((n-1)/N)+1 th column.

%% Calculate the coupling from the heirarchy

%general up and down coupling terms
coup_com = zeros(size(nn,1));  %coefficient will be the sqrt of number it couples to (upward)
%coup_p1_save= zeros(size(nn,1),size(nn,1),N);  %takes up too much memory

coup_com_save{N}= sparse(coup_com);

%must also have the term that accounts for the fact that down coupling
%terms are different except for the matsubara terms
%coup_m1 = coup_com ;
coup_m1_acom = coup_com ;  %anti commutator part
coup_acom_save= coup_com_save;

rng_j =0; %full_rng = 1:totpoles; 

for j = 1:N
    if ~isempty(cc_com{j}) %no bath at this mode for some reason
    rng_j = rng_j(end)+1:length(cc_com{j})+rng_j(end); 
    rng_rest = [1:rng_j(1)-1 , rng_j(end)+1:totpoles];
    cc = cc_com{j};  ccI = cc_acom{j};
for k =1:sum(numwithn(1:end-1))
    
    currenttier = tierofnn(k);
    tierabove = currenttier+1;

    tierpickr = abs(tierofnn-tierabove)<eps(10);

    nnsec = nn(k,rng_j); %section of interest
    nnsec2 = nn(k,rng_rest); %rest
    %which obviously must still match, but must be element of j section
    %changed

    temp0 = repmat(nnsec,numwithn(tierabove+1),1); 
    %make same size as nn in tier above in dimension 1
    temp02 = repmat(nnsec2,numwithn(tierabove+1),1); 
    
    temp = temp0 - nn(tierpickr,rng_j); 
    %take away elements in the tier above

    temp4 = temp02 - nn(tierpickr,rng_rest);
    
    temp2 = sum(temp~=0,2) ==1; %only one element +/- 1 diff
    
    temp3 = sum(temp,2) < -1/2 & sum(temp,2) > -3/2; %ones with a -1 net 
    temp4 = all(temp4==0,2); %match rest of the terms
    
    comb_lg = temp2 & temp3 & temp4; %full condition required
    tierpickr(tierpickr) = comb_lg ; 
    %use to put elements in the right row position in coupling matrix,
    %comb_lg is displaced by the number of tiers below
    
    if any(comb_lg)
        
    [~,temp44] = find(temp(comb_lg,:)); %find row position, i.e. which coefficient
    temp7 = sqrt(sqrt(abs(cc(temp44)).^2+abs(ccI(temp44)).^2));
    coup_com(k,tierpickr)=temp7.*sqrt(1+nnsec(temp44));
    
    %elements with minus are in the positions which are the transpose of this 
    %due to the rescaling the coefficients should be the same up and down
    %for the matsubara terms and the commutator part of the drude terms
    %now second down coupling operator
    
    coup_com(tierpickr,k) = cc(temp44).*sqrt(1+nnsec(temp44))./temp7;%...
       % ./ sqrt(sqrt(abs(cc(j,temp44)).^2+abs(ccI(j,temp44)).^2));
     %  if any(cc_acom{j}(temp44))
    coup_m1_acom(tierpickr,k) = ccI(temp44).*sqrt(1+nnsec(temp44))./temp7;%...
      %  ./ sqrt(sqrt(abs(cc(j,temp44)).^2+abs(ccI(j,temp44)).^2));    
     %  end
    end

end
    %Note I am scaling things w.r.t. (abs(a)^2+abs(b)^2)^(1/4)
    %for a term with operator a V^O + b V^X
    
    end

coup_com_save{j} = sparse(coup_com); coup_com = 0*coup_com;
coup_acom_save{j}= sparse(coup_m1_acom); coup_m1_acom=coup_com;
  
end
%down coupling from real frequency and residue modes is just  coup_p1.'


% calculate the sum_{j=1}^N sum_{k=0}^Kappa n_{jk} v_{ik}
if size(nn,1)==1
  const_factor = 0; 
else
    %reshape(vv.',numel(vv),1)
    %([vv(1,:).';vv(2,:).';vv(3,:).';vv(4,:).'])
    
const_factor = nn*totfreq;%([vv(1,:).';vv(2,:).']);
%const_factor = const_factor -imag(sum(const_factor))/length(const_factor);
end

%%  Calculate the Uber operator that propogates the entire thing
t_indep_prop_op = sparse(length(rho_vec),length(rho_vec));
%t_dep_prop_op = -1i*(LE_fwd+LE_bak); %not full size of HEOM yet
   %plus_term = total_prop_op;  minus_term = total_prop_op; 

        for j = 1:N
    
    Qjsec = zeros(N);  Qjsec(j,j) = 1;
    Qjsec = kron(Qjsec,eye(viblvls));
    eyeQ = eye(length(Qjsec)); 
    Qj = -1i*sparse(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    %keep factor of -1i from -1i tilde(V)^X term
    Uj=  sparse(kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ)); %anti commutator
    %factor of -1i cancels with +1i in front of imag part
    % plus_term = plus_term + kron(sparse(coup_p1_save{j}),Qj);
    %minus_term = minus_term + kron(sparse(coup_m1_acom_save{j}),Uj) ...
    %              +  kron(sparse(coup_m1_save{j}),Qj);
    
    t_indep_prop_op = t_indep_prop_op + kron(sparse(coup_com_save{j}),Qj);
    t_indep_prop_op = t_indep_prop_op + kron(sparse(coup_acom_save{j}),Uj);
        end

    t_indep_prop_op = t_indep_prop_op + kron(eye(length(rho_vec)/numel(Htot)),sparse(L0-Q));
    temp = kron(const_factor,ones(numel(Htot),1));
    %identical constant factor 
    t_indep_prop_op = t_indep_prop_op - sparse(1:length(temp),1:length(temp),temp);
    
%% Propogate in time and solve system

%loop over euler angles until integral over the angles converges,
easet = [0,0,0]; %first angle to take



t_dep_prop_op = subs(-1i*(LE_fwd+LE_bak),{a,b,c},{easet(1),easet(2),easet(3)});

topass{1} = t_indep_prop_op; %total prop for whole Heirarchy
topass{2} = TT;   %projector to new basis
topass{3} = sum(numwithn);  %number of tiers
topass{4} = t_dep_prop_op; 

rhs_dif_eq(1,topass); %pass parameters to function

toc

clearvars -except tend basis_proj convfact rho_vec nn use_reduced_mat ...
            rho_0 numpoints saveuptotier save_file_name numwithn

%no point keeping two copies really, save memory etc

tic
options = odeset('OutputFcn',@outputfun);
%call with number of points, 
%in general one can call the output function with 
% outputfun([numpoints,timetosave1,...timetosaveN],rho_vec,'file.mat')
% to get it to save N different sections numpoints long to a file, useful
% for long intergrations that may have to be killed.
%options.refine = 1; %less output
if length(numpoints)>1
   save(save_file_name,'convfact') %needed to convert time to ps   
end


    
       outputfun(numpoints,rho_vec(1:(sum(numwithn(1:(saveuptotier+1)))...
           *numel(rho_0))),save_file_name);

    %[Time_units,rho_vec] = ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    %[~,~] = ode23(@rhs_dif_eq,[0,tend],rho_vec,options);
    [time_units,rho_vec] =outputfun(numpoints,rho_vec,'get_data');

    non_redundant_points = [true;time_units(2:end)~=0]; 
    time_units = time_units(non_redundant_points)/convfact;
    rho_vec = rho_vec(non_redundant_points,:);
toc
%end
nntot = nn;
save(save_file_name,'nntot','-append')
clear nntot
nn = nn(1:sum(numwithn(1:(saveuptotier+1))),:);
if length(numpoints)>1
   save(save_file_name,'nn','-append')
   
end
%trim off nn that isn't part of the saved component


end
%%
function drho = rhs_dif_eq(t,rho_vc) %Couples to same numwithn value

persistent t_indep_prop light_int_op Tproj_fn Tinv_fn padmat

        if isempty(t) %pass empty t to get persistent vars back
            drho = t_indep_prop;
                        clear t_indep_prop 
            return
        end
        if iscell(rho_vc)
                        
            t_indep_prop = rho_vc{1}; %propogation matrix
            TT = rho_vc{2}; %Proj
            Tproj_fn = matlabfunction(TT);
            Tinv_fn = matlabfunction(TT^-1);
            hier_len = rho_vc{3};
            padmat = sparse(1:hier_len,1:hier_len,1);
            
            light_int_op = matlabfunction(rho_vc{4}); 
%already in interaction basis, make sure only dependence is time not angle
            
            drho =[];

            return           
        end
        
        Tproj = kron(padmat,Tproj_fn(t));
        Tinv = kron(padmat,Tinv_fn(t));
        
        L_light = kron(padmat,light_int_op(t));
        
        drho = full(L_light*rho_vc);
        rho_vc = Tinv*rho_vc;
        
        drho = drho + full(Tproj*(t_indep_prop*rho_vc));
        
end
%to get rho_vec from rho_vec red 
%             temp = reshape(tril(true(N)) , 1,N^2);
%            temp2 = repmat(temp,1,size(nn,1));
%       rho_vec(temp2) = rho_vec_red;  
% extra density matrix as normal and then just add the conjugate of matrix
% minus the leading diagonal

function [status,othervar] = outputfun(t,rho_v,flag)
persistent filetosave whentosave cnt numpoints saved_timee saved_rho Tpoints lastpoint %wtfmatlab
        %This function processes the output of the ODE

status= 0; 
if ~isempty(t) && strcmp(flag,'')
    %size(saved_rho,2)
   % size(rho_v)
    tt = t(end); 
    rr = rho_v(1:size(saved_rho,2),end);
end
        if strcmp(flag,'') %empty flag = standard call during time-ev                
  
            
            oldlastpoint = lastpoint;
            while tt>=Tpoints(lastpoint+1) %if time is past the spacing interval save it

                lastpoint = lastpoint+1;
                %if so many intervals are picked and the solver takes big
                %steps the saved things will have lots of zeros, I can't be
                %arsed to feed back into the options with refine for shit I
                %will never need.
            end
            
            
            if  oldlastpoint~=lastpoint

                saved_timee(lastpoint + 1) = tt;
                saved_rho(lastpoint + 1,:) = rr;     
            end
            
        elseif strcmp(flag,'init')
         % wtfmatlab = 1:10
                
            if isempty(whentosave)
            Tpoints = [linspace(t(1),t(2),numpoints),inf];
            else
                %function will save to disk as it goes, dumping all
                %previous data
             Tpoints = [linspace(t(1),whentosave(1),numpoints),inf];
             cnt = 1;
  
            end
            saved_timee(1) = t(1);
            saved_rho(1,:) = rho_v(1:size(saved_rho,2));

            lastpoint=1;
            return
            %when time is larger than this save the value to the var
        elseif strcmp(flag,'get_data')
            
            status = saved_timee;
            othervar = saved_rho;
               %clearvars -except status othervar
               return
        elseif ~strcmp(flag,'done') && ~isempty(flag)
         %   most initial step of all, pass savefilename as flag

            if length(t) == 1 %past a single number t equal to the 
                %number of points you wish to save
            numpoints = t;
            saved_timee = zeros(numpoints,1);
            %size(saved_timee)
            saved_rho = zeros(numpoints,length(rho_v));
            whentosave = [];
  
            else %past a vector t indicating how many points after which  
                %to save the output to disk under a filename given by flag
                numpoints = t(1);
                saved_timee = zeros(numpoints,1);
                saved_rho = zeros(numpoints,length(rho_v));
                filetosave = flag;
                whentosave = t(2:end);
                
            end
            return

        elseif strcmp(flag,'done') %finish integration
                    clearvars -except saved_rho saved_time
                      status = 1;
        return          
        end
        
        
    if ~isempty(whentosave)
%out  = flag
        if tt >= whentosave(1)
     
            if cnt==1
           
            non_redundant_points = [true;saved_timee(2:lastpoint)~=0];  %#ok<NASGU>
            else
            non_redundant_points = saved_timee~=0;     %#ok<NASGU>
            end
eval([strcat('saved_time',num2str(cnt)) '= saved_timee(non_redundant_points);']);
eval([strcat('saved_rho',num2str(cnt))  '= saved_rho(non_redundant_points,:);']);
        
            saved_timee = 0*saved_timee;
            saved_rho   = 0*saved_rho;

            save(filetosave,strcat('saved_time',num2str(cnt)),...
                strcat('saved_rho',num2str(cnt)),'-append') ;
            %toc
            if length(whentosave)>1  %shouldn't not be the case                        
            Tpoints =  [linspace(whentosave(1),whentosave(2),length(saved_timee)),inf]; 
            whentosave = whentosave(2:end); %get rid of first value
            else
            Tpoints =  [whentosave(1),inf];    %no more gets saved
                status = 1; %stop the integration
            end
            lastpoint = 1;  cnt = cnt+1;
        end
    end
end