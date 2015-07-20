function [Time_units,rho_vec,nn,basis_proj,total_prop_op]=...
            dimer_elec_vib_HOM(E1,V,omegavib,viblvls,...
            coupg,tendps,gamma_dru,lambda_dru,Kappa,Kap2,Temp,rho_0)
%This code is used to model a dimer system of two chromophores with the
%Heirachial orders method.
% Unit setup
%All terms in inverse cm, conversion to sensible units achieved by
%E = hc / lamba => 6.62606957 *10^{-33} 299792458 *10 * "E in code" gives
%
tic
if nargin == 0  %use a select set nargin
%%
 E1 = 1042; %Energy of initial state (compared to excitation in other chromopore)
 V = 92; %Coupling
 omegavib = 1111; %Frequency of vibrations / level spacing,
% %assumed to be the same in each protein, choose (omega_1+omega_2)/2
 viblvls = 10; %Number of vibrational levels to include, integer
 coupg = 267.1; %coupling term between electronic and vibrational DOF
 tendps = 2; %end time in pico seconds
 gamma_dru = [100;100]; %drude decay constant
 lambda_dru = [1;1]; % weighting of distrubtion
 % %J_j(omega) = 2 lamda_j gamma_j /hbar * (omega/(omega^2+lambda_j^2)
 Kappa = 1; %truncation parameter satisfying Kappa >>omega_0 beta hbar/2pi
 %beyond this e^(-v_(k>kappa) t)*v_k ~delta(t) is made
 Kap2 = 2; %second truncation parameter, there are Kappa + 1 explicitly 
 %treated exponentials, this number is the maximum order the sum of the
 %number of powers can take, e.g. n_1 + n_2 +... + n_Kappa <= Kap2
 Temp = 300; %temperature in Kelvin
 include_truncation_correction = true; %for debugging#
 use_reduced_mat =1; %include only linearly independent elements of rho
 %i.e. for 2 X 2 rho: rho_vec_red rho_11, rho_21, p_22
 %tested to work

 saveonlyrho00 = true;

end
N=2; %number of sites, set to two for dimer code.
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}
% to run the code not as a function just uncomment and use each module
% Energy in Joules,  extra factor of 10 is because thats speed in cm/s
convfact = 2*pi * light_speed *length_unit  *10^(-12); %2 pi * "c" in cm/s * 1ps 
tend =convfact*tendps; %end time in units you get by putting everything in wavenumbers

beta = (2 * pi * hbar * light_speed * length_unit)/ ( Temp * boltz_const);

%% Hamiltonian terms considered explicitly
%Excitonic Hamilonian
Hel = zeros(2); 

Hel(1,1) = E1; Hel(2,2) = 0;  %On diagonal terms
Hel(1,2) = V;  Hel(2,1) = V;  % Coupling terms from dipole dipole interaction

%Vibrational Hamilonian, relative degrees of freedom only
if viblvls > 1 %zero/1 means no quantised vibration considered
tmp = ones(viblvls,1);
tmp(1) = 0;
for k=2:viblvls
    tmp(k) = tmp(k-1) + omegavib;
end
Hvib = diag(tmp); clear tmp 
%Note this is H_vib = omega_vib * (b_r^dag b_r) not the COM excitation
notquitezero = kron(eye(size(Hel)),diag(E1*(-eps*(viblvls-1):2*eps:eps*(viblvls-1))));
%weighted in order to just keep the vibrational exited states ahead, this
%is quite a lazy solution I could just reorder afterwards
[basis_proj,~] = eig( kron(Hel,eye(viblvls))+notquitezero ); 
%projector to delocalised basis

Htot = kron(Hel,eye(viblvls))+kron(eye(2),Hvib); %total Hamiltonian
%which will be explicitly considered QM

%Hamilonian coupling electronic and vibrational levels
%Hexvib = -g/sqrt(2)  * (n_1 - n_2)*(b^dag_rel + b_rel)

Hexvib = diag(-coupg/sqrt(2) * sqrt(1:viblvls-1),1);
Hexvib = Hexvib+Hexvib.';
Hexvib = kron([-1,0;0,1],Hexvib);

Htot = Htot + Hexvib;  %Full Hamiltonian with electric + vib + coupling
%in chromophore basis
%Htot = basis_proj'*Htot*basis_proj; %project to exciton basis
else
    [basis_proj,~] = eig(Hel); 
    Htot = Hel;
   %rho Htot = basis_proj'*Hel*basis_proj;
end

%coefficients used in the HOM sum
cc = zeros(N,Kappa+1); vv = cc;
vv(:,1) = gamma_dru;  %Drude decay constant

vv(:,2:end) = (2*pi/beta)*kron(ones(N,1),1:Kappa); %These are the Matsuraba frequencies

cc(:,1) = gamma_dru.*lambda_dru.*(cot(beta.*gamma_dru/2)-1i);

cc(:,2:end) = (4/beta) * repmat(gamma_dru.*lambda_dru,1,Kappa)...
                .*   (         vv(:,2:Kappa+1)./ ...
                (vv(:,2:Kappa+1).^2 -repmat(gamma_dru.^2,1,Kappa))   );
            
%% Construct Ltot, the total Louville for self coupling density matricies

% For N X N A and B we have the following
% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho

L = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));

Q = L*0;

if include_truncation_correction
for j = 1:N
    
    Qj = zeros(N);
    Qj(j,j) = 1;  Qj = kron(Qj,eye(viblvls)); eyeQ = eye(length(Qj));
      
    Q = Q + (lambda_dru(j).*(2./(beta*gamma_dru(j))-cot(beta*gamma_dru(j)/2))...
        - dot(cc(j,2:end),1./vv(j,2:end)))...
        .*(kron(Qj.',eyeQ) + kron(eyeQ,Qj) - 2*kron(Qj.',Qj)); 
    %reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
    %[Q_i,[Q_i,rho_n]] = (1 X Q_i + Q_i^T X 1 -2 Q_i^T X Q_i ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
end
end

Ltot = (L-Q);

%% Construct density matrix and higher orders 

rho = zeros(size(Htot));  % zeroth order in Heirarchy, ACTUAL density matrix

if Kap2>0 %else trivial

numwithn = zeros(1,Kap2+1);
for kk = 0:Kap2
numwithn(kk+1) = nchoosek((Kappa+1)*N-1+kk,kk); 
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+1)) 
end

nn = zeros(sum(numwithn),N*(Kappa+1)); 
count1=0;
for k =2:Kap2+1

        count2 = count1 + numwithn(k-1); 
        
        tmp3 = zeros(length(count1+1:count2)*N*(Kappa+1),N*(Kappa+1));
       for j = count1+1:count2
        tmp = repmat(nn(j,:),N*(Kappa+1),1); %select element from tier 
        tmp2 = tmp+eye(N*(Kappa+1)); %add one to ever position
        tmp3(1+(j-(count1+1))*N*(Kappa+1):(j-count1)*N*(Kappa+1) ,:) = tmp2  ;     
       end
       %remove duplicates
       jj=1;
       while jj < size(tmp3,1)-1
           lg = any(tmp3(jj+1:end,:)  ~= repmat(tmp3(jj,:),size(tmp3,1)-jj,1),2);
           lg = [true(jj,1);lg]; %#ok<AGROW>
           tmp3 = tmp3(lg,:);
           jj=jj+1;
       end
       count1 = count2;
       nn(count1+1:count1 + numwithn(k),:) = tmp3;
    
end
else
   nn=0; numwithn=1; %trivial case
end
%% Set initial condition
if ~exist('rho_0','var') %no initial condition, have to create standard one
    if viblvls > 1

        vibpop = 1./(exp((1:viblvls).*beta*omegavib)-1); %BE stats
        vibpop = vibpop/sum(vibpop); %normalise, sum not sum squared
        rho_0 = full(sparse(1:length(rho),1:length(rho),[zeros(1,length(rho)/2),vibpop]));
        %1 excitation in exiton basis (excited state) and random vib pops
        rho_0 = basis_proj*rho_0*basis_proj'; %excited state in chromophore basis
    else
        rho_0 = zeros(N);
        rho_0(end) = 1; %excited state in exciton basis (as energy ordered)
        rho_0 = basis_proj*rho_0*basis_proj'; %excited state in chromophore basis
    end
end
rho_vec = zeros(size(nn,1)*numel(rho_0),1);
%initial condition
rho_vec(1:numel(rho_0)) = reshape(rho_0,numel(rho_0),1);
clear rhoharch
% Can reshape to zeros(size(rho_heir,1)*numel(rho),1); 
%Column vector of ever possible element of every order of density matrix.
% Ordered as you would get if you called the density matrix by one number
% e.g. if rho is an N X N matrix, rho(n) is the [n-1,mod(N)]+1th row and 
%the floor((n-1)/N)+1 th column.

%% Calculate the coupling from the heirarchy

%general up coupling
coup_p1 = zeros(size(nn,1));  %coefficient will be the sqrt of number it couples to (upward)
coup_p1_save= {};
%must also have the term that accounts for the fact that down coupling to 
% the k=0 states are different
coup_m1 = coup_p1 ; 
coup_m1_acom = coup_p1 ;  %anti commutator part
coup_m1_save= {}; coup_m1_acom_save= {};
for j = 1:N
for k =1:sum(numwithn(1:end))
    
    nnsec = nn(k,1+(j-1)*(Kappa+1):j*(Kappa+1)); %section
    nnsec2 = [nn(k,1:(j-1)*(Kappa+1)),nn(k,j*(Kappa+1)+1:end)]; %rest
    
    temp0 = repmat(nnsec,size(nn,1),1); %make same size as nn
    temp02 = repmat(nnsec2,size(nn,1),1); 
    
    temp = temp0 - nn(:,1+(j-1)*(Kappa+1):j*(Kappa+1)); %take away rest
    temp4 = temp02 - [nn(:,1:(j-1)*(Kappa+1)),nn(:,j*(Kappa+1)+1:end)];
    
    temp2 = sum(temp~=0,2) ==1; %only one element +/- 1 diff
    
    temp3 = sum(temp,2) < -1/2 & sum(temp,2) > -3/2; %ones with a -1 net 
    temp4 = all(temp4==0,2); %match rest of the terms
    
    comb_lg = temp2 & temp3 & temp4; %full condition required
    if any(comb_lg)
    [~,temp44] = find(temp(comb_lg,:)); %find row position
    
    temp7 = sqrt(abs(cc(j,temp44)));
    coup_p1(k,comb_lg)=temp7.*sqrt(1+nnsec(temp44));
    end
    %now second down coupling operator

    temp33 = sum(temp,2) > 1/2 & sum(temp,2) < 3/2; %ones with a + net about
    comb_lg = temp2 & temp33 & temp4;
    [~,temp44] = find(temp(comb_lg,:)); %find row position
    if any(comb_lg)
        
    temp7 = real(cc(j,temp44)) ./ sqrt(abs(cc(j,temp44)));
    temp8 = 1i*imag(cc(j,temp44)) ./ sqrt(abs(cc(j,temp44)));
    
    coup_m1(k,comb_lg)= temp7.*sqrt(nnsec(temp44));
    coup_m1_acom(k,comb_lg) = temp8.*sqrt(nnsec(temp44));
    end
    %elements with minus are the transpose of this 
    
end
coup_p1_save{j} = coup_p1; coup_p1 = 0*coup_p1;
coup_m1_save{j}= coup_m1;  coup_m1 = coup_p1;
coup_m1_acom_save{j}= coup_m1_acom; coup_m1_acom=coup_p1;
  
end
%standard down coupling is just  coup_p1.'


% calculate the sum_{j=1}^N sum_{k=0}^Kappa n_{jk} v_{ik}
if size(nn,1)==1
  const_factor = 0; 
else
const_factor = nn*reshape(vv.',numel(vv),1);%([vv(1,:).';vv(2,:).']);
end




%%  Calculate the Uber operator that propogates the entire thing
total_prop_op = sparse(length(rho_vec),length(rho_vec));
   plus_term = total_prop_op;  minus_term = total_prop_op; 

        for j = 1:N
    
    Qjsec = zeros(N);  Qjsec(j,j) = 1;  Qjsec = kron(Qjsec,eye(viblvls)); 
    eyeQ = eye(length(Qjsec)); 
    Qj = -1i*sparse(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    Uj=  -1i*sparse(kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ)); %anti commutator
    %Uj=  +1i*sparse(kron(Qjsec,eyeQ) + kron(eyeQ,Qjsec.'));

     plus_term = plus_term + kron(sparse(coup_p1_save{j}),sparse(Qj));
    minus_term = minus_term + kron(sparse(coup_m1_acom_save{j}),sparse(Uj)) ...
                  +  kron(sparse(coup_m1_save{j}),sparse(Qj));
    
    total_prop_op = total_prop_op + kron(sparse(coup_p1_save{j}+coup_m1_save{j}),sparse(Qj));
    total_prop_op = total_prop_op + kron(sparse(coup_m1_acom_save{j}),sparse(Uj));
        end

    total_prop_op = total_prop_op + kron(eye(length(rho_vec)/numel(Htot)),sparse(Ltot));

             temp = kron(const_factor,ones(numel(Htot),1));
    total_prop_op = total_prop_op - sparse(1:length(temp),1:length(temp),temp);

%% Propogate in time and solve system
topass{1} = total_prop_op; 
if  use_reduced_mat
    topass{2} = length(Htot);   topass{3} = size(nn,1);   
rhs_dif_eq2(1,topass); %pass parameters to function    
else
rhs_dif_eq(1,topass); %pass parameters to function
end
clearvars -except tend basis_proj convfact rho_vec nn use_reduced_mat rho_0 saveonlyrho00 

%no point keeping two copies really, save memory etc

%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
%[T,Y] = ode45(@rhs_dif_eq,tspan,initcond,options);
%split over several time intervals so as not to fill memory in general
toc

% if 1==0 %this really doesn't work for some reason
% tintervals = [linspace(0,tend*(9/10),10);linspace(tend*(1/10),tend,10)];
% rho_vec_init = rho_vec;
% points_per_interval = 40;
% Time_units = zeros(points_per_interval*10,1);
% rho_vec = zeros(size(rho_vec,1),length(Time_units));
% size(rho_vec )
% for k = 1:10
% [time,temp] = ode45(@rhs_dif_eq,tintervals(:,k),rho_vec_init);
% temp = temp.';
% %save about 40 points in time, 
% if length(time) > points_per_interval
% ptosave = 1+round((0:points_per_interval-1)*length(time)/(points_per_interval));
% 
% Time_units((k-1)*points_per_interval+1:k*points_per_interval) = time(ptosave);
% rho_vec(:,(k-1)*points_per_interval+1:k*points_per_interval) =temp(:,ptosave);
% else
% Time_units((k-1)*points_per_interval+(1:length(time))) = time;   
% rho_vec(:,(k-1)*points_per_interval+(1:length(time))) = temp;    
% end
% 
% rho_vec_init = temp(:,end);
%  
% clear temp time
% end
% 
% Time_units = Time_units/convfact; %convert to ps from cm^-1
% else
tic
options = odeset('OutputFcn',@outputfun);
%call with number of points, 
%in general one can call the output function with 
% outputfun([numpoints,timetosave1,...timetosaveN],rho_vec,'file.mat')
% to get it to save N different sections numpoints long to a file, useful
% for long intergrations that may have to be killed.
options.refine = 1; %less output
if  use_reduced_mat
                 temp = reshape(tril(true(size(rho_0))) , 1 , numel(rho_0));
                temp2 = repmat(temp,1,size(nn,1));
 
       rho_vec = rho_vec(temp2);  %reduced rho_vec
       if saveonlyrho00
    outputfun(200,rho_vec(1:sum(temp)),'num_points');
       else
    outputfun(200,rho_vec,'num_points');       
       end
    %[Time_units,rho_vec] = ode45(@rhs_dif_eq2,[0,tend],rho_vec,options);
    [~,~] = ode45(@rhs_dif_eq2,[0,tend],rho_vec,options);
    [Time_units,rho_vec] =outputfun(200,rho_vec,'get_data'); %also cleans
else
    outputfun(200,rho_vec,'num_points');
    %[Time_units,rho_vec] = ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    [~,~] = ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    [Time_units,rho_vec] =outputfun(200,rho_vec,'get_data');
end

    non_redundant_points = [true;Time_units(2:end)~=0]; 
    Time_units = Time_units(non_redundant_points)/convfact;
    rho_vec = rho_vec(non_redundant_points,:);
toc
%end
total_prop_op = rhs_dif_eq2([],[]);
end
%%
function drho = rhs_dif_eq(t,rho_vc) %Couples to same numwithn value

persistent total_prop

        if isempty(t) %pass empty t to get persistent vars back
            drho = total_prop;
                        clear total_prop 
            return
        end
        if iscell(rho_vc)
                        
            total_prop = rho_vc{1}; %propogation matrix

            drho =[];

            return           
        end
 
        drho = full(total_prop*rho_vc);
end

function drho = rhs_dif_eq2(t,rho_vc_red) %Couples to same numwithn value

%rho_vc_red is the reduced vector (length = numstates * (N^2+N)/2 ) 
%from taking only the independent vectors of each density matrix 


persistent total_prop_1  total_prop_2 pickr 

        if isempty(t) %pass empty t to get persistent vars back
            drho{1} = total_prop_1;
            drho{2} = total_prop_2;
                        clear total_prop_1 total_prop_2  pickr 
            return
        end
        if iscell(rho_vc_red)
                        
            total_prop = rho_vc_red{1}; %propogation matrix
            N = rho_vc_red{2}; %size of Htot
            NN = rho_vc_red{3}; %total number of density matricies
            temp = reshape(tril(true(N)) , 1,N^2);
            temp2 = repmat(temp,1,NN);  %picks out L.I. elements in lower d
            total_prop_1 = total_prop(temp2,temp2);
            %total_1 doesn't include anything that multiples the conjugates
%              temp3 = repmat(temp',1,NN);
%              total_prop_2 = total_prop(temp2,temp3);
%              temp4 = reshape(eye(N)~=0 , 1,N^2); %diagonal eles
%              temp5 =  repmat(temp4,1,NN);
%              total_prop_2(:,temp5(temp3)) =0;
%             
            total_prop_2 = total_prop(temp2,~temp2);
            %element that multiples the conjugate terms
            
 % must construct matrix that maps the off diagonal positions to
 % the order they need to appear in the matrix multiplication,
 % the flattening does not preserve this order... must map to the order in
 % which the non LI elements (upper half) appear to the order in which the
 % elements on the lower diagonal (on which they depend) will appear.
            
            count1 =1;  count2 = 0;
            test = zeros(N);
     for k = 1:(N-1)
     count2 = count2 + N-k;
     test(1+k:N,k) =count1:count2;
     count1= count2+1; 
     end
%     for k = 2:N
%     count2 = count2 + k-1;
%     test(1:k-1,k) =count1:count2;
%     count1= count2+1; 
%     end
    test2 = reshape(test',1,N^2);  
    [~,order] = sort(test2(test2~=0));

            pickr = reshape(tril(true(N),-1) , 1,N^2);
            pickr = pickr(temp); 
            %pic elements of rho_vec_red which are on diagonal
            pickr = repmat(pickr,1,NN);
            %prop 2 multiplies, i.e. off diag element
            %size(total_prop_1)
           % size(total_prop_2)
            drho =[];

                neworder = zeros(sum(pickr),1);
                neworder(1:count2) = order; 
    for k = 2:NN
        neworder(1+((k-1)*count2):(k*count2))  = ...
            neworder(1+((k-2)*count2):((k-1)*count2)) + max(order);
    end
    total_prop_2 = total_prop_2(:,neworder);

            return           
        end
        drho = full(total_prop_1*rho_vc_red);

        drho = drho + full(total_prop_2*conj(rho_vc_red(pickr)));
       % test = full(total_prop_3*conj(rho_vc_red(pickr)));
      %  test2 = full(total_prop_2*conj(rho_vc_red));
      %  test3 = norm(test-test2)
      %  drho = drho + full(total_prop_2*conj(rho_vc_red));
end
%to get rho_vec from rho_vec red 
%             temp = reshape(tril(true(N)) , 1,N^2);
%            temp2 = repmat(temp,1,size(nn,1));
%       rho_vec(temp2) = rho_vec_red;  
% extra density matrix as normal and then just add the conjugate of matrix
% minus the leading diagonal

function [status,othervar] = outputfun(t,rho_v,flag)

persistent filetosave whentosave cnt numpoints saved_timee saved_rho Tpoints lastpoint %wtfmatlab
        %This function controls the output from the ODE solver to stop it
        %saving so much fucking data!  

status= 0;  
if ~isempty(t)
    tt = t(end); rr = rho_v(1:size(saved_rho,2),end);
end
        if strcmp(flag,'') %empty flag = standard call                 
  
            
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
            Tpoints = [linspace(t(1),t(2),numpoints-1),inf];
            else
                %function will save to disk as it goes, dumping all
                %previous data
             Tpoints = [linspace(t(1),whentosave(1),numpoints-1),inf];
             cnt = 1;
  
            end
            saved_timee(1) = t(1);
            saved_rho(1,:) = rho_v(1:size(saved_rho,2));

            lastpoint=1;
            %when time is larger than this save the value to the var
        elseif strcmp(flag,'get_data')
            
            status = saved_timee;
            othervar = saved_rho;
               clearvars -except status othervar
               return
        elseif ~strcmp(flag,'done') && ~isempty(flag)


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
            

        elseif strcmp(flag,'done') %finish integration
           status = 1;
           size(saved_timee)
           size(saved_rho)
            clearvars -except saved_rho saved_time
            return
            
        end
    if ~isempty(whentosave)
out  = flag
        if t >= whentosave(1)
            eval([strcat('saved_time',num2str(cnt)) '= saved_time']);
            eval([strcat('saved_rho',num2str(cnt))  '= saved_rho']);
        
            save(filetosave,strcat('saved_time',num2str(cnt)),...
                strcat('saved_time',num2str(cnt)),'-append') ;
            if length(whentosave)>1  %shouldn't not be the case                        
            Tpoints =  [linspace(whentosave(1),whentosave(2),length(saved_timee)+1),inf]; 
            whentosave = whentosave(2:end); %get rid of first value
            else
            Tpoints =  [whentosave(1),inf];    %no more gets saved
                status = 1; %stop the integration
            end
            lastpoint = 1;  cnt = cnt+1;
                
        end
    end
    %wtfmatlab
%size(saved_timee)
end