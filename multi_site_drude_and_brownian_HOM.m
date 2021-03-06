function [Time_units,rho_vec,nn,total_prop_op]=multi_site_drude_and_brownian_HOM...
            (Htot,QQ,rho_0,cc1,cc2R,cc2I,vv1,vv2,Kap1,Kap2,viblvls,...
            numpoints,tendps,Temp,saveonlyrho00,use_reduced_mat,save_file_name)
        %NOT WORKING PROPERLY YET
%This code is used to model a dimer system of N chromophores with the
%Heirachial orders method.
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
 N = size(cc1,1); 
 Kappa = size(cc1,2);
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
% Temp e.g= 300; %temperature in Kelvin
% rho_0 is some initial condition
% saveonlyrho00 is a flag which makes the solver save only rho00
if ~exist('saveonlyrho00','var')
    saveonlyrho00 =true;
end
if ~exist('save_file_name','var')
   save_file_name  = 'saved_data.mat'; 
end
 %use_reduced_mat =1; %include only linearly independent elements of rho
 %i.e. for 2 X 2 rho: rho_vec_red rho_11, rho_21, p_22
 %tested to work
 
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); %units J K^{-1}
% to run the code not as a function just uncomment and use each module
% Energy in Joules,  extra factor of 10 is because thats speed in cm/s
convfact = 2*pi * light_speed *length_unit  *10^(-12); %2 pi * "c" in cm/s * 1ps 
tend =convfact*tendps; %end time in units you get by putting everything in wavenumbers

beta = (2 * pi * hbar * light_speed * length_unit)/ ( Temp * boltz_const);

tic

% for j = 1:N
% lg = abs(imag(vv2(j,:)))>eps(abs(vv2(j,:))); %find elements with imaginary
% 
% temp = vv2(j,lg);
% 
% for k = 1:size(temp,2) %find poles with opposite frequency
%     
%     [temp2,temp3] = min(abs(temp(k) + temp));
%     if abs(temp2) < eps(temp2)
%         pair(j,k) = temp3;
%     end       
% end
% end        
%% Construct Ltot, the total Louville for self coupling density matricies

% For N X N A and B we have the following
% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho

 LL = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));
% QQ as input should be the constants in the correction 
if size(QQ) == size(LL) %assume matrix given that is Q
    Q = QQ; clear QQ
else
if numel(QQ) == 1 % same for every site
    QQ = repmat(QQ,1,N*viblvls);
end
 Q = LL*0;
% 
% %if include_truncation_correction
for j = 1:N
    
    Qj = zeros(N);
    Qj(j,j) = 1;  Qj = kron(Qj,eye(viblvls)); eyeQ = eye(length(Qj));
      
    Q = Q + QQ(j).*(kron(Qj.',eyeQ) + kron(eyeQ,Qj) - 2*kron(Qj.',Qj)); 
   
    %reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
    %[Q_i,[Q_i,rho_n]] = (1 X Q_i + Q_i^T X 1 -2 Q_i^T X Q_i ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
    
    % .*(kron(Vj.',eyeQ) + kron(eyeQ,Vj)); 
end
% %end
end
Ltot = (LL-Q);
%in principle Q can also include high frequency decaying components from
% poles in J(omega) 

%% Construct density matrix and higher orders 


if isempty(Kap1)
    Kap1 = inf; %no truncation based on Kap1
end
    
if Kap2>0  %else trivial

numwithn = zeros(1,Kap2+1); %size(cc2,2) is number of poles in J(omega)
for kk = 0:Kap2
numwithn(kk+1) = nchoosek((Kappa+size(cc2R,2))*N-1+kk,kk); 
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+size(cc2,2))) 
end
cc = [cc1, cc2R]; ccI = [cc1*0, cc2I];
vv = [vv1, vv2];
nn = zeros(sum(numwithn),N*size(cc,2)); %if Kap1 <inf then less elements may be present
totfreq = real(reshape(vv.',numel(vv),1));
count1=0; tmp3 = 0;
for k =2:Kap2+1

        count2 = count1 + size(tmp3,1); 
        
        tmp3 = zeros(length(count1+1:count2)*size(nn,2),size(nn,2));
       for j = count1+1:count2
        tmp = repmat(nn(j,:),size(nn,2),1); %select element from tier below
        tmp2 = tmp+eye(size(nn,2)); %add one to every position
        tmp3(1+(j-(count1+1))*size(nn,2):(j-count1)*size(nn,2) ,:) = tmp2  ;     
       end
       %remove any elements of tmp3 which fail to satisfy the high
       %frequency cut off condition
       lg = false(size(tmp3,1),1);
       for j =1:length(count1+1:count2)*size(nn,2)
           lg(j) = tmp3(j,:)*totfreq < Kap1;
       end
       tmp3 = tmp3(lg,:);
       %now remove duplicate entries
       jj=1;
       while jj < size(tmp3,1)-1
           lg = any(tmp3(jj+1:end,:)  ~= repmat(tmp3(jj,:),size(tmp3,1)-jj,1),2);
           lg = [true(jj,1);lg]; %#ok<AGROW>
           tmp3 = tmp3(lg,:);
           jj=jj+1;
       end
       count1 = count2;
       nn(count1+1:count1 + size(tmp3,1),:) = tmp3;
        numwithn(k) = size(tmp3,1); %set actual number included
end
%would be good to include a complicated truncation condition such that
%weakly coupled modes could be taken only to lower orders

%clean any nn values which are all zero besides the first
nn = nn([true;any(nn(2:end,:)~=0,2)],:);

else
   nn=0; numwithn=1; %trivial case
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
    
    nnsec = nn(k,1+(j-1)*(Kappa+size(cc2R,2)):j*(Kappa+size(cc2R,2))); %section
    nnsec2 = [nn(k,1:(j-1)*(Kappa+size(cc2R,2))),nn(k,j*(Kappa+size(cc2R,2))+1:end)]; %rest
    
    temp0 = repmat(nnsec,size(nn,1),1); %make same size as nn
    temp02 = repmat(nnsec2,size(nn,1),1); 
    
    temp = temp0 - nn(:,1+(j-1)*(Kappa+size(cc2R,2)):j*(Kappa+size(cc2R,2))); %take away rest
    temp4 = temp02 - [nn(:,1:(j-1)*(Kappa+size(cc2R,2))),nn(:,j*(Kappa+size(cc2R,2))+1:end)];
    
    temp2 = sum(temp~=0,2) ==1; %only one element +/- 1 diff
    
    temp3 = sum(temp,2) < -1/2 & sum(temp,2) > -3/2; %ones with a -1 net 
    temp4 = all(temp4==0,2); %match rest of the terms
    
    comb_lg = temp2 & temp3 & temp4; %full condition required
    if any(comb_lg)
    [~,temp44] = find(temp(comb_lg,:)); %find row position
    
    temp7 = sqrt(sqrt(abs(cc(j,temp44)).^2+abs(ccI(j,temp44)).^2));
    coup_p1(k,comb_lg)=temp7.*sqrt(1+nnsec(temp44));
    end
    %elements with minus are the transpose of this 
    %now second down coupling operator

    temp33 = sum(temp,2) > 1/2 & sum(temp,2) < 3/2; %ones with a + net about
    comb_lg = temp2 & temp33 & temp4;
    [~,temp44] = find(temp(comb_lg,:)); %find row position
    if any(comb_lg)
        
    temp7 = cc(j,temp44) ./ sqrt(sqrt(abs(cc(j,temp44)).^2+abs(ccI(j,temp44)).^2));
    temp8 = ccI(j,temp44) ./ sqrt(sqrt(abs(cc(j,temp44)).^2+abs(ccI(j,temp44)).^2));
    
    coup_m1(k,comb_lg)= temp7.*sqrt(nnsec(temp44));
    coup_m1_acom(k,comb_lg) = temp8.*sqrt(nnsec(temp44));
    end
    
    
end
coup_p1_save{j} = coup_p1; coup_p1 = 0*coup_p1;
coup_m1_save{j}= coup_m1;  coup_m1 = coup_p1;
coup_m1_acom_save{j}= coup_m1_acom; coup_m1_acom=coup_p1;
  
end
%down coupling from real frequency and residue modes is just  coup_p1.'


% calculate the sum_{j=1}^N sum_{k=0}^Kappa n_{jk} v_{ik}
if size(nn,1)==1
  const_factor = 0; 
else
    %reshape(vv.',numel(vv),1)
    %([vv(1,:).';vv(2,:).';vv(3,:).';vv(4,:).'])
    
const_factor = nn*reshape(vv.',numel(vv),1);%([vv(1,:).';vv(2,:).']);
%const_factor = const_factor -imag(sum(const_factor))/length(const_factor);
end

%%  Calculate the Uber operator that propogates the entire thing
total_prop_op = sparse(length(rho_vec),length(rho_vec));
   %plus_term = total_prop_op;  minus_term = total_prop_op; 

        for j = 1:N
    
    Qjsec = zeros(N);  Qjsec(j,j) = 1;  Qjsec = kron(Qjsec,eye(viblvls)); 
    eyeQ = eye(length(Qjsec)); 
    Qj = -1i*sparse(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    %keep factor of -1i from -1i tilde(V)^X term
    Uj=  -1i*sparse(kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ)); %anti commutator
    % factor of -1i cancels with the factor here anyway
    % plus_term = plus_term + kron(sparse(coup_p1_save{j}),Qj);
    %minus_term = minus_term + kron(sparse(coup_m1_acom_save{j}),Uj) ...
    %              +  kron(sparse(coup_m1_save{j}),Qj);
    
    total_prop_op = total_prop_op + kron(sparse(coup_p1_save{j}+coup_m1_save{j}),Qj);
    total_prop_op = total_prop_op + kron(sparse(coup_m1_acom_save{j}),Uj);
        end

    total_prop_op = total_prop_op + kron(eye(length(rho_vec)/numel(Htot)),sparse(Ltot));
    temp = kron(const_factor,ones(numel(Htot),1));
    %identical constant factor 
    total_prop_op = total_prop_op - sparse(1:length(temp),1:length(temp),temp);
%% Propogate in time and solve system
topass{1} = total_prop_op; 
if  use_reduced_mat
    topass{2} = length(Htot);   topass{3} = size(nn,1);   
rhs_dif_eq2(1,topass); %pass parameters to function    
else
rhs_dif_eq(1,topass); %pass parameters to function
end
clearvars -except tend basis_proj convfact rho_vec nn use_reduced_mat ...
            rho_0 numpoints saveonlyrho00 save_file_name

%no point keeping two copies really, save memory etc

toc
tic
options = odeset('OutputFcn',@outputfun);
%call with number of points, 
%in general one can call the output function with 
% outputfun([numpoints,timetosave1,...timetosaveN],rho_vec,'file.mat')
% to get it to save N different sections numpoints long to a file, useful
% for long intergrations that may have to be killed.
options.refine = 1; %less output
if length(numpoints)>1
   save(save_file_name,'convfact') %needed to convert time to ps 
end
if  use_reduced_mat
    
       temp = reshape(tril(true(size(rho_0))) , 1 , numel(rho_0));
       temp2 = repmat(temp,1,size(nn,1));
 
       rho_vec = rho_vec(temp2);  %reduced rho_vec
       
       if saveonlyrho00
    outputfun(numpoints,rho_vec(1:sum(temp)),save_file_name);  %save only denmat
       else
    outputfun(numpoints,rho_vec,save_file_name);
       end
    %[Time_units,rho_vec] = ode45(@rhs_dif_eq2,[0,tend],rho_vec,options);
    
    [~,~] = ode45(@rhs_dif_eq2,[0,tend],rho_vec,options);
    [Time_units,rho_vec] =outputfun(numpoints,rho_vec,'get_data'); %also cleans
else
       if saveonlyrho00
    outputfun(numpoints,rho_vec(1:numel(rho_0)),save_file_name);  %save only denmat
       else
    outputfun(numpoints,rho_vec,save_file_name);
       end

    %[Time_units,rho_vec] = ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    [~,~] = ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    [Time_units,rho_vec] =outputfun(numpoints,rho_vec,'get_data');
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
        %outtest = [t,max(max(abs(rho_vc)))]
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
         %   most initial step of all

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