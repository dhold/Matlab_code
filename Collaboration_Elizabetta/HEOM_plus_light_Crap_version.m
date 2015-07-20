function [Time_units,rho_vec,nn,total_prop_op]=HEOM_plus_light_Crap_version...
            (Htot,QQ,rho_0,cc_com,cc_acom,vv,Kap1,Kap2,viblvls,...
            E_t, omega,mu, wave_vec_trunc , use_RWA ,...
            numpoints,tendps,saveuptotier,use_reduced_mat,save_file_name)

%In this version cc and vv and cell arrays
%All terms in inverse cm, conversion to sensible units achieved by
%E = hc / lamba => 6.62606957 *10^{-33} 299792458 *10 * "E in code" gives
%
%% Inputs
 % Htot is the total N site + ground state system and
 % whatever vibrations are treated quantum mechanically.
 % QQ is the extra convergence parameter made up of all the terms of the 
 % form int_0^t d tau exp(-v (t-tau)) * ( c_1 V^X(tau) + c_2 V^O(tau) )
 %  which are evaluated via exp(-v (t-tau)) ~ delta(t-tau)/v, as v->inf
 % and are thus just markovian terms.
 % cc_com = coefficients with exp(-vv*t) prefact contributing to V^x
 % cc_acom = coefficients with exp(-vv*t) prefact contributing to V^o
 % i.e. the anti commutator part
 % each of these is paired with a frequency vv
 % N is the total number of states
 N = length(Htot); 
% E_t time dependent electric field envelope, SYMBOLIC in t, e.g. exp(-(t-1)^2/10)
% omega is the frequency
% mu dipole coupling operator, in general this is a tensor quantity
% wave_vec_trunc: how many additional density matricies with factors of
% exp(1i*j*dot(k,r)) in front to include.  this sets the maximum j.  Note
% this is related to the order of time dep pert theory considered, chosing
% wave_vec_trunc = n is equiv to 2nth order TDPT
% use_RWA: logical, whether rotating wave approx used
%V_int = -(mu{1}*E_t(1) + mu{2}*E_t(2) + mu{3}*E_t(3));
V_int = -(mu*E_t);
%V_int_p = -(mu{1}*E_t(1) + mu{2}*E_t(2) + mu{3}*E_t(3))*exp(1i*omega*t);
%V_int_m = -(mu{1}*E_t(1) + mu{2}*E_t(2) + mu{3}*E_t(3))*exp(-1i*omega*t);
%interaction terms in the hamiltonian split into fwd and backward waves

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
tend =convfact*tendps; %end time in units you get by putting everything in wavenumbers
totpoles = sum(cellfun(@length,cc_com));
tic


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
    %straight to the function if it has already been calculated before
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
   nn=0; numwithn=1; %trivial case, DO NOT ENTER Kap1 as nn and also Kap2 =0;
end
nntot = nn; %this is usefuly to have
save(save_file_name,'nntot') 
clear nntot
toc
tic
%% Set initial condition as lioville space
rho_vec = zeros(size(nn,1)*numel(rho_0),1);
%initial condition
rho_vec(1:numel(rho_0)) = reshape(rho_0,numel(rho_0),1);

% Can reshape to zeros(size(rho_heir,1)*numel(rho),1); 
% Column vector of ever possible element of every order of density matrix.
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
%non_hermit = abs(imag(const_factor))>eps(max(imag(const_factor)));
%const_factor = const_factor -imag(sum(const_factor))/length(const_factor);
end

%%  Calculate the Uber operator that propogates the entire thing
      %this whole section will need lots of testing
tim = sym('t');
U_tev = expm(-1i*tim*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot)))));
%rowpos,colpos,expression needed to construct sparse operator

bigK = length(coup_com_save{1});
    cnt=0;
    rowpos = zeros(bigK^2*N,1); colpos=rowpos; 
    valu = sym(colpos); %preallocate as symbolic, might help? 
    
    cnt2=0; %this count relates to the number of 
    L2 = -1i*(kron(eye(length(Htot)),V_int)-kron(V_int.',eye(length(Htot))));
    L2 = U_tev*L2*(U_tev'); %interaction picture
    
    rowpos2 = zeros(bigK^2*N,1); colpos2=rowpos2; 
    valu2 = sym(colpos); %preallocate as symbolic, might help?     
     
    
for j = 1:N    
tic
    Qjsec = zeros(N);  Qjsec(j,j) = 1;
    Qjsec = kron(Qjsec,eye(viblvls)); eyeQ = eye(length(Qjsec)); 
    
    trunc_opj = QQ(j,1).*(kron(Qjsec.',eyeQ) + kron(eyeQ,Qjsec) - 2*kron(Qjsec.',Qjsec))...
                + QQ(j,2) .*(kron(Qjsec.',eyeQ) - kron(eyeQ,Qjsec));
    trunc_opj = U_tev*trunc_opj*(U_tev'); %interaction picture       
    
    Qj = -1i*(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    Qj = U_tev*Qj*(U_tev'); %interaction picture
    %keep factor of -1i from -1i tilde(V)^X term
    Uj =  (kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ));     
    Uj = U_tev*Uj*(U_tev'); %note U_tev is symbolic
    

    for k = 1:bigK
        for kk = 1:bigK
            if coup_com_save{j}(k,kk)~=0  || coup_acom_save{j}(k,kk)~=0 
                %non zero element is coupling
                for j1 = 1:N^2
                    for j2 = 1:N^2
                        if Uj(j1,j2) ~=0 
                            %will give symbolic for non zero entries
                            %which is interpreted as "true" by if
                cnt = cnt+1;
              rowpos(cnt) = (kk-1)*N^2+j2;
              colpos(cnt) = (k-1)*N^2+j1;
              valu(cnt) = Qj(j1,j2)*coup_com_save{j}(k,kk);
              valu(cnt) = valu(cnt) + Uj(j1,j2)*coup_acom_save{j}(k,kk);
              
                        end
                    end
                end
            end
            if k ==kk %on diagonal 
                for j1 = 1:N^2
                    for j2 = 1:N^2
                        if trunc_opj(j1,j2) ~=0 

                cnt = cnt+1;
              rowpos(cnt) = (kk-1)*N^2+j2;
              colpos(cnt) = (k-1)*N^2+j1;
              valu(cnt) = trunc_opj(j1,j2);
              
                        end
                        if j1  == j2 %decay factor appears
                cnt = cnt+1;
              rowpos(cnt) = (kk-1)*N^2+j2;
              colpos(cnt) = (k-1)*N^2+j1;
              valu(cnt) = -const_factor(k);                           
                        end
                        if j==1 %only need to do once
                            if L2(j1,j2) ~=0 
              %calculate elements of super op for light-matter-int
                                cnt2 = cnt2 + 1;
              rowpos2(cnt2) = (kk-1)*N^2+j2;
              colpos2(cnt2) = (k-1)*N^2+j1;
              valu2(cnt2) = L2(j1,j2);
                            end
                        end    
                           
                    end
                end
            end
        end
    end
    toc
end
%trim only ones that were included
rowpos = rowpos(1:cnt); colpos = colpos(1:cnt); valu = valu(1:cnt);
rowpos2 = rowpos2(1:cnt2); colpos2 = colpos2(1:cnt2); valu2 = valu2(1:cnt2);
%% collate duplicates, there may be some in sup op but not light_mat_op
bignum = bigK*N^2+1;
lg = sum([rowpos,colpos*bignum],2);
[lg2,lg3] = sort(lg); 
lg4 = diff(lg3)==0; %identical position entries
lg5 = true(size(rowpos)); %entires that will be included in final
jlast = -1; 

for j = find(lg4)
%note lg2 maps j back to original positions 
    if j-jlast == 1 %part of same block
        jlast = j; 
        valu(lg2(jblk_start)) = valu(lg2(jblk_start)) + valu(lg2(j+1));
        lg5(lg2(j+1)) = false; %this element is not needed
    else %start of new block
        if jlast >0 %not for the very first
           lg5(lg2(jlast+1)) = false;
        end
        jblk_start = j; jlast = j;
        valu(lg2(j)) = valu(lg2(j)) + valu(lg2(j+1)); %collate with next value
    end
end

rowpos = rowpos(lg5); colpos = colpos(lg5); valu = valu(lg5);

% to get the super operator sub in for time value and do
% sup_op = sparse(rowpos,colpos,double(sub(valu,time,t)))
% for the interaction with light one must take  
% L_plus = sparse(rowpos2,colpos2,double(sub(exp(-1i*omega*time)*valu,time,t)))
% L_minus = sparse(rowpos2,colpos2,double(sub(exp(+1i*omega*time)*valu,time,t)))
% the reason for plus and minus seeming counter intuitive is these raise /
% lower the index related to the factor of exp(i*l*dot(k,r))

%% Propogate in time and solve system
topass{1} = {rowpos,colpos,valu}; 
topass{2} = {rowpos2,colpos2,valu2};  topass{3} = omega; 
topass{4} = wave_vec_trunc; topass{5} = tim; %pass symbolic variable
if  use_reduced_mat
    topass{6} = length(Htot);   topass{7} = size(nn,1);  
rhs_dif_eq2(1,topass); %pass parameters to function    
else
rhs_dif_eq(1,topass); %pass parameters to function
end
toc

clearvars -except tend basis_proj convfact rho_vec nn use_reduced_mat wave_vec_trunc...
            rho_0 numpoints saveuptotier save_file_name non_hermit numwithn

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

if  use_reduced_mat
    warning('havent got this to work with brownian modes yet...')
       temp = reshape(tril(true(size(rho_0))) , 1 , numel(rho_0));
       temp2 = repmat(temp,1,size(nn,1));

       rho_vec = rho_vec(temp2);  %reduced rho_vec
       
    outputfun(numpoints,rho_vec(1:sum(numwithn(1:saveuptotier+1))*sum(temp)...
        ),save_file_name);  %save only denmat, only works with no nonhermitian

    %[Time_units,rho_vec] = ode45(@rhs_dif_eq2,[0,tend],rho_vec,options);
    rho_vec = [rho_vec; repmat(rho_vec*0,2*wave_vec_trunc,1)];
       %add extra tiers on the bottom relating to terms with spatial
       %dependence, note none of these are saved
    ode45(@rhs_dif_eq2,[0,tend],rho_vec,options);

    [Time_units,rho_vec] =outputfun(numpoints,rho_vec,'get_data'); %also cleans
else
    
       outputfun(numpoints,rho_vec(1:(sum(numwithn(1:(saveuptotier+1)))...
           *numel(rho_0))),save_file_name);

       rho_vec = [rho_vec; repmat(rho_vec*0,2*wave_vec_trunc,1)];
       %add extra tiers on the bottom relating to terms with spatial
       %dependence, note none of these are saved
       
    %[Time_units,rho_vec] = ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    %[~,~] = ode23(@rhs_dif_eq,[0,tend],rho_vec,options);
    [Time_units,rho_vec] =outputfun(numpoints,rho_vec,'get_data');
end

    non_redundant_points = [true;Time_units(2:end)~=0]; 
    Time_units = Time_units(non_redundant_points)/convfact;
    rho_vec = rho_vec(non_redundant_points,:);
toc
%end

nn = nn(1:sum(numwithn(1:(saveuptotier+1))),:);
if length(numpoints)>1
   save(save_file_name,'nn','-append')
end
%trim off nn that isn't part of the saved component

if nargout ==4
total_prop_op = rhs_dif_eq2([],[]);
end
end
%%
function drho = rhs_dif_eq(t,rho_vc) %Couples to same numwithn value

persistent tim a1 a2 a3 b1 b2 b3  omega num_wvc rng rngex coup_plus coup_minus

        if isempty(t) %pass empty t to get persistent vars back
            total_prop = sparse(a1,a2,double(subs(a3,tim,0)));
            drho = total_prop;
                        clear total_prop p_plus p_minus
            return
        end
        if iscell(rho_vc)
            tim = rho_vc{5}; %symbolic
            a1 = rho_vc{1}{1}; %propogation matrix
            a2 = rho_vc{1}{2}; %propogation matrix
            a3 = matlabFunction(rho_vc{1}{3},'vars',tim); %propogation matrix   
            
            
            b1 = rho_vc{2}{1}; %propogation matrix of light
            b2 = rho_vc{2}{2}; %propogation matrix of light
            omega = rho_vc{3}; %frequency
          
            b3 = matlabFunction(rho_vc{2}{3},'vars',tim); %forward dir
           
            
            num_wvc = rho_vc{4};
            total_prop = sparse(a1,a2,double(subs(a3,tim,0)));
            %this just speeds up what sections couple etc
            rng = 1:length(total_prop);  rngex(1) = 0;
            %arranged like 0,1,-1,2,-2 ....
            tmp = zeros(2*num_wvc+1,1);
            for j = 2:(2*num_wvc-1)
                if tmp(j-1)>0
                 tmp(j) = -tmp(j-1);
                else
                    if j==2
                        tmp(j)=1;
                    else
                 tmp(j) = tmp(j-2) + 1;
                    end
                end
            end
            coup_plus = false(2*num_wvc+1);  coup_minus = coup_plus';
            for j = 1:(2*num_wvc-1)
                for jj = 1:(2*num_wvc-1)
                    if tmp(j)==tmp(jj)+1
                coup_plus(j,jj) = true;  
                    elseif tmp(j)==tmp(jj)-1
                coup_minus(j,jj) = true;          
                    end
                end
            end
               
            
            for k =2:(2*num_wvc+1)
                rngex(k) = rngex(k-1) + length(total_prop);
            end
            
            
            drho =[];

            return           
        end
        %outtest = [t,max(max(abs(rho_vc)))]
        %construct propogation matricies at these times
        tic

        tmp =a3(t); %takes too long... dammit
        tmp2=b3(t);
        
    total_prop = sparse(a1,a2,tmp); 
    pp_plus  = sparse(b1,b2,exp(-1i*omega*t)*tmp2);
    pp_minus = sparse(b1,b2,exp(+1i*omega*t)*tmp2);
        drho = zeros(rngex(end) + rng(end),1);
    toc
    tic
        for k = 1:num_wvc
            rng1 = rng+rngex(k);
            drho(rng1) = drho(rng1) + full(total_prop*rho_vc(rng1));
             %rngex(coup_plus(k,:))
            % rng + rngex(coup_plus(k,:))
            tmp =  rngex(coup_plus(k,:));
            if ~isempty(tmp)
            drho(rng1) = drho(rng1) + full(...
                pp_plus*rho_vc(rng +tmp))  ;
            end
            tmp =  rngex(coup_minus(k,:));
            if ~isempty(tmp)
            drho(rng1) = drho(rng1) + full(...
                pp_minus*rho_vc(rng + rngex(coup_minus(k,:))) ) ;   
            end
        end
        toc
end

function drho = rhs_dif_eq_RWA(t,rho_vc) %Couples to same numwithn value
warning('this function isnt implimented yet')
persistent a1 a2 a3 b1 b2 b3 omega num_wvc rng rngex coup_plus coup_minus

        if isempty(tt) %pass empty t to get persistent vars back
            total_prop = sparse(a1,a2,double(subs(a3,tim,0)));
            drho = total_prop;
                        clear total_prop p_plus p_minus
            return
        end
        if iscell(rho_vc)
            a1 = rho_vc{1}{1}; %propogation matrix
            a2 = rho_vc{1}{2}; %propogation matrix
            a3 = rho_vc{1}{3}; %propogation matrix                        
            b1 = rho_vc{2}{1}; %propogation matrix of light
            b2 = rho_vc{2}{2}; %propogation matrix of light
            b3 = rho_vc{2}{3}; %propogation matrix of light
            omega = rho_vc{3}; %frequency
            num_wvc = rho_vc{4};
            total_prop = sparse(a1,a2,double(subs(a3,tim,0)));
            %this just speeds up what sections couple etc
            rng = 1:length(total_prop);  rngex(1) = 0;
            coup_plus = false(num_wvc+1);
            for j = (1:(num_wvc-1))
                coup_plus(j,j+1) = true;  
            end
                coup_minus = coup_plus';
            %arranged like 0,1,2 .... only fwd prop taken
            for k =2:(num_wvc+1)
                rngex(k) = rngex(k-1) + length(total_prop);
            end
            
            
            drho =[];

            return           
        end
        %outtest = [t,max(max(abs(rho_vc)))]
        %construct propogation matricies at these times
    total_prop = sparse(a1,a2,double(subs(a3,tim,0)));
    pp_plus  = sparse(b1,b2,exp(-1i*omega*tt)*double(sub(b3,time,t)));
    pp_minus = sparse(b1,b2,exp(1i*omega*tt)*double(sub(b3,time,t)));
        drho = zeros(rngex(end) + rng(end),1);

        for k = 1:num_wvc
            rng1 = rng+rngex(k);
            drho(rng1) = drho(rng1) + full(total_prop*rho_vc(rng1));
            drho(rng1) = drho(rng1) + full(...
                pp_plus*rho_vc(rng + rngex(coup_plus(k,:))) ) ;
            drho(rng1) = drho(rng1) + full(...
                pp_minus*rho_vc(rng + rngex(coup_minus(k,:))) ) ;            
        end
end


function drho = rhs_dif_eq2(t,rho_vc_red) 

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