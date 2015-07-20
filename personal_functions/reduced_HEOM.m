function [out1,t_range1,out2,t_range2,savehrarch] = reduced_HEOM(input_state,t_end,...
                    Htot,QQ,cc_com,cc_acom,vv,Kap1,Kap2,numpoints,save_to_tier)

%used when input state is of the form
%input_state = rho0  = (0,a,b,c,d;A,0,0,0;B,0,0,0;C,0,0,0,0;D,0,0,0,0) etc, and that the 
% Hamiltonian does not mix the ground state into excited states
%mainly for use in linear response theory in spectroscopy
% 
% This set of equations is generally much simpler as there is no coupling 
% between anything but the state and the same element up and down the
% Heirarchy and to other states within the same heirarchy if the different
% sites interact (or the problem is totally seperable anyway)

%Check input is of correct form
if numel(input_state) == length(input_state)
   input_state = reshape(input_state,sqrt(numel(input_state)),sqrt(numel(input_state))); 
end
test = all(all(input_state(2:end,2:end)==0));
if ~test
    warning('input state is not of the correct form, may give funny results')
end
if ~exist('save_to_tier','var')
    save_to_tier = Kap2;
end
Htot = Htot(2:end,2:end); %reduce size
N = length(Htot);
totpoles = sum(cellfun(@length,cc_com)); %only give cc_com as bath on sites

%% Calculate HEOM stuff 
if ~iscell(Kap1) 

if isempty(Kap1)
    Kap1 = inf; %no truncation based on Kap1
end
    
if Kap2>0  %else trivial

numwithn = zeros(1,Kap2+1); %size(cc2,2) is number of poles in J(omega)
tot_poles = sum(cellfun(@length,cc_com));
for kk = 0:Kap2
numwithn(kk+1) = nchoosek(tot_poles-1+kk,kk);
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
%clean any nn values which are all zero besides the first
nn = nn([true;any(nn(2:end,:)~=0,2)],:);
tierofnn = sum(nn,2);
for k  =2:Kap2+1    
    numwithn(k) = sum(tierofnn == (k-1)); %set actual number included
end
else
   nn=0; numwithn=1; %trivial case, DON'T ENTER Kap1 as nn and also Kap2 =0;
end

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

end


else

 nn = Kap1{1}; coup_com_save = Kap1{2};    coup_acom_save = Kap1{3};   
 const_factor = Kap1{4};  

     savehrarch =[];    
end
%%  Calculate the Uber operator that propogates the entire thing

if length(Kap1) ~= 6  %else prop operators have already been passed
set1 = input_state(2:end,1); %first column
set2 = input_state(1,2:end); %first row


for j = 1:N 
    Q1(j) = QQ(j,1) + QQ(j,2);
    Q2(j) = QQ(j,1) - QQ(j,2); %this part is like a commutator and so
    %this picks up a negative element for the upper row which it acts on
end

rho_vec = zeros(length(set1)*size(nn,1),1);

decay_term = -sparse(1:length(rho_vec),1:length(rho_vec),...
    kron(const_factor,ones(length(Htot),1)));
decay_term1 = decay_term-sparse(1:length(rho_vec),1:length(rho_vec),...
                kron(ones(size(nn,1),1),Q1));
decay_term2 = decay_term-sparse(1:length(rho_vec),1:length(rho_vec),...
                kron(ones(size(nn,1),1),Q2));            

spseye = sparse(1:size(nn,1),1:size(nn,1),ones(size(nn,1),1));
total_prop_op1 = (kron(spseye ,sparse(-1i*Htot)));
total_prop_op2 = total_prop_op1 + decay_term2;
total_prop_op1 = total_prop_op1 + decay_term1;

        for j = 1:N
    
    Qjsec = zeros(N);  Qjsec(j,j) = 1;
    %Qjsec = kron(Qjsec,eye(viblvls));
    
    total_prop_op1 = total_prop_op1 -1i*kron(sparse(coup_com_save{j}),Qjsec)...
                       +kron(sparse(coup_acom_save{j}),Qjsec) ;
    total_prop_op2 = total_prop_op2 -1i*kron(sparse(coup_com_save{j}),Qjsec)...
                       -kron(sparse(coup_acom_save{j}),Qjsec) ;
        end

        options = odeset('OutputFcn',@outputfun);
        len_extra = sum(sum(nn,2)<= save_to_tier); 
        set1_ext = [set1;zeros(length(set1)*(len_extra-1),1)];
        outputfun(numpoints,set1_ext,'notsavingnayway'); 
        %set number of points to take and how many to save
        % I am not saving any of the auxilliary density matricies here
if nargout == 5
savehrarch = {nn,coup_com_save,coup_acom_save,const_factor,total_prop_op1,total_prop_op2};
if isempty(t_end) %just used to calculate HEOM stuff
    out1=[];t_range1=[];out2=[];t_range2=[];
   return 
end
%pass this back so this section can be skipped for further calculations
end
else
    total_prop_op1 = Kap1{5}; total_prop_op2 = Kap1{6}; 
end
 rho_vec1 = rho_vec;  rho_vec1(1:length(set1)) = set1;
 if norm(rho_vec1)~=0
  rhs_dif_eq(0,{total_prop_op1});   %pass to function     
ode45(@rhs_dif_eq,[0,t_end],rho_vec1,options);
[t_range1,out1]  =  outputfun(numpoints,set1,'get_data');

%trim values not output, happens when num points is high and ode45 takes
%big steps

lg1 = t_range1~=0; lg1(1) =true; %first is zero and this is sensible
t_range1=t_range1(lg1); out1 = out1(lg1,:);
 else
     t_range1=0; out1 = 0;     
 end
  rho_vec2 = rho_vec;  
  rho_vec2(1:length(set2)) = set2;
 if norm(rho_vec2)~=0 
rhs_dif_eq(0,{total_prop_op2});
ode45(@rhs_dif_eq,[0,t_end],rho_vec2,options);
[t_range2,out2]  =  outputfun(numpoints,set2,'get_data');

lg2 = t_range2~=0; lg2(1) =true; %first is zero and this is sensible
t_range2=t_range2(lg2); out2 = out2(lg2,:);
 else
t_range2=0; out2 = 0;     
 end
 



end

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

function [status,othervar] = outputfun(t,rho_v,flag)
persistent filetosave whentosave cnt numpoints saved_timee saved_rho Tpoints lastpoint %wtfmatlab
        %This function controls the output from the ODE solver to stop it
        %saving so much fucking data!  
status= 0; 
if ~isempty(t) && strcmp(flag,'')
    %size(saved_rho,2)
   % size(rho_v)
    tt = t(end); 
    %rr = rho_v(1:size(saved_rho,2),end);
end
        if strcmp(flag,'') %empty flag = standard call during time-ev                
  
            if tt >=Tpoints(lastpoint+1) %only loop if last point is
          for lp = 1:length(t)
              
              tt = t(lp);
              rr = rho_v(1:size(saved_rho,2),lp);
            oldlastpoint = lastpoint;
            while tt>=Tpoints(lastpoint+1) %if time is past the spacing interval save it

                lastpoint = lastpoint+1;
                %if so many intervals are picked and the solver takes big
                %steps the saved things will have lots of zeros. Trim these
            end
            
            
            if  oldlastpoint~=lastpoint

                saved_timee(lastpoint + 1) = tt;
                saved_rho(lastpoint + 1,:) = rr;     
            end
          end
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
