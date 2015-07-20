function [Time_units,rho_vec]=HEOM_exmanifold_solver...
            (Htot,Qtrunc,nn,H_prop_op,rho_0,numpoints,tend,saveuptotier)

%This code is used to calculate the evolution of the excited states
%dynamics in a multipigment complex, it requires parameters passed from
%another sub program
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

    % gsinc = length(Htot) - N; %one if ground state is included
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

% Construct Ltot, the total Louville for self coupling density matricies

eye_mat = sparse(1:length(Htot),1:length(Htot),ones(1,length(Htot)));
Ltot = sparse(-1i*(kron(eye_mat,Htot)-kron(Htot.',eye_mat))-Qtrunc);
%in principle Q can also include high frequency decaying components from
% poles in J(omega) 

 %reshape input into Liouville space
rho_vec = reshape(rho_0,numel(rho_0),1);
rho_vec = [rho_vec;zeros(length(rho_vec)*(length(nn)-1),1) ];
%  Combine into the Uber operator that propogates the entire thing

eye_mat = sparse(1:(length(rho_vec)/numel(Htot)),1:...
    (length(rho_vec)/numel(Htot)),ones(1,length(rho_vec)/numel(Htot)));
    H_prop_op = H_prop_op + kron(eye_mat,Ltot);
numwithn = sum(nn,2);
    
%% Propogate in time and solve system
topass{1} = H_prop_op; 

rhs_dif_eq(1,topass); %pass parameters to function

clearvars -except tend basis_proj convfact rho_vec nn use_reduced_mat ...
            rho_0 numpoints saveuptotier numwithn

%no point keeping two copies really, save memory etc

options = odeset('OutputFcn',@outputfun);
%call with number of points, 
%in general one can call the output function with 
% outputfun([numpoints,timetosave1,...timetosaveN],rho_vec,'file.mat')
% to get it to save N different sections numpoints long to a file, useful
% for long intergrations that may have to be killed.
%options.refine = 1; %less output
    
outputfun(numpoints,rho_vec(1:(sum(numwithn(1:(saveuptotier+1)))...
           *numel(rho_0))),'init_str');

    %[Time_units,rho_vec] = ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    ode45(@rhs_dif_eq,[0,tend],rho_vec,options);
    %[~,~] = ode23(@rhs_dif_eq,[0,tend],rho_vec,options);
    [Time_units,rho_vec] =outputfun(numpoints,rho_vec,'get_data');
    
    non_redundant_points = [true;Time_units(2:end)~=0]; 
    Time_units = Time_units(non_redundant_points);
    rho_vec = rho_vec(non_redundant_points,:);

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