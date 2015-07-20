function [status,othervar] = output_fun(t,rho,flag)
% This function is used as an output control for ODE solver based on
% density matrix "rho" evolution.  Rho is a COLUMN VECTOR i.e. Louiville 
% space representation.  
% Call this function with a flag 'set_vars' or something similar passing
% rho a vector of length number of density elements you intend to save and
% t either a single number, the number of points to save over the time
% range or [numpoints, whentosave(1),...,whentosave(n)]
persistent filetosave whentosave cnt numpoints saved_time saved_rho Tpoints lastpoint
        %This function controls the output from the ODE solver to stop it
        %saving so much data!  

status= 0;  

        if strcmp(flag,'') %empty flag = standard call                 
  
        for j = 1:length(t) %multiple points output in some solvers
            lgvar = false;
            while t(j)>=Tpoints(lastpoint+1) %if time is past the spacing interval save it

                lastpoint = lastpoint+1; lgvar = true; % point to save
            end
            if lgvar
                saved_time(lastpoint) = t(j);
                saved_rho(lastpoint,:) = rho(1:size(saved_rho,2),j);
            end
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
            saved_time(1) = t(1);
            saved_rho(1,:) = rho(1:size(saved_rho,2));

            lastpoint=1;
            %when time is larger than this save the value to the var
        elseif strcmp(flag,'get_data')
            
            status = saved_time;
            othervar = saved_rho;
               clearvars -except status othervar
               return
        elseif ~strcmp(flag,'done') && ~isempty(flag)


            if length(t) == 1 %past a single number t equal to the 
                %number of points you wish to save
            numpoints = t;
            saved_time = zeros(numpoints,1);
            saved_rho = zeros(numpoints,length(rho));
            whentosave = [];
  
            else %past a vector t indicating how many points after which  
                %to save the output to disk under a filename given by flag
                numpoints = t(1);
                saved_time = zeros(numpoints,1);
                saved_rho = zeros(numpoints,length(rho));
                filetosave = flag;
                whentosave = t(2:end);

            end
            

        elseif strcmp(flag,'done') %finish integration
           status = 1;
            clearvars -except saved_rho saved_time
            return
            
        end
    if ~isempty(whentosave)

        if t(end) >= whentosave(1)
            
            if cnt==1
            non_redundant_points = [true;saved_time(2:lastpoint)~=0]; 
            %if num points is quite high and the function solver takes big
            %timesteps then some points won't have a value for them
            else
            non_redundant_points = saved_time~=0;     %#ok<*NASGU>
            end

            eval([strcat('saved_time',num2str(cnt)) '= saved_time(non_redundant_points);']);
            eval([strcat('saved_rho',num2str(cnt))  '= saved_rho(non_redundant_points,:);']);
            saved_time  = 0*saved_time;
            saved_rho   = 0*saved_rho;

            save(filetosave,strcat('saved_time',num2str(cnt)),...
                strcat('saved_time',num2str(cnt)),'-append') ;
            if length(whentosave)>1  %should always be the case                        
            Tpoints =  [linspace(whentosave(1),whentosave(2),length(saved_time)+1),inf]; 
            whentosave = whentosave(2:end); %get rid of first value
            else
            Tpoints =  [whentosave(1),inf];    %no more gets saved
            warning('This shouldnt happen unless you only specified one point to save at')
                status = 1; %stop the integration
            end
            lastpoint = 1;  cnt = cnt+1;
                
        end
    end
end