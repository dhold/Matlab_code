function [status,othervar] = output_DE_fun(t,rho_v,flag)
persistent filetosave whentosave cnt numpoints save_time saved_rho Tpoints lastpoint %wtfmatlab
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

                save_time(lastpoint + 1) = tt;
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
            save_time(1) = t(1);
            
            saved_rho(1,:) = rho_v(1:size(saved_rho,2));
            
            lastpoint=1;
            return
            %when time is larger than this save the value to the var
        elseif strcmp(flag,'get_data') || strcmp(flag,'get_data+clean')
            
            status = save_time;
            othervar = saved_rho;
            if  strcmp(flag,'get_data+clean')
          % remove bad points with no data included in them or no gap
    tmp = [true; diff(status)~=0] & [true; status(2:end)~=0];
    status = status(tmp); othervar = othervar(tmp,:);
            end
            
               clearvars -except status othervar
               return
        elseif ~strcmp(flag,'done') && ~isempty(flag)
         %   most initial step of all, pass savefilename as flag
         %clear persistent 
%clear filetosave whentosave cnt numpoints save_time saved_rho Tpoints lastpoint
            if length(t) == 1 %past a single number t equal to the 
                %number of points you wish to save
            numpoints = t;
            save_time = zeros(numpoints,1);
            %size(save_time)
            saved_rho = zeros(numpoints,length(rho_v));
            whentosave = [];
  
            else %past a vector t indicating how many points after which  
                %to save the output to disk under a filename given by flag
                numpoints = t(1);
                save_time = zeros(numpoints,1);
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
           
            non_redundant_points = [true;save_time(2:lastpoint)~=0];  %#ok<NASGU>
            else
            non_redundant_points = save_time~=0;     %#ok<NASGU>
            end
eval([strcat('saved_time',num2str(cnt)) '= save_time(non_redundant_points);']);
eval([strcat('saved_rho',num2str(cnt))  '= saved_rho(non_redundant_points,:);']);
        
            save_time = 0*save_time;
            saved_rho   = 0*saved_rho;

            save(filetosave,strcat('saved_time',num2str(cnt)),...
                strcat('saved_rho',num2str(cnt)),'-append') ;
            %toc
            if length(whentosave)>1  %shouldn't not be the case                        
                     Tpoints =  [linspace(whentosave(1),whentosave(2),length(save_time)),inf]; 
            whentosave = whentosave(2:end); %get rid of first value
            else
            Tpoints =  [whentosave(1),inf];    %no more gets saved
                status = 1; %stop the integration
            end
            lastpoint = 1;  cnt = cnt+1;
        end
    end
end
