function PP_sig = PP_from_PE(S_PE,E_pump,om_pump,om_1,E_probe,om_probe,om_3)

%S_PE = S_PE(om_1,T,om_3, aux_dim) , if it isn't in this order then permute
%dimensions so it is
sz = size(S_PE); 
if sz(1) ~= length(om_1)
    error('S_PE is not equal to size of om_1 in dimension 1')
end

if isempty(E_probe) %assume frequency resolved detection

   if sum(abs(diff(om_1,2))) < eps*length(om_1)
       use_simpson = true;
           n=numel(om_1)-1; h=(om_1(end)-om_1(1))/n;
   else
       use_simpson = false;
   end
   PP_sig  = zeros([length(om_pump),sz(2:end)]); 
        
        for k=1:length(om_pump)
            %{
            if ishandle(E_pump)
            interp_fn = @(x) repmat(reshape(E_pump(x-om_pump(k)),[length(x),1])...
                         ,[1,sz(2:end)]) .*interp1(om_1,S_PE,x-om_pump(k));    
            else %assume this is a width "tau"
            interp_fn = @(x) repmat(reshape(exp((x-om_pump(k)).^2.*E_pump^2)...
                ,[length(x),1]),[1,sz(2:end)]).*interp1(om_1,S_PE,x-om_pump(k));
            end
            PP_sig(k,:) = integral(interp_fn,om_1(1),om_1(end));
            %}
            %couldn't get ^^ that shit to work
            if strcmp(class(E_pump),'function_handle')
                interp_E = E_pump(om_1-om_pump(k));    
            else %assume this is a width "tau"
                interp_E =exp(-(om_1-om_pump(k)).^2.*E_pump^2);
                norm_fct = trapz(om_1,interp_E);
                interp_E = interp_E/norm_fct;
            end
            interp_E = reshape(interp_E,length(interp_E),1);
            interp_E = repmat(interp_E,[1,sz(2:end)]).*S_PE;
            if use_simpson
                
    PP_sig(k,:) =  h/3*(interp_E(1,:,:,:)+2*sum(interp_E(3:2:end-2,:,:,:))+...
                4*sum(sum(interp_E(2:2:end,:,:,:)))+interp_E(end,:,:,:));
            else
                PP_sig(k,:) = trapz(om_1,interp_E);
            end
        end 

else %assume integrated over all frequency detection, need probe pulse shape
    warning('not actually finished yet as not required, complete before use')      
   PP_sig  = zeros([length(om_pump),sz(2),length(om_probe),sz(4:end)]); 
        
        for k=1:length(om_pump)
            for j=1:length(om_probe)
            if strcmp(class(E_pump),'function_handle')
            interp_fn = @(x) E_pump(x-om_pump(k)).*E_probe(x-om_probe(j))...
                        .*interp1(om_1,S_PE,x-om_pump(k));    
            else
            interp_fn = @(x) exp((x-om_pump(k)).^2.*E_pump^2).*interp1(om_1,S_PE,x-om_pump(k));
            end
            PP_sig(k,:) = trapz(interp_fn,om_1(1),om_1(end));
            end
        end    
    
    
    
end