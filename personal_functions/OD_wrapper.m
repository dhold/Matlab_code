function out_fun = OD_wrapper(t_range,de_fn,init_state,len_to_save,solver,tol)

%this is a wrapper for both the ODsolver and output function
%trange is time range to solve for
%init_state is initial value
%len_to_save is used when some of the output of the ODE solver is truncated
%for whatever reason (typically HEOM stuff)
%solver is either ode45 or ode23 (default 45)
%de_fn is function handle of the rhs of the ode
%tol is a 2 component vector with the absolute and relative tolerances



numpoints = length(t_range)*5;
if ~exist('len_to_save','var')
    len_to_save = length(init_state);
elseif isempty(len_to_save)
    len_to_save = length(init_state);
end

if length(t_range) == 1 && t_range(1) ==0 %no point if there is only 1 point
    out_fun = init_state(1:len_to_save); 
    if size(out_fun,2) ~= 1; out_fun = out_fun.';end %make it a column vector
    return
end

if ~exist('solver','var')
    solver = 'ode45';
end
     
if exist('tol','var')
options = odeset('OutputFcn',@output_DE_fun,'AbsTol',tol(1),'RelTol',tol(2));      
else
options = odeset('OutputFcn',@output_DE_fun);   
end

output_DE_fun(numpoints,init_state(1:len_to_save),'notsavingnayway');
 
if strcmp(solver,'ode45') %choose solver method

    ode45(de_fn,[t_range(1),t_range(end)],init_state,options);
    
elseif strcmp(solver,'ode23')
    
    ode23(de_fn,[t_range(1),t_range(end)],init_state,options);
end
        
[t_out,fun_out]  = output_DE_fun(numpoints,init_state(1:len_to_save),'get_data');
%this will have an excess of values which are at zero, trim all but the
%first

if size(fun_out,1)~=length(t_out)
    fun_out = fun_out.';
end

lg = t_out ~=0; lg(1) = true;

out_fun = interp1(t_out(lg),fun_out(lg,:),t_range,'pchip'); %interpolate
%this interpolation may cause issues if the phase rotates fast compared to
%the point spacing, linear would ignore phase oscillations
%diagnostic test, uncomment to test
%out_fun_test = interp1(t_out(lg),fun_out(lg,:),t_range,'linear');
%diff_test = abs(out_fun-out_fun_test).^2;
%sum(diff_test(:))./sum(abs(out_fun(:).^2))

if size(fun_out,2) == 1
    out_fun = out_fun.'; %annoying quirk of interp1 that it flips this
end

end