function out_fun = OD_wrapper2(t_range,de_fn,init_state,filter_fun,...
                              numpoints,len_to_save,solver,tol) %optional arguments

%this is a wrapper for both the ODsolver and output function
%trange is time range to solve for
%init_state is initial value
%len_to_save is used when some of the output of the ODE solver is truncated
%for whatever reason (typically HEOM stuff)
%solver is either ode45 or ode23 (default 45)
%de_fn is function handle of the rhs of the ode
%tol is a 2 component vector with the absolute and relative tolerances

if ~exist('numpoints','var')
numpoints = length(t_range)*5;
elseif isempty(numpoints) || numpoints == 0
numpoints = length(t_range)*5;    
end

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
%this may have an excess of values which are at zero due to the way the
%output function deals with being asked for an excess of time points, 
%trim all but the first

if size(fun_out,1)~=length(t_out)
    fun_out = fun_out.';
end

lg = t_out ~=0; lg(1) = true;
fun_out =fun_out(lg,:); t_out = t_out(lg); 

t_long = linspace(t_range(1),t_range(end),2*floor(length(t_out)/2)+1); %longer range
dt = t_long(2) - t_long(1); tL = length(t_long);

fun_out2 = interp1(t_out,fun_out,t_long,'pchip'); %interpolate to constant freq
LL=size(fun_out,2); 
if LL == 1
    fun_out2 = fun_out2.'; %annoying quirk of interp1, keep as column
end
Nfft = 2^nextpow2(2*tL-1); %next fast fourier transform
 %pad with zeros then take fft and rearrange w/ 0 freq centre
fun_out2 = fftshift(fft([fun_out2;zeros(Nfft-tL,LL)],[],1),1);
om_rng = pi/dt*linspace(-1,1,Nfft); %effective om_rng for longer time rng
%enumerate the filter_function over this time_range
filt_f = filter_fun(om_rng); filt_f = reshape(filt_f,length(om_rng),1);
 %multiply this by every component of the solution
fun_out2 = fun_out2.*repmat(filt_f,[1,LL]);
%inverse fft to get back in time space
phase_corr = repmat(exp(1i*pi*(linspace(0,Nfft-1,Nfft).')),[1,LL]);
fun_out2 = phase_corr.*ifft(fun_out2,[],1);  
fun_out2 = fun_out2(1:tL,:); %trim extra crap from using FFT, also post ringing

%interpolate to the ACTUAL range desired
out_fun = interp1(t_long,fun_out2,t_range,'pchip'); 

if size(fun_out,2) == 1
    out_fun = out_fun.'; %again correct for this shit
end

end