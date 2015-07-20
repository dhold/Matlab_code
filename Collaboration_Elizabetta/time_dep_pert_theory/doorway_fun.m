function  [rho_neq_g,rho_neq_e,rho_neq_fg] = doorway_fun(...
          des_freq_u,H_exciton,R_red_sc,supop_ge ,V_ge,V_eg,sz2,V_fe) 


%% Calculate to second order the density matrix before the probe beam hits
%A density matrix associated to each path in Louiville space for
%averaging purposes
if nargout ==3 && nargin < 8
    error('require V_fe to be passed to calculate the ground excited coherences')
end

mid_freq = mean( des_freq_u );
om_s = -1i*mid_freq ; %~frequency of oscillation to shift 
N = length(V_ge);  %number of sites

t_step = 2*pi/(des_freq_u(end)-des_freq_u(1)); %set by range
t_end = 2*pi/(des_freq_u(2)-des_freq_u(1)); %set by frequency spacing
t1_range = 0:t_step:t_end ;
fft_sh = repmat(exp(om_s*t1_range),length(supop_ge),1);

%length of the ground + 1st excited manifold, range of time_seperations and
%number of transitions
rho_neq_g = zeros([sz2,sz2,length(t1_range),N,N]); %this ends up in the gs manifold
rho_neq_e = zeros([sz2*N,sz2*N,length(t1_range),N,N]); %end in ex-state manifold

if nargin == 8 %include %gs-2nd ex-state coherences
rho_neq_fg = zeros([sz2,N*(N-1)/2*sz2,length(t1_range),N,N]); 
end

for e1 = 1:N %first interation
        
    op1 = V_eg{e1}; 
    temp1 = op1*rho_0; temp1 = temp1(sz2+1:sz2*(N+1),1:sz2);
    temp1  = reshape(temp1,(numel(temp1)),1);    
        
    Gamma_s = real(R_red_sc(1,1+e1,1,1+e1)); %~coherence decay rate

    %additional rescaling om_eg om_ge beyond om_s, om_s is simply the shift
    %of the actual frequency range from the Fourier transform, this will
    %need to be shifted
    
    om_shift = H_exciton(1+(1+e1)*sz2,1+(1+e1)*sz2)-H_exciton(1,1); 
    om_eg =  om_shift;

    num_scaling_eg = exp((1i*om_eg-Gamma_s)*t1_range); 
    
    %approx scale out oscillations and decay from electronic DOF for easier
    %numerical solutions, could generalise for different rates for each
    %element
    
    supop_eg_scaled = supop_eg  + (1i*om_eg+Gamma_s)*sparse(eye(length(supop_eg)));
     de_fn_eg_scaled  = @(t,v) supop_eg_scaled*v;   
     
    output_DE_fun(length(t1_range),temp1,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);

    ode45(de_fn_eg_scaled,[0,t_end],temp1,options);

    
    [tmp1,tmp2]  =  output_DE_fun(length(t1_range),temp1,'get_data+clean');

    %clean gets rid of empty points with no data, now interpolate to
    %correct time range and reintroduce scaling function
    tmp3 = (interp1(tmp1,tmp2,t1_range).').*repmat(num_scaling_eg,length(temp1),1); 
    %flip matrix so time is in dimension two
    
    %Fourier transform into frequency space, range will be around om_char(e1)
    rho_eg = fftshift(fft(tmp3.*fft_sh,[],2),2)/length(t1_range);
        %reshape this to a matrix with 3rd dimension time and introduce
        %prefactor of 1i, 
    rho_eg = 1i*reshape(rho_eg,sz2*N,sz2,length(t1_range));
    %now the other part of the first order matrix is the HC of this
    
    %rho^(1) = [0,rho_eg'; rho_eg,0];
    
    for e2 = 1:N %loop over second operator
        
        op_2 = V_ge{e2}; %select elements
        op_2 = op_2(1:sz2,sz2+1:sz2*(N+1));%take only the important section
  
%  rho_g^(2) = [i(V_ge*rho_eg - rho_eg'*V_eg)]= iV_ge*rho_eg + h.c .       
        
        tmp1 = 1i*mtimesx(rho_eg,op_2);
 rho_neq_g(:,:,:,e1,e2) = tmp1 +conj(permute(tmp1,[2,1,3]));  
 
 %  rho_e^(2) = [i(V_eg*rho_eg' - rho_eg*V_eg)]= -i*rho_eg*V_eg + h.c .  
 
        tmp1 = -1i*mtimesx(op_2,'C',rho_eg); %V_eg = (V_ge)' obviously
        %in mtimesx the C means conjugate the first term
 
 rho_neq_e (:,:,:,e1,e2)  = tmp1 +conj(permute(tmp1,[2,1,3]));         
             
            %prop these in time for selected frequencies to get the density
            %matrix at any time step.  For shorter pulse seperations this
            %will not be true  
            
    end
      if nargin == 8
    for f2 = 1:N*(N-1)/2 %loop over e-f transitions
    % rho_fg =  i V_fe rho_eg
 
    op_3 = V_fe{f2}; %select elements
    tmp1 = 1i*mtimesx(op_3,rho_eg);
    rho_neq_fg (:,:,:,e1,e2) = tmp1 +conj(permute(tmp1,[2,1,3]));    
    % rho_neq_gf is just the complex conjugate to this
    end
      end
end