function  [rho_neq_g,rho_neq_e,rho_neq_fg] = doorway_fun_basic(rho_0,...
          t1_range,E_u_t,om_u_range,supop ,V_ge,V_eg,N,sz2,V_fe) 
%This version assumes only transitions between excitons with no
%mixing between them that can occur due to coupling to vibrations. 
%t1_range is given in cm^-1 ,E_u_t is Efields integrated aong t'-inf to inf
%calculates the doorway function in segments, but does not propogate it in
%the delay time range before the probe pulse
%sz2 is number of vib states
%this currently includes no line broadening functions
%% Calculate to second order the density matrix before the probe beam hits
%A density matrix associated to each path in Louiville space for
%averaging purposes
if nargout ==3 && nargin < 12
    error('require V_fe to be passed to calculate the ground excited coherences')
end

%calculate int_{-infinity}^{infinity} dt' E_u(t') E_u(t'-t1)
%this is analytic for Gaussian pulses, if not taken as far as infinity it
%would be exp(-t1^2/4tau^2)*erfc((t1-2*t)/(2*tau))/2
if isfloat(E_u_t) %assume a width has been passed
    E_u_fct =  exp(-t1_range.^2/4/E_u_t^2);
else
   E_u_fct  =  E_u_t(t1_range); 
end
E_u_fct = reshape(E_u_fct,1,1,length(t1_range));

if iscell(supop) %pass line broadening functions as first part
    g_t = supop{1}; %g_t should be N by length(t1_range) or symbolic
    c_nk = supop{2};
    supop = supop{3};
    use_lbf = true;
else
    use_lbf = false;
end

sz1 = sqrt(length(supop))/sz2; %total size of elect Hamiltonian
%length of the ground + 1st excited manifold, range of time_seperations and
%number of transitions
rho_neq_g = zeros([sz2,sz2,length(om_u_range),N]); %this ends up in the gs manifold
rho_neq_e = zeros([sz2,sz2,length(om_u_range),N]); %end in ex-state manifold

if nargin == 12 && nargout ==3   %include %gs-2nd ex-state coherences, not tested or finished
rho_neq_fg = zeros([sz2,N*(N-1)/2*sz2,length(t1_range),N,N]); 
end
t_end = t1_range(end);

for e1 = 1:N %first interation is between ground and Nth exciton state
 
    if use_lbf
      
            if size(g_t,1)==1 %nosite dep baths
           g_ex = sum(c_nk(:,e1).*c_nk(:,e1).*c_nk(:,e1).*c_nk(:,e1))*g_t;      
            else
          g_ex = (c_nk(:,e1).*c_nk(:,e1).*c_nk(:,e1).*c_nk(:,e1)).'*g_t;
            end
          if issym(g_t)    
          g_fun = matlabFunction(g_ex);
          else
        lin_interp_g(t1_range,g_ex);      
        g_fun = @(t) lin_interp_g(t); %probably a better way to do this..
          end
    end    
    
    
    op1 = V_eg{e1}; 
    temp1 = op1*rho_0; temp1 = temp1(sz2*e1+1:sz2*(e1+1),1:sz2);
    temp1  = reshape(temp1,(numel(temp1)),1);    
        
%     Gamma_s = real(R_red_sc(1,1+e1,1,1+e1)); %~coherence decay rate (RF)
%     om_shift = H_exciton(1+(1+e1)*sz2,1+(1+e1)*sz2)-H_exciton(1,1); %approx freq
    %num_scaling_eg = exp((1i*om_eg-Gamma_s)*t1_range);     
    %approx scale out oscillations and decay from electronic DOF for easier
    %numerical solutions, could generalise for different rates for each
    %element  in hindsight it is better just to look at the operator
    
    %select only the elements mapping from ground to exciton state e1's
    %coherences, assumes weakish mixing from vibrations between excitons
    tmpeg  = zeros(sz1*sz2);tmpeg(sz2*e1+1:sz2*(e1+1),1:sz2)=1; %lower diag
   % [ab,bb,sb] = find(tmpeg);
  
	tmpeg = logical(reshape(tmpeg,sz2^2*sz1^2,1)); 
    supop_red = supop(tmpeg,tmpeg);
    %scale out the worst of this from ever element
    Gamma_s = -real(supop_red(1,1)); om_eg = -imag(supop_red(1,1));

    supop_eg_scaled = supop_red  + (1i*om_eg+Gamma_s)*sparse(eye(sz2^2));
    if use_lbf
     de_fn_eg_scaled  = @(t,v) supop_eg_scaled*v - g_fun(t).*v;   
    else
     de_fn_eg_scaled  = @(t,v) supop_eg_scaled*v;
    end 
    output_DE_fun(length(t1_range),temp1,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn_eg_scaled,[0,t_end],temp1,options);

    
    [tmp1,tmp2]  =  output_DE_fun(length(t1_range),temp1,'get_data+clean');

    %clean gets rid of empty points with no data, now interpolate to
    %correct time range but don't reintroduce scaling function yet
    rho_eg = interp1(tmp1,tmp2,t1_range).'; 
    %flip matrix so time is in dimension two
    rho_eg = reshape(rho_eg,sz2,sz2,length(t1_range));
%     if e1==1
%        plot(t1_range,reshape(rho_eg,sz2^2,length(t1_range))) 
%     end
    %rho^(1) = [0,rho_eg'; rho_eg,0];
        op_2 = V_ge{e1}; %select elements
        op_2 = op_2(1:sz2,sz2*e1+1:sz2*(e1+1));%take only the important section
  
%  rho_g^(2) = [i(V_ge*rho_eg - rho_eg'*V_eg)]= iV_ge*rho_eg + h.c .       
        for lp = 1:length(om_u_range)
               
        tmp1 = -1i*mtimesx(rho_eg,op_2).*...
                repmat(exp((1i*(om_eg-om_u_range(lp))-Gamma_s).*reshape...
                (t1_range,1,1,length(t1_range))).*E_u_fct,sz2,sz2,1);

 rho_neq_g(:,:,lp,e1)  = trapz(t1_range,tmp1+conj(permute(tmp1,[2,1,3])),3);   
         end
 %  rho_e^(2) = [i(V_eg*rho_eg' - rho_eg*V_eg)]= -i*rho_eg*V_eg + h.c .  
 
        %V_eg = (V_ge)' obviously
        %in mtimesx the C means conjugate the first term
 for lp = 1:length(om_u_range)
     
     tmp1 = 1i*mtimesx(op_2,'C',rho_eg).*... %include extra factor from rescale and E
             repmat(exp((1i*(om_eg-om_u_range(lp))-Gamma_s)...
              .*reshape(t1_range,1,1,length(t1_range))).*E_u_fct,sz2,sz2); 
 rho_neq_e(:,:,lp,e1)  = trapz(t1_range,(tmp1 + conj(permute(tmp1,[2,1,3]))),3);         
 end            
            %prop these in time for selected frequencies to get the density
            %matrix at any time step.  For shorter pulse seperations this
            %will not be true          
            
      if nargin == 12 && nargout ==3 
    for f2 = 1:N*(N-1)/2 %loop over e-f transitions
    % rho_fg =  i V_fe rho_eg
 
    op_3 = V_fe{f2}; %select elements
    tmp1 = 1i*mtimesx(op_3,rho_eg);
    rho_neq_fg (:,:,:,e1,e2) = tmp1 +conj(permute(tmp1,[2,1,3]));    
    % rho_neq_gf is just the complex conjugate to this
    end
      end
end
    function out = lin_interp_g(t,setparams)
        persistent gvals trng
        
        if nargin ==2
           gvals = setparams; trng = t;
            out=[];
            return
        end
        
        tt = find(trng>t,1,'first'); %first value greated
        t2 = trng(tt);
        if abs(t2-t)>10*eps(t2)
            t1 = trng(tt-1); 
        out = gvals(:,tt-1)*(t2-t)  + gvals(:,tt)*(t-t1);
        out = out/(t2-t1);
        else
        out = gvals(:,tt);   
        end

      