function  [rho_neq_g_t,rho_neq_e_t,rho_neq_fg_t] = doorway_prop(...
          om_rng,E_u_w,om_u_range,N,sz2,t_sep_rng,...
          supop_g,supop_e,rho_neq_g,rho_neq_e,supop_fg,rho_neq_fg) 
%This function calculates the propogation of the doorway matrix with the
%pump probe delay time tau
%  N = number of sites
% sz2 = size of vibrational manifold
% E_u_w = (@(w) fn handle) electric field envelope from pump E field in freq space
% om_u_range = range of carrier frequencies considered
% om_rng = freq range for which initial states are calculated
% rho_neq_g  %this ends up in the gs manifold
% rho_neq_e %end in ex-state manifold
% t_sep_rng are the time points for the seperation tau

de_fn_gg = @(t,v) supop_g*v;
de_fn_ee = @(t,v) supop_e*v;
om_rng = reshape(om_rng,1,length(om_rng)); %make horizontal

rho_neq_g_t = zeros(length(supop_g),length(t_sep_rng),length(om_u_rng),N,N);
rho_neq_e_t = zeros(length(supop_e),length(t_sep_rng),length(om_u_rng),N,N);
if narout==3
rho_neq_fg_t = zeros(sz2^2*N*(N-1)/2,length(t_sep_rng),length(om_u_rng),N,N);
end

for j = 1:length(om_u_range) %loop over carrier envelope freq
        
    freq_rng1 = om_rng-om_u_range(j);
   % freq_rng2 = om_rng+om_u_range(j); %the bit usually lost in RWA
    
    E_fct1 = abs(E_u_w(freq_rng1)).^2; 
   % E_fct2 = abs(E_u_w(freq_rng2)).^2; 
for e1 = 1:N
    for e2 = 1:N
%     om_shift = H_exciton(1+(1+e1)*sz2,1+(1+e1)*sz2)-H_exciton(1,1); 
%     om_eg =  om_shift;
% 
%     num_scaling_eg = exp((1i*om_eg-Gamma_s)*t1_range); 
%     
%     %approx scale out oscillations and decay from electronic DOF for easier
%     %numerical solutions, could generalise for different rates for each
%     %element
%     
%     supop_eg_scaled = supop_eg  + (1i*om_eg+Gamma_s)*sparse(eye(length(supop_eg)));
%      de_fn_eg_scaled  = @(t,v) supop_eg_scaled*v;   
     
%first ground state
    tmpg = reshape(rho_neq_g(:,:,:,e1,e2),length(supop_g),length(om_rng));
    tempg = trapz(om_rng,tmpg.*repmat(E_fct1,length(supop_g),1),2);

    output_DE_fun(length(t1_range),tempg,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-12);
    ode45(de_fn_gg,[0,t_sep_rng(end)],tempg,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_rng),tempg,'get_data+clean');

    %clean gets rid of empty points with no data, now interpolate to
    %correct time range and reintroduce scaling function
    tmp3 = (interp1(tmp1,tmp2,t_sep_rng).'); 
    rho_neq_g_t(:,:,j,e1,e2) = tmp3;

        %now for the excited state
    tmpe = reshape(rho_neq_e(:,:,:,e1,e2),length(supop_e),length(om_rng));
    tempe = trapz(om_rng,tmpe.*repmat(E_fct1,length(supop_e),1),2);

    output_DE_fun(length(t1_range),tempe,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-12);
    ode45(de_fn_ee,[0,t_sep_rng(end)],tempe,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_rng),tempe,'get_data+clean');
    tmp3 = (interp1(tmp1,tmp2,t_sep_rng).');    
    rho_neq_e_t(:,:,j,e1,e2) = tmp3;    
    

      if narout==3
        for f2 = 1:N*(N-1)/2
    tmpfg = reshape(rho_neq_fg(:,:,:,e1,e2),length(supop_fg),length(om_rng));
    tempfg = trapz(om_rng,tmpfg.*repmat(E_fct1,length(supop_fg),1),2);

    output_DE_fun(length(t1_range),tempfg,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-12);
    ode45(de_fn_fg,[0,t_sep_rng(end)],tempfg,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_rng),tempfg,'get_data+clean');
    tmp3 = (interp1(tmp1,tmp2,t_sep_rng).');    
    rho_neq_fg_t(:,:,j,e1,e2) = tmp3;                
        end
      end
    end
end
end