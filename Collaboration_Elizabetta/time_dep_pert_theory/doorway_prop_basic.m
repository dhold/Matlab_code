function [rho_t_g,rho_t_e] = doorway_prop_basic(rho_neq_g,rho_neq_e,...
                t_sep_range,supop_gg,supop_ee,sz2,N,tol,pop_only) 

        if nargin <8
            tol = 1e-10; %solver absolute tolerance
        end
        if nargin <9
            pop_only = false; %consider coherences between states
        end        
 rho_t_g = zeros([sz2,sz2,length(t_sep_range),size(rho_neq_g,3),size(rho_neq_g,4)]);
 rho_t_e = zeros([sz2*N,sz2*N,length(t_sep_range),size(rho_neq_e,3),size(rho_neq_e,4)]); 
 %excited state manifold has to be larger due to the exciton mixing that
 %occurs in the waiting time.  Initial exciton states spread

  de_gg = @(t,v) supop_gg*v;
 if pop_only 
     %just consider the population dynamics in ex state, assume coherences are zero
     if length(supop_ee)~=N*sz2 %else it have already been reduced
     tmpe  = eye(N*sz2); 
	tmpe = logical(reshape(tmpe,N^2*sz2^2,1));    
    supop_ee = supop_ee(tmpe,tmpe);        
     end
 end
  de_ee = @(t,v) supop_ee*v;
  
 t_end = t_sep_range(end);
 
for j=1:size(rho_neq_g,3) %loop over pump frequencies
        
    for j2 = 1:size(rho_neq_g,4) %loop over each of the possible initial interactions
    
        temp1 = rho_neq_g(:,:,j,j2); %initial condition
        trace_val1 = trace(temp1);
        if trace_val1~=0 %otherwise it doesn't contribute
        temp1 = reshape(temp1,numel(temp1),1)/trace_val1;
    output_DE_fun(length(t_sep_range),temp1,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',tol);
   % tic
    ode45(de_gg,[0,t_end],temp1,options);
   % toc

    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_range),temp1,'get_data+clean');
    rho_gg = interp1(tmp1,tmp2,t_sep_range).'; 
    rho_gg = reshape(rho_gg,sz2,sz2,length(t_sep_range));

        rho_t_g(:,:,:,j,j2) = trace_val1*rho_gg;
        
        %now consider excited state manifold dynamics (much slower)
     temp2a = rho_neq_e(:,:,j,j2); %initial condition
     trace_val2 = trace(temp2a); %set trace to unity
     %expand to the full size of the 1st ex state manifold
     temp2b = zeros(N*sz2,N*sz2);
     temp2b(1+(j2-1)*sz2:j2*sz2,1+(j2-1)*sz2:j2*sz2) = temp2a/trace_val2;
     if ~pop_only
     temp2b = reshape(temp2b,numel(temp2b),1);
     else
        temp2b = diag(temp2b);    
     end
     
    output_DE_fun(length(t_sep_range),temp2b,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',tol);
 %   tic
    ode45(de_ee,[0,t_end],temp2b,options);     
 %   toc
     
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_range),temp2b,'get_data+clean');
    rho_ee = interp1(tmp1,tmp2,t_sep_range).'; 
    if ~pop_only
    rho_ee = reshape(rho_ee,sz2*N,sz2*N,length(t_sep_range));     
    rho_t_e(:,:,:,j,j2) = trace_val2*rho_ee;
    else
        for k = 1:sz2*N
     rho_t_e(k,k,:,j,j2)   = trace_val2*rho_ee(k,:);
        end
    end
        end
    end
end