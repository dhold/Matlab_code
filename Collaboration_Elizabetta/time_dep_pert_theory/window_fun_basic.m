function  [V_gg_t,V_ee_t,V_eeff_t,V_fg_t] = window_fun_basic(...
          t3_range,supop ,V_ge,V_eg,N,sz2,sz1,inc_double) 
%%  Calculate the window function
% Assume the pulse is short enough that no time evolution occurs during the
% interaction.
% Do not perform integral over t1 as this function is useful for
% calculating linear spectra as well
if ~exist('inc_double','var')
    inc_double = true;
end
%length of the ground + 1st excited manifold, range of time_seperations and
%number of transitions
V_gg_t = zeros([sz2,sz2,length(t3_range),N]); %ge->gg
V_ee_t = zeros([sz2,sz2,length(t3_range),N]); %ge->ee
if inc_double
V_eeff_t = zeros([sz2,sz2,length(t3_range),N,N*(N-1)/2]); %ef->ee
end

if nargout ==4   %include %gs-2nd ex-state coherences, not tested or finished
V_fg_t = zeros([sz2,N*(N-1)/2*sz2,length(t1_range),N,N]); 
end
t_end = t3_range(end);

use_simple_exp =false; 
if iscell(supop) %pass line broadening functions as first part
    g_t = supop{1}; %g_t should be N by length(t1_range) or symbolic
    lam_tot = supop{2};
     %make same size as g_t
    c_nk = supop{3};
    c_nm_f = supop{4};
    if length(supop)==6
       use_simple_exp =true;    
       lam_tot = lam_tot.*repmat(t3_range,size(lam_tot,1),1);
    %assumes g_t is the exponential comp of soln not derivative!
    else
    lam_tot = repmat(lam_tot,1,length(t3_range));
    end
    supop = supop{5}; %remaining operator, write over other mat
    use_lbf = true;
else
    use_lbf = false;
end

for e4 = 1:N

      if use_simple_exp %don't pass a sym g, g_t should not be the derivative
          %assumes sz2=1
            if size(g_t,1)==1 %nosite dep baths
           g_ex = sum(c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4))*g_t;      
            else
          g_ex = (c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4)).'*g_t;
            end    
       tmpge  = zeros(sz1);tmpge(1,e4+1)=1;
	tmpge = logical(reshape(tmpge,sz1^2,1)); 
    trans_freq = imag(full(supop(tmpge,tmpge)));        

   V_gg_t(1,1,:,e4) = 1i*exp(-1i*trans_freq*t3_range-g_ex);  
   V_gg_t(1,1,:,e4) = V_gg_t(1,1,:,e4) + conj(permute(V_gg_t(1,1,:,e4),[2,1,3]));
        if size(g_t,1)==1 %assume all sites have same bath
       g_ex = sum(c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4))...
                        *(-2i*lam_tot+conj(g_t));                  
        else
        g_ex = (c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4)).'...
                        *(-2i*lam_tot +conj(g_t));   
        end          
   V_ee_t(1,1,:,e4) =  -1i*exp(-1i*trans_freq*t3_range-g_ex);       
   V_ee_t(1,1,:,e4) = V_ee_t(1,1,:,e4) + conj(permute(V_ee_t(1,1,:,e4),[2,1,3]));
        
      else  %solve the differential equation numerically
    if use_lbf     
            if size(g_t,1)==1 %nosite dep baths
           g_ex = sum(c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4))*g_t;      
            else
          g_ex = (c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4)).'*g_t;
            end
          if isa(g_t,'sym')    
          g_fun = matlabFunction(g_ex);
          else
        lin_interp_g(t3_range,g_ex);      
        g_fun = @(t) lin_interp_g(t); %probably a better way to do this..
          end
    end        

    init_state = V_ge{e4}; init_state = init_state(1:sz2,sz2*e4+1:sz2*(e4+1));
    init_state  = reshape(init_state,(numel(init_state)),1).';    %select elements

    tmpge  = zeros(sz1*sz2);tmpge(sz2*e4+1:sz2*(e4+1),1:sz2)=1;%tmpge(1:sz2,sz2*e4+1:sz2*(e4+1))=1;
	tmpge = logical(reshape(tmpge,sz2^2*sz1^2,1)); 
    supop_red = supop(tmpge,tmpge);
    if ~use_lbf
    Gamma_s = -real(supop_red(1,1)); om_eg = -imag(supop_red(1,1));
    else %also scale out Markovian contributions
    Gamma_s = -real(supop_red(1,1))+real(g_fun(t3_range(end)));
    om_eg = -imag(supop_red(1,1))+imag(g_fun(t3_range(end)));
    end

    num_scaling_ge = exp((-1i*om_eg-Gamma_s)*t3_range); 
    supop_scaled = supop_red  + (1i*om_eg+Gamma_s)*eye(sz2^2);

    if use_lbf
     de_fn = @(t,v) mtimesx(v,'T',supop_scaled).' -g_fun(t).*v;   
    else
     de_fn = @(t,v) mtimesx(v,'T',supop_scaled).';   
    end 
        
    %ode45 only takes column so need two transposes
    
    output_DE_fun(length(t3_range),init_state,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn ,[t3_range(1),t_end],init_state,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),init_state ,'get_data+clean');
    
    if sz2 == 1
    tmp3 = interp1(tmp1,tmp2,t3_range) ;    
    tmp3 = tmp3  .* num_scaling_ge; 
    else
    tmp3 = (interp1(tmp1,tmp2,t3_range).') ;    
    tmp3 = tmp3  .* repmat(num_scaling_ge,size(tmp2,2),1);
    end
      
     op_3 = V_eg{e4};    op_3 = op_3(sz2*e4+1:sz2*(e4+1),1:sz2); 
     tmp3 = reshape(full(tmp3),sz2,sz2,length(t3_range));

        V_gg_t(:,:,:,e4) = -1i*mtimesx(tmp3,op_3);  
        V_gg_t(:,:,:,e4) =  V_gg_t(:,:,:,e4) + conj(permute(V_gg_t(:,:,:,e4),[2,1,3]));    
        
 %Next calculate elements ending up in the excited state   
    
    if use_lbf
        
        if size(g_t,1)==1 %assume all sites have same bath
       g_ex = sum(c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4))...
                        *(-2i*lam_tot +conj(g_t));                  
        else
        g_ex = (c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4)).'...
                        *(-2i*lam_tot +conj(g_t));   
        end
        if isa(g_t,'sym')
            g_fun = matlabFunction(g_ex);
        else
            lin_interp_g(t3_range,g_ex);      
            g_fun = @(t) lin_interp_g(t); %probably a better way to do this..
        end
    end    
    
    %temp1 = V_eg{e4}; temp1 = temp1(sz2*e4+1:sz2*(e4+1),1:sz2);
    init_state = V_ge{e4}; init_state = init_state(1:sz2,sz2*e4+1:sz2*(e4+1));
    init_state  = reshape(init_state,(numel(init_state)),1).';    %select elements
    %temp11 = eye(sz2);  %temp11 = reshape(temp11,(numel(temp11)),1).'; 
    %isequal(temp11,temp1)
  
    %select only the elements mapping from ground to exciton state e4's
    %coherences, assumes weakish mixing from vibrations between excitons
    tmpge  = zeros(sz1*sz2);tmpge(sz2*e4+1:sz2*(e4+1),1:sz2)=1;
	tmpge = logical(reshape(tmpge,sz2^2*sz1^2,1)); 
    supop_red = supop(tmpge,tmpge);
    %scale out the worst of this from ever element
    if ~use_lbf
    Gamma_s = -real(supop_red(1,1)); om_eg = -imag(supop_red(1,1));
    else %also scale out Markovian contributions
    Gamma_s = -real(supop_red(1,1))+real(g_fun(t3_range(end)));
    om_eg = -imag(supop_red(1,1))+imag(g_fun(t3_range(end)));
    end

    num_scaling_ge = exp((-1i*om_eg-Gamma_s)*t3_range); 
    supop_scaled = supop_red  + (1i*om_eg+Gamma_s)*eye(sz2^2);

    if use_lbf
     de_fn = @(t,v) mtimesx(v,'T',supop_scaled).' -g_fun(t).*v;   
    else
     de_fn = @(t,v) mtimesx(v,'T',supop_scaled).';   
    end 
        
    %ode45 only takes column so need two transposes
    
    output_DE_fun(length(t3_range),init_state,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn ,[t3_range(1),t_end],init_state,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),init_state ,'get_data+clean');

    if sz2 == 1
    tmp3 = interp1(tmp1,tmp2,t3_range) ;    
    tmp3 = tmp3  .* num_scaling_ge; 
    else
    tmp3 = (interp1(tmp1,tmp2,t3_range).') ;    
    tmp3 = tmp3  .* repmat(num_scaling_ge,size(tmp2,2),1);
    end  

     op_3 = V_eg{e4};    op_3 = op_3(sz2*e4+1:sz2*(e4+1),1:sz2); 
     tmp3 = reshape(full(tmp3),sz2,sz2,length(t3_range));

        V_ee_t(:,:,:,e4) = +1i*mtimesx(op_3,tmp3);       
        V_ee_t(:,:,:,e4) = V_ee_t(:,:,:,e4) + conj(permute(V_ee_t(:,:,:,e4),[2,1,3]));
      end
 
if inc_double
    for f = 1:N*(N-1)/2 %transitions from single ex state e4 to double 
        %ex manifold state f.   
    %These are in principle of the same form as those of ground - single ex
      
    if use_lbf
       if size(g_t,1)==1 %assume all sites have same bath
       tmp1 = sum(c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4));     
       tmp = 0; tmp2 = 0;
       for m=2:N
            for n = 1:m-1
       tmp = tmp + (c_nm_f(n,m,f)*c_nk(n,e4))^2+ (c_nm_f(n,m,f)*c_nk(m,e4))^2;
                for mm = 2:N
                    for nn = 1:mm-1
                        if n==nn || n==mm 
                         tmp2 = tmp2 + (c_nm_f(n,m,f)*c_nm_f(nn,mm,f))^2;                       
                        end
                        if  m==mm || m==nn
                         tmp2 = tmp2 + (c_nm_f(n,m,f)*c_nm_f(nn,mm,f))^2;                       
                        end                                               
                    end
                end       
            end
       end
       
       if tmp1+tmp2-2*tmp <0 %if this is negative grows exponentially 
           f 
       end
%(tmp1+tmp2-2*tmp)*g_t  +2i*(tmp1-tmp)*lam_tot
       g_ex = tmp1*(g_t+2i*lam_tot) -2*tmp*(g_t+1i*lam_tot) +tmp2*g_t;      

        else %each site different bath
 
        g_ex = (c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4).*c_nk(:,e4)).'*(g_t+2i*lam_tot);
 
            tmp = zeros(size(t3_range)); tmp2 = zeros(size(t3_range));
            for m=2:N
                for n = 1:m-1

            tmp = tmp + (c_nm_f(n,m,f)*c_nk(n,e4))^2*(g_t(n,:)+1i*lam_tot(n,:));
            tmp = tmp + (c_nm_f(n,m,f)*c_nk(m,e4))^2*(g_t(m,:)+1i*lam_tot(m,:));

                for mm = 2:N
                    for nn = 1:mm-1
                        if n==nn || n==mm
                         tmp2 = tmp2 + (c_nm_f(n,m,f)*c_nm_f(nn,mm,f))^2*g_t(n,:); 
                        end 
                        if m==mm || m==nn
                         tmp2 = tmp2 + (c_nm_f(n,m,f)*c_nm_f(nn,mm,f))^2*g_t(m,:);          
                        end
                    end
                end
                end
            end
   
            g_ex = g_ex - 2*tmp+tmp2;
       end
        if isa(g_t,'sym')   
            g_fun = matlabFunction(g_ex);
        else
            lin_interp_g(t3_range,g_ex);      
            g_fun = @(t) lin_interp_g(t); %probably a better way to do this..
        end
    end    

    if use_simple_exp
        
    tmpef = zeros(sz1);  tmpef(N+f+1,e4+1)=1;%tmpef(e4+1,N+f+1)=1;
	tmpef = logical(reshape(tmpef,sz1^2,1));         
        trans_freq = imag(full(supop(tmpef,tmpef)));
        
        V_eeff_t(:,:,:,e4,f) = +1i*exp(-1i*trans_freq*t3_range-g_ex);        
        V_eeff_t(:,:,:,e4,f) = V_eeff_t(:,:,:,e4,f) + conj(permute(...
                               V_eeff_t(:,:,:,e4,f),[2,1,3]));         
    else
    
    Vef_L = eye(sz2);
    Vef_L = reshape( Vef_L ,(numel( Vef_L )),1).';    %select elements
    
    tmpef = zeros(sz1*sz2);
    tmpef(sz2*e4+1:sz2*(e4+1),sz2*N+(sz2*f+1:sz2*(f+1)))=1;
	tmpef = logical(reshape(tmpef,sz2^2*sz1^2,1)); 
    %select elements which contribute to evolution
    supop_red = supop(tmpef,tmpef); 
    %scale out the worst of this from ever element
    
     if ~use_lbf
    Gamma_f = -min(real(full(supop_red(:)))); om_ef = -imag(full(supop_red(1,1)));  
    else %also scale out Markovian contributions, I.E. late time
    Gamma_f = -min(real(full(supop_red(:))))+real(g_fun(t3_range(end)));
    om_ef = -imag(supop_red(1,1))+imag(g_fun(t3_range(end)));
    end   
    
  %  [Gamma_f,om_ef]
    num_scaling_ef = exp((-1i*om_ef-Gamma_f)*t3_range); 

   supop_ef_scaled = supop_red + (1i*om_ef+Gamma_f)*sparse(eye(length(supop_red))); 
    if use_lbf
     de_fn_bk_ef = @(t,v) mtimesx(v,'T',supop_ef_scaled).' -g_fun(t).*v;   
    else
      de_fn_bk_ef = @(t,v) mtimesx(v,'T',supop_ef_scaled).';    
    end        
        
        
    output_DE_fun(length(t3_range),Vef_L,'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-15);
    ode45(de_fn_bk_ef ,[t3_range(1),t_end],Vef_L,options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t3_range),Vef_L,'get_data+clean'); 
    if sz2 == 1
    tmp4 = interp1(tmp1,tmp2,t3_range,'pchip')  ;    
    tmp4 = tmp4.*num_scaling_ef;      
    else
    tmp4 = (interp1(tmp1,tmp2,t3_range,'pchip').')  ;    
    tmp4 = tmp4.*repmat(num_scaling_ef,size(tmp2,2),1);
    end
    tmp4 = reshape(full(tmp4),sz2,sz2,length(t3_range));

        V_eeff_t(:,:,:,e4,f) = +1i*mtimesx(tmp4,eye(sz2),'C');       
        V_eeff_t(:,:,:,e4,f) = V_eeff_t(:,:,:,e4,f) + ...
                            conj(permute(V_eeff_t(:,:,:,e4,f),[2,1,3]));    
    end
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

        tt = find(trng>=t,1,'first'); %first value greated
        t2 = trng(tt);
        if abs(t2-t)>10*eps(t2)
            t1 = trng(tt-1); 
        out = gvals(:,tt-1)*(t2-t)  + gvals(:,tt)*(t-t1);
        out = out/(t2-t1);
        else
        if isempty(tt)
            out=0; return
        end
        out = gvals(:,tt);
        end