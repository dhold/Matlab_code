function  [D_t,W_t,Wf_t,D0,Wgg,W1ee,W2ee]= spec_fun_simple(...
          t1_range,t3_range,om_r_rng,om_u_rng,E_r_t,E_u_t,E1,E2,...
          g_t,g_t_door,lam,c_nk,c_nm_f,N,inc_double,use_fft) 
%%  Calculate the window and doorway function
% Assumes simple format where linebroadening functions (in site basis) are
% given and uses this to calculate the elements of the window function
% No explicit vibrations can be included in the Hamiltonian
% E1,E2 are the energies in the exciton basis, g_t the linbroadening functions
% lam are the sum of all reorg energies in the site basis
if ~exist('inc_double','var')
    inc_double = true;
end
if ~exist('use_fft','var')
    use_fft = true;
end 
 
if isa(E_r_t,'function_handle') 
    E_r_fct = E_r_t(t3_range);
elseif length(E_r_t) == length(t3_range)
    E_r_fct = E_r_t;
elseif ~isempty(E_r_t) %assume a width has been passed
    E_r_fct =  exp(-t3_range.^2/4/E_r_t^2);
end
if isa(E_u_t,'function_handle') 
    E_u_fct = E_u_t(t1_range);
elseif length(E_u_t) == length(t1_range)
    E_u_fct = E_u_t;
elseif ~isempty(E_u_t) %assume a width has been passed
    E_u_fct =  exp(-t1_range.^2/4/E_u_t^2);
end
    
D_t = zeros([length(t1_range),N]);
D_t2 = zeros([length(t3_range),N]);

W_t = zeros([length(t3_range),N]); 
if inc_double
Wf_t = zeros([length(t3_range),N,N*(N-1)/2]);
else
    Wf_t = [];
end
if size(g_t,1)==1
    sm_bth = true;
else
    sm_bth = false;
end
ex_basis_coeff(c_nk,c_nm_f); %pass coefficients to sub function


if use_fft
        t3_r = t3_range(end)-t3_range(1);
        dt3 = t3_range(2)-t3_range(1);
        t1_r = t1_range(end)-t1_range(1);
        dt1 = t1_range(2)-t1_range(1);      
        om_rng_1 = pi*(-1/dt1:2/t1_r:1/dt1);
        om_rng_3 = pi*(-1/dt3:2/t3_r:1/dt3);
end     
    W1ee = zeros(length(om_r_rng),N);
    W2ee = zeros(length(om_r_rng),N,N*(N-1)/2);
    Wgg = zeros(length(om_r_rng),N);
    D0 = zeros(length(om_u_rng),N);
for j = 1:N
 %first calculate the transitions from the ground state
 deltaom = E1(j);
  
 coeff1 = ex_basis_coeff(j,j,j,j,1,sm_bth).';
 
  g_ex_door = coeff1*g_t_door;
  
  D_t(:,j) = exp(-1i*deltaom*t1_range - g_ex_door);
  if use_fft
  tmp = fftshift(ifft(D_t(:,j).'.*E_u_fct));
  D0(:,j) = 2*dt1.*length(t1_range).*real(interp1(om_rng_1,tmp,om_u_rng));
  end
  g_ex = coeff1*g_t;
  
  D_t2(:,j) = exp(-1i*deltaom*t3_range - g_ex);  %need a sep one over t3_range
  if use_fft
  tmp = fftshift(ifft(D_t2(:,j).'.*E_r_fct));
  Wgg(:,j) = 2*dt3.*length(t3_range).*real(interp1(om_rng_3,tmp,om_r_rng));
  end
  
 g_ex = conj(coeff1*g_t);
 E_shift = 2*(coeff1*lam);
                    
  W_t(:,j) = exp(-1i*(deltaom-E_shift)*t3_range - g_ex);
  if use_fft
  tmp = fftshift(ifft(W_t(:,j).'.*E_r_fct));
  W1ee(:,j) = 2*dt3.*length(t3_range)*real(interp1(om_rng_3,tmp,om_r_rng));
  end
  
        
if inc_double
    for f = 1:N*(N-1)/2 %transitions from single ex state e4 to double 
       
        %calculate linbroadening fn
        
        coeff2 = ex_basis_coeff(j,j,f,f,2,sm_bth).';
        coeff3 = ex_basis_coeff(f,f,f,f,3,sm_bth).';
        
 g_ex = (coeff1 + coeff3-2*coeff2)*g_t;
 E_shift = 2*((coeff1-coeff2)*lam);
 %could just add this shift to deltaom
 %extra_shift = 2i*((coeff2-coeff1)*lam);
 
        deltaom = E2(f) - E1(j);
        
        Wf_t(:,j,f) = exp(-1i*(deltaom-E_shift)*t3_range - g_ex);
  if use_fft
  tmp = fftshift(ifft(Wf_t(:,j,f).'.*E_r_fct));
  W2ee(:,j,f) = 2*dt3.*length(t3_range)*real(interp1(om_rng_3,tmp,om_r_rng));
  end

    end
end
end

if~use_fft
 if~isempty(om_r_rng)
    W1ee = zeros(length(om_r_rng),N);
     if inc_double
    W2ee = zeros(length(om_r_rng),N,N*(N-1)/2);
     end
    Wgg = zeros(length(om_r_rng),N);
    if sum(abs(diff(t3_range,2)))<eps(length(t3_range))
        use_simp =true; %equal spacing, use simpsons rule
    else
        use_simp = false;
    end
    for lp = 1:length(om_r_rng)
        
        om_r = om_r_rng(lp);
        tmp = repmat(reshape(E_r_fct.*exp(1i*om_r*t3_range),length(t3_range),1),1,N);
        if use_simp
           Wgg(lp,:) = 2*real(simpsons(D_t2.*tmp,t3_range(1),t3_range(end)));
           W1ee(lp,:) = 2*real(simpsons(W_t.*tmp,t3_range(1),t3_range(end)));
        else
        Wgg(lp,:) = trapz(t3_range,D_t2.*tmp); 
        Wgg(lp,:) = Wgg(lp,:) + conj(Wgg(lp,:));
        W1ee(lp,:) = trapz(t3_range,W_t.*tmp); 
        W1ee(lp,:) = W1ee(lp,:) + conj(W1ee(lp,:));  
        end
 if inc_double
    for f = 1:N*(N-1)/2     
        if use_simp
        W2ee(lp,:,f) = 2*real(simpsons(Wf_t(:,:,f).*tmp,t3_range(1),t3_range(end)));    
        else
        W2ee(lp,:,f) = trapz(t3_range,Wf_t(:,:,f).*tmp);  
        W2ee(lp,:,f) = W2ee(lp,:,f) + conj(W2ee(lp,:,f));
        end
    end
 end
 
    end  
 end
if~isempty(om_u_rng)
 
    D0 = zeros(length(om_u_rng),N);
    
    for lp = 1:length(om_u_rng)
        
        om_u = om_u_rng(lp);
        D0(lp,:) = trapz(t1_range,D_t.*repmat(reshape(E_u_fct.*...
               exp(1i*om_u*t1_range),length(t1_range),1),1,N));    
        D0(lp,:) = D0(lp,:) + conj(D0(lp,:));
    end
    
    
end
end


function coeff = ex_basis_coeff(k1,k2,k3,k4,type,same_bath)
%calculates coefficients for exciton basis, give 6th argument as true iff
%the baths are the same on all sites
persistent c_nk c_nm_f N

if nargin == 2
    c_nk = k1; c_nm_f = k2; N = size(c_nk,1); coeff =[]; return
end
coeff = zeros(N,1);
pickr = eye(N); %picks that element, does nothing if all sites same bath
    if type == 1 %k k' k'' k'''
       coeff = c_nk(:,k1).*c_nk(:,k2).*c_nk(:,k3).*c_nk(:,k4);         
    elseif type ==2 %kk qq function  
        f = k3; check_v = k1==k2 && k4==k3;
        if ~check_v
           error('you must pass k1==k2 and k3==k4') 
        end
       for m=2:N
            for n = 1:m-1
       coeff = coeff + (c_nm_f(n,m,f)*c_nk(n,k1))^2*pickr(:,n)...
                    + (c_nm_f(n,m,f)*c_nk(m,k1))^2*pickr(:,m);
            end
       end
    elseif type ==3 %qqqq
        f = k1; check_v = k1==k2 && k1==k3 && k1==k4;
        if ~check_v
           error('only works for all inputs same') 
        end
       for m=2:N
            for n = 1:m-1
                for mm = 2:N
                    for nn = 1:mm-1
                        if n==nn || n==mm 
                         coeff = coeff + (c_nm_f(n,m,f)*c_nm_f(nn,mm,f))^2*pickr(:,n);                       
                        end
                        if  m==mm || m==nn
                         coeff = coeff + (c_nm_f(n,m,f)*c_nm_f(nn,mm,f))^2*pickr(:,m);                       
                        end                                               
                    end
                end       
            end
       end
       
    end
    if nargin==6
        if same_bath
       coeff = sum(coeff); 
        end
    end