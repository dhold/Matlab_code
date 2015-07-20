function  [D_markov,D_t_markov,W_gg_markov,W_kk_markov,W_qq_markov,R_red_markov,fock_rep]...
    = doorway_fun_simple(H_single,H_double,tau,om_u_rng,om_r_rng,t_sep_rng,params) 
% Calculates the doorway function assuming a diagonal matrix for propogation
%uses line shapes from Zhang et al 1998, sums over all transitions which
%are assumed to be vertical
[coeff1,H_single_ex]  = eig(H_single); %1 ex manifold of H_exciton
[coeff2,H_double_ex]  = eig(H_double); %ground and double excited
%H_single is site ordered, 1:N are sites 1:N
%H_double is groundstate a double excited states ordered
%[1,2],[1,3],...,[1,N], [2,3],[2,4],...,[2,N],..,[N-1,N]
N = length(H_single);
fock_rep = zeros(length(H_single)+length(H_double),N);
fock_rep(2:N+1,1:N) = eye(N); cnt1 = 1; cnt2 =2;
for q = 1:length(H_double)-1 
    fock_rep(N+1+q,cnt1) = 1; fock_rep(N+1+q,cnt2) = 1;
    if cnt2 == N
        cnt1 = cnt1+1; cnt2=1;
    else
        cnt2 = cnt2 +1;
    end          
end
N2 = length(H_double)-1;
fock_rep = logical(fock_rep);

%H_tot = blkdiag(H_double(1,1),H_single,H_double(2:end,2:end));
%H_ex = blkdiag(H_double_ex(1,1),H_single_ex,H_double_ex(2:end,2:end));


B = params{1}; %beta  = 1/k_B T 
gam_rng= params{2}; om0= params{3}; lam_rng = params{4}; %underdamped
gam_dru= params{5}; lam_dru= params{6}; %overdamped on each site
lambda_tot = cellfun(@sum, lam_dru) + cellfun(@sum,lam_rng); %size N
lambda_tot = reshape(lambda_tot,N,1);
%create line broadening functions
%syms lambda om_0 gamma real 
%syms beta t v_n real

% [g1,g2,gmat] = line_broadening_fn(lambda,om_0,gamma, beta);
% g_ud = matlabFunction(g1+g2,'vars',{lambda,om_0,gamma,beta,t});
% g_ud_mat = matlabFunction(gmatsu,'vars',{lambda,om_0,gamma,beta,t,v_n});
% 
% [G1,G2,Gmat] = line_broadening_fn(lambda,[],gamma, beta);
% g_od = matlabFunction(G1+G2,'vars',{lambda,gamma,beta,t});
% g_od_mat = matlabFunction(gmatsu,'vars',{lambda,gamma,beta,t,v_n});

%Electric field parameters
tau_u = tau(1); tau_r = tau(2);
%E_u = @(t) exp(-t.^2/2/tau_u^2)./sqrt(sqrt(pi)*tau_u);
%E_u_w = @(w) exp(-w.^2*tau_u^2/2).*sqrt(tau_u)/sqrt(sqrt(pi));
%E_r = @(t) exp(-t.^2/2/tau_r^2)./sqrt(sqrt(pi)*tau_r);
%E_r_w = @(w) exp(-w.^2*tau_r^2/2).*sqrt(tau_r/sqrt(pi));
%Calculate line broadening functions for each site
g_markov = zeros(N,1);  Gam_site = g_markov; 
for j = 1:N %site loop  
    g_markov(j) = line_broad_fn_markov(B,gam_rng{j},om0{j},lam_rng{j},gam_dru{j},lam_dru{j});
    Gam_site(j) = real(g_markov(j)); %damping due to bath coupling
  %  lambda_site(j) = -imag(g_markov(j)); %frequency shift due to bath
end %ft of this is just 1/(g-i om) = (g+iom)/(g^2+om^2)
%lambda_site ==lambda_tot %just to test these are indeed, the same

% g_fn = @(g,k1,k2,k3,k4) (coeff1(k1,1:N).*coeff1(k2,1:N).*coeff1(k3,1:N)...
%                     .*coeff1(k4,1:N))*g(1:N,:);
%                 %sum_{n<m},n' c^q_{nm} c^k_n'(delta_nn'+delta_mn')g_n'(t)
% g_fn2 = @(g,k,q,j) (coeff2(1+q,j).*coeff1(k,fock_rep(j,:))).^2*g(fock_rep(j,:),:);                
%                 %need to call for all j and sum 
%                 
% g_fn3 = @(g,q,j) (coeff2(1+q,j).*coeff1(k,fock_rep(j,:))).^2*g(fock_rep(j,:),:);        
%pass these to functions, could use global variables but I really dont know
%how to do this...
g_kkkk(coeff1); 
g_kkqq(coeff1,coeff2,fock_rep);
g_qqqq(coeff2,fock_rep);

%% Calculate Doorway and window functions

D_markov = zeros(N,length(om_u_rng)); 
W_gg_markov = zeros(N,length(om_r_rng));
Gamma_ex = zeros(N,1); lambda_ex = Gamma_ex;
for j=1:N
    
        lambda_ex(j) = g_kkkk(lambda_tot,j,j,j,j); %-imag part
        Gamma_ex(j) = g_kkkk(Gam_site,j,j,j,j); %real part of g
        for lp = 1:length(om_u_rng)
            
            om_u = om_u_rng(lp);
        deltaom = (H_single_ex(j,j)-om_u); %diff between energy and carrier wavelength               
                    
D_markov(j,lp) =  4*sqrt(pi)*real(Faddeeva_w(tau_u*(deltaom - lambda_ex(j) + 1i*Gamma_ex(j))));            
                  
        end
        for lp = 1:length(om_r_rng)
            om_r = om_r_rng(lp);
        deltaom = (H_single_ex(j,j)-om_r); %diff between energy and carrier wavelength
                    % - 1i*(deltaom+om_int_range) term vanished when cc taken
       %this can be evaluated analytically to be a voigt profile             
W_gg_markov(j,lp) = 4*sqrt(pi)*real(Faddeeva_w(tau_r*(lambda_ex(j)-deltaom +  1i*Gamma_ex(j))));
        end
end
%hole in the ground state is just - the sum of this, should be down
%dimension 1 so this is fine
%D_hole = -sum(D_full); D_hole_mark = -sum(D_markov);

%now calculate window functions for the stim-emission and ES-abs conts
%W_kk_markov = zeros(N,length(om_r_rng)); 
W_kk_markov = W_gg_markov; %within the Markovian approx these are the same
W_qq_markov = zeros(N,N2,length(om_r_rng));

%calculate linebroadening fns (g-matricies) in exciton basis
 lambda_kq = zeros(N,N2);  lambda_qq = zeros(N2,1); 
Gamma_kq = zeros(N,N2);  Gamma_qq = zeros(N2,1);
for j = 1:N  
    for q = 1:N2
        if j==1
        Gamma_qq(q) = g_qqqq(Gam_site,q);
        lambda_qq(q) = g_qqqq(lambda_tot,q);        
        end
        Gamma_kq(j,q) = g_kkqq(Gam_site,j,q);
        lambda_kq(j,q) = g_kkqq(lambda_tot,j,q);
    end
end

for lp = 1:length(om_r_rng)
    om_u = om_r_rng(lp);
    for j = 1:N
        for q = 1:N2 %Note H_double also includes ground state
    deltaom_q = (H_double_ex(q+1,q+1)-H_single_ex(j,j)-om_u);          
    
    W_qq_markov(j,q,lp) = 4*sqrt(pi)*real(Faddeeva_w(tau_r*(...
        lambda_qq(q)-lambda_ex(j)-deltaom_q +...
        1i*(Gamma_ex(j)+Gamma_qq(q)-2*Gamma_kq(j,q)))));        
 
        end
    end

end

%% Calculate evolution of the doorway wavepacket via modified redfield

R_red_markov = zeros(N,N); %pretty sure this stays as zero since there is
%no mixing from the bath between exciton states in this shitty
%approximation
for k =1:N
    for j = 1:N
        if k ~=j
            deltaom = (H_single_ex(j,j)-H_single_ex(k,k)); 
            deltaom = deltaom-lambda_ex(k)+lambda_ex(j);
            
            Gamma_fc = (Gamma_ex(j)+Gamma_ex(k)-2*g_kkkk(Gam_site,k,k,j,j));        
 
            tmp1 = g_kkkk(Gam_site,j,k,j,j) - g_kkkk(Gam_site,k,j,k,k)+...
                   1i*(g_kkkk(lambda_tot,j,k,j,j) - g_kkkk(lambda_tot,k,j,k,k)) ;
            tmp2 = g_kkkk(Gam_site,j,j,k,j) - g_kkkk(Gam_site,k,k,j,k)+...
                   1i*(g_kkkk(lambda_tot,j,j,k,j) - g_kkkk(lambda_tot,k,k,j,k)) ;               
                          
        R_red_markov(k,j) = -2*real(tmp1*tmp2/(Gamma_fc+1i*deltaom));        
            
        end
    end
%R_{kkkk} = -sum_{k'} R_{kk,k'k'}    
 R_red_markov(k,k) = -sum(R_red_markov(k,:));   
end

%% solve the ODE d/dt D_k(t) = -sum(k,k') R_{k,k'} D_{k'}(t)

    de_fn = @(t,v) -mtimesx(R_red_markov,v);       
    %t_sep_rng is range of sep of pulses
    D_t_markov = zeros(length(t_sep_rng),N,length(om_u_rng));
    %ode45 only takes column so can't do left acting directly
  for j = 1:length(om_u_rng)  
    output_DE_fun(length(t_sep_rng),D_markov(:,j),'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-10);
    ode45(de_fn ,[0,t_sep_rng(end)],D_markov(:,j),options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_rng),D_markov(:,j),'get_data+clean');        
    tmp3 = interp1(tmp1,tmp2,t_sep_rng);
    D_t_markov(:,:,j) = tmp3;
  end
  
end
%these functions can be used to calculate exciton participation averages of
%either g (the lineshape broadening) or lambda (the reorganisation energy)
function  out = g_kkkk(g,k1,k2,k3,k4) %single excited exciton states

persistent CC

if nargin == 1
    CC = g; out=[];return
end
out =(CC(k1,:).*CC(k2,:).*CC(k3,:).*CC(k4,:))*g;
end
function out = g_kkqq(g,k,q) %double + single excited exciton states

persistent CC CCq occ N

if nargin == 1
    CC = g{1}; CCq = g{2};  occ = g{3}; N = size(occ,2); out=[];return
end
    out = zeros(size(g(1,:)));
for qq = 1:length(CCq)-1
    out = out + CCq(1+q,1+qq)^2.*(CC(k,occ(N+1 + qq,:)).^2*g(occ(N+1+qq,:),:));
end
end
function out = g_qqqq(g,q) %just double exciton states

 persistent CCq occ

if nargin == 1
    CCq = g{1}; occ = g{2};
    double_ex_lg = sum(double(occ),2)==2;
    occ = occ(double_ex_lg,:); %only need double exciton
    out=[];
    return
end
    
    out = zeros(size(g(1,:)));
for qq = 1:size(occ,1)
    occ1 = occ(qq,:); occ1 = repmat(occ1,size(occ,1));
    occ2 = occ & occ1; %coincidences in occupation
    rng = 1:size(occ,1); rng = rng(any(occ2,2));
    for qqlp = length(rng)
        qq2 = rng(qqlp);
    out = out + CCq(1+q,1+qq)^2.*(CCq(1+q,1+occ2(qq2,:)).^2*g(occ2(qq2,:),:));
    end
end
end
