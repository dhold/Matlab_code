function  [D_full,D_t_full,W_gg_full,W_kk_full,W_qq_full,...
        D_markov,W_gg_markov,W_kk_markov,W_qq_markov,R_red,D_markov_test,fock_rep] = doorway_fun_simple(...
          H_single,H_double,tau,om_u_rng,om_r_rng,t_sep_rng,params) 
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
E_u = @(t) exp(-t.^2/2/tau_u^2)./sqrt(sqrt(pi)*tau_u);
E_u_w = @(w) exp(-w.^2*tau_u^2/2).*sqrt(tau_u)/sqrt(sqrt(pi));
E_r = @(t) exp(-t.^2/2/tau_r^2)./sqrt(sqrt(pi)*tau_r);
E_r_w = @(w) exp(-w.^2*tau_r^2/2).*sqrt(tau_r/sqrt(pi));
%Calculate line broadening functions for each site
g_markov = zeros(N,1);  Gam_site = g_markov; lambda_site = Gam_site;
for j = 1:N %site loop  
    g_markov(j) = line_broad_fn_markov(B,gam_rng{j},om0{j},lam_rng{j},gam_dru{j},lam_dru{j});
    Gam_site(j) = real(g_markov(j)); %damping due to bath coupling
  %  lambda_site(j) = -imag(g_markov(j)); %frequency shift due to bath
end %ft of this is just 1/(g-i om) = (g+iom)/(g^2+om^2)
%lambda_site ==lambda_tot %just to test these are indeed, the same

%take om_range to be -5/tau to 5/tau, spaced in units relivant to both
om_space = min((1e-4)*mean(Gam_site),(1e-3)/max(tau));
om_max = min(20*mean(Gam_site),5/max(tau));
om_int_range = -om_max:om_space:om_max;
%choose a timerange that maps to this one
t_range = 0:pi/om_max:2*pi/om_space;
g_full = zeros(N,length(t_range)); %not Markovian approx

for j=1:N
    g_full(j,:) = line_broad_fn_full(B,gam_rng{j},om0{j},lam_rng{j},gam_dru{j},lam_dru{j},t_range);
    %ft not analytic due to horrible things.  By "full" I am still neglecting
    %some non markovian parts of the Matsubara terms
end

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

%calculate in frequency space
E_u_prefct = E_u_w(om_int_range).*E_u_w(-om_int_range); %conjugate for general
E_r_prefct = E_r_w(om_int_range).*E_r_w(-om_int_range); %conjugate for general

D_markov = zeros(N,length(om_u_rng)); D_full = D_markov;
W_gg_markov = zeros(N,length(om_r_rng));W_gg_full = W_gg_markov;
g_mark_ex = zeros(N,1); g_ex = zeros(length(t_range),N);
for j=1:N
    
        g_mark_ex(j) = g_kkkk(g_markov,j,j,j,j);
        g_ex(:,j) = g_kkkk(g_full,j,j,j,j);
        for lp = 1:length(om_u_rng)
            om_u = om_u_rng(lp);
        deltaom = (H_single_ex(j,j)-om_u); %diff between energy and carrier wavelength
               
D_markov_test(j,lp) = trapz(om_int_range,E_u_prefct./...
                        (g_mark_ex(j) + 1i*(deltaom-om_int_range)));
D_markov_test(j,lp) = D_markov_test(j,lp)+conj(D_markov_test(j,lp));
                    
D_markov(j,lp) =  4*sqrt(pi)*real(Faddeeva_w(tau_r*(deltaom + 1i*g_mark_ex(j))));            
                    
  tmp = exp(-1i*deltaom*t_range-g_ex(:,j).');
  tmp2 = fft(tmp)/length(tmp); tmp3 = fftshift(tmp2); %shift ft to have zero freq centre
                  
D_full(j,lp) = 2*real(trapz(om_int_range,E_u_prefct.*tmp3));                    
        end
        for lp = 1:length(om_r_rng)
            om_r = om_r_rng(lp);
        deltaom = (H_single_ex(j,j)-om_r); %diff between energy and carrier wavelength
        deltaom_2 = deltaom + imag(g_mark_ex(j)); %shift from bath       
W_gg_markov_test(j,lp) = 2*trapz(om_int_range,E_u_prefct.*(real(g_mark_ex(j)))./...
                        (real(g_mark_ex(j)).^2 + (deltaom_2+om_int_range).^2));
                    % - 1i*(deltaom+om_int_range) term vanished when cc taken
       %this can be evaluated analytically to be a voigt profile             
W_gg_markov(j,lp) = 4*sqrt(pi)*real(Faddeeva_w(tau_r*(-deltaom_2 + 1i*real(g_mark_ex(j)))));
                    
  tmp = exp(-1i*deltaom*t_range-g_ex(:,j).');
  tmp2 = ifft(tmp)/length(tmp); tmp3 = fftshift(tmp2); 
  %ifft because I want fft at minus omega
                    
W_gg_full(j,lp) = 2*real(trapz(om_int_range,E_r_prefct.*tmp3));               
            
        end
end
%hole in the ground state is just - the sum of this, should be down
%dimension 1 so this is fine
%D_hole = -sum(D_full); D_hole_mark = -sum(D_markov);

%now calculate window functions for the stim-emission and ES-abs conts
W_kk_markov = zeros(N,length(om_r_rng)); W_kk_full = W_kk_markov;
W_qq_markov = zeros(N,N2,length(om_r_rng)); W_qq_full = W_qq_markov;

%calculate linebroadening fns (g-matricies) in exciton basis
lambda_ex = zeros(N,1); lambda_ex_kq = zeros(N,N2); 
g_ex_kq_mark = zeros(N,N2);  g_ex_kq = zeros(length(t_range),N,N2);
g_ex_qq_mark = zeros(N2,1); g_ex_qq = zeros(length(t_range),N2);
for j = 1:N
    lambda_ex(j) = g_kkkk(-imag(g_markov),j,j,j,j); 
        %this is also the imaginary component of derive of g_full as t->inf
    for q = 1:N2
        if j==1
        g_ex_qq_mark(q) = g_qqqq(g_markov,q);
        g_ex_qq(:,q) = g_qqqq(g_full,q);
        end
        g_ex_kq_mark(j,q) = g_kkqq(g_markov,j,q);
        g_ex_kq(:,j,q) = g_kkqq(g_full,j,q);
        lambda_ex_kq(j,q) = g_kkqq(-imag(g_markov),j,q);
    end
end

for lp = 1:length(om_r_rng)
    om_u = om_r_rng(lp);
    for j = 1:N
     deltaom = (H_single_ex(j,j)-om_u); %diff between energy and carrier wavelength  
     deltaom2 = (H_single_ex(j,j)-om_u)-2*lambda_ex(j)-imag(g_mark_ex(j));  %shifted
     %this is the same as + imag(g_mark_ex(j))
   W_kk_markov(j,lp)  = 4*sqrt(pi)*real(Faddeeva_w(tau_r*(-deltaom2 +...
                       -1i*real(g_mark_ex(j)))));
                                       
  tmp = exp(1i*(2*lambda_ex(j)-deltaom)*t_range-conj(g_ex(:,j).'));
  tmp2 = ifft(tmp)/length(tmp); tmp3 = fftshift(tmp2); %                   
     W_kk_full(j,lp) = 2*real(trapz(om_int_range,E_r_prefct.*tmp3));
     
        for q = 1:N2 %H_double also includes ground state
    deltaom_q = (H_double_ex(q+1,q+1)-H_single_ex(j,j)-om_u);          
    deltaom_q2 = deltaom_q -2*(lambda_ex_kq(j,q) - lambda_ex(j)) + ...
                    imag(g_mark_ex(j) +g_ex_qq_mark(q) -2*g_ex_kq_mark(j,q));
    
    W_qq_markov(j,q,lp) = 4*sqrt(pi)*real(Faddeeva_w(tau_r*(-deltaom_q2 +...
        1i*real(imag(g_mark_ex(j) +g_ex_qq_mark(q) -2*g_ex_kq_mark(j,q))))));
    
                    
  tmp = exp(-1i*(deltaom_q+2*lambda_ex_kq(j,q)-2*lambda_ex(j))*t_range...
            -g_ex(:,j).' -g_ex_qq(:,q).' +2*g_ex_kq(:,j,q).');
  tmp2 = ifft(tmp)/length(tmp); tmp3 = fftshift(tmp2); %shift ft to have zero freq centre

 W_qq_full(j,q,lp) = 2*real(trapz(om_int_range,E_r_prefct.*tmp3));   
 
        end
    end

end

%% Calculate evolution of the doorway wavepacket via modified redfield

R_red = zeros(N,N);

for k =1:N
    for j = 1:N
        if k ~=j
        tmp1 = g_kkkk(g_full,k,k,j,j); tmp2 = g_kkkk(lambda_tot,k,k,j,j);
        deltaom = (H_single_ex(j,j)-H_single_ex(k,k)); 
        
        tmp3 = exp(-1i*(deltaom+2*tmp2-2*lambda_ex(j))*t_range...
            -g_ex(:,j).' -g_ex(:,k).' +2*tmp1);  
        
        tmp4 = diff(g_kkkk(g_full,j,k,j,j)-g_kkkk(g_full,j,k,k,k));
        tmp4 = [tmp4,(3*tmp4(end)-tmp4(end-1))/2]+2i*g_kkkk(lambda_tot,j,k,j,j); 
        %interp last point and add on renorm E
        tmp5 = diff(g_kkkk(g_full,j,j,k,j)-g_kkkk(g_full,k,k,k,j));
        tmp5 = [tmp5,(3*tmp5(end)-tmp5(end-1))/2]+2i*g_kkkk(lambda_tot,j,j,k,j); 
        
        tmp6 = diff(g_kkkk(g_full,j,k,k,j),2);
        tmp6 = [(3*tmp6(1)-tmp6(2))/2,tmp6,(3*tmp6(end)-tmp6(end-1))/2];   
        %interp last 2 points
        %finally integrate to obtain the modified redfield tensor
        R_red(k,j) = -2*real(trapz(t_range,tmp3.*(tmp6-tmp4.*tmp5)));        
            
        end
    end
%R_{kkkk} = -sum_{k'} R_{kk,k'k'}    
 R_red(k,k) = -sum(R_red(k,:));   
end



R_red_markov = zeros(N,N);

for k =1:N
    for j = 1:N
        
        tmp1 = g_kkkk(g_markov,k,k,j,j); tmp2 = g_kkkk(lambda_tot,k,k,j,j);
        deltaom = (H_single_ex(j,j)-H_single_ex(k,k)); 
        tmp3 = exp(-1i*(deltaom+2*tmp2-2*lambda_ex_kk(j))*t_range...
            -g_ex(:,j) -g_ex_qq(:,k) +2*tmp1);  
        
        tmp4 = diff(g_kkkk(g_full,j,k,j,j)-g_kkkk(g_full,j,k,k,k));
        tmp4 = [tmp4,(3*tmp4(end)-tmp4(end-1))/2]+2i*g_kkkk(lambda_tot,j,k,j,j); 
        %interp last point and add on renorm E
        tmp5 = diff(g_kkkk(g_full,j,j,k,j)-g_kkkk(g_full,k,k,k,j));
        tmp5 = [tmp5,(3*tmp5(end)-tmp5(end-1))/2]+2i*g_kkkk(lambda_tot,j,j,k,j); 
        
        tmp6 = diff(diff(g_kkkk(g_full,j,k,k,j)));
        tmp6 = [tmp6,(3*tmp6(end)-tmp6(end-1))/2];
        tmp6 = [tmp6,(3*tmp6(end)-tmp6(end-1))/2];       
        %interp last 2 points
        
        %finally integrate to obtain the modified redfield tensor
        R_red_markov(k,j) = -2*real(trapz(t_range,tmp3.*(tmp6-tmp4.*tm5)));
        
    end
end

%% solve the ODE d/dt D_k(t) = -sum(k,k') R_{k,k'} D_{k'}(t)

    de_fn = @(t,v) -mtimesx(R_red,v);       
    %t_sep_rng is range of sep of pulses
    D_t_full = zeros(length(t_sep_rng),N,length(om_u_rng));
    %ode45 only takes column so can't do left acting directly
  for j = 1:length(om_u_rng)  
    output_DE_fun(length(t_sep_rng),D_full(:,j),'notsavingnayway'); 
    options = odeset('OutputFcn',@output_DE_fun,'AbsTol',1e-10);
    ode45(de_fn ,[0,t_sep_rng(end)],D_full(:,j),options);
    
    [tmp1,tmp2]  =  output_DE_fun(length(t_sep_rng),D_full(:,j),'get_data+clean');        
    tmp3 = interp1(tmp1,tmp2,t_sep_rng);
    D_t_full(:,:,j) = tmp3;
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
