%% 
Temperature=300;
[convfact, B,speed_unit]= inv_cm_unit_sys(Temperature);
lambda2 = 0*lambda{1}; n=1;
%g_markov   = line_broad_fn_markov(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n});
g_markov  = line_broad_fn_markov(B,35,100,0,600,100);
Gamma = real(g_markov); lam_tot = -imag(g_markov);
g_mat = zeros(N,N,N,N); lam_mat  = zeros(N,N,N,N); Gam_mat = lam_mat;

%PP = M_prj(2:N+1,2:N+1); c_nk = PP; ex_E = diag(H_exciton(2:N+1,2:N+1));

VV = 50; E1 = 0; E2 = 600;
Htest = [E1,VV;VV,E2]; theta_val = 1/2*atan(2*VV/(E1-E2));
Proj = [cos(theta_val),sin(theta_val);-sin(theta_val),cos(theta_val)];
f1 = 2*Gamma*(cos(theta_val)^2-sin(theta_val)^2)^2;
f2 = 2*cos(theta_val)^2*sin(theta_val)^2;

[c_nk,ex_E] =eig(Htest);  ex_E = diag(ex_E); N=2;
%Gamma
%    lam_tot
pred_val =2*Gamma*f1^2*f2/((ex_E(2)-ex_E(1))^2+f1^2);

for k =1:N
    for k2 = 1:N
        for k3= 1:N
            for k4=1:N
                %fct(k,k2,k3,k4) = sum(c_nk(:,k).*c_nk(:,k2).*c_nk(:,k3).*c_nk(:,k4));
                tmp = sum(c_nk(:,k).*c_nk(:,k2).*c_nk(:,k3).*c_nk(:,k4));
                g_mat(k,k2,k3,k4)  = tmp*g_markov;
                lam_mat(k,k2,k3,k4)   = tmp*lam_tot;
                 Gam_mat(k,k2,k3,k4)  =  tmp*Gamma ;
            end
        end
    end
end

R_markov = zeros(N); R_markov2 = zeros(N);
for k = 1:N
    for kk = 1:N %k'
        if k~=kk
    dom = ex_E(k)-ex_E(kk);
     exp_fct_imag = -1i*dom + 1i*(lam_mat(k,k,k,k) - lam_mat(kk,kk,kk,kk)) ;
     exp_fct_real =  -Gam_mat(k,k,k,k) - Gam_mat(kk,kk,kk,kk) + 2*Gam_mat(kk,kk,k,k);
            
        
       tmp =  Gam_mat(kk,k,kk,kk) -  Gam_mat(kk,k,k,k) ...
             +1i*(lam_mat(kk,k,kk,kk) + lam_mat(kk,k,k,k) );
       
     tmp = tmp^2/(exp_fct_imag+exp_fct_real); %int is analytic
      R_markov(k,kk) = -2*real(tmp);     
       
            
       tmp2 = g_mat(kk,k,kk,kk) -  g_mat(kk,k,k,k) + 2i*lam_mat(kk,k,kk,kk);
       
       tmp3 = g_mat(kk,kk,k,kk) -  g_mat(k,k,k,kk) + 2i*lam_mat(kk,kk,k,kk);

       
       exp_arg = -1i*dom + 2i*(lam_mat(kk,kk,k,k)-lam_mat(kk,kk,kk,kk)) ...
          -g_mat(k,k,k,k) - g_mat(kk,kk,kk,kk) + 2*g_mat(kk,kk,k,k);
        R_markov2(k,kk) = -2*real( tmp3*tmp2/(exp_arg));

        end
    end
end
    R_markov = R_markov  - diag(sum(R_markov,1));
        R_markov2  = R_markov2  - diag(sum(R_markov2,1));
%%
% 
% g_broad   = line_broad_fn_full(B,[],[],[],600,100,t_int_rng,nmax);   
% d_g_broad = diff(g_broad)./diff(t_int_rng);
% g_deriv  = line_broad_fn_deriv(B,[],[],[],600,100,t_int_rng,nmax);
% d_g_deriv = diff(g_deriv)./diff(t_int_rng);
% g_sec_der = line_broad_fn_sec_der(B,[],[],[],600,100,t_int_rng,nmax);
tol  = 1e-8;
n=8;
t_int_rng = [linspace(Gamma_min*1e-8,1e-3*Gamma_min,n*1000),...
             linspace(1e-3*Gamma_min,2*Gamma_min,n*10000),...
            linspace(2*Gamma_min,20*Gamma_min,n*10000),...
            linspace(20*Gamma_min,150*Gamma_min,n*10000)];
 t_int_rng = t_int_rng([true,diff(t_int_rng)~=0]);

tic
 [g_broad,g_mat,n_req] = line_broad_fn_full(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng,tol); 
 toc
 % g_broad_1 = diff(g_broad)./diff(t_int_rng);
 %   g_broad_2 = diff(g_broad_1)./diff(t_int_rng(2:end));
 tic
 [g_deriv,gmat_1]  = line_broad_fn_deriv(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng,tol);
 toc
 
 tic
 [g_sec_der,gmat_2] = line_broad_fn_sec_der(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng,tol);
 toc
 
 g_broad = g_broad + g_mat;   g_deriv = g_deriv+gmat_1;  g_sec_der = g_sec_der+gmat_2;
 % g_broad  =  g_broad /pi; g_deriv = g_deriv/pi;  g_sec_der = g_sec_der/pi;
 %%
N=8;
PP = M_prj(2:N+1,2:N+1); ex_E = diag(H_exciton(2:N+1,2:N+1));


c_nk = PP'; %as I have defined it c_nk = PP(k,n)
g_markov  = line_broad_fn_markov(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n});


% crap_to_pass_to_fn ={g_broad,g_deriv,g_sec_der,t_int_rng,beta,...
%                     gam{1},om_0{1},lambda{1},gam_dru{1},lam_dru{1} };
%  mod_redfield_calc4(crap_to_pass_to_fn); clear crap_to_pass_to_fn
%             tic
% R_mod_red_test =  mod_redfield_calc4(c_nk,ex_E,1e-4);
%              toc

% figure
% plot(t_int_rng,g_broad )
% hold on
% plot(t_int_rng,g_markov.*t_int_rng,'r')
% figure
% plot(t_int_rng,imag(g_broad))
% hold on
% plot(t_int_rng,imag(g_markov).*t_int_rng,'r')

%R_mod_red2 =  mod_redfield_calc2(g_markov.*t_int_rng/pi,g_markov.*t_int_rng.^0/pi ...
%            ,g_sec_der*0,lam_tot, PP ,ex_E,t_int_rng);
                 % R_mod_red  = R_mod_red  - diag(sum(R_mod_red,1));
%R_mod_red =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,...
%                      c_nk ,ex_E,t_int_rng);                  
R_mod_red =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,B,...
                     c_nk ,ex_E,t_int_rng,1e-2);       
%           tic
% R_mod_red =  mod_redfield_calc3(B,gam{n},om_0{n},lambda{n},gam_dru{n},...
%                  lam_dru{n},c_nk',ex_E,1e-6,1e-2);
%              toc
       %%
       lg = R_mod_red>0;lg2 = logical(eye(size(R_mod_red)));
       if any(any(lg &~lg2))
           warning('positive elements in Redfield of diag removed')
       R_mod_red(R_mod_red>0) = 0;
    R_mod_red  = R_mod_red  - diag(sum(R_mod_red,1)); %balance condition     
       end
de_ee = @(t,v) -R_mod_red*v;
%de_ee = @(t,v) -R_reduced*v;
%de_ee = @(t,v) -saved_val*v;

init_state = zeros(N);
init_state(6,6) = 1;
init_state = PP'*init_state*PP;
init_state = diag(init_state);
[a,b]=ode45(de_ee,[0,10*convfact],init_state); 

% figure
% plot(a/convfact,b)
% hold on
% plot(a/convfact,sum(b,2),'--')

bb = zeros(8,8,length(a));
for k = 1:8
bb(k,k,:) = b(:,k);
end
bb = mtimesx(PP,bb);
bb = mtimesx(bb,PP,'C');
bbb = reshape(bb,64,length(a));

thermal_state = expm(-B*H_exciton(2:N+1,2:N+1));
thermal_state = thermal_state/trace(thermal_state);
thermal_state_site_basis = mtimesx(PP,thermal_state);
thermal_state_site_basis = mtimesx(thermal_state_site_basis,PP,'C');
thermal_state_site_basis  = reshape(thermal_state_site_basis,64,1);
thermal_state_site_basis  = repmat(thermal_state_site_basis,1,length(a));

tmp = eye(N); tmp = reshape(logical(tmp),N^2,1);

[~,alex_order]=sort(diag(H0(2:N+1,2:N+1)));
figure1 = figure;
axes1 = axes('Parent',figure1);
plot1 = plot(a/convfact,bbb(tmp,:),'Parent',axes1);

for k = 1:N
set(plot1(alex_order(k)),'LineWidth',2,'DisplayName',strcat('Site',num2str(k)));
end
% for k = 1:N
% set(plot1(alex_order(k)+N),'LineWidth',2,'DisplayName',strcat('Site',num2str(k),'Thermal pop'));
% end
legend(axes1,'show');
hold on
plot(a/convfact,thermal_state_site_basis(tmp,:),'Parent',axes1,'Linestyle','--')




%% This part was to specifically test my results vs Richard
%either give half width at half max or sd
sys =4; %choose  system
[H_site,mu,R,lambda,gam,om_0,lam_dru, gam_dru,om_vib,displ,pdm,sd_shift] =sys_chooser(sys);
N = length(H_site); sz1 = N+1+N*(N-1)/2;
use_markov = true; inc_double_ex = true; use_mod_red  =true;
site_dep_bath = false; %sites have same bath / diff baths

Hwhm_u_fs = 75; %half width at half maximum, pump
tau_u = (Hwhm_u_fs/sqrt(2*log(2)))/1000*convfact; %pumpp pulse SD in inverse cm

points_to_save = 10000;
t1_range = linspace(0,10*tau_u,points_to_save);
 n=1; %bath at first site same as all others
 g_cnst = line_broad_fn_markov(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n});
Gamma_min = 1/real(g_cnst);
t_int_rng = [linspace(0,B/(10*pi),2000),...
             linspace(B/(10*pi),2*Gamma_min,10000),...
            linspace(2*Gamma_min,20*Gamma_min,10000),...
            linspace(20*Gamma_min,150*Gamma_min,10000)];
 t_int_rng = t_int_rng([true,diff(t_int_rng)~=0]);
 

 tol = 1e-11;
 g_broad   = line_broad_fn_full(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);    
g_deriv  = line_broad_fn_deriv(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);
g_sec_der = line_broad_fn_sec_der(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t_int_rng, tol);
g_sec_der(1) = (3*g_sec_der(2)-g_sec_der(3))/2;
lam_tot = sum(lambda{n}) + sum(lam_dru{n});
g_t_door  = line_broad_fn_full(B,gam{n},om_0{n},lambda{n},gam_dru{n},lam_dru{n},t1_range, tol);
%%
num_realisations = 10; %choose decent number of realisations. 
site_shift = randn(num_realisations,N); %randn has variance 1, assumes indep
for j =  1:N
  site_shift(:,j) = site_shift(:,j)*sd_shift(j); %change variance
end
delta_w = sum(site_shift,2)/N; %mean shift to excited state manifold energy
site_shift_centred = site_shift - repmat(delta_w,1,N); 
use_HEOM = true;
if use_HEOM
    Kappa=1; QQ_topass = zeros(N,2);
    for j = 1:N
    % [cc1,cc2R,cc2I,vv1,vv2,QQ] = coeffients_from_brownian_new...
    %(lambda{j},gam{j},om_0{j},Temp,Kappa,lam_dru{j},gam_dru{j});
[cc1,cc2R,cc2I,vv1,vv2,QQ] = coeffients_from_brownian_new...
    ([],[],[],Temperature,Kappa,lam_dru{j},gam_dru{j});    
  cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
 cc_acom{j}= [cc2I,cc1*0]; 
 
 QQ_topass(j,:) = [QQ,0];
    end
    kap2 = 2;
   % [ Qtrunc,H_prop_op,nn]=HEOM_propogator(QQ_topass,cc_com,cc_acom,vv,inf,kap2,1);
    [ Qtrunc,H_prop_op,nn,nn_alt]=HEOM_propogator2(QQ_topass,cc_com,cc_acom,vv,kap2,1);
end
%%
D0_save = zeros(num_realisations,N);
D0_save_ex = D0_save;  D0_save_site = zeros(num_realisations,N,N);
t_sep_rng_ps = linspace(0,30,1000); t_sep_rng = t_sep_rng_ps*convfact;
D_actual_save = zeros(N,length(t_sep_rng),num_realisations);
for real_lp=1:num_realisations
H_site_shifted = H_site + diag(site_shift(real_lp,:));
[Prj,H_ex] = eig(H_site_shifted); E1 = diag(H_ex);
[D_t,~,~,D0,~,~,~] = spec_fun_simple(t1_range,[0,0.5,1],...
          1.0e+04*1.9802,1.0e+04*1.9802,tau_u,tau_u,E1,ones(N*(N-1)/2),...
          [1,0.5,0.1],g_t_door,lam_tot,Prj',ones(N*(N-1)/2),N,false,false) ;

       c_nk = Prj';
  mu_ex = zeros(3,N); %transition dipoles between exciton states   
  tmp = zeros(N,1);
  for k = 1:N
    mu_ex(:,k) = mu.'*c_nk(:,k);   
    tmp(k) = dot(mu_ex(:,k),mu_ex(:,k))*D0(k)/3;
  end      

      D0_save_ex(real_lp,:) = tmp;
tmp = mtimesx(Prj,diag(tmp) ); %proj to site basis
tmp = mtimesx(tmp,Prj,'C');       
      D0_save_site(real_lp,:,:) = tmp;
  
      tic
if use_HEOM
   % for lp  =1:size(D0,1)
   
    rho_0 = squeeze(D0_save_site(real_lp,:,:)); 
    tr_save = trace(rho_0); rho_0 = rho_0 / tr_save;
  [Time_units,rho_vec]=HEOM_exmanifold_solver(H_site_shifted,...
            Qtrunc,nn,H_prop_op,rho_0,length(t_sep_rng),t_sep_rng(end),1);
       D_actual = interp1( Time_units,rho_vec,t_sep_rng);
      D_actual = reshape(D_actual.',N,N,size(D_actual,1));  
      D_actual = mtimesx(Prj,'C',D_actual);
      D_actual = mtimesx(D_actual,Prj);
     D_actual = reshape(D_actual,N^2,length(t_sep_rng));
   % end
else      
   R_mod_red =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,B,...
           c_nk,E1,t_int_rng,2e-3);
    t_for_red(real_lp) = toc;
    lg = R_mod_red>0 &~logical(eye(N));
    if any(any(lg))
        warning('positive elements found and set to zero')
        R_mod_red(lg) = 0;
    end
    
       ode_toler = 1e-12;
    de_fn = @(t,v) -R_mod_red*v;
     options = odeset('AbsTol', ode_toler);
     t_end_ps = 10; t_end = t_end_ps*convfact;
     
    % for jj = 1:N %loop over different exciton dipole transitions
   %      init_state = D0*0; init_state(jj) = D0(jj);
   
    [t_rng,D_tev] =  ode45(de_fn,[0,t_end],D0_save_ex(real_lp,:),options);   

       D_tev = interp1(t_rng,D_tev,t_sep_rng);

     D_actual = zeros(N,N,length(t_sep_rng));  
     for j=1:N %loop over final exciton occupations
D_actual(j,j,:) = D_tev(:,j) ;
     end
    % end
D_actual = mtimesx(Prj,D_actual ); %proj to site basis
D_actual = mtimesx(D_actual,Prj,'C');  

D_actual = reshape(D_actual,N^2,length(t_sep_rng_ps));
end
D_actual_save(:,:,real_lp) = D_actual(1:N+1:end,:);
end

%%
save('data_for_Richard.mat','site_shift','D0_save_ex','D0_save_site')
save('data_for_Richard.csv','site_shift','D0_save_ex','D0_save_site')

%%
test = sum(D_actual_save,3); %test=test/sum(test(:,1));
figure
plot(t_sep_rng_ps,test/num_realisations)

test2 = sqrt(sum(D_actual_save.^2-repmat...
            ((test/num_realisations).^2,1,1,num_realisations),3)); 
        legend show
        
hold on
plot(t_sep_rng_ps,test/num_realisations+3*test2/num_realisations^(3/2),'--')
plot(t_sep_rng_ps,test/num_realisations-3*test2/num_realisations^(3/2),'--')

%% test late just using normal states

    sys =4; %choose  system
[H_site,mu,R,lambda,gam,om_0,lam_dru, gam_dru,om_vib,displ,pdm,sd_shift] =sys_chooser(sys);
N = length(H_site); sz1 = N+1+N*(N-1)/2;
num_test = 10000; %choose decent number of realisations. 
site_shift = randn(num_test,N); %randn has variance 1, assumes indep
for j =  1:N
  site_shift(:,j) = site_shift(:,j)*sd_shift(j); %change variance
end
state_at_early_time = zeros(size(H_site));
state_at_late_time =state_at_early_time;

Temperature=77;
[convfact, B,speed_unit]= inv_cm_unit_sys(Temperature);

Hwhm_u_fs = 75; %half width at half maximum, pump
tau_u = (Hwhm_u_fs/sqrt(2*log(2)))/1000*convfact; %pumpp pulse SD in inverse cm
points_to_save = 50000;
t1_range = linspace(0,15*tau_u,points_to_save);

lam_tot = sum(lambda{1}) + sum(lam_dru{1});
g_t_door  = line_broad_fn_full(B,gam{1},om_0{1},lambda{1},gam_dru{1},lam_dru{1},t1_range, tol);

for some_lp = 1:num_test
 H_site_shifted = H_site + diag(site_shift(some_lp,:));
[Prj,H_ex] = eig(H_site_shifted); E1 = diag(H_ex);
%Prj*H_ex*Prj'=H_site_shifted
therm_state =diag(exp(-B*(E1-E1(1)))); 
therm_state = therm_state/trace(therm_state);
[D_t,~,~,D0,~,~,~] = spec_fun_simple(t1_range,[0,0.5,1],...
          1.0e+04*1.9802,1.0e+04*1.9802,tau_u,tau_u,E1,ones(N*(N-1)/2),...
          [1,0.5,0.1],g_t_door,lam_tot,Prj',ones(N*(N-1)/2),N,false,false) ;
      
  c_nk = Prj';
  mu_ex = zeros(3,N); %transition dipoles between exciton states   
  tmp = zeros(N,1);
  for k = 1:N
    mu_ex(:,k) = mu.'*c_nk(:,k);   
    tmp(k) = dot(mu_ex(:,k),mu_ex(:,k))*D0(k)/3;
  end           
  
state_at_late_time = state_at_late_time + sum(tmp)*Prj* therm_state *Prj';
state_at_early_time = state_at_early_time+ Prj* diag(tmp)*Prj';
end


% Create multiple lines using matrix input to bar
figure
plot(diag(state_at_early_time)/trace(state_at_early_time))
hold on
bar(diag(state_at_late_time)/trace(state_at_late_time),'r')
xlabel('Site','FontSize',14);
% Create ylabel
ylabel('Population','FontSize',14);
