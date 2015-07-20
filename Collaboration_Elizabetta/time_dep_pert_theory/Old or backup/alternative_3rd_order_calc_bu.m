% precompute dimer rot averages note R12 points along z
Rsp = norm(R12);
%|kr| and |ku| are just 2 pi/ lambda, or k = omega/c -> 2*pi*omega in cm^-1 
%kr = 2*pi*om_r*kpr; %third beam to interact
kr = 2*pi*om3*kpr;
ku = 2*pi*om_u*kpu; %first/second beam to interact
%this approximates the wavevector by the average wavevector
%in reality each frequency is associated with a given wavevector
%could have k1 = 2*pi*om_1*kpu;  k2 = 2*pi*om_2*kpu; 
%k_rp = -k1+k2+k3; k_nr = k1-k2+k3;

av_rp = sym(zeros(3,2,2,2,2)); %first interaction is conjugated
%av_nr = sym(zeros(3,2,2,2,2)); %second interaction is conjugated
av_rp_FO = sym(zeros(3,2,2,2,2));  av_nr_FO = sym(zeros(3,2,2,2,2));
%first order terms in k1,k2,k3

for lp = 1:3 %dimension
for ell = 1:2
    for e3 = 1:2
        for e2 = 1:2
            for e1 = 1:2
deltaRmat = -[R(ell,:)- R(e1,:);R(ell,:)- R(e2,:);R(ell,:)- R(e3,:)].';          
%mumat = [mu(ell,:);mu(e3,:);mu(e2,:);mu(e1,:)].';  
mumat = [mu(e1,:);mu(e2,:);mu(e3,:);mu(ell,:)].';  

av_rp(lp,ell,e3,e2,e1) = tensor_average(mumat,1,lp);  
for lp2=1:3
    
av_rp_FO(lp,ell,e3,e2,e1) = av_rp_FO(lp,ell,e3,e2,e1)+...
                    tensor_average([mumat,deltaRmat(:,1)],lp,lp2)*(-ku(lp2))...
                     + tensor_average([mumat,deltaRmat(:,2)],lp,lp2)*ku(lp2)...
                     + tensor_average([mumat,deltaRmat(:,3)],lp,lp2)*kr(lp2); 
                 
av_nr_FO(lp,ell,e3,e2,e1) = av_nr_FO(lp,ell,e3,e2,e1)+...
                    tensor_average([mumat,deltaRmat(:,1)],lp,lp2)*ku(lp2)...
                     + tensor_average([mumat,deltaRmat(:,2)],lp,lp2)*(-ku(lp2))...
                     + tensor_average([mumat,deltaRmat(:,3)],lp,lp2)*kr(lp2);                 
                 
end                
            end
        end
    end
end
end

av_nr  = av_rp; %to zeroth order these are the same

%% Calculate all contributions via redfield / heom in SITE basis

% 4 different F-man diagrams contribute
% R / L for acting to the left (usual) or right (commutator term)
% 2Stim emmission
% mu*_e1_R mu_e2_L mu*_e2_L mu_e1_R %rephasing
% mu*_e1_R mu_e2_L mu_e1_R mu*_e2_L %nonrephasing
%
% 2 GS bleaching
% mu*_e2_R mu_e2_R mu_e1_L mu^*_e1_L %nonrephasing
% mu*_e2_R mu_e2_R mu*_e1_R mu_e1_R %rephasing
S_nr = sym([0;0;0]); S_rp = S_nr;
S_rp_FO = sym([0;0;0]); S_nr_FO  = S_rp_FO ;
rho00 = reshape(rho0,3,3);
clear mu_op
for e1 = 1:2
    temp = zeros(3);
    temp(1,e1+1) = 1; 
    mu_op_1{e1} = temp; %minus op
    mu_op_2{e1} = temp'; %plus operator 
end


for ell = 1:2
    for e3 = 1:2
        for e2 = 1:2
            for e1 = 1:2
     
 tmp_rp = ex_basis*reshape(red_propop1*reshape(ex_basis'*(rho00*mu_op_1{e1})*ex_basis,9,1),3,3)*ex_basis';
 
 tmp_rp1 = ex_basis*reshape(red_propop2*reshape(ex_basis'*(mu_op_2{e2}*tmp_rp)*ex_basis,9,1),3,3)*ex_basis';
 tmp_rp1 = ex_basis*reshape(red_propop3*reshape(ex_basis'*(tmp_rp1*mu_op_2{e3})*ex_basis,9,1),3,3)*ex_basis';
 %tmp_rp1 = trace(mu_op_1{ell}*tmp_rp1);
 
 tmp_rp2 = ex_basis*reshape(red_propop2*reshape(ex_basis'*(tmp_rp*mu_op_2{e2})*ex_basis,9,1),3,3)*ex_basis';
 tmp_rp2 = ex_basis*reshape(red_propop3*reshape(ex_basis'*(mu_op_2{e3}*tmp_rp2)*ex_basis,9,1),3,3)*ex_basis';
 %tmp_rp2 = trace(mu_op_1{ell}*tmp_rp2); 
 
 tmp = trace(mu_op_1{ell}*(tmp_rp1+tmp_rp2));
 
S_rp = S_rp + tmp*av_rp(:,ell,e3,e2,e1); 
S_rp_FO = S_rp_FO + tmp*av_rp_FO(:,ell,e3,e2,e1); 

 tmp_rp = ex_basis*reshape(red_propop1*reshape(ex_basis'*(mu_op_2{e1}*rho00)*ex_basis,9,1),3,3)*ex_basis';
 
 tmp_rp1 = ex_basis*reshape(red_propop2*reshape(ex_basis'*(tmp_rp*mu_op_1{e2})*ex_basis,9,1),3,3)*ex_basis';
 tmp_rp1 = ex_basis*reshape(red_propop3*reshape(ex_basis'*(tmp_rp1*mu_op_2{e3})*ex_basis,9,1),3,3)*ex_basis';
 %tmp_rp1 = trace(mu_op_1{ell}*tmp_rp1);
 
 tmp_rp2 = ex_basis*reshape(red_propop2*reshape(ex_basis'*(mu_op_1{e2}*tmp_rp)*ex_basis,9,1),3,3)*ex_basis';
 tmp_rp2 = ex_basis*reshape(red_propop3*reshape(ex_basis'*(mu_op_2{e3}*tmp_rp2)*ex_basis,9,1),3,3)*ex_basis';
 %tmp_rp2 = trace(mu_op_1{ell}*tmp_rp2); 
 
 tmp = trace(mu_op_1{ell}*(tmp_rp1+tmp_rp2));
 
S_nr = S_nr + tmp*av_rp(:,ell,e3,e2,e1); 
S_nr_FO = S_nr_FO + tmp*av_rp_FO(:,ell,e3,e2,e1); 
       
            end
        end
    end
end
% S_nr should just point along x and S_nr_FO along y
%% Take fourier transforms and such like

S_rp_ft = sym([0;0]); S_nr_ft =sym([0;0]);
 
   S_rp_ft(1) = fourier(subs(S_rp(1),om_r,om3)*heaviside(t1),t1,om1);
   S_rp_ft(1) = fourier(S_rp_ft(1)*heaviside(t2),t2,om1+om2);
   S_rp_ft(1) = fourier(S_rp_ft(1)*heaviside(t3),t3,om1+om2+om3);
%could just assume the t3 dep is dirac delta at t3-t_sep and reduce to subs
   
   S_rp_ft(2) = fourier(subs(S_rp_FO(2),om_r,om3)*heaviside(t1),t1,om1);
   S_rp_ft(2) = fourier(S_rp_ft(2)*heaviside(t2),t2,om1+om2);
   S_rp_ft(2) = fourier(S_rp_ft(2)*heaviside(t3),t3,om1+om2+om3);
   
   S_nr_ft(1) = fourier(subs(S_nr(1),om_r,om3)*heaviside(t1),t1,om1);
   S_nr_ft(1) = fourier(S_nr_ft(1)*heaviside(t2),t2,om1+om2);
   S_nr_ft(1) = fourier(S_nr_ft(1)*heaviside(t3),t3,om1+om2+om3);   
   
   S_nr_ft(2) = fourier(subs(S_nr_FO(2),om_r,om3)*heaviside(t1),t1,om1);
   S_nr_ft(2) = fourier(S_nr_ft(2)*heaviside(t2),t2,om1+om2);
   S_nr_ft(2) = fourier(S_nr_ft(2)*heaviside(t3),t3,om1+om2+om3);
   
   %% Simplify expressions
   
 %  S3_rp_ftx = simplify(S_rp_ft(1));
 %   S3_nr_ftx = simplify(S_rp_ft(1));
 
   %%  Compute signal assuming short probe pulse (many frequencies)
   syms om real
 
    S3_rp_ftx = subs(S_rp_ft(1),{om3},{om-om2-om1});
    S3_nr_ftx = subs(S_nr_ft(1),{om3},{om-om2-om1});
    S3_rp_fty = subs(S_rp_ft(2),{om3},{om-om2-om1});
    S3_nr_fty = subs(S_nr_ft(2),{om3},{om-om2-om1});
    tdelay_fc = exp(-1i*t_sep*(om-om2-om1));
    
   om_u_range = linspace(0.6e3,2e3,11);
   %range of pump wavelengths
   numpoints = 300;
   om_plot_rng = linspace(0.6e3,2e3,numpoints);
   tsep_range = linspace(0.05,10)*convfact;
   quad_order= 1;
   tau = tau_r;
  [omm, weig] = GaussHermite(quad_order,eps(20));
  omm = omm*sqrt(2)/tau; %rescale to appropriate size
  %[omm1,omm2] = meshgrid(omm,omm); [weig1,weig2] = meshgrid(weig,weig);
  prefact = 1/sqrt(pi)/tau; 
  %total prefct after normalisation factor of Gaussian wave packets and
  %rescale
  Ix = zeros(length(om_plot_rng),length(om_u_range),length(tsep_range));Iy=Ix;
  
  sym_to_fn('temp_fn_rp_x.m',vpa(S3_rp_ftx*tdelay_fc,12),[om,om1,om2,om_u,t_sep])
  sym_to_fn('temp_fn_nr_x.m',vpa(S3_nr_ftx*tdelay_fc,12),[om,om1,om2,om_u,t_sep])
  sym_to_fn('temp_fn_rp_y.m',vpa(S3_rp_fty*tdelay_fc,12),[om,om1,om2,om_u,t_sep])
  sym_to_fn('temp_fn_nr_y.m',vpa(S3_nr_fty*tdelay_fc,12),[om,om1,om2,om_u,t_sep])

  for lp = 1:length(om_u_range)
      omu = om_u_range(lp);
      for lp1 = 1:length(tsep_range)
          tsp = tsep_range(lp1);
          temp1 = 0; temp2 = 0;
      for lp2 = 1:quad_order
          for lp3 = 1:quad_order
              
          temp1 = temp1 + weig(lp2)*weig(lp3)*(temp_fn_rp_x(...
             om_plot_rng,omm(lp2)-omu,omm(lp3)+omu,omu,tsp) +...
             temp_fn_nr_x(om_plot_rng,omm(lp2)+omu,omm(lp3)-omu,omu,tsp)).';
          
          temp2 = temp2 + weig(lp2)*weig(lp3)*(temp_fn_rp_y(...
             om_plot_rng,omm(lp2)-omu,omm(lp3)+omu,omu,tsp) +...
             temp_fn_nr_y(om_plot_rng,omm(lp2)+omu,omm(lp3)-omu,omu,tsp)).';          
 
%           temp1 = temp1 + weig(lp2)*weig(lp3)*(temp_fn_rp_x(...
%              om_plot_rng,omm(lp2)+omu,omm(lp3)-omu,omu,tsp) +...
%              temp_fn_nr_x(om_plot_rng,omm(lp2)-omu,omm(lp3)+omu,omu,tsp)).';
%           
%           temp2 = temp2 + weig(lp2)*weig(lp3)*(temp_fn_rp_y(...
%              om_plot_rng,omm(lp2)+omu,omm(lp3)-omu,omu,tsp) +...
%              temp_fn_nr_y(om_plot_rng,omm(lp2)-omu,omm(lp3)+omu,omu,tsp)).';              
          
          end
      end
      Ix(:,lp,lp1) = temp1.*(om_plot_rng.'); %things I plot have this prefct
      Iy(:,lp,lp1) = temp2 .*(om_plot_rng.');     
      end
  end
   Ix = Ix*prefact;     Iy = Iy*prefact; 
   
   Dalpha_J  = -(8*pi^2)*real(Ix); %P(omega) = i S(omega)
Deta_J  = (16*pi^2)*real(Iy); 
Ddelta_J  = -(8*pi^2)*imag(Iy); 

%these third order contributions are the shifts in these values from the
%pump beam



%% plot parameters at given ranges 
lam_plot_rng = 10^7./om_plot_rng ;

figure
plot(lam_plot_rng,Dalpha_J(:,floor(end/2),floor(end/2)),'LineWidth',2)
xlabel('Wavelength, \lambda (nm)');
ylabel('Absorption change, \Delta \alpha (a.u.)');

figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'YColor',[0 0 1]);
box(axes1,'on');
hold(axes1,'all');
plot(lam_plot_rng,Deta_J(:,floor(end/2),floor(end/2)),'Parent',axes1,'LineWidth',2);
xlabel('Wavelength, \lambda  (nm)');
ylabel('Circular dichorism (relative to absolute absorption)','Color',[0 0 1]);
axes2 = axes('Parent',figure1,'YAxisLocation','right','YColor',[0 0.5 0],...
    'ColorOrder',[0 0.5 0;1 0 0;0 0.75 0.75;0.75 0 0.75;0.75 0.75 0;0.25 0.25 0.25;0 0 1],...
    'Color','none');
hold(axes2,'all');
plot(lam_plot_rng,Ddelta_J(:,floor(end/2),floor(end/2)),'Parent',axes2,'LineWidth',2);
ylabel('Optical rotation (units of wavevector)','VerticalAlignment','cap',...
    'Color',[0 0.5 0]);

 %% Compute signal
 
 % assuming tau (pulse width) is the same for pump & probe
 % P(omega) =  %pi^(-3/2)*tau^(-1/2)*exp(-tau^2/2 *(om^2+om_r^2+2 om_u^2-2*om*om_r)  ) 
 %           * int d om_1 int d om_2 exp(-(om_1^2 + om_2^2))*...           
 %           exp(i*tsep/4 *(4*om_1/tau+2*om_2/tau +3*om_r + om  +/-  om_u))...
 %          S3_ft^{+/-} ( (om +/- om_u)/4+3*om_r/4 + (om_1+om_2/2)/tau,
 %          om_2/tau - (om_r-om -/+ om_u)/2, (om_1 - om_2/2)/tau +
 %          (om_r-om+/-om_u)/4 )
 % rf is + in +/-


 om_r_range = linspace(0.6e3,2e3,7); om_u_range = [1,1.2,1.4]*1e3;
 if tau_r == tau_u
     tau = tau_r;
 else
     warning('you wrote this code to work assuming tau_r == tau_u, taking average for tau_u = tau but not strictly correct')
    tau = tau_u;
 end
 prefact = exp(-tau^2/2 *(om^2+om_r^2+2*om_u^2-2*om*om_r));
 
  [omm, weig] = GaussHermite(11,eps(20));
  [omm1,omm2] = meshgrid(omm,omm); [weig1,weig2] = meshgrid(weig,weig);
 

for lp1 = 1:length( om_r_range)
    for lp2 = 1:length( om_u_range)
        
        
    end
end