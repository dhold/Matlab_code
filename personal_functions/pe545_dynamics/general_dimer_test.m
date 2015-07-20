%% Consider explicit vibrational mode
 V = 400; %Coupling electronic dipole
  coupg = 267.1; %coupling term between electronic and vibrational DOF
  inc_mode_in_ham = true;
H_el = [0,V;V,0];   omegavib = 1000; 
  if inc_mode_in_ham
  
%Frequency of vibrations / level spacing,
 viblvls = 5; %Number of vibrational levels to include, integer
 viblvls = 1 %remove mode by decommenting
 % %J_j(omega) = 2 lamda_j gamma_j /hbar * (omega/(omega^2+lambda_j^2)


Hvib = diag(0:viblvls-1)*omegavib; %vibrational Hamiltonian
%subtract zero point energy
%relative mode only, com mode trivially uncouples
Hexvib = diag(-coupg/sqrt(2) * sqrt(1:viblvls-1),1);
Hexvib = Hexvib+Hexvib.';
Hexvib = kron([-1,0;0,1],Hexvib);

 %sites have same energy

Htot = kron(H_el,eye(size(Hvib))) + kron(eye(size(H_el)),Hvib) + Hexvib;

[ projm,H_ex] = eig(H_el);
proj_full = kron(projm,eye(size(Hvib)));
H_ex_full = proj_full'*Htot *proj_full;
  else
      viblvls = 1;
   Htot = H_el;
   [ projm,H_ex] = eig(H_el);
   proj_full = projm; H_ex_full = H_ex;
      
  end

%% Set up dimer parameters 
Temp = 300; %temperature in Kelvin
[convfact, B,speed_unit]= inv_cm_unit_sys(Temp); N=2;
 
 gam_dru = [600;600]; %drude decay constant / cut off frequency
 lam_dru = [60;60]; % renorm energy 
   Kappa = 3; %truncation parameter satisfying Kappa >>omega_0 beta hbar/2pi
 %beyond this e^(-v_(k>kappa) t)*v_k ~delta(t) is made
 
 if ~inc_mode_in_ham
     
     lam_bro = omegavib/sqrt(coupg)*[1,1];
     gam_bro = 5*[1,1]; %damping otherwise not included
     om_0 =omegavib*[1,1];
 else
     lam_bro=[]; gam_bro = []; om_0=[];
 end
 
 QQ = zeros(length(gam_dru),2); cnt=0;
 cc = zeros(1,length(gam_dru) + Kappa);
 clear cc_com cc_acom vv nn
 for j=1:length(gam_dru)
   
      if ~inc_mode_in_ham %add mode to spec density
 [cc1,cc2R,cc2I,vv1,vv2,QQ(j,1),tmp] = coeffients_from_brownian_new(...
    lam_bro(j),gam_bro(j),om_0(j),Temp,Kappa,lam_dru(j),gam_dru(j));
      else %no underdamped mode in spec density
  [cc1,cc2R,cc2I,vv1,vv2,QQ(j,1),tmp] = coeffients_from_brownian_new(...
   [],[],[],Temp,Kappa,lam_dru(j),gam_dru(j));         
      end
 % QQ is in general length two for each point as I could (potentially)
 % included high frequency decay terms from poles in J(omega) which would
 % have anti commutator parts
 cc(cnt+1:cnt+length(tmp)) = tmp; cnt =cnt+length(tmp);
 cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
 cc_acom{j}= [cc2I,cc1*0];
 
 end

%%
aa  = sqrt(0.5); bb = sqrt(0.5); %coeffs of init wavefn
%populations in initial state are these mod squared
rho_0 = [aa;bb] *  [aa,bb];
rho_vib = diag(1./(1+exp(B*(1/2:viblvls-1/2)* omegavib)));
rho_vib  = rho_vib ./ trace(rho_vib);
rho_0 = kron(rho_0,rho_vib);
%project to site basis
rho_0 = proj_full'*rho_0*proj_full;

 tendps = 1; %end time in pico seconds 
 Kap2 = 4; %second truncation parameter, there are Kappa + 1 explicitly 
 %treated exponentials, this number is the maximum order the sum of the
 %number of powers can take, e.g. n_1 + n_2 +... + n_Kappa <= Kap2
 save_to_tier = 0; numpoints= 10000; use_red_mat = false;
 [Time_units,rho_vec,nn,total_prop_op]=multi_site_drude_and_brownian_HOM_3...
            (Htot,QQ,rho_0,cc_com,cc_acom,vv,inf,Kap2, viblvls,...
            numpoints,tendps, save_to_tier,use_red_mat);
        
 
        %% Perform partial trace over vibrational states

 if use_red_mat
             temp = reshape(tril(true(N*viblvls)) , 1,(N*viblvls)^2);
           rho_vec_full = zeros(length(Time_units),size(nn,1)*(N*viblvls)^2);
            temp2 = repmat(temp,1,size(nn,1)); 
       rho_vec_full(:,temp2) = rho_vec; 
       rho_vec = rho_vec_full; clear rho_vec_full temp temp2
     
 end
        
        if viblvls == 1 %no vibrational state to consider
            rho_red = reshape(rho_vec,size(rho_vec,1),...
                        sqrt(size(rho_vec,2)),sqrt(size(rho_vec,2)));
            rho_red = permute(rho_red,[3,2,1]);        
            rho_red =  mtimesx(proj_full,'C',rho_red); 
            rho_red =  mtimesx(rho_red,proj_full); 
  purity = real(diagsum(mtimesx(rho_red,rho_red),1,2));  
  
  entropy_vn = purity*0;
  for lp = 1:length(Time_units)
      
      entropy_vn(lp) = -trace(rho_red(:,:,lp)*logm(rho_red(:,:,lp))); 
  end
            rho_red = permute(rho_red,[3,2,1]);   
        
        else
    rho_red = zeros(2,2,length(Time_units));
    rho_00 = reshape(rho_vec,size(rho_vec,1),...
                        sqrt(size(rho_vec,2)),sqrt(size(rho_vec,2)));
    rho_00 = permute(rho_00,[3,2,1]);
    rho_00_ex =  mtimesx(proj_full,'C',rho_00);
    rho_00_ex =  mtimesx(rho_00_ex,proj_full);
     entropy_vn = zeros(size(Time_units));
  for lp = 1:length(Time_units)
   
rho_red(:,:,lp) = TrX(rho_00_ex(:,:,lp) ,2,[N,viblvls]);
entropy_vn(lp) = -trace(rho_red(:,:,lp)*logm(rho_red(:,:,lp)));    
  end
  rho_red = permute(rho_red,[3,2,1]);
  purity = real(diagsum(mtimesx(rho_00_ex,rho_00_ex),1,2));
        end
      trace_test = real(diagsum(rho_red,1,2));    
 entropy_vn(isnan(entropy_vn))=0; %0*log(0) = 0
      %% Fourier transform to see frequencies and stuff
   
      time_smooth = linspace(0,max(Time_units),numpoints);
  entropy_vn_smooth = interp1(Time_units,entropy_vn-entropy_vn(end),time_smooth)  ;  
  coherence_smooth = interp1(Time_units,rho_red(:,1,2) ,time_smooth);   
  
[om_rng,sig] = ezfft(time_smooth*convfact,entropy_vn_smooth); 
if 1==0   
figure
plot(om_rng,sig)

sig2 = fft(rho_red(:,1,2));
      
figure
plot(om_rng,abs(sig2(1:length(om_rng))))
end
        %% plot stuff
        
figure
plot(Time_units,abs(rho_red(:,2,2)),'r')
hold on
plot(Time_units,abs(rho_red(:,1,1)),'g')
plot(Time_units,abs(rho_red(:,2,1))) % coherence
xlabel('Time (ps)')
ylabel('|\rho(X,Y)|')

figure
plot(Time_units,entropy_vn)
xlabel('Time (ps)')
ylabel('-Tr(\rho log(\rho))')
% figure
% plot(Time_units,purity)
% xlabel('Time (ps)')
% ylabel('Tr(\rho^2)')