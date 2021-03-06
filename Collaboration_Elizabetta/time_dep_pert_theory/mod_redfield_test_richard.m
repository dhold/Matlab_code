%%
 Hamfn = @(delta_E,V) [delta_E/2, V;V, -delta_E/2]; %hamiltonian
 thetafn = @(delta_E,V)  atan(2*V/(delta_E))/2;
 c_nkfn = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)];
 ex_Efn = @(delta_E,V) [1,-1]*sqrt(delta_E^2+4*V^2);
damping =3;
 mode_params = [ 97.,0.02396,damping ;...
                          138.,0.02881,damping ;...
                          213.,0.03002,damping ;...
                          260.,0.02669,damping ;...
                          298.,0.02669,damping ;...
                          342.,0.06035,damping ;...
                          388.,0.02487,damping ;...
                          425.,0.01486,damping ;...
                          518.,0.03942,damping ;...
                          546.,0.00269,damping ;...
                          573.,0.00849,damping ;...
                          585.,0.00303,damping ;...
                          604.,0.00194,damping ;...
                          700.,0.00197,damping ;...
                          722.,0.00394,damping ;...
                          742.,0.03942,damping ;...
                          752.,0.02578,damping ;...
                          795.,0.00485,damping ;...
                          916.,0.02123,damping ;...
                          986.,0.01031,damping ;...
                          995.,0.02274,damping ;...
                          1052.,0.01213,damping ;...
                          1069.,0.00636,damping ;...
                          1110.,0.01122,damping ;...
                          1143.,0.04094,damping ;...
                          1181.,0.01759,damping ;...
                          1190.,0.00667,damping ;...
                          1208.,0.01850,damping ;...
                          1216.,0.01759,damping ;...
                          1235.,0.00697,damping ;...
                          1252.,0.00636,damping ;...
                          1260.,0.00636,damping ;...
                          1286.,0.00454,damping ;...
                          1304.,0.00576,damping ;...
                          1322.,0.03032,damping ;...
                          1338.,0.00394,damping ;...
                          1354.,0.00576,damping ;...
                          1382.,0.00667,damping ;...
                          1439.,0.00667,damping ;...
                          1487.,0.00788,damping ;...
                          1524.,0.00636,damping ;...
                          1537.,0.02183,damping ;...
                          1553.,0.00909,damping ;...
                          1573.,0.00454,damping ;...
                          1580.,0.00454,damping ;...
                          1612.,0.00454,damping ;...
                          1645.,0.00363,damping ;...
                          1673.,0.00097,damping ];
    
lam_ud = mode_params(:,1).*mode_params(:,2);
 reorg_energy = 37; % wavenumbers
cutoff_freq = 30;
       
        
delta_E_values = [0:10:90,linspace(100,2000,100)]; %wavenumbers
coupling_values = [225, 100, 55];
temperature = 77; % Kelvin
[convfact, B,speed_unit]= inv_cm_unit_sys(temperature);

%% Calculate line broadening functions

    g_cnst = line_broad_fn_markov(B, mode_params(:,3),mode_params(:,1),lam_ud,...
                                cutoff_freq,reorg_energy);

Gamma_min = 1/min(real(g_cnst));
t_int_rng = [linspace(eps(Gamma_min),1.5*Gamma_min,150000),...
            Gamma_min*exp(linspace(log(1.5),log(100),150000))];
 t_int_rng = t_int_rng([true,diff(t_int_rng)~=0]);
%t_int_rng = linspace(0,0.1,1000);
 tol = 1e-13;
 g_broad   = line_broad_fn_full(B, mode_params(:,3),mode_params(:,1),lam_ud,...
                                cutoff_freq,reorg_energy,t_int_rng, tol);    
g_deriv  = line_broad_fn_deriv(B, mode_params(:,3),mode_params(:,1),lam_ud,...
                                cutoff_freq,reorg_energy,t_int_rng, tol);
g_sec_der = line_broad_fn_sec_der(B, mode_params(:,3),mode_params(:,1),lam_ud,...
                                cutoff_freq,reorg_energy,t_int_rng, tol);
lam_tot = sum(lam_ud)+sum(reorg_energy);
%% Calculate modified Redfield for each one
ratesave  = zeros(length(coupling_values),length(delta_E_values));
ratesave2  = zeros(length(coupling_values),length(delta_E_values));
for lp1=  1:length(coupling_values)
    V = coupling_values(lp1);
    for lp2 = 1:length(delta_E_values)
            delta_E = delta_E_values(lp2);
            
            Ham = Hamfn(delta_E,V);
            [a,b]= eig(Ham); ex_E = diag(b); c_nk = a.';
            
% theta = atan(2*V/(delta_E))/2;
% c_nk = [cos(theta),-sin(theta);sin(theta),cos(theta)];
% ex_E =  [1,-1]*sqrt(delta_E^2+4*V^2);            
R_reduced =  mod_redfield_calc2(g_broad,g_deriv,g_sec_der,lam_tot,...
                     c_nk,ex_E,t_int_rng);

% R_reduced =  mod_redfield_calc3(B, mode_params(:,3),mode_params(:,1),...
%     lam_ud, cutoff_freq,reorg_energy,c_nk,ex_E);

ratesave(lp1,lp2) = R_reduced(1,2);
ratesave2(lp1,lp2) = R_reduced(2,1);
    end
end
%%
figure
plot(delta_E_values,-ratesave*convfact)
ylabel('rate ps^{-1}')
xlabel('\Delta E')
% figure
% plot(delta_E_values,-ratesave2*convfact)