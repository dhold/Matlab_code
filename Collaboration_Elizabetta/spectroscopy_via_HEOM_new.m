
reorg_E = {[],25,25,25,25,25,25,25,25}; 
drude_gam = {[],35,35,35,35,35,35,35,35};
Temp_kel = 77;
light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); 
%%
dipole_vec = ...
[ -0.823072771604003,	-0.5246284653639453,    -0.21752284012024087;...
-0.21108138821716999,	-0.9529594677893477,	-0.2175152875063172;...
0.5246284653639459, 	-0.8230727716040029,	-0.21752284012024128;...
0.9529594677893479, 	-0.211081388217168, 	-0.21751528750631696;...
0.8230727716040026, 	0.524628465363946,      -0.2175228401202412;...
0.21108138821716746,	0.9529594677893481,     -0.2175152875063173;...
-0.5246284653639459,	0.8230727716040029,     -0.21752284012024128;...
-0.9529594677893481,	0.21108138821716843,	-0.2175152875063174];
%magnitude of dipole_vec isn't important here, just constant of
%proportionality, as it is not used to calculate the coupling between sites
R_disp = [72.062,-12.081,90.095; 72.913	,9.924,90.122;...
19.538	,11.081	,90.095; 18.687,-10.924	,90.122;...
57.881	,25.262	,90.095;  35.876,26.113	,90.122;...
33.719,	-26.262	,90.095; 55.724	,-27.113,90.122];
%%
om_0 = 12590; om_0j = om_0 + cellfun(@sum,reorg_E); om_0j(1) = 0;

Ham = diag(om_0j); N = length(Ham);
%%%%
av_freq = mean(om_0j(2:end));
freq_scale = [0,ones(1,N-1)*av_freq];
Ham = Ham - diag(freq_scale);
%%%%%
MM = -20;  V = MM*(diag(ones(N-2,1),1) + diag(ones(N-2,1),-1));
%periodic boundary conditions mean that first and last are also coupled
%not that M(1,:) relates to the ground state only
V(1,end) = MM;  V(end,1) = MM;

Ham(2:end,2:end) = Ham(2:end,2:end) + V;

rho_0 = zeros(size(Ham)); rho_0(1,1) = 1;

mu_x = Ham*0;  mu_x(:,1)=  [0;dipole_vec(:,1)]; mu_x = mu_x+mu_x';
mu_y = Ham*0;  mu_y(:,1)=  [0;dipole_vec(:,2)]; mu_y = mu_y+mu_y';
mu_z = Ham*0;  mu_z(:,1)=  [0;dipole_vec(:,3)]; mu_z = mu_z+mu_z';

rho_0x = mu_x * rho_0;   rho_0y = mu_y * rho_0;  rho_0z = mu_z * rho_0;  

 Kappa = 0; 
 Kap2 = 6; %max level truncation
  Kap1= inf; %max frequency truncation
%%
N2 = sum(cellfun(@length,reorg_E));
 QQ = zeros(N,2); cnt=0;
 cc = zeros(1,sum(N2 + Kappa));
 for j=1:N
 [cc1,cc2R,cc2I,vv1,vv2,QQ(j,1),tmp] = ...
    coeffients_from_brownian_new([],[],[],Temp_kel,Kappa,reorg_E{j},drude_gam{j});
 % QQ is in general length two for each point as I could (potentially)
 % included high frequency decay terms from poles in J(omega) which would
 % have anti commutator parts
 cc(cnt+1:cnt+length(tmp)) = tmp; cnt =cnt+length(tmp);
 cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];
  cc_acom{j}= [cc2I,cc1*0];
 end


 tendps =2; 
 convfact = 2 * pi * light_speed * length_unit * 10^(-12);
 numpoints =[50, linspace(0,tendps,150).*convfact]; 
 saveuptotier =0;
 
 use_reduced_mat = false;
 vibstates = 1; flename_x = 'saved_data_specx.mat';
 flename_y = 'saved_data_specy.mat';
 flename_z = 'saved_data_specz.mat';
%%
 [Time_units_x,rho_vec_x,nn]=multi_site_drude_and_brownian_HOM_3...
            (Ham,QQ,rho_0x,cc_com,cc_acom,vv,Kap1,Kap2,vibstates,...
   numpoints,tendps,saveuptotier,use_reduced_mat,flename_x);
%%
 [Time_units_y,rho_vec_y]=multi_site_drude_and_brownian_HOM_3...
            (Ham,QQ,rho_0y,cc_com,cc_acom,vv,Kap1,Kap2,vibstates,...
   numpoints,tendps,saveuptotier,use_reduced_mat,flename_y);
%%
 [Time_units_z,rho_vec_z]=multi_site_drude_and_brownian_HOM_3...
            (Ham,QQ,rho_0z,cc_com,cc_acom,vv,Kap1,Kap2,vibstates,...
   numpoints,tendps,saveuptotier,use_reduced_mat,flename_z);
%%
if length(numpoints)>1
tmp = open(flename); convfact = tmp.convfact; cnt=1;
brklp = false; Time_units=[];  rho_vec=[];
for k=1:1000
    try
        Y=strcat('saved_time',num2str(cnt)); YY=strcat('saved_rho',num2str(cnt));
        temp = tmp.(Y);  temp2 = tmp.(YY);

        Time_units = [Time_units;temp]; %#ok<*AGROW>
        rho_vec = [rho_vec;temp2];
        cnt=cnt+1;
    catch ME
        %ME
       brklp  = true; 
    end
    if brklp
        break
    end
end
end
%Time_units = Time_units/convfact;

if use_reduced_mat
             temp = reshape(tril(true(N*vibstates)) , 1,(N*vibstates)^2);
           rho_vec_full = zeros(length(Time_units),size(nn,1)*(N*vibstates)^2);
            temp2 = repmat(temp,1,size(nn,1)); %only nn to tier saved now
       rho_vec_full(:,temp2) = rho_vec; 
       rho_vec = rho_vec_full; clear rho_vec_full
end
%%

rho00 = zeros(N*vibstates,N*vibstates,length(Time_units));

clear test pop test2 
for k = 1:length(Time_units)
   
    rho00(:,:,k) = reshape(rho_vec(k,1:(N*vibstates)^2),N*vibstates,N*vibstates);
    if use_reduced_mat
    rho00(:,:,k) = rho00(:,:,k) + rho00(:,:,k)'-diag(diag(rho00(:,:,k))); %add cc
    end

end

C_t = zeros(size(Time_units));
for k = 1:length(Time_units)

    temp = TrX(rho00(:,:,k),2,[N,vibstates]);
    %C_t(k) = trace(mu_x * (temp.*exp(-1i*Time_units(k)*diag(freq_scale))));
    %temp(1,1) = temp(1,1) * exp(1i*Time_units(k)*freq_scale);
    C_t(k) = trace(mu_x * temp);
end

%% interpolate to get constant time spacing

wanted_samp_freq = av_freq*4;
LL = ceil(Time_units(end)*wanted_samp_freq)+1;

wanted_t_spacing = linspace(Time_units(1),Time_units(end),LL);
   C_t_interp = interp1(Time_units, C_t,wanted_t_spacing,'spline');

samp_F = 1/(wanted_t_spacing(2) - wanted_t_spacing(1));
NFFT = 2^nextpow2(LL);
%add back in frequency
I_om = 2*real(fft( exp(1i*wanted_t_spacing*av_freq).*C_t_interp,NFFT)/LL);

figure
plot( samp_F/2*linspace(0,1,NFFT/2 +1) ,I_om(1:NFFT/2+1))
xlabel('Freq (cm^{-1})|')
ylabel('|I(\omega)')
