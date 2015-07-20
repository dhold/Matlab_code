%%

light_speed = 299792458;  length_unit = 100;
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23);  
    convfact = 2 * pi * light_speed * length_unit * 10^(-12);
    
%%   
if 1==1
flename = {'test_specxR1.mat','test_specyR1.mat','test_speczR1.mat',...
     'test_specxL1.mat','test_specyL1.mat','test_speczL1.mat'};
reorg_E = {[],25,25,25,25,25,25,25,25}; 
drude_gam = {[],35,35,35,35,35,35,35,35};
Temp_kel = 77;

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
R_disp = R_disp- repmat(mean(R_disp),8,1);
R_disp = R_disp*10^(-8); % convert from angstrom to inverse cm
N =length(reorg_E);
om_0 = 12590; om_0j = om_0 + cellfun(@sum,reorg_E); om_0j(1) = 0;

Ham = diag(om_0j); MM = -20;  V = MM*(diag(ones(N-2,1),1) + diag(ones(N-2,1),-1));
%periodic boundary conditions mean that first and last are also coupled
%not that M(1,:) relates to the ground state only
V(1,end) = MM;  V(end,1) = MM;
Ham(2:end,2:end) = Ham(2:end,2:end) + V;

Kappa = 0; Kap1 = inf; Kap2 = 6; tendps = 3;
 numpoints =[400, linspace(0,tendps,20).*convfact]; 
 
else
flename = {'test_specxR.mat','test_specyR.mat','test_speczR.mat',...
     'test_specxL.mat','test_specyL.mat','test_speczL.mat'};

reorg_E = {[],20,20}; 
drude_gam = {[],53,53};
Temp_kel = 300;
dipole_vec = [1,0,0;1,0,0];
R_disp = [0,0,0;0,0,0]; %not important here, but in general is

 Kappa = 0; 
 Kap2 = 10; %max level truncation
  Kap1= inf; %max frequency truncation
  
om_0 = 10000; om_0j = cellfun(@sum,reorg_E)+[0, om_0 , om_0+300]; 
Ham = diag(om_0j); N = length(Ham);
MM = 100;  V = MM*(diag(ones(N-2,1),1) + diag(ones(N-2,1),-1));
V(1,end) = MM;  V(end,1) = MM;
Ham(2:end,2:end) = Ham(2:end,2:end) + V;    
end
 %%
 
[V_mat,nn,av_freq,freq_scale] = spectroscopy_via_HEOM_fn(flename, Ham,reorg_E, ...
                drude_gam, Temp_kel,dipole_vec,R_disp,Kappa,Kap1,Kap2,tendps, numpoints);


%% this version considers a more general case and included CD
 clear I_om
for kk = 1:6

        tmp = open(flename{kk});

 convfact = tmp.convfact; cnt=1;
 
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
 
%Time_units = Time_units/convfact;

if use_reduced_mat
             temp = reshape(tril(true(N*vibstates)) , 1,(N*vibstates)^2);
           rho_vec_full = zeros(length(Time_units),size(nn,1)*(N*vibstates)^2);
            temp2 = repmat(temp,1,size(nn,1)); %only nn to tier saved now
       rho_vec_full(:,temp2) = rho_vec; 
       rho_vec = rho_vec_full; clear rho_vec_full
end
clear test pop test2


rho00 = zeros(N*vibstates,N*vibstates,length(Time_units));
for k = 1:length(Time_units)
    rho00(:,:,k) = reshape(rho_vec(k,1:(N*vibstates)^2),N*vibstates,N*vibstates);
    if use_reduced_mat
    rho00(:,:,k) = rho00(:,:,k) + rho00(:,:,k)'-diag(diag(rho00(:,:,k))); %add cc
    end
end

  C_t = zeros(size(Time_units));
  tmp2 = V_mat(:,:,kk);
  
for k = 1:length(Time_units)
    if vibstates > 1
    temp = TrX(rho00(:,:,k),2,[N,vibstates]);
    else
        temp = rho00(:,:,k);
    end
    %CC_t{kk}(k) = trace(tmp2 * temp);
    %tsave{kk}(k) = Time_units(k);
    C_t(k) = trace(tmp2 * temp);

end

%% interpolate to get constant time spacing

wanted_samp_freq = av_freq;
LL = ceil(Time_units(end)*wanted_samp_freq)+1;

wanted_t_spacing = linspace(Time_units(1),Time_units(end),LL);
%should be the same for all of them
   C_t_interp = interp1(Time_units, C_t,wanted_t_spacing,'spline');
% hold on
 %plot(wanted_t_spacing ,log(abs(C_t_interp)))
samp_F = 1/(wanted_t_spacing(2) - wanted_t_spacing(1));
NFFT = 2^nextpow2(LL);
%add back in frequency
I_om{kk} = 2*real(fft(exp(1i*wanted_t_spacing*av_freq).*C_t_interp,NFFT)/LL);
%I_om{kk} = 2*real(fft(exp(1i*wanted_t_spacing*av_freq/10).* C_t_interp,NFFT)/LL);
end
%%
I_om_x = (I_om{1}+I_om{4})/2; I_om_y = (I_om{2}+I_om{5})/2; I_om_z = (I_om{3}+I_om{6})/2;
CD_x = I_om{1}-I_om{4}; CD_y = I_om{2}-I_om{5}; CD_z = I_om{3}-I_om{6};

I_om_av = (I_om_x + I_om_y +I_om_z)/3;

CD_val = (CD_x+CD_y+CD_z)/3;

figure

freqrng = 2*pi*samp_F/2*linspace(0,1,NFFT/2 +1);
freqrng_hz = freqrng.*length_unit*light_speed;
wlrng_nm = 10^7./freqrng;
%%
plot(freqrng,I_om_av (1:NFFT/2+1),'k')
xlabel('\omega (cm^{-1})|')
ylabel('I(\omega)')

 hold on 
 plot(freqrng ,[I_om_x(1:NFFT/2+1);...
                I_om_y(1:NFFT/2+1);I_om_z(1:NFFT/2+1)])
            
figure

plot(freqrng ,CD_val (1:NFFT/2+1),'k')
xlabel('\omega (cm^{-1})|')
ylabel('CD(\omega)')

 hold on 
 plot(freqrng ,[CD_x(1:NFFT/2+1);...
                CD_y(1:NFFT/2+1);CD_z(1:NFFT/2+1)])
   %%
   figure
plot(wlrng_nm ,I_om_av (1:NFFT/2+1),'k')
xlabel('\lambda (nm)')
ylabel('Signal (au)')

 hold on 
 plot(wlrng_nm  ,[I_om_x(1:NFFT/2+1);...
                I_om_y(1:NFFT/2+1);I_om_z(1:NFFT/2+1)])
            
figure

plot(wlrng_nm  ,CD_val (1:NFFT/2+1),'k')
xlabel('\lambda (nm)')
ylabel('Signal (au)')

 hold on 
 plot(wlrng_nm  ,[CD_x(1:NFFT/2+1);...
                CD_y(1:NFFT/2+1);CD_z(1:NFFT/2+1)])   
   
   
  %%
figure
plot(freqrng ,[I_om{1};I_om{2};I_om{3};I_om{4};I_om{5};I_om{6}])
xlabel('\omega (cm^{-1})|')
ylabel('|CD(\omega)')          