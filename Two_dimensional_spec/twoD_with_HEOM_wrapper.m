% This code calculates 2D spectroscopy using the HEOM in the exciton basis
% for an excitonic system with a spectral density of underdamped and
% overdamped brownian oscillator modes on each side.

%clear old response functions and perm dipole shifts
clear R1 R2 R3 R4 R5 R6 pdm_shift pdm
if ~exist('params_for_wrapper','var')
    params_for_wrapper = [];
end
if ~exist('t_verbose','var')
    t_verbose = true;
end
%unless explicitly given assume all options are false
    pop_only = false; %don't consider transitions to coherences
    conv_w_gaussian = false; %deal with static disorder in mean energy
    freq_diff = []; %select things that beat at this frequency
    beam_width = []; % time width in inverse cm to use for windowing
    
if exist(params_for_wrapper,'file')==2  
    %if this variable is present and it is a string corresponding to a file
    %it uses this.  Clear the name if this isn't desired and you wish to
    %enter them into the wrapper manually
    
max_tier = 3; Kappa = 4; num_realisations = 1; %default values if nothing loaded

    
 load(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file',...
     'beam_param_set','om1_rng','t2_range_fs','om3_rng','t1_coh','om2_coh','om3_coh',...
    'max_tier','Kappa','num_realisations',...
    'pop_only','conv_w_gaussian','freq_diff','beam_width')

    calc_rp = true;  calc_nr = true;  calc_coh = false; %not saved 
    [t_scale, B,speed_unit]= inv_cm_unit_sys(Temperature_in_Kelvin);
    t2_range = t2_range_fs*t_scale/1e3;   %in inverse cm 
    
    max_tier_final = max_tier;
else %enter them into this format into the wrapper
    for notactuallyaloop=1
chooser = 2; %choose which system stuff to use, this is just a lazy setup 

%% System related parameters
%enter the name of the file where the parameters are stored
if chooser==0
system_parameters_file = 'Parameters/dimer_system_params.mat';
Temperature_in_Kelvin = 300; %self explainatory
elseif chooser  == 1
system_parameters_file = 'Parameters/Adv_HEOM_for_Efficient_eval_params.mat'; 
Temperature_in_Kelvin = 77; %self explainatory    
else %to reproduce for alexandra etc etc
system_parameters_file = 'Parameters/plenio_paper_parameters.mat'; 
Temperature_in_Kelvin = 77; %self explainatory
end
%calculate scaling from ps->inv cm, thermo beta in cm^(-1) 
[t_scale, B,speed_unit]= inv_cm_unit_sys(Temperature_in_Kelvin);
 %multiply time in ps by t_scale to get it in inverse cm
 %t_in_ps * t_scale == t_in_inverse_cm
 
num_realisations = 1; %realisations of static disorder
%% HEOM truncation parameters

max_tier = 3; %max_tier to include
Kappa = 4; %number of matsubara frequencies to include
max_tier_final = max_tier; %max tier to save for final thing.
%max_tier is used for the calculates of doorway and window parts but not
%the final function

%% Beam params, wavevector and polarization

k1 = [1/sqrt(2),0,1/sqrt(2)]; pol_1 = [0,1,0];  mag_1 = cross(k1,pol_1);
k2 = [1/sqrt(2),0,-1/sqrt(2)]; pol_2 = [0,1,0]; mag_2 = cross(k2,pol_2); 
k3 = [0,0,1]; pol_3 = [0,1,0]; mag_3 = cross(k3,pol_2); 
%also specify the signal used for Heterodyne detection (will be conjugated)
 pol_HET = [0,1,0];  pol_HET2 = [1,0,0]; 

%beam_param_set = cat(3,[k1;k2;k3;pol_1;pol_2;pol_3;pol_HET],...
%                 [k1;k2;k3;pol_1;pol_2;pol_3;pol_HET2]);
beam_param_set = cat(3,[pol_1;pol_2;pol_3;pol_HET;k1;k2;k3],...
                 [pol_1;pol_2;pol_3;pol_HET2;k1;k2;k3]);                 
 %put all parameters into a single set 
             
calc_rp = true;  calc_nr = true;  calc_coh = false; 

%% Choose the range of frequencies to tabulate for RP/NRP
%Calculate R(om_3,t_2,om_1) = iint_dt_3 dt_1 R(t_3,t_2,t_1)
if chooser ==0
 om1_rng = linspace(10^7./(850),10^7./(730),100); %range in cm^(-1)
t2_range_fs = (0:2:60);  %delay time t_2 in fs
om3_rng = linspace(10^7./(920),10^7./(680),1000); 
t2_range = t2_range_fs*t_scale/1e3;   %in inverse cm
elseif chooser == 1

om1_rng = linspace(10000-800,10000+800,100); %range in cm^(-1)
t2_range_fs = [0,100,200,300,400,500];  %delay time t_2 in fs
om3_rng = linspace(10000-800,10000+800,1000); 
t2_range = t2_range_fs*t_scale/1e3;    %in inverse cm
else
    
om1_rng = linspace(1.21e4,1.27e4,50); %range in cm^(-1)

t2_range = linspace(0,0.1586,100);
t2_range_fs = t2_range/t_scale*1e3;  %delay time t_2 in fs
om3_rng = linspace(1.21e4,1.27e4,300) ;
    
end
%{
% dom = om1_rng(2)-om1_rng(1);
% om_end = om1_rng(end)-om1_rng(1);
% t_ef = 2*pi*(0:1/om_end:1/dom);
%take long well spaced time intervals when calculating the coherences,
%possible as I am using FFT
t1_range = (0:1e-4:2) /t_scale; 
t3_range = (0:1e-4:2) /t_scale;
%}
%% If required choose range for the coherence pathway

t1_coh = 40; %pass single number to get the program to choose a range 
%based on markovian decay with this number of points
om2_coh= linspace(10^7./(1000),10^7./(700),100); %range 
om3_coh= linspace(10^7./(1000),10^7./(700),100); %range 
    end
end
%% Pass these parameters to the function which calculates quantities 
%independent of static disorder

[coup_com_save,coup_acom_save,const_factor,QQ_topass,nn,fock_space_rep,...
    H_site,site_shift,g_cnst,lam_tot,V_n,av_2,av_3,av2_3,av_4,av_5,av2_5,chk] = ...
    TwoD_with_HEOM_prerun(system_parameters_file,Temperature_in_Kelvin,...
    num_realisations,max_tier,Kappa,beam_param_set) ;

if conv_w_gaussian
%subtract mean difference so average excitation is the same, just include
%the extra shift at the end.
site_shift = site_shift - repmat(sum(site_shift,2),1,length(sd_mat))/2; 
end
%also find pdm_shift is an N X N cell array with the
%shift of the excited state with sites j and k excited from the PDM
load(system_parameters_file,'pdm_shift');
pdm_shift{1,2} =0; pdm_shift{2,1} =0;
%[LASTMSG, LASTID] = lastwarn; %turn off warnings that this variable hasn't loaded
%warning('off',LASTID)
if ~exist('pdm_shift','var') 
   load(system_parameters_file,'pdm');
   if exist('pdm','var') %perm dipole moments known but not shifts to states
        load(system_parameters_file,'R');
        for j =1:N
            for kk = 1:N
pdm_shift{kk,j} = dot(pdm(kk,:),pdm(j,:))-3*dot(pdm(kk,:),R(j,:))...
              *dot(pdm(j,:),R(kk,:))/norm(R(j,:)-R(kk,:))^2;
pdm_shift{kk,j} = pdm_shift/norm(R(j,:)-R(kk,:))^3  ;   
            end
        end       %loop over states
   else
       pdm_shift=[]; %no perm dipole moment
   end
end

%% Now calculate the response functions

sdlp_rng = 1:num_realisations; % in general one can do this multiple 
%times until convergence is achieved

%{
[R1,R2,R3,R4,R5,R6,~,~,lin_spec] = twoD_with_HEOM_main(...
    sdlp_rng,coup_com_save,coup_acom_save,const_factor,...
             QQ_topass,nn,max_tier_final,fock_space_rep,...
            H_site, site_shift,V_n,av_2,av_3,av2_3,av_4,av_5,av2_5,...
            calc_rp,calc_nr,t1_range,t2_range,t3_range,om1_rng,om3_rng,...
            calc_coh,t1_coh,om2_coh,om3_coh);
%}
%Version two doesn't include the possibility of calculating the coherence
%contribution and uses an alternative method calculating things in
%frequency space rather than time
[R1,R2,R3,R4,R5,R6,lin_spec,tmp_gg_t,tmp_ee_t]=twoD_with_HEOM_main2(...
        sdlp_rng,coup_com_save,coup_acom_save,const_factor,QQ_topass...
        ,nn,max_tier_final,fock_space_rep,...
       H_site,site_shift,pdm_shift,V_n,av_2,av_3,av2_3,av_4,av_5,av2_5,  ...
 calc_rp,calc_nr,om1_rng,t2_range,om3_rng,pop_only,t_verbose,freq_diff,beam_width);
    
%% Convolve with Gaussian if required

% Mean frequency shift Delta om_eg = sum_k Delta om_kg is also distributed
% Gaussianly, and is simple to include
N=length(H_site);
sd_var = sqrt(dot(sd_mat,sd_mat))/N; %sigma of normal dist of mean freq shift

%check there is actually static D otherwise no point
if sd_var == 0;  conv_w_gaussian = false;    end

if conv_w_gaussian 
load(system_parameters_file,'sd_mat');

zz = linspace(-10*sd_var,10*sd_var,4000);

Conv_dist = exp(-zz.^2/2/sd_var^2)/sqrt(2*pi*sd_var^2);
%Need to calculate 
%S'(om_3,t,om_1) = int_{-inf}^{inf} dw S(om_3+w,t,om_1+w)*P(w)
%%
%for frequency space rephasing signal
Conv_dist = full(sparse(1:length(Conv_dist),1:length(Conv_dist),Conv_dist));%diag(Conv_dist);

R1 = convn(R1,Conv_dist,'same');
R2 = convn(R2,Conv_dist,'same');
R3 = convn(R3,Conv_dist,'same');
R4 = convn(R4,Conv_dist,'same');
R5 = convn(R5,Conv_dist,'same');
R6 = convn(R6,Conv_dist,'same');

%% Inverse FFT for each value of t2, multiply by xi(t1 +/- t3) fft back

%can pass 
if iscell(om1_rng); t1_rng = om1_rng{2};om1_rng = om1_rng{1};  end;
if iscell(om3_rng); t3_rng = om3_rng{2};om3_rng = om3_rng{1};  end;

test1 = diff(om1_rng,2); test2 = diff(om3_rng,2); 
if any(abs(test1)>eps(max(om1_rng))) ||  any(abs(test2)>eps(max(om3_rng)))
warning('om1_rng and om3_rng are not equally spaced, this will cause problems')
end

N1 = length(om1_rng); N3 = length(om3_rng);
dt1 = 2*pi/(om1_rng(end)-om1_rng(1)); t1 = 0:dt1:(N1-1)*dt1;
dt3 = 2*pi/(om3_rng(end)-om3_rng(1)); t3 = 0:dt3:(N3-1)*dt3;

for lp =1:6 %signals
    for lp2 = 1:(size(av_4,5)+size(av_5,6)) %beam configurations
    Rbroad = zeros(length(om1_rng),length(t2_rng),length(om3_rng));
    eval(strcat('Rtmp=R',num2str(lp),'(:,:,:,',num2str(lp2),');')) 

        for t2lp = 1:length(t2_range)
            Rtmp2 = squeeze(Rtmp(1:length(om3_rng),t2lp,1:length(om1_rng)));
            Rtmp2 = ifft2(Rtmp2); %take inverse fft
            %correct for the fact that om1/3 don't start from zero
            expfact = exp(1i*(t1+t3)) 
            'unfinished'
        end
        
eval(strcat('R',num2str(lp),'(:,:,:,',num2str(lp2),')=Rbroad;'))   %replace with the correct thing        
    end
end

%% Shitty method
%{
om3m = min(om3_rng); om3M = max(om3_rng);
om1m = min(om1_rng); om1M = max(om1_rng);

PP = @(w) exp(-w.^2/2/sd_var^2)/sqrt(2*pi*sd_var^2);
rng1 = 1:length(t2_range_fs);

for lp =1:6
eval(strcat('Rbroad=R',num2str(lp),';')) 
if~isempty(Rbroad) %don't bother if I haven't set this
    Rtmp = Rbroad;  rng2 = 1:size(Rtmp,4);
    for lp1 = 1:length(om1_rng); om1 = om1_rng(lp);
        for lp3 = 1:length(om3_rng);  om3 = om1_rng(lp);

wmin = max(om3-om3m,om1-om1m); wmax =  min(om3M-om3,om1M-om1);
%range of w which I can use given the range of om3 and om1

norm_fct = integral(PP,wmin,wmax);
%effective normalisation due to trunction, should be ~ 1 in most cases and
%greater than a half
if length(rng2)>1
FF = @(w) PP(w).*interpn(om3_rng,rng1,om1_rng,rng2,Rtmp,om3+w,rng1,om1+w,rng2);
else
FF = @(x) PP(x).*interpn(om3_rng,rng1,om1_rng,Rtmp,om3+x,rng1,om1+x);    
end
Rbroad(lp3,:,lp1,:) = integral(FF,wmin,wmax,'RelTol',1e-6,'ArrayValued',true)/norm_fct; %#ok<SAGROW>
        end
    end
eval(strcat('R',num2str(lp),'=Rbroad;'))   %replace with the correct thing
end
end
clear Rtmp Rbroad
%}
end