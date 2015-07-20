
 
for kk = 1:3
    if kk ==1
        tmp = open(flename_x);
    elseif kk==2
        tmp = open(flename_y);
    else
        tmp = open(flename_z);
    end
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
  
    if kk==1
    C_t(k) = trace(mu_x * temp);
    elseif kk==2
    C_t(k) = trace(mu_y * temp);    
    else
    C_t(k) = trace(mu_z * temp);    
    end
end

%% interpolate to get constant time spacing

wanted_samp_freq = av_freq*4;
LL = ceil(Time_units(end)*wanted_samp_freq)+1;

wanted_t_spacing = linspace(Time_units(1),Time_units(end),LL);
%should be the same for all of them
   C_t_interp = interp1(Time_units, C_t,wanted_t_spacing,'spline');

samp_F = 1/(wanted_t_spacing(2) - wanted_t_spacing(1));
NFFT = 2^nextpow2(LL);
%add back in frequency
I_om = 2*real(fft( exp(1i*wanted_t_spacing*av_freq).*C_t_interp,NFFT)/LL);
    if kk==1
    C_tx = C_t; I_om_x = I_om;
    elseif kk==2
    C_ty = C_t; I_om_y = I_om;   
    else
    C_tz = C_t; I_om_z = I_om;   
    end
end

I_om_av = (I_om_x+I_om_y+I_om_z)/3;
figure

plot( 2*pi*samp_F/2*linspace(0,1,NFFT/2 +1) ,I_om_av (1:NFFT/2+1))
xlabel('Freq (cm^{-1})|')
ylabel('|I(\omega)')

hold on 
plot( 2*pi*samp_F/2*linspace(0,1,NFFT/2 +1) ,I_om_x (1:NFFT/2+1),'g')
plot( 2*pi*samp_F/2*linspace(0,1,NFFT/2 +1) ,I_om_y (1:NFFT/2+1),'r')
plot( 2*pi*samp_F/2*linspace(0,1,NFFT/2 +1) ,I_om_z (1:NFFT/2+1),'k')