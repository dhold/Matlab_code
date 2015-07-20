
if length(numpoints)>1
tmp = open(save_file_name); convfact = tmp.convfact; cnt=1;
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
       brklp  =true; 
    end
    if brklp
        break
    end
end
end
Time_units = Time_units/convfact;
%%
%Script to analyse data produced by the dimer_elec_vib_HOM thingy

clear rho00 test pop test2 popdiag
for k = 1:length(Time_units) 
    rho00{k} = reshape(rho_vec(k,1:(N*viblvls)^2),N*viblvls,N*viblvls);
    rho00{k} = basis_proj'*rho00{k}*basis_proj;
    %projected to exciton basis
end
cnt=0; polarization = zeros(length(Time_units),3); %3d vector
mufull = zeros(N,N,3);
for k = 1:N-1
mufull(1,k+1,:) = mu(k,:);
end
mufull = mufull + permute(mufull,[2,1,3]);

    Efield = subs(E_t_1,[Env_1, omega_1],[exp(-(tt-t01)^2/tau1^2), om1])+...
        subs(E_t_2,[Env_2,omega_2],[exp(-(tt-t02)^2/tau2^2),om2]);   
    Efield = double(subs(Efield,tt,Time_units*convfact)); %convert back to cm^-1
    
for k = 1:length(Time_units)

    temp = TrX(rho00{k},2,[N,viblvls]);

    
    polarization(k,1) = trace(mufull(:,:,1)*temp);
    polarization(k,2) = trace(mufull(:,:,2)*temp);
    polarization(k,3) = trace(mufull(:,:,3)*temp);
   for j = 1:numel(temp) 
       
        pop(k,j) = temp(j);
    %population of exciton states and coherences
   end
       popdiag(k,1:length(temp)) = diag(temp);
end

%%  Calculate density matrix

figure
hold on
plot(Time_units, abs(popdiag))
plot(Time_units, abs(pop(:,2:3)),'--') %ground excited coherences
%plot(Time_units, abs(pop(:,6)),'.-')
%include also the electric field amplitudes
plot(Time_units,Efield(:,3)/max(Efield(:,3)),'k')

xlabel('Time (ps)');
ylabel('Abs value of density element');


%% Polarization calculation

figure
plot(Time_units,polarization)

signal = trapz(Time_units,real(dot(polarization.',Efield.').^2));

figure
plot(Time_units,real(dot(polarization.',Efield.')).^2)