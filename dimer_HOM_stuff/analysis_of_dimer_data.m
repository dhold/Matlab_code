%Script to analyse data produced by the dimer_elec_vib_HOM thingy
N=2; viblvls=3;
clear rho00 test pop test2
for k = 1:length(Time_units)
   
    rho00{k} = basis_proj'*reshape(rho_vec(k,1:(N*viblvls)^2),N*viblvls,N*viblvls)*basis_proj;
    test(k) = sum(sum(abs(rho00{k}-rho00{k}'))); %should be zero or nearly
    %projected to exciton basis
end
 for k = 1:length(Time_units)
    
     rho10{k} = reshape(rho_vec(k,(N*viblvls)^2+1:2*(N*viblvls)^2),N*viblvls,N*viblvls);
     test2(k) = sum(sum(abs(rho10{k}-rho10{k}')));
     %projected to exciton basis
 end


for k = 1:length(Time_units)
    temp = TrX(rho00{k},2,[N,viblvls]);
   for j = 1:numel(temp) 
    pop(k,j) = temp(j);
    %population of exciton states and coherences
   end
end

figure
plot(Time_units, abs(pop))
% Create xlabel
xlabel('Time (ps)');

% Create ylabel
ylabel('Abs value of density element');
% figure
% plot(Time_units, abs(test))
% % Create xlabel
% xlabel('Time (ps)');
% 
% % Create ylabel
% ylabel('Non hermitivity (bad)');

% figure
% plot(Time_units,[test;test2])