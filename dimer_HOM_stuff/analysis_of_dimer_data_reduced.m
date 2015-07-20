N=2; viblvls=10;

%%
% use to analyse the output of the reduced vector (only LI elements)
%to get rho_vec from rho_vec red 
         rho_vec_full = zeros(length(Time_units),size(nn,1)*(N*viblvls)^2);
         %rho_vec_full = zeros(length(Time_units),(N*viblvls)^2);
             temp = reshape(tril(true(N*viblvls)) , 1,(N*viblvls)^2);
            temp2 = repmat(temp,1,size(nn,1));
            %temp2 = temp;
       rho_vec_full(:,temp2) = rho_vec; 
       rho_vec = rho_vec_full; clear rho_vec_full
% extra density matrix as normal and then just add the conjugate of matrix
% minus the leading diagonal

%%
%Script to analyse data produced by the dimer_elec_vib_HOM thingy

clear rho00 test pop test2
for k = 1:length(Time_units)
   
    rho00{k} = reshape(rho_vec(k,1:(N*viblvls)^2),N*viblvls,N*viblvls);
    rho00{k} = rho00{k} + rho00{k}'-diag(diag(rho00{k})); %add cc
    rho00{k} = basis_proj'*rho00{k}*basis_proj;
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
