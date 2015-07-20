clear rho pop
tmp =sqrt(length(rho_vc(1,:))); viblvls = 6; 
for k = 1:length(time_units)
   
    rho{k} = reshape(rho_vc(k,:),tmp,tmp);
    rho{k} = basis_proj'*rho{k}*basis_proj;
    %projected to exciton basis
end

for k = 1:length(time_units)
    temp = TrX(rho{k},2,[2,viblvls+1]);
   for j = 1:numel(temp) 
    pop(k,j) = temp(j);
    %population of exciton states and coherences
   end
end

figure
plot(time_units, abs(pop))
% Create xlabel
xlabel('Time (ps)');

% Create ylabel
ylabel('Abs value of density element');
