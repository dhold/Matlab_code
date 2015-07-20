N=2; viblvls=3;
clear rho00 test pop pop2 pop2_intp
for k = 1:length(Time_units)
   
    rho00{k} = basis_proj'*reshape(rho_vec(k,1:(N*viblvls)^2),N*viblvls,N*viblvls)*basis_proj;

    temp = TrX(rho00{k},2,[N,viblvls]);
   for j = 1:numel(temp) 
    pop(k,j) = abs(temp(j));
    %population of exciton states and coherences
   end
    
   
end
for k = 1:length(Time_units2)
   
    rho002{k} = basis_proj'*reshape(rho_vec2(k,1:(N*viblvls)^2),N*viblvls,N*viblvls)*basis_proj;

    temp = TrX(rho002{k},2,[N,viblvls]);
   for j = 1:numel(temp) 
    pop2(k,j) = abs(temp(j));
    %population of exciton states and coherences
   end
end
for j = 1:numel(temp)

        %pop2_intp(:,j) = pchip(Time_units2,pop2(:,j),Time_units);
        pop2_intp(:,j) = interp1(Time_units2,pop2(:,j),Time_units);
        
end
figure
plot(Time_units,abs(pop-pop2_intp))