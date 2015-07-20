    
viblvls = {[],4,4}; %number of modes to include on each site

for j = 1:length(om_0)
coupg{j} = lambda{j}.^2./om_0{j};  Hvib_sec{j} = 0; Hvib_couple_sec{j} = 0;
if ~isempty(om_0{j})
    
%     om_0_set = om_0{j}; 
%     for k = 1:length(om_0_set)
%         Hvib_sec{j,k} = diag(0:viblvls{j}(k))*om_0_set(k);
%         tmp = diag(coupg{j}(k)*sqrt(1:viblvls{j}(k)),1);
%         Hvib_couple_sec{j,k} =  tmp+tmp.';
%     end
    %tensor prouct sum these together  kron(ones(length(H_{j,k}
 
    %assume there is only one mode per site for now

        Hvib_sec{j} = diag(1/2 +(0:viblvls{j}))*om_0{j};
       
        tmp = diag(coupg{j}*sqrt(1:viblvls{j}),1);
        Hvib_couple_sec{j} =  tmp+tmp.';   
    
end     
end

%general vibrational Hamiltonian
szmat = cellfun(@length,Hvib_sec); szmat2 = cellfun(@length,Hvib_couple_sec);
H_vib = kron(Hvib_sec{1},eye(prod(szmat(2:end))));

temp = zeros(size(H0)); temp(1,1) = 1;
H_vib_couple = kron(temp ,kron(Hvib_couple_sec{1},eye(prod(szmat2(2:end)))));


for j = 2:length(om_0)  
    H_vib = H_vib + kron(kron(eye(prod(szmat(1:j-1))),Hvib_sec{j})...
                        ,eye(prod(szmat(j+1:end))));   
                    temp = zeros(size(H0)); temp(j,j) = 1;
    H_vib_couple = H_vib_couple + kron(temp,kron(kron(eye(prod(szmat2(1:j-1))),...
        Hvib_couple_sec{j}) ,eye(prod(szmat2(j+1:end)))));    
end

H0full = kron(H0,eye(size(H_vib))) + kron(eye(size(H0)),H_vib) +...
        H_vib_couple ;