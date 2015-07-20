save_file_name ='saved_data_test4.mat';
debug_md = true;
%use_reduced_mat=0
%% If loading from a save file run this one
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
% use to analyse the output of the reduced vector (only LI elements)
%to get rho_vec from rho_vec red 
         if use_reduced_mat
             temp = reshape(tril(true(N*vibstates)) , 1,(N*vibstates)^2);
        % if saveonlyrho00
        % rho_vec_full = zeros(length(Time_units),(N*vibstates)^2);
        %    temp2 = temp;
        % else
           rho_vec_full = zeros(length(Time_units),size(nn,1)*(N*vibstates)^2);
            temp2 = repmat(temp,1,size(nn,1)); %only nn to tier saved now
        % end
       rho_vec_full(:,temp2) = rho_vec; 
       rho_vec = rho_vec_full; clear rho_vec_full
         end
% extra density matrix as normal and then just add the conjugate of matrix
% minus the leading diagonal

%%
%Script to analyse data produced by the dimer_elec_vib_HOM thingy
rho00 = zeros(N*vibstates,N*vibstates,length(Time_units));

clear test pop test2
for k = 1:length(Time_units)
   
    rho00(:,:,k) = reshape(rho_vec(k,1:(N*vibstates)^2),N*vibstates,N*vibstates);
    if use_reduced_mat
    rho00(:,:,k) = rho00(:,:,k) + rho00(:,:,k)'-diag(diag(rho00(:,:,k))); %add cc
    end
    rho00(:,:,k) = basis_proj'*rho00(:,:,k)*basis_proj;
    %projected to exciton basis
end
tmp = TrX(rho00(:,:,1),2,[N,vibstates]);
    pop = zeros(length(Time_units),length(tmp));
    coherences = zeros(length(Time_units),(numel(tmp)-length(tmp))/2);
for k = 1:length(Time_units)

    temp = TrX(rho00(:,:,k),2,[N,vibstates]);
    if debug_md
    test(k) = trace(abs(rho00(:,:,k))); %#ok<*UNRCH>
    test2(k) = trace(rho00(:,:,k)^2);
    end
   for j = 1:numel(temp) 
    pop(k,:) = diag(temp);
    tmp = tril(temp,-1);
    coherences(k,:) = reshape(tmp(tmp~=0),(numel(tmp)-length(tmp))/2,1);
    %population of exciton states and coherences
   end
end
%%used for testing
if debug_md
figure
plot(Time_units,log(abs(1-test)))
figure
plot(Time_units,real(test2-1)) %must be -ve
end
%figure
hold on
plot(Time_units, abs(pop))
% Create xlabel
xlabel('Time (ps)');

% Create ylabel
ylabel('Abs value of density element');
%% might take longer

if saveuptotier>=1
    
    %recursively generate coeffmat
% coeffmat = zeros(n+1);
% coeffmat(1,1) = 1;
% %recursion relation c(n + 1,j + 1) = -c(n,j) + (n-1) A_0 c(n-1,j+1)
% A_0 = sum(cc,2); %in general length N vector
% coeffmat(2,2) = -1;  coeffmat(3,3) = 1;
% coeffmat(3,1) = A_0;
% 
% for nlp = 2:n
%     for jlp = 1:n
%             
%         coeffmat(nlp+1,jlp+1) = -coeffmat(nlp,jlp) + (nlp-1)*A0*coeffmat(nlp-1,jlp+1);
%         
%     end
% end
A_0 = sum(cc,2); %in general length N vector
for lp = 1:length(A_0)
coeffmat{lp} = coeffmat_gen(saveuptotier,saveuptotier,A_0(lp));
end
    
    %undoall rescaling
    scalefct = repmat(sqrt(abs(horzcat(cc_acom{:})).^2 + abs(horzcat(cc_com{:})).^2),size(nn,1),1);
    scalefct = sqrt(scalefct).^nn;
    scalefct = prod(sqrt(factorial(nn)).*scalefct,2);
    
    tier1 = sum(nn,2)==1;
    
    X1 = zeros(length(Time_units),N);
for k = 1:length(Time_units)
   %for j =1:N
       temp2 = zeros(N,1); jjrng = find(tier1);
       for j =1:N
           N2 = size([cc2R,cc1],2);
           jjrng2 = jjrng(1+(j-1)*N2:N2*j);
       for jj = jjrng2(1):jjrng2(end)
           rng = 1+(jj-1)*(N*vibstates)^2:jj*(N*vibstates)^2;
            temp= reshape(rho_vec(k,rng),N*vibstates,N*vibstates);
%     if use_reduced_mat
%     temp =temp + temp'-diag(diag(temp)); %add cc
%     end
            temp2(j) = temp2(j) - trace(temp)*scalefct(jj);
       end
       end
   %end
   X1(k,:) = temp2 ;
end
end
%% Further debug
if size(nn,1)>1 && debug_md  %check growth of higher order components
   nnlvl =1; %order
   maxnlvl = max(max(nn)); prevpos = (N*vibstates)^2;
    while nnlvl <= maxnlvl
        clear svd
       for j = 1:sum(sum(nn,2)==nnlvl)
           pos = ((j-1)*(N*vibstates)^2+1+prevpos):j*(N*vibstates)^2+prevpos;
        for k = 1:length(Time_units)
    temp = reshape(rho_vec(k,pos),N*vibstates,N*vibstates);
    if use_reduced_mat
    temp = temp + temp'-diag(diag(temp)); %add cc
    end
    svd(k,j) = trace(abs(temp));
    temp2 = temp-diag(diag(temp));
    svd2(k,j) = sum(sum(abs(temp2-temp2')));
        end
       end
       growth_check{nnlvl} = svd;
       ishermitian_check{nnlvl} = svd2;
        nnlvl = nnlvl +1;
        prevpos = max(pos);
        
    end
    
end    
    
