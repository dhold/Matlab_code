function [total_prop_op,nn]=HEOM_op_gen(Htot,fock_space_rep,...
                        QQ,cc_com,cc_acom,vv,Kap1,Kap2,viblvls)
%%
 % Htot is total system Hamiltonian
 % fock_space_rep shows the locations of excitations in the system the
 % Hamiltonian relates to, this code can in theory deal with many
 % excitations in the system
 % QQ is the extra convergence parameter made up of all the terms of the 
 % form int_0^t d tau exp(-v (t-tau)) * ( c_1 V^X(tau) + c_2 V^O(tau) )
 %  which are evaluated via exp(-v (t-tau)) ~ delta(t-tau)/v, as v->inf
 % and are thus just markovian terms.
 % cc_com cc_acom and vv are coefficients of exponential decay terms in the
 % bath correlation function these should be N (number of site) long 
 % 
 % Kap1 is max decay frequency truncation
 % Kap2 is max heirarchy order
 % viblvls is number of viblvls included explicitly

 N = length(cc_com); %pass only the baths on each site
 N2 = length(fock_space_rep); %total number, should be 1+N+N*(N-1)/2
 
totpoles = sum(cellfun(@length,cc_com));

%% Construct Ltot, the total Louville for self coupling density matricies

% For N X N A and B we have the following
% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho

 LL = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));
% QQ as input should be the constants in the correction 
if size(QQ) == size(LL) %assume matrix given that is Q
    Q = QQ; clear QQ
else
if size(QQ,1) == 1 %assume same for every site
    QQ = repmat(QQ,N*viblvls,1);
end
 Q = LL*0;
% 
% %if include_truncation_correction
for j = 1:N
    
    Qj = zeros(N2); lg = fock_space_rep(:,j); %states with excitation on site j
    Qj(lg,lg) = 1;   Qj = diag(diag(Qj)); %only want diagonal ones
    Qj = kron(Qj,eye(viblvls));   eyeQ = eye(length(Qj));
    
      %note Qj^2 = Qj so I don't need to worry about this

    Q = Q + QQ(j,1).*(kron(Qj.',eyeQ) + kron(eyeQ,Qj) - 2*kron(Qj.',Qj)); 
    Q = Q  + QQ(j,2) .*(kron(Qj.',eyeQ) - kron(eyeQ,Qj)); 
    
    %reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
    %[Q_i,[Q_i,rho_n]] = (1 X Q_i*Q_i + Q_i^T*Q_i^T X 1 - 2 Q_i^T X Q_i ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
    
    % .*(kron(Vj.',eyeQ) + kron(eyeQ,Vj)); 
end
% %end
end
Ltot = (LL-Q);
%in principle Q can also include high frequency decaying components from
% poles in J(omega) 

%% Construct labels for auxilliary density matrix

if isempty(Kap1)
    Kap1 = inf; %no truncation based on Kap1
end
    
if Kap2>0  %else trivial

numwithn = zeros(1,Kap2+1); %size(cc2,2) is number of poles in J(omega)
tot_poles = sum(cellfun(@length,cc_com));
for kk = 0:Kap2
numwithn(kk+1) = nchoosek(tot_poles-1+kk,kk);
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+size(cc2,2))) 
end

if length(Kap1)>1 %this is an overloading to allow you to just pass nn 
    %straight to the function if it has already been calculated
   nn = Kap1; totfreq = horzcat(vv{:}).';
else
    
nn = zeros(sum(numwithn),tot_poles); %if Kap1 <inf then less elements may be present
totfreq = horzcat(vv{:}).';
count1=0; tmp3 = 0;
for k =2:Kap2+1

        count2 = count1 + size(tmp3,1); 
        
        tmp3 = zeros(length(count1+1:count2)*size(nn,2),size(nn,2));
       for j = count1+1:count2
        tmp = repmat(nn(j,:),size(nn,2),1); %select element from tier below
        tmp2 = tmp+eye(size(nn,2)); %add one to every position
        tmp3(1+(j-(count1+1))*size(nn,2):(j-count1)*size(nn,2) ,:) = tmp2;       
       end
       %remove any elements of tmp3 which fail to satisfy the high
       %frequency cut off condition
       lg = true(size(tmp3,1),1);
       if isfinite(Kap1)
            for j =1:length(count1+1:count2)*size(nn,2)
                lg(j) = tmp3(j,:)*real(totfreq) < Kap1;
            end
       end
       tmp3 = tmp3(lg,:);
       %now remove duplicate entries
       jj=1;
       while jj < size(tmp3,1)-1
           lg = any(tmp3(jj+1:end,:)  ~= repmat(tmp3(jj,:),size(tmp3,1)-jj,1),2);
           tmp3 = tmp3([true(jj,1);lg],:);
           jj=jj+1;
       end
       count1 = count2;
       nn(count1+1:count1 + size(tmp3,1),:) = tmp3;
end
end

%would be good to include a complicated truncation condition such that
%weakly coupled modes could be taken only to lower orders

%clean any nn values which are all zero besides the first
nn = nn([true;any(nn(2:end,:)~=0,2)],:);
tierofnn = sum(nn,2);
for k  =2:Kap2+1    
    numwithn(k) = sum(tierofnn == (k-1)); %set actual number included
end
else
   nn=0; numwithn=1; %trivial case, DON'T ENTER Kap1 as nn and also Kap2 =0;
end

%% Calculate the coupling from the heirarchy

%general up and down coupling terms
coup_com = zeros(size(nn,1));  %coefficient will be the sqrt of number it couples to (upward)
%coup_p1_save= zeros(size(nn,1),size(nn,1),N);  %takes up too much memory

coup_com_save{N}= sparse(coup_com);

%must also have the term that accounts for the fact that down coupling
%terms are different except for the matsubara terms
%coup_m1 = coup_com ;
coup_m1_acom = coup_com ;  %anti commutator part
coup_acom_save= coup_com_save;

rng_j =0; %full_rng = 1:totpoles; 

for j = 1:N
    if ~isempty(cc_com{j}) %no bath at this mode for some reason
    rng_j = rng_j(end)+1:length(cc_com{j})+rng_j(end); 
    rng_rest = [1:rng_j(1)-1 , rng_j(end)+1:totpoles];
    cc = cc_com{j};  ccI = cc_acom{j};
for k =1:sum(numwithn(1:end-1))
    
    currenttier = tierofnn(k);
    tierabove = currenttier+1;

    tierpickr = abs(tierofnn-tierabove)<eps(10);

    nnsec = nn(k,rng_j); %section of interest
    nnsec2 = nn(k,rng_rest); %rest
    %which obviously must still match, but must be element of j section
    %changed

    temp0 = repmat(nnsec,numwithn(tierabove+1),1); 
    %make same size as nn in tier above in dimension 1
    temp02 = repmat(nnsec2,numwithn(tierabove+1),1); 
    
    temp = temp0 - nn(tierpickr,rng_j); 
    %take away elements in the tier above

    temp4 = temp02 - nn(tierpickr,rng_rest);
    
    temp2 = sum(temp~=0,2) ==1; %only one element +/- 1 diff
    
    temp3 = sum(temp,2) < -1/2 & sum(temp,2) > -3/2; %ones with a -1 net 
    temp4 = all(temp4==0,2); %match rest of the terms
    
    comb_lg = temp2 & temp3 & temp4; %full condition required
    tierpickr(tierpickr) = comb_lg ; 
    %use to put elements in the right row position in coupling matrix,
    %comb_lg is displaced by the number of tiers below
    
    if any(comb_lg)
        
    [~,temp44] = find(temp(comb_lg,:)); %find row position, i.e. which coefficient
    temp7 = sqrt(sqrt(abs(cc(temp44)).^2+abs(ccI(temp44)).^2));
    coup_com(k,tierpickr)=temp7.*sqrt(1+nnsec(temp44));
    
    %elements with minus are in the positions which are the transpose of this 
    %due to the rescaling the coefficients should be the same up and down
    %for the matsubara terms and the commutator part of the drude terms
    %now second down coupling operator
    
    coup_com(tierpickr,k) = cc(temp44).*sqrt(1+nnsec(temp44))./temp7;%...
       % ./ sqrt(sqrt(abs(cc(j,temp44)).^2+abs(ccI(j,temp44)).^2));
     %  if any(cc_acom{j}(temp44))
    coup_m1_acom(tierpickr,k) = ccI(temp44).*sqrt(1+nnsec(temp44))./temp7;%...
      %  ./ sqrt(sqrt(abs(cc(j,temp44)).^2+abs(ccI(j,temp44)).^2));    
     %  end
    end

end
    %Note I am scaling things w.r.t. (abs(a)^2+abs(b)^2)^(1/4)
    %for a term with operator a V^O + b V^X
    
    end

coup_com_save{j} = sparse(coup_com); coup_com = 0*coup_com;
coup_acom_save{j}= sparse(coup_m1_acom); coup_m1_acom=coup_com;
  
end
%down coupling from real frequency and residue modes is just  coup_p1.'


% calculate the sum_{j=1}^N sum_{k=0}^Kappa n_{jk} v_{ik}
if size(nn,1)==1
  const_factor = 0; 
else
    %reshape(vv.',numel(vv),1)
    %([vv(1,:).';vv(2,:).';vv(3,:).';vv(4,:).'])
    
const_factor = nn*totfreq;%([vv(1,:).';vv(2,:).']);
end

%%  Calculate the Uber operator that propogates the entire thing
total_prop_op = sparse(numel(Htot)*size(nn,1),numel(Htot)*size(nn,1));
   %plus_term = total_prop_op;  minus_term = total_prop_op; 

        for j = 1:N
            
    Qjsec = zeros(N2); lg = fock_space_rep(:,j); %states with excitation on site j
    Qjsec(lg,lg) = 1;   Qjsec = diag(diag(Qjsec)); %only want diagonal ones
    Qjsec = kron(Qjsec,eye(viblvls));
    eyeQ = eye(length(Qjsec)); 
    Qj = -1i*sparse(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    %keep factor of -1i from -1i tilde(V)^X term
    Uj =  sparse(kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ)); %anti commutator
    %factor of -1i cancels with +1i in front of imag part
    % plus_term = plus_term + kron(sparse(coup_p1_save{j}),Qj);
    %minus_term = minus_term + kron(sparse(coup_m1_acom_save{j}),Uj) ...
    %              +  kron(sparse(coup_m1_save{j}),Qj);
    
    total_prop_op = total_prop_op + kron(sparse(coup_com_save{j}),Qj);
    total_prop_op = total_prop_op + kron(sparse(coup_acom_save{j}),Uj);
        end

    total_prop_op = total_prop_op + kron(eye(size(nn,1)),sparse(Ltot));
    temp = kron(const_factor,ones(numel(Htot),1));
    %identical constant factor 
    total_prop_op = total_prop_op - sparse(1:length(temp),1:length(temp),temp);
 