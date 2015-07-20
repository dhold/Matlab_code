function [ Qtrunc,H_prop_op,nn,nn_alt]=HEOM_propogator2...
            (QQ,cc_com,cc_acom,vv,max_tier,viblvls)
%Calculates the propogator for the Heirarchy and also truncation
%corrections
% Calculates seperately the propogator for the ground state manifold, the
% coherences between the single and double excited states, the excited
% state manifold and the coherences between single and double excited
% states


% Unit setup
%All terms in inverse cm, conversion to sensible units achieved by
%E = hc / lamba => 6.62606957 *10^{-33} 299792458 *10 * "E in code" gives
%
%if nargin == 0  %use a select set nargin
%%
 % LL is the total Louivillian for the system of N chromophores and 
 % whatever vibrations are treated quantum mechanically.
 % QQ is the extra convergence parameter made up of all the terms of the 
 % form int_0^t d tau exp(-v (t-tau)) * ( c_1 V^X(tau) + c_2 V^O(tau) )
 %  which are evaluated via exp(-v (t-tau)) ~ delta(t-tau)/v, as v->inf
 % and are thus just markovian terms.
 % cc1 = residues at all the poles of 1/(1-e^{-beta omega}) 
 % cc2R = residues at all the poles of J(omega), contribution to V^x
 % cc2I = residues at all the poles of J(omega), contribution to V^o
 % i.e. the anti commutator part
 % each of these is paired with a frequency
 % vv1 = position of poles of 1/(1-e^{-beta omega})  times -1i
 % vv2 = position of poles of J(omega)  times -1i
 % N is the number of sites / Chromophores
 N = length(cc_com); 

% max_tier e.g.= 4; %truncation parameter, there are Kappa + poles in J(om)
%explicitly treated exponentials, this number is the maximum order the sum of the
 %number of powers can take, e.g. n_1 + n_2 +... + n_end <= max_tier
 % numpoints is the number of time points to save of the solved equation
 % for all the auxilary density matricies

totpoles = sum(cellfun(@length,cc_com));


%% Construct trunc correction

% For N X N A and B we have the following
% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho

if size(QQ,1) == 1 %assume same for every site
    QQ = repmat(QQ,N*viblvls,1);
end
 Qtrunc = zeros(N^2*viblvls^2);

for j = 1:N
    
    Qj = zeros(N);
    Qj(j,j) = 1;  Qj = kron(Qj,eye(viblvls));eyeQ = eye(length(Qj));
    
      %note Qj^2 = Qj so I don't need to worry about this

    Qtrunc = Qtrunc + QQ(j,1).*(kron(Qj.',eyeQ) + kron(eyeQ,Qj) - 2*kron(Qj.',Qj)); 
    Qtrunc = Qtrunc  + QQ(j,2) .*(kron(Qj.',eyeQ) - kron(eyeQ,Qj)); 
    
    %reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
    %[Q_i,[Q_i,rho_n]] = (1 X Q_i*Q_i + Q_i^T*Q_i^T X 1 - 2 Q_i^T X Q_i ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
    
    % .*(kron(Vj.',eyeQ) + kron(eyeQ,Vj)); 
end



%% Construct labels for auxilliary density matrix


cnt = 0; sz = zeros(N,1);
for lp = 1:N
    cc = cc_com{lp};  ccI = cc_acom{lp}; sz(lp) = length(cc);
    scale_fct{lp} = sqrt(sqrt(abs(cc).^2+abs(ccI).^2));  
    cnt2 = cnt+sz(lp);
    coup_s_j{lp} = sparse(1:cnt2-cnt,cnt+1:cnt2,ones(cnt2-cnt,1),cnt2-cnt,totpoles);
    coup_s_alt{lp} = (cnt+1:cnt2)';
    cnt = cnt2;
end

hash_fn = @(state) sum(state.^(repmat(1:size(state,2),size(state,1),1)),2);
% index_fn = (1:max_tier).^(; %hash fn->unique identifier for states
% index_fn  = reshape(index_fn ,length(index_fn ),1);

numwithn = zeros(max_tier+1,1);
for kk = 0:max_tier
numwithn(kk+1) = nchoosek(totpoles-1+kk,kk);
%irreducible nth power rep of sym group U(N*(Kappa+size(cc2,2))) 
end
nn = zeros(sum(numwithn),totpoles); %list of indexes for states
nn_alt = zeros(sum(numwithn),max_tier);%alternative form

tot_couplings = ceil(sum(numwithn)*totpoles/(N-1));

if max_tier>=1 %enumerate first tier and couplings
    lst = 1; %last value
    for j = 1:N %calculate elements coupling to tier above (those coupling down from below)
 
    rng = lst+1:lst+sz(j); %rng to give these
    nn(rng,:) = coup_s_j{j}; % enumerate
    nn_alt(rng,1) = coup_s_alt{j};       
    %sqrt(n+1)=1 here as n = 0 for all

coup_com_save{j} = sparse(1,rng,scale_fct{j},size(nn,1),size(nn,1),tot_couplings)+...
                    sparse(rng,1,cc_com{j}./scale_fct{j},size(nn,1),size(nn,1),tot_couplings);
coup_acom_save{j} =  sparse(rng,1,cc_acom{j}./scale_fct{j},size(nn,1),size(nn,1),tot_couplings);
    lst = lst + sz(j);
    end
    nn_current = nn(2:lst,:); max_tier_below = lst; max_tier_2below = 1;
    nn_alt_current = nn_alt(2:lst,:);
    
for k = 2:max_tier
   
    tier_below = nn_current;  tier_below_alt = nn_alt_current;
    nn_current = zeros(numwithn(k+1),totpoles); %current tier 
    nn_alt_current = zeros(numwithn(k+1),k);
    
    lst = 0; %set last value to zero
    tmp3 = zeros(size(nn_current,1),1); %will contain hashfn indexing
    for lp = 1:N %site loop index
       try 
           coup_com = zeros(numwithn(k+1)); coup_acom = coup_com;
           use_sparse = false;
       catch
          coup_com = sparse([],[],[],size(nn,1),size(nn,1),numwithn(k+1)*totpoles);
          coup_acom = coup_com; 
          use_sparse = true; %memory issues force the use of sparse matricies
       end
        for j =1:size(tier_below,1)
            %enumerate with all up couplings
  
            tmp = repmat(tier_below(j,:),size(coup_s_j{lp},1),1) + coup_s_j{lp};
            %some of these will map to a state which already has a label
            %unless this is the first
            
            %do same in alternative form for the hash_fn
            tmp2 = [repmat(tier_below_alt(j,1:k-1),size(coup_s_j{lp},1),1),coup_s_alt{lp}];
            tmp2 = sort(tmp2,2);
  
            tmphash = hash_fn(tmp2);          
   
            [ismem,pos_found] = ismember(tmphash,tmp3(1:lst)); %find elements of tmp2 in tmp3

      %add elements that are not already in the set to nn and hash index     
            tmp = tmp(~ismem,:); tmp2 = tmp2(~ismem,:);

            rng = lst+1:lst+size(tmp,1);
              tmp3(rng) = tmphash(~ismem);          
             nn_current(rng,:) = tmp;   
             nn_alt_current(rng,1:k) = tmp2;
      
      %calculate couplings from the new elements       
             rng2 = max_tier_below + rng; j2 = max_tier_2below+j;
                                  
             
             if use_sparse
coup_com = coup_com + sparse(j2,rng2,scale_fct{lp}(~ismem),size(nn,1),size(nn,1))+...
           sparse(rng2,j2,cc_com{lp}(~ismem)./scale_fct{lp}(~ismem),size(nn,1),size(nn,1));
coup_acom =  coup_acom + sparse(rng2,j2,cc_acom{lp}(~ismem)./scale_fct{lp}(~ismem),size(nn,1),size(nn,1));       
             else
coup_com(j2,rng2) = scale_fct{lp}(~ismem);
coup_com(rng2,j2) = cc_com{lp}(~ismem)./scale_fct{lp}(~ismem);
coup_acom(rng2,j2) = cc_acom{lp}(~ismem)./scale_fct{lp}(~ismem);                        
             end
             if any(ismem)
                 rng3 = max_tier_below+pos_found(ismem);                
             if use_sparse
coup_com = coup_com + sparse(j2,rng3,scale_fct{lp}(ismem),size(nn,1),size(nn,1))+...
           sparse(rng3,j2,cc_com{lp}(ismem)./scale_fct{lp}(ismem),size(nn,1),size(nn,1));
coup_acom = coup_acom+ sparse(rng3,j2,cc_acom{lp}(ismem)./scale_fct{lp}(ismem),size(nn,1),size(nn,1));                                  
            else
coup_com(j2,rng3) = scale_fct{lp}(ismem);
coup_com(rng3,j2) = cc_com{lp}(ismem)./scale_fct{lp}(ismem);
coup_acom(rng3,j2) = cc_acom{lp}(ismem)./scale_fct{lp}(ismem);                        
             end            
             
             end
             
lst = lst +size(tmp,1);   %change last elements          
             
        end
                     if use_sparse
coup_com_save{lp} = coup_com_save{lp} + sparse(coup_com);
coup_acom_save{lp} =  coup_acom_save{lp} + sparse(coup_acom);    
                     else
                         
  [ii,jj,ss] = find(coup_com);    S = sparse(ii,jj,ss,size(nn,1),size(nn,1));
coup_com_save{lp} = coup_com_save{lp} + sparse(S);
  [ii,jj,ss] = find(coup_acom);    S = sparse(ii,jj,ss,size(nn,1),size(nn,1));
coup_acom_save{lp} =  coup_acom_save{lp} + sparse(S);                             
                     end
    end
    max_tier_below = max_tier_below + lst;
    max_tier_2below  = max_tier_below;
    tot_rng = sum(numwithn(1:k))+1:sum(numwithn(1:k))+numwithn(k+1);
    nn(tot_rng,:) = nn_current; 
    nn_alt(tot_rng,:) = [nn_alt_current,zeros(length(tot_rng),max_tier-k)];
end
    
    
    
end

%down coupling from real frequency and residue modes is just  coup_p1.'

% calculate the sum_{j=1}^N sum_{k=0}^Kappa n_{jk} v_{ik}
if size(nn,1)==1
  const_factor = 0; 
else
    %reshape(vv.',numel(vv),1)
    %([vv(1,:).';vv(2,:).';vv(3,:).';vv(4,:).'])
totfreq = horzcat(vv{:}).';    
const_factor = nn*totfreq;%([vv(1,:).';vv(2,:).']);
end

%%  Calculate the Uber operator that propogates the entire thing

   %plus_term = total_prop_op;  minus_term = total_prop_op; 

        for j = 1:N

    Qjsec = zeros(N);  Qjsec(j,j) = 1;
    Qjsec = kron(Qjsec,eye(viblvls));
    eyeQ = eye(length(Qjsec)); 
    Qj = -1i*sparse(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    %keep factor of -1i from -1i tilde(V)^X term
    Uj=  sparse(kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ)); %anti commutator
    %factor of -1i cancels with +1i in front of imag part
    % plus_term = plus_term + kron(sparse(coup_p1_save{j}),Qj);
    %minus_term = minus_term + kron(sparse(coup_m1_acom_save{j}),Uj) ...
    %              +  kron(sparse(coup_m1_save{j}),Qj);
    
    if j==1
     H_prop_op =  kron(sparse(coup_com_save{j}),Qj);   
    else
    H_prop_op = H_prop_op + kron(sparse(coup_com_save{j}),Qj);
    end
    H_prop_op = H_prop_op + kron(sparse(coup_acom_save{j}),Uj);
        end

    %H_prop_op = H_prop_op + kron(eye(length(rho_vec)/numel(Htot)),sparse(Ltot));
    temp = kron(const_factor,ones(N^2*viblvls^2,1));
    H_prop_op = H_prop_op - sparse(1:length(temp),1:length(temp),temp);
    
    