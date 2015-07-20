function [ Qtrunc_set,H_prop_op_set,nn,nn_alt, sec_ops,Qtrunc_set2,...
            H_prop_op_set2,H_prop_full,debug]...
        =HEOM_propogator_sep(QQ,cc_com,cc_acom,vv,max_tier,occ_mat,viblvls)
% Calculates the propogator for the Heirarchy and also truncation
% corrections
% Calculates seperately the propogator for the ground state manifold, the
% coherences between the single and double excited states, the excited
% state manifold and the coherences between single and double excited
% states
if nargout >=6 %also gives two further matricies which are the same things
% acting on an operator (to the left)
calc_2 = true;
else
    calc_2 = false;
end

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
%occ_mat is the matrix of a_j^dag a_j, should be diagonal.  If you want to
%work in a different basis this could in principle still work.
% max_tier e.g.= 4; %truncation parameter, there are Kappa + poles in J(om)
%explicitly treated exponentials, this number is the maximum order the sum of the
 %number of powers can take, e.g. n_1 + n_2 +... + n_end <= max_tier
 % numpoints is the number of time points to save of the solved equation
 % for all the auxilary density matricies

totpoles = sum(cellfun(@length,cc_com));
if viblvls == 0
    viblvls = 1; %set to unity in order to avoid bugs
end
%% Operators which pick sections of particular operators
sz_full = 1+N+N*(N-1)/2; sz_full = sz_full * viblvls;

    tmpg = zeros(sz_full); tmpg(1:viblvls,1:viblvls)=1;
	tmpg = logical(reshape(tmpg,sz_full^2,1)); %picks this out

%reduced operators acting only on ground excited coherences, these don't
%mix to p_gg or p_ee'

    tmpge  = zeros(sz_full); tmpge(1:viblvls,viblvls+1:viblvls*(N+1))=1; %upper diag
	tmpge = logical(reshape(tmpge,sz_full^2,1));
    
    tmpeg  = zeros(sz_full);tmpeg(viblvls+1:viblvls*(N+1),1:viblvls)=1; %lower diag
	tmpeg = logical(reshape(tmpeg,sz_full^2,1));

%reduced operators acting only on 1st ex state manifold, 

    tmpe  = zeros(sz_full); tmpe(viblvls+1:viblvls*(N+1),viblvls+1:viblvls*(N+1))=1;
	tmpe = logical(reshape(tmpe,sz_full^2,1));
    
%reduced operator acting on gs-double ex coherence

    tmpgf  = zeros(sz_full); tmpgf(1:viblvls,viblvls*(N+1)+1:end)=1;
	tmpgf = logical(reshape(tmpgf,sz_full^2,1));
    
    tmpfg  = zeros(sz_full); tmpfg(viblvls*(N+1)+1:end,1:viblvls)=1;
	tmpfg = logical(reshape(tmpfg,sz_full^2,1));    
    
%reduced operators acting only on 1st-2nd ex state manifold, 

    tmpef  = zeros(sz_full); tmpef(viblvls+1:viblvls*(N+1),viblvls*(N+1)+1:end)=1; %upper diag
	tmpef = logical(reshape(tmpef,sz_full^2,1));
    
    tmpfe  = zeros(sz_full); tmpfe(viblvls*(N+1)+1:end,viblvls+1:viblvls*(N+1))=1; %lower diag
	tmpfe = logical(reshape(tmpfe,sz_full^2,1));


%% Construct trunc correction

% For N X N A and B we have the following
% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho

if size(QQ,1) == 1 %assume same for every site
    QQ = repmat(QQ,N*viblvls,1);
end

   Qtrunc = zeros(size(occ_mat,1)^2*viblvls^2);  
   if calc_2
   Qtrunc2 = zeros(size(occ_mat,1)^2*viblvls^2);
   end

for j = 1:N
    
    tmp = occ_mat(:,:,j);  %pick site
    Qj = kron(tmp,eye(viblvls)); %add vib lvls
   eyeQ = sparse(1:length(Qj),1:length(Qj),ones(length(Qj),1));

    Qtrunc = Qtrunc + QQ(j,1).*(kron(Qj.',eyeQ) + kron(eyeQ,Qj) - 2*kron(Qj.',Qj)); 
    Qtrunc = Qtrunc + QQ(j,2).*(kron(Qj.',eyeQ) - kron(eyeQ,Qj));    
% Qtrunc = Qtrunc + QQ(j,1).*(kron(Qj.',eyeQ) - kron(eyeQ,Qj))^2; 
% Qtrunc = Qtrunc + QQ(j,2).*(kron(Qj.',eyeQ) - kron(eyeQ,Qj))*(kron(Qj.',eyeQ) + kron(eyeQ,Qj));     
      if calc_2
    Qtrunc2 = Qtrunc2 + QQ(j,1).*(kron(Qj.',eyeQ) + kron(eyeQ,Qj) - 2*kron(Qj.',Qj)); 
    Qtrunc2 = Qtrunc2 - QQ(j,2).*(kron(Qj.',eyeQ) - kron(eyeQ,Qj));    %commutator has -sgn
% Qtrunc2 = Qtrunc2 + QQ(j,1).*(kron(Qj.',eyeQ) - kron(eyeQ,Qj))^2; 
 %Qtrunc2 = Qtrunc2 + QQ(j,2).*(kron(Qj.',eyeQ) - kron(eyeQ,Qj))*(kron(Qj.',eyeQ) + kron(eyeQ,Qj));     
      end
end
Qtrunc = sparse(Qtrunc);
%Qtrunc_00 is trivially zero
Qtrunc_01 =  Qtrunc(tmpeg,tmpeg); %N*viblvls^2
Qtrunc_10 =  Qtrunc(tmpge,tmpge);
Qtrunc_02 =  Qtrunc(tmpfg,tmpfg); %N*(N-1)/2*viblvls^2
Qtrunc_20 =  Qtrunc(tmpgf,tmpgf);
Qtrunc_11 =  Qtrunc(tmpe,tmpe); %N^2*viblvls^2
Qtrunc_12 =  Qtrunc(tmpfe,tmpfe); %(N*N*(N-1)/2)^2*viblvls^2

Qtrunc_set = {Qtrunc_01,Qtrunc_02,Qtrunc_11,Qtrunc_12,Qtrunc_10,Qtrunc_20...
             'Qtrunc_01,Qtrunc_02,Qtrunc_11,Qtrunc_12,Q_trunc_10,Q_trunc_20'}; 
         %also string saying what ops are

clear Qtrunc_01 Qtrunc_02 Qtrunc_11 Qtrunc_12 Qtrunc

if calc_2

    Qtrunc2 = sparse(Qtrunc2);

Qtrunc_01 =  Qtrunc2(tmpge,tmpge); %N*viblvls^2
Qtrunc_02 =  Qtrunc2(tmpgf,tmpgf); %N*(N-1)/2*viblvls^2
Qtrunc_11 =  Qtrunc2(tmpe,tmpe); %N^2*viblvls^2
Qtrunc_12 =  Qtrunc2(tmpef,tmpef); %(N*N*(N-1)/2)^2*viblvls^2

Qtrunc_set2 = {Qtrunc_01,Qtrunc_02,Qtrunc_11,Qtrunc_12,...
             'Qtrunc_10,Qtrunc_20,Qtrunc_11,Qtrunc_21'}; 
         %also string saying what ops are

clear Qtrunc_01 Qtrunc_02 Qtrunc_11 Qtrunc_12 Qtrunc2
    
end
%% Construct labels for auxilliary density matrix


cnt = 0; sz = zeros(N,1);
for lp = 1:N
    cc = cc_com{lp}; sz(lp) = length(cc); ccI = cc_acom{lp}; 
    %scale_fct{lp} = sqrt(abs(cc));   %#ok<*AGROW> %use other factor to
    %compare with the old code
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
    
for k = 2:max_tier %this part does not seem to calculate the coupling correctly
   
    tier_below = nn_current;  tier_below_alt = nn_alt_current;
    nn_current = zeros(numwithn(k+1),totpoles); %current tier 
    nn_alt_current = zeros(numwithn(k+1),k);
    
    lst = 0; %set last value to zero
    tmp3 = zeros(size(nn_current,1),1); %will contain hashfn indexing
    for lp = 1:N %site loop index

           coup_com = zeros(sum(numwithn(1:k+1))); coup_acom = coup_com;
           sc1 = scale_fct{lp};
           
        for j =1:size(tier_below,1)
            %enumerate with all up couplings
            
            sc2 = sqrt(1+tier_below(j,coup_s_alt{lp})); %n scale factor
             
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
             rng2 = max_tier_below + rng; %new elements in tier above
             j2 = max_tier_2below+j;  %element in this tier              
 
coup_com(j2,rng2) = sc1(~ismem).*sc2(~ismem); %sqrt(c_k *(n+1))
%also downcouples back, with an anticommutator part as well
coup_com(rng2,j2) = cc_com{lp}(~ismem)./sc1(~ismem).*sc2(~ismem);
coup_acom(rng2,j2) = cc_acom{lp}(~ismem)./sc1(~ismem).*sc2(~ismem);                       
%n factor is same for the down coupling as n is one larger
             if any(ismem) %one that were members already
                 
                 rng3 = max_tier_below+pos_found(ismem); 

coup_com(j2,rng3) = sc1(ismem).*sc2(ismem);
coup_com(rng3,j2) = cc_com{lp}(ismem)./sc1(ismem).*sc2(ismem);
coup_acom(rng3,j2) = cc_acom{lp}(ismem)./sc1(ismem).*sc2(ismem);                      

             end
             
lst = lst +size(tmp,1);   %change last elements          
             
        end

                         
  [ii,jj,ss] = find(coup_com);    S = sparse(ii,jj,ss,size(nn,1),size(nn,1));
coup_com_save{lp} = coup_com_save{lp} + S;
  [ii,jj,ss] = find(coup_acom);    S = sparse(ii,jj,ss,size(nn,1),size(nn,1));
coup_acom_save{lp} =  coup_acom_save{lp} + S;                             

    end
    max_tier_below = max_tier_below + lst;
    max_tier_2below  = max_tier_below;
    tot_rng = sum(numwithn(1:k))+1:sum(numwithn(1:k))+numwithn(k+1);
    nn(tot_rng,:) = nn_current; 
    nn_alt(tot_rng,:) = [nn_alt_current,zeros(length(tot_rng),max_tier-k)];
end 
else
    for j=1:N
        coup_com_save{j} = 0; coup_acom_save{j} =0;
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

if nargout >=9 
    debug = {coup_com_save,coup_acom_save};
end

%%  Calculate the Uber operator that propogates the entire thing

   %plus_term = total_prop_op;  minus_term = total_prop_op; 

        for j = 1:N
    tmp = occ_mat(:,:,j);  %pick site
    Qjsec = sparse(kron(tmp,eye(viblvls))); %add vib lvls
    eyeQ = sparse(1:length(Qjsec),1:length(Qjsec),ones(length(Qjsec),1));
   
    Qj = -1i*sparse(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    %keep factor of -1i from -1i tilde(V)^X term
    Uj = sparse(kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ)); %anti commutator
    %factor of -1i cancels with +1i in front of imag part assumed in coeffs
    %from brownian new function
    % plus_term = plus_term + kron(sparse(coup_p1_save{j}),Qj);
    %minus_term = minus_term + kron(sparse(coup_m1_acom_save{j}),Uj) ...
    %              +  kron(sparse(coup_m1_save{j}),Qj);
    
    if j==1
     H_prop_op =  kron(sparse(coup_com_save{j}),Qj);   
    else
    H_prop_op = H_prop_op + kron(sparse(coup_com_save{j}),Qj);
    end
    H_prop_op = H_prop_op + kron(sparse(coup_acom_save{j}),Uj);
    if calc_2 %acting to the left, commutator reverses sign
    if j==1
     H_prop_op2 =  -kron(sparse(coup_com_save{j}),Qj);   
    else
    H_prop_op2 = H_prop_op2 - kron(sparse(coup_com_save{j}),Qj);
    end
    H_prop_op2 = H_prop_op2 + kron(sparse(coup_acom_save{j}),Uj);        
        
    end
        end

    %H_prop_op = H_prop_op + kron(eye(length(rho_vec)/numel(Htot)),sparse(Ltot));

    temp = kron(const_factor,ones(sz_full^2,1));

    H_prop_op = H_prop_op - sparse(1:length(temp),1:length(temp),temp);
    %take out the correct bits of this    
    tmp = repmat(tmpeg,[size(nn,1),1]); tmp2 = repmat(tmpge,[size(nn,1),1]);
    H_prop_op_01 =  H_prop_op(tmp ,tmp); %N*viblvls^2
    H_prop_op_10 =  H_prop_op(tmp2 ,tmp2);
    tmp = repmat(tmpfg,[size(nn,1),1]); tmp2 = repmat(tmpfg,[size(nn,1),1]);
    H_prop_op_02 =  H_prop_op(tmp ,tmp); %N*(N-1)/2*viblvls^2
    H_prop_op_20 =  H_prop_op(tmp2 ,tmp2);     
    tmp = repmat(tmpe,[size(nn,1),1]);
    H_prop_op_11 =  H_prop_op(tmp ,tmp); %N^2*viblvls^2
    tmp = repmat(tmpfe,[size(nn,1),1]);
    H_prop_op_12 =  H_prop_op(tmp ,tmp); %(N*N*(N-1)/2)^2*viblvls^2
    if nargout >= 8 
        H_prop_full  = H_prop_op; %output everything
    end
    clear H_prop_op %just in case out of mem
    H_prop_op_set = {H_prop_op_01,H_prop_op_02,H_prop_op_11,H_prop_op_12,...
                    H_prop_op_10,H_prop_op_20,...
                    'H_prop_op_01,H_prop_op_02,H_prop_op_11,H_prop_op_12,H_prop_op_10,H_prop_op_20'};
    sec_ops = {tmpg,tmpge,tmpeg,tmpe,tmpgf,tmpfg,tmpef,tmpfe,...
                'tmpg,tmpge,tmpeg,tmpe,tmpgf,tmpfg,tmpef,tmpfe'};
if calc_2
 
    H_prop_op = H_prop_op2 - sparse(1:length(temp),1:length(temp),temp);
    clear H_prop_op2 %don't need this anyway
    %take out the correct bits of this     tmpge
    tmp = repmat(tmpge,[size(nn,1),1]);
    H_prop_op_01 =  H_prop_op(tmp ,tmp); %N*viblvls^2
    tmp = repmat(tmpgf,[size(nn,1),1]);
    H_prop_op_02 =  H_prop_op(tmp ,tmp); %N*(N-1)/2*viblvls^2
    tmp = repmat(tmpe,[size(nn,1),1]);
    H_prop_op_11 =  H_prop_op(tmp ,tmp); %N^2*viblvls^2
    tmp = repmat(tmpef,[size(nn,1),1]);
    H_prop_op_12 =  H_prop_op(tmp ,tmp); %(N*N*(N-1)/2)^2*viblvls^2
    clear H_prop_op %just in case out of mem
    H_prop_op_set2 = {H_prop_op_01,H_prop_op_02,H_prop_op_11,H_prop_op_12,...
                    'H_prop_op_10,H_prop_op_20,H_prop_op_11,H_prop_op_21'};   
    
end
    