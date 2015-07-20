function [H_ex_vib,fock_space_rep,mu_ex,H_exciton,indiv_op] = ...
    generate_ex_vib_ham(om_dip,om_vib,numvib,displ,mu,pdm_n_pos) 
% Generates the total Hamiltonian for a system with Nsites with
% N2 vibrational modes with frequencies om_vib(1,...,N2), this includes for
% each of these num_vib (vector, can be diff for each vib) modes and disp
% is the equilibrium displacement, length (N2, 1 + N + N(N-1)/2)
% pdm and pos contains for the perminant dipole moment and position (a 6
% by N vector) for the excited states
% sum_m E_m B^dag_m B_m + sum_{m neq n} J_{nm} B^dag_m B_n + ...
%               1/2*sum_{m neq n} K_{nm} B^dag_m B^dag_n B_m B_n
% mu are the new resulting dipole moments in the new basis,
% mu_ex =  N +[ N(N-1)/2] by 3
% Fock space rep is the representation of the excitation on each
% Chromophore in Hilbert space
mu_ex_sym =[];
if nargin <= 5
    %only 1 level exciton Ham
    N = length(om_dip);
    dip_size = 1+N;
    flag = false; 
    %flag is commonly used as a colourmap in matlab, don't use this in future
    pdm_n_pos =[];
else
    N = length(om_dip);
    dip_size = 1+N+N*(N-1)/2; flag = true;
    if isempty(pdm_n_pos) %don't know perm dipole moments so dont include
        pdm_n_pos = ones(6,N);
        pdm_n_pos(1:3,:) = 0;
        pdm_n_pos(4:6,:) = cumsum(pdm_n_pos(4:6,:),2);
    end
end
% Firstgenerate two level excition hamiltonian

H_dip = zeros(dip_size); cnt = 0; 
fock_space_rep = false(dip_size,N); 
fock_space_rep(2:N+1,:) = logical(diag(ones(N,1)));

if isempty(mu)
    mu = zeros(N,3); mu(:,1)  =1;
end

if length(om_dip) == length(H_dip)-1;
 H_dip(2:end,2:end) = om_dip;   %can pass om_dip as whole site basis ham
 for j = 1:N
     rng = j+1:N;
     for k = 1:length(rng)
         cnt = cnt+1; kk = rng(k);
         fock_space_rep(N+1+cnt,j) = 1; fock_space_rep(N+1+cnt,kk) = 1;
     end
 end
else
    
H_dip(2:N+1,2:N+1) = om_dip;


if flag %include 2 exciton manifold
    for j = 1:N
        rng = j+1:N; %cnt2 = cnt;
    for k = 1:length(rng)
        cnt = cnt+1;  kk = rng(k);
        state_shift = dot(pdm_n_pos(1:3,kk),pdm_n_pos(1:3,j))-...
                        3*dot(pdm_n_pos(1:3,kk),pdm_n_pos(4:6,j))*...
                        dot(pdm_n_pos(1:3,j),pdm_n_pos(4:6,kk))/...
                         norm(pdm_n_pos(4:6,j)-pdm_n_pos(4:6,kk))^2;
        state_shift = state_shift/norm(pdm_n_pos(4:6,j)-pdm_n_pos(4:6,kk))^3  ;
        %this is the shift from the perminant dipole moments of the excited
        %states, if this is not known just set all pdm(1:3,:) to zero,rest
        %to one
       H_dip(N+1+cnt,N+1+cnt)= om_dip(j,j)+om_dip(kk,kk)+state_shift;
       fock_space_rep(N+1+cnt,j) = 1; fock_space_rep(N+1+cnt,kk) = 1; 
       %also J_nm B^dag_m B_n+c.c. can mix state |n;l> to |m;l> etc
       %i.e. all double exciton state with an excitation on oone other site
%        H_dip(N+1+cnt,N+1+cnt2+rng(k+1:end)) = om_dip(j,rng(k+1:end));
%        H_dip(N+1+cnt,N+1+cnt2+rng(1:k-1)) = om_dip(j,rng(1:k-1));
%        H_dip(N+1+cnt2+rng(k+1:end),N+1+cnt) = conj(om_dip(j,rng(k+1:end)));
%        H_dip(N+1+cnt2+rng(1:k-1),N+1+cnt) = conj(om_dip(j,rng(1:k-1)));
    end
    end
    
    lg =  sum(double(fock_space_rep),2)==2; tmp = 1:N;
    rnj = 1:size(fock_space_rep,1); rnj=rnj(lg); 
    for k = 1:length(rnj)
        lg1 = fock_space_rep(rnj(k),:); s1 = find(lg1);
        for kk = k+1:length(rnj)
            lg2 = fock_space_rep(rnj(kk),:); 
            lg3 = lg2 & lg1; %shared element
            if any(lg3)
                s2 = find(lg2); s3 = s2(s2~=tmp(lg3)); s4 = s1(s1~=tmp(lg3));
                H_dip(rnj(k),rnj(kk)) = om_dip(s3,s4);
                H_dip(rnj(kk),rnj(k)) =H_dip(rnj(k),rnj(kk))';
            end            
        end
    end
end
end

   [proj_sm,~]=eig(H_dip);

%generate the transition dipole moments in the exciton basis
mu_sym = sym('mu',[N,1]); 
mu_ge = zeros(N,3); mu_ge_sym = sym(zeros(N,1));
for k = 1:N
   for ell = 1:N
    mu_ge(k,:) = mu_ge(k,:) + proj_sm(k+1,ell+1)*mu(ell,:);
    mu_ge_sym (k) = mu_ge_sym(k,:) + proj_sm(k+1,ell+1)*mu_sym(ell);
   end
end

if flag 
    mu_ef = zeros(N*(N-1)/2,N,3); %mu_ef(a, b,dim) = transition dipole 
    mu_ef_sym = sym(zeros(N*(N-1)/2,N));
    %from state b to double excited state with label a, dim is xyz 
    cq_nm = zeros(N*(N-1)/2 ,N); %participation of each site in double exciton state 
for f = 1:N*(N-1)/2
    for f2 = 1:N*(N-1)/2 
    cq_nm(f,:)= cq_nm(f,:) + proj_sm(1+N+f,1+N+f2)*double(fock_space_rep(1+N+f2,:)); 
    %effective occupation of each site in the double exciton
    end
    
    for k = 1:N  
    mu_ef(f, k,:) = cq_nm(f,k)*( proj_sm(k+1,2:N+1) * mu(1:N,:) );
    mu_ef_sym(f, k) = cq_nm(f,k)*( proj_sm(k+1,2:N+1) * mu_sym(1:N) );
    end
end  
mu_ex = {mu_ge,mu_ef};
mu_ex_sym = {mu_ge_sym,mu_ef_sym};
else
    mu_ex = mu_ge;
end
% Next generate vibrational hamiltonian and displacement op term
%note this assumes the modes found are the same on all sites
if isempty(om_vib) %no vibrations included
H_ex_vib = H_dip;
H_vib = 0; H_coup = [];
else
N2 = length(om_vib);
if any(numvib==0)
    warning('numvib should not be set to zero, this is will be set to unity')
end
numvib(numvib==0) = 1; %don't put zeros in this

H_vib = zeros(prod(numvib)); %do not set any numvib to zero
% vib_fock = zeros(length(H_vib),length(numvib));
% cnt=0;
% for j=1:N2
%     vib_fock(cnt+1:cnt+numvib(j),j) = 0:numvib(j)-1;
%     
% end
H_coup = zeros(prod(numvib),prod(numvib),dip_size  );

for j = 1:N2 %loop over number of vibrations
 
    %tmp = kron(eye(prod(numvib(1:j-1))),...
    %    diag(1/2:1:numvib(j)-1/2)*om_vib(j));
    %tmp =  kron(tmp , eye(prod(numvib(j+1:N2))));
    tmp = kron(diag(1/2:(numvib(j)-1/2))*om_vib(j),eye(prod(numvib(1:j-1))));
    tmp = kron(eye(prod(numvib(j+1:N2))),tmp);    

    H_vib = H_vib + tmp;
    
    for ell = 1:dip_size %loop over all included states
        % om(j)* (Delta(j,ell1,ell2)^2/2 - Delta(j,ell1,ell2)(b_j+b_j^dagger)/sqrt(2)) 
        %=(om(j)*Delta(j,ell1,ell2)/sqrt(2))*( (Delta(j,ell1,ell2)/sqrt(2)) 
        %-(b_j+b_j^dagger) )
        bj = diag(sqrt(1:numvib(j)-1),1);
    tmp = kron((displ(j,ell)*eye(numvib(j))/2 -(bj+bj')/sqrt(2))...
           *om_vib(j)*displ(j,ell),eye(prod(numvib(1:j-1))));
    tmp =  kron(eye(prod(numvib(j+1:N2))),tmp );        
    H_coup(:,:,ell) = H_coup(:,:,ell) + tmp;    
    end
    
    
end
%combine together

H_ex_vib = kron(H_dip,eye(size(H_vib))) +  kron(eye(size(H_dip)),H_vib);

   for ell = 1:dip_size 
        tmp = zeros(size(H_dip)); tmp(ell,ell) = 1;
       H_ex_vib = H_ex_vib + kron(tmp,H_coup(:,:,ell));
       %if values of H_dip are those from spectroscopy then don't take the
       %diagonal part of this matrix
   end
    
end
if nargout >=4
   [proj,H_ex]=eig(H_dip);

       proj = kron(proj,eye(size(H_vib)));
       H_exciton = proj'*H_ex_vib*proj;
end

if nargout >=5
   indiv_op{1} =  H_dip;
   indiv_op{2} =  H_vib;
   indiv_op{3} =  H_coup;
   indiv_op{4} =  H_ex;
   indiv_op{5} =  proj;
   indiv_op{6} = mu_ex_sym;
   indiv_op{7} = proj_sm(2:N+1,2:N+1);
   indiv_op{8} = proj_sm(N+2:end,N+2:end);
end
