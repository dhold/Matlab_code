function [H_el,H_vib,H_e,H_f,e_1,e_2,fock_space_rep,M_p,M_pfull,mu_ex,mag_ex,H_coup] = ...
    generate_ex_vib_ham2(om_dip,om_vib,numvib,displ,mu,mdip,R,pdm) 
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
% Note that om_dip should include shifts from reorganisation energy etc
if nargin <= 5
    %only 1 level exciton Ham
    N = length(om_dip);
    dip_size = 1+N;
    flag = false; 
    %flag is commonly used as a colourmap in matlab, don't use this in future
    pdm =[];
else
    N = length(om_dip);
    dip_size = 1+N+N*(N-1)/2; flag = true;
end
% Firstgenerate two level excition hamiltonian


fock_space_rep = false(dip_size,N); 
fock_space_rep(2:N+1,:) = logical(diag(ones(N,1)));

H_e = om_dip;
H_f = zeros(N*(N-1)/2); cnt = 0; 

if flag %include 2 exciton manifold
    for j = 1:N
        rng = j+1:N; %cnt2 = cnt;
    for k = 1:length(rng)
        cnt = cnt+1;  kk = rng(k);
        
        if iscell(pdm) %pass cell to just give coefficients
            state_shift = pdm{j,kk}; 
        elseif ~isempty(pdm)
        state_shift = dot(pdm(kk,:),pdm(j,:))-...
                        3*dot(pdm(kk,:),R(j,:))*dot(pdm(j,:),R(kk,:))/...
                         norm(R(j,:)-R(kk,:))^2;
        state_shift = state_shift/norm(R(j,:)-R(kk,:))^3  ;
        else
           state_shift = 0; % assume all zero
        end
        %this is the shift from the perminant dipole moments of the excited
        %states, if this is not known/small just set all pdm(1:3,:) to zero
       H_f(cnt,cnt)= om_dip(j,j)+om_dip(kk,kk)+state_shift;
       fock_space_rep(N+1+cnt,j) = 1; fock_space_rep(N+1+cnt,kk) = 1; 

    end
    end
       %also J_nm B^dag_m B_n+c.c. can mix state |n;l> to |m;l> etc
       %i.e. all double exciton state with an excitation on oone other site    
    fock_f = fock_space_rep(N+2:end,:); tmp = 1:N;
    for k = 1:size(fock_f,1)
        lg1 = fock_f(k,:); s1 = find(lg1);
        for kk = k+1:size(fock_f,1)
            lg2 = fock_f(kk,:); 
            lg3 = lg2 & lg1; %shared element
            if any(lg3)
                s2 = find(lg2); s3 = s2(s2~=tmp(lg3)); s4 = s1(s1~=tmp(lg3));
                H_f(k,kk) = om_dip(s3,s4);
                H_f(kk,k) = H_f(k,kk)';
            end            
        end
    end
end
   [M_1,E_1]=eig(H_e); E_1 = diag(E_1); [M_2,E_2]=eig(H_f); E_2=diag(E_2);
   M_p = {M_1,M_2}; H_el = blkdiag(0,H_e,H_f);
   M_tot = blkdiag(1,M_1,M_2);

%%generate the transition dipole moments in the exciton basis

if ~isempty(mu) %don't generate if no dipole moments passed
    % expand magnetic dipole moments to include the contribution from the 
    % intra chromophore quadrupole moment
    mdip2 = +1i*pi*cross(R,mu); 
   

mu_ge = zeros(N,3);  mag_ge = mu_ge;
for k = 1:N
   for ell = 1:N
    mu_ge(k,:) = mu_ge(k,:) + M_tot(k+1,ell+1)*mu(ell,:);
    mag_ge(k,:) = mag_ge(k,:) + M_tot(k+1,ell+1)*(mdip(ell,:)+ mdip2(j,:)*E_1(j));    
    %approximates lambda to resonant value, f*lambda = c -> f/c = 1/lambda
   end
end

if flag 
    mu_ef = zeros(N*(N-1)/2,N,3); %mu_ef(a, b,dim) = transition dipole 
    mag_ef = zeros(N*(N-1)/2,N,3);
    %from state b to double excited state with label a, dim is xyz 
    cq_nm = zeros(N*(N-1)/2 ,N); %participation of each site in double exciton state 
for f = 1:N*(N-1)/2
    for f2 = 1:N*(N-1)/2 
    cq_nm(f,:)= cq_nm(f,:) + M_tot(1+N+f,1+N+f2)*double(fock_space_rep(1+N+f2,:)); 
    %effective occupation of each site in the double exciton
    end
    
    for k = 1:N  
    mu_ef(f, k,:) = cq_nm(f,k)*( M_tot(k+1,2:N+1)*mu );
    mag_ef(f,k,:) = cq_nm(f,k)*( M_tot(k+1,2:N+1)*...
                                  ( mdip+ mdip2*(E_2(f)-E_1(k))) );
    end
end  
mu_ex = {mu_ge,mu_ef};
mag_ex = {mag_ge,mag_ef};
else
    mu_ex = mu_ge; mag_ex  = mag_ge;
end
else
   mu_ex = []; mag_ex = []; 
end
%% Next generate vibrational hamiltonian and displacement op term

if all(cellfun(@isempty,om_vib)) %no vibrations included
    H_vib = 0; H_coup = [];
else
    
    tot_vib_each_site = cellfun(@prod,numvib);
H_vib = zeros(prod(tot_vib_each_site));
H_coup = zeros(length(H_vib),length(H_vib),dip_size-1  );

for j = 1:N
    N2 = length(om_vib{j}); nvib = numvib{j};
    H_vib_sec{j} = zeros(prod(nvib)); %vibrations on site j
    disp = displ{j}; %dislacement of modes on site j in excited state
    for j2 = 1:N2 %loop over number of vibrations on site N
 
        %pad each section with identity matricies featuring the other
        %vibrations
    tmp = kron(diag(1/2:(nvib(j2)-1/2))*om_vib{j}(j2),eye(prod(nvib(1:j2-1))));
    tmp = kron(eye(prod(nvib(j2+1:N2))),tmp);    
    %collect together
    H_vib_sec{j} = H_vib_sec{j} + tmp;
    
    for ell = 1:dip_size-1 %loop over all included states except GS
        % om(j)* (Delta(j,ell1,ell2)^2/2 - Delta(j,ell1,ell2)(b_j+b_j^dagger)/sqrt(2)) 
        %=(om(j)*Delta(j,ell1,ell2)/sqrt(2))*( (Delta(j,ell1,ell2)/sqrt(2)) 
        %       -(b_j+b_j^dagger) )
        if fock_space_rep(ell+1,j)~=0
        bj = diag(sqrt(1:nvib(j2)-1),1);
        tmp = (disp(j2)*eye(nvib(j2))/2 -(bj+bj')/sqrt(2)).*om_vib{j}(j2)*disp(j2);
        %pad out to the full Hilbert space
        
        tmp = kron(tmp,eye(prod(tot_vib_each_site(1:j-1))*prod(nvib(1:j2-1))));
        tmp = kron(eye(prod(tot_vib_each_site(j+1:end))*prod(nvib(j2+1:end))),tmp );       
        H_coup(:,:,ell) = H_coup(:,:,ell) + tmp;    
        end
    end
    end
    %pad this section out to the full size oh H_vib
       tmp = kron(H_vib_sec{j},eye(prod(tot_vib_each_site(1:j-1))));
     H_vib = H_vib +  kron(eye(prod(tot_vib_each_site(j+1:end))),tmp);    
end
%combine together

%H_ex_vib = kron(H_dip,eye(size(H_vib))) +  kron(eye(size(H_dip)),H_vib);
H_e = kron(H_e,eye(size(H_vib))) +  kron(eye(size(H_e)),H_vib);
H_f = kron(H_f,eye(size(H_vib))) +  kron(eye(size(H_f)),H_vib);

   for ell = 1:N
        tmp = zeros(N); tmp(ell,ell) = 1;
       H_e = H_e + kron(tmp,H_coup(:,:,ell));
       %if values of H_dip are those from spectroscopy then don't take the
       %diagonal part of this matrix
   end
   for ell = N+1:dip_size-1
        tmp = zeros(N*(N-1)/2); tmp(ell-N,ell-N) = 1;
       H_f = H_f + kron(tmp,H_coup(:,:,ell));
   end   
    
end
   [M_1,e_1]=eig(H_e); e_1 = diag(e_1); [M_2,e_2]=eig(H_f); e_2=diag(e_2);
   M_pfull  = {M_1,M_2}; 

