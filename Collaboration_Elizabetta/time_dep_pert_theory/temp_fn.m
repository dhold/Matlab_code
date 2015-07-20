
om_dip = H_site;
om_vib = [om_0{:}];
numvib = [4,4];
mu = ones(1,2);

    %only 1 level exciton Ham
    N = length(om_dip);
    dip_size = 1+N;
    flag = false;  pdm_n_pos =[];

% Firstgenerate two level excition hamiltonian

H_dip = zeros(dip_size); cnt = 0; mu_full = zeros(dip_size);
fock_space_rep = false(dip_size,N); 
fock_space_rep(2:N+1,:) = logical(diag(ones(N,1)));

mu_full(1,2:N+1) = mu;

if length(om_dip) == length(H_dip)-1;
 H_dip(2:end,2:end) = om_dip;   %can pass om_dip as whole site basis ham
 for j = 1:N
     rng = j+1:N;
     for k = 1:length(rng)
         cnt = cnt+1; kk = rng(k);
         mu_full(j,N+1+cnt) = mu(kk);
         fock_space_rep(N+1+cnt,j) = 1; fock_space_rep(N+1+cnt,kk) = 1;
     end
 end
    
H_dip(2:N+1,2:N+1) = om_dip;

% Next generate vibrational hamiltonian and displacement op term

N2 = length(om_vib);

H_vib = zeros(prod(numvib)); %do not set any numvib to zero

H_coup = zeros(prod(numvib),prod(numvib),dip_size  );

for j = 1:N2
    
    tmp = kron(eye(prod(numvib(1:j-1))),...
        diag(1/2:1:numvib(j)-1/2)*om_vib(j));
    tmp =  kron(tmp , eye(prod(numvib(j+1:N2))));

    H_vib = H_vib + tmp;
    
    for ell = 1:dip_size 
        %(om(j)*Delta(j,ell1,ell2)/sqrt(2))*( (Delta(j,ell1,ell2)/sqrt(2)) 
        %-(b_j+b_j^dagger) )
    tmp = kron(eye(prod(numvib(1:j-1))),(displ(j,ell)/sqrt(2)*eye(numvib(j))...
        - diag(sqrt(1:numvib(j)-1),1) - diag(sqrt(1:numvib(j)-1),-1))...
        *om_vib(j)*displ(j,ell)/sqrt(2));
    tmp =  kron(tmp , eye(prod(numvib(j+1:N2))));        
    H_coup(:,:,ell) = H_coup(:,:,ell) + tmp;    
    end
    
    
end
%combine together

H_ex_vib = kron(H_dip,eye(size(H_vib))) +  kron(H_vib,eye(size(H_dip)));

   for ell = 1:dip_size 
        tmp = zeros(size(H_dip)); tmp(ell,ell) = 1;
       H_ex_vib = H_ex_vib + kron(H_coup(:,:,ell),tmp);
       %if values of H_dip are those from spectroscopy then don't take the
       %diagonal part of this matrix
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
   indiv_op{5} = proj;
end
end