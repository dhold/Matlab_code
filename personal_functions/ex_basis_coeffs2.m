function [mu_ex,m_ex,V_ge,V_ef]=ex_basis_coeffs(M_g,M_e,M_f,fock_rep1,fock_rep2,mu,m_intr,R,vib_pad)

% M_{g,e,f} are projectors from site basis to the exciton basis
% mu and m_intr are the chromophore electric and magnetic transition dipole
% moments and R are the positions of the chromophores
% if vib_pad is given will pad with an appropriate number of viblvls.

Ng = length(M_g); N = length(M_e); N2 = length(M_f);

for j=1:N
    
    
       
    for g = 1:Ng
        
        V_ge(:,:,j,g) = kron(fock_rep1(j,:),II)*M_e;  
        
    end
    
    coup_coeff = fock_rep2(:,logical(fock_rep1(j,:)));
    
    for f = 1:N*(N-1)/2
         V_ef(:,:,j,f) = kron(fock_rep2(f,:),M_e)*M_f; 
         V01_tmp(rng1+j*length(M_g),rng1) = eye(length(M_g));
    end
    
end


II = eye(vib_pad); %identity operator, equal to number of viblvls
M_g2 = kron(M_g,II); M_e2 = kron(M_e,II);  M_f2 = kron(M_f,II);

V_ge = zeros(length(M_g2),length(M_e2),N,Ng); rng1 = 1:N; rng2 = 1:N2;
V_ef = zeros(length(M_e2),length(M_f2),N,N2);

alpha_n = mtimesx(MM,'C',mtimesx(V_n,MM));
Cg  = squeeze(alpha_n(1,2:N+1,:)) ;  Cf = alpha_n(2:N+1,N+2:end,:);

for j=1:N

    coupto = fock_rep2(:,logical(fock_rep1(j,:))) == 1;
    for g = 1:Ng
        
        V_ge(:,:,j,g) = kron(fock_rep1(j,:),M_g2)*M_e2;  
        
    end
    
    for f = 1:N*(N-1)/2
         V_ef(:,:,j,f) = kron(fock_rep2(f,:),M_e2)*M_f2; 
    end
    
end