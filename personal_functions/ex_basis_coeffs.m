function [mu_ex1,mu_ex2,m_ex,V_ge,V_ef]=ex_basis_coeffs(M_e,M_f,fock_rep1,...
                    fock_rep2,mu,m_intr,R,vib_pad)

% M_{e,f} are projectors from site basis to the exciton basis
% mu and m_intr are the chromophore electric and magnetic transition dipole
% moments and R are the positions of the chromophores
% fock_rep1 and fock_rep2 are the excitons occupied in the 
% if vib_pad is given will pad with an appropriate number of viblvls.

N = length(M_e); N2 = length(M_f);

    C_g = M_e'*fock_rep1;
    C_f = zeros(N2,N,N);
    for j=1:N
        tmp = fock_rep2; tmp(:,[1:j-1,j+1:N]) = 0;
    C_f(:,:,j) = M_f'*tmp*M_e;
    end
mu_ex1 = C_g*mu; mu_ex2 = C_f*mu;

II = eye(vib_pad); %identity operator, equal to number of viblvls
M_g2 = kron(M_g,II); M_e2 = kron(M_e,II);  M_f2 = kron(M_f,II);

V_ge = zeros(length(M_g2),length(M_e2),N,Ng); 
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