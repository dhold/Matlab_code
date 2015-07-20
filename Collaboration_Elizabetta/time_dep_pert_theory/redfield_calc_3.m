function [R_red,R_red_op,Ctilde, Gam_RF,Q_j] = redfield_calc_3(H_single,H_double...
            ,fock_space_rep,beta,gam_dru,lam_dru,gamma,lambda,om_0)
%Calculate redfield tensor
%R_{abcd} = delta_{ab} Gam_{a,e,e,c}(om_ce) + delta_{ac} Gam^*_{b,e,e,d}
%           - Gam_{d,b,a,c}(om_ca) - Gam^*_{c,a,b,d}(om_db)
%sum over e inplicit
%diagonalise the first excited state manifold into exciton states
% assumes no significant chance of transition gs
%H_single is the first excited state basis Hamiltonian
%H_double should also include the ground state, which doublely excited
%states mix into via perm dipole moments...

N = length(H_single);
[ex_st1,Hex] = eig(H_single);  %Hamiltonian is block diagonal in this way
[ex_st2,Hex_2]=eig(H_double); %ground state can mix into double

[ex_st,Hexfull] = eig(blkdiag(H_double(1,1),H_single,H_double(2:end,2:end)));
Hexfull = diag(Hexfull); LL = length(Hexfull);

Ctilde_re = zeros(N,LL,LL);

for a=1:LL
    for b = 1:LL

            deltaom = Hexfull(a)-Hexfull(b);
%note this assumes non correlated baths, I could include correlated baths
%with Ctilde(lp1,lp2,a,b) but I don't need to yet
n_fc = 1/(1-exp(-beta*deltaom));
            for lp = 1:N
                    
gam = gam_dru{lp}; lam = lam_dru{lp}; 
Ctilde_re(lp,a,b) = sum(2*gam.*lam.*deltaom./(deltaom^2+gam.^2));
 gam = gamma{lp}; lam = lambda{lp}; om0 = om_0{lp}; 
Ctilde_re(lp,a,b) = Ctilde_re(lp,a,b) + sum(2*gam.*om0.^2.*lam.*ww./((om0.^2-ww^2).^2+ww^2*gam.^2));                              
Ctilde_re(lp,a,b) = Ctilde_re(lp,a,b)*n_fc;
            end
    end
end


 %Gam_RF is tensor which makes up the components of the redfield tensor
 % The symbolic component from int_0^inf C_j(t) exp(iwt) is held seperately


 Gam_RF = zeros(LL,LL,LL,LL,length(Hex));
        
 %find mixings <a | Q_j |b>   
 fock_space_rep  = double(fock_space_rep); %convert from logical
 manifold_chk = round(sum(fock_space_rep,2));
 Q_j = zeros(LL,LL,N); %mixes only within
 %the same excitation manifold
 %tic
 for lp = 1:N
     for lp2 = 2:LL
         occ_L = ex_st(:,lp2).*fock_space_rep(:,lp); %occupation of this state on bra
        for lp3 = 2:LL
            occ_R = ex_st(:,lp3).*fock_space_rep(:,lp); %occupation of this state on ket
            if manifold_chk(lp2) == manifold_chk(lp3) %else no mixing
                    Q_j(lp2,lp3,lp) =  occ_L'*occ_R;
            end
        end
     end
 end
     
Ctilde = permute(Ctilde,[3,2,1]);
tmp = permute(Ctilde,[2,1,3]).*Q_j;
%tic
for a = 1:LL
    for b = 1:LL
         tmp2 = tmp.*repmat(Q_j(a,b,:),LL,LL,1);
         Gam_RF(a,b,:,:,:)=tmp2;
    end
end
Gam_RF = sum(Gam_RF,5); %take the sum over sites


R_red = -(permute(Gam_RF,[3,2,4,1])+conj(permute(Gam_RF,[2,3,1,4])));  
 
tmp = diagsum(Gam_RF,2,3); 
 II = eye(LL); %rep of delta functions

R_red_tmp1 = zeros(size(R_red));R_red_tmp2 = zeros(size(R_red));

for a=1:LL
    %b=a; 
    for c = 1:LL
       % d=c;
         R_red_tmp1(a,:,c,:) = tmp(a,c).*II; 
         R_red_tmp2(:,a,:,c) = conj(tmp(a,c)).*II;          
    end
end

R_red = R_red+R_red_tmp1+R_red_tmp2;
end