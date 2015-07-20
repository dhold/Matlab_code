function [R_red,R_red_op] = redfield_calc_old(H_site,beta,gam_dru,lam_dru,gamma,lambda,om_0,use_markov)
%Calculate redfield tensor
% This version does not include the ground state coherences
%R_{abcd} = delta_{ab} Gam_{a,e,e,c}(om_ce) + delta_{ad} Gam^*_{b,e,e,d}
%           - Gam_{d,b,a,c}(om_ca) - Gam^*_{c,a,b,d}(om_db)
%sum over e inplicit
%diagonalise the first excited state manifold into exciton states
% assumes no significant chance of transition gs
if nargin == 7
    use_markov = true;
end
[ex_st,Hex]=eig(H_site  );
Hex = diag(Hex);
% assume on site baths so
% Gam_{abcd}(om) = sum_j Q_{j,ab} Q_{j,cd} * int_0^t dt' e^{i om t') C_j(t)
% no C_{ij}(t) cross terms etc, also take Markovian t->inf
%
%Calculate Fourier Laplace transform of the site bath time correlations
Ctild = sym(zeros(length(Hex),1));

%markovian calculation so int_0^inf C_j(t) exp(iwt) taken
syms deltaom
%nom = 1/(exp(beta*deltaom)-1);
%tic
for lp = 1:length(Hex)
    
    tmp = gam_dru{lp+1}; tmp2 = lam_dru{lp+1}; 
    for k =1:length(tmp)
              if tmp(k) ~= 0
        tmp3 = sym(0);
        if use_markov
        tmp3  = tmp3 + (cot(beta*tmp(k)/2)-1i)/(1i*deltaom-tmp(k));
        vn = 2*pi*(1:20)/beta; %take 20 matsubara, probably more than enough
         %Ctild(lp) = Ctild(lp) + (nom+1)*4*tmp(k)*tmp2(k)*om/(om^2+tmp2(k)^2);
         %that's the real part
        tmp3  = tmp3 + 4*vn./beta./(vn.^2-tmp(k)^2)./(1i*deltaom-vn);
        Ctild(lp) =Ctild(lp) + tmp3 *(tmp(k)*tmp2(k));
        end
              else
                  error('you havent written non markov yet')
              end
    end
    
    tmp = gamma{lp+1}; tmp2 = lambda{lp+1};     tmp3 = om_0{lp+1}; 
        for k =1:length(tmp)
            if tmp(k) ~= 0
                if use_markov
           tmp4 =sym(0);
            eta = sqrt(tmp3(k)^2-tmp(k)^2/4); pm = [+1,-1];
            phi = tmp(k)/2 + pm*1i*eta;
            Phi = pm*(1i/2/eta) .* (cot(phi.*beta/2)-1i);
       tmp4 = tmp4 + Phi(1)/(1i*deltaom-phi(1)) + Phi(2)/(1i*deltaom-phi(2));
       vn = 2*pi*(1:20)/beta; %take 20 matsubara, probably more than enough
       tmp4 = tmp4 + 4*tmp(k)/beta*...
           sum(vn./((tmp3(k)^2+vn.^2).^2-tmp(k)^2*vn.^2)./(1i*deltaom-vn));    
       Ctild(lp) = Ctild(lp) + tmp4*tmp2(k)*tmp3(k)^2; %all have a lambda*om_0^2 prefac
                else
                end
            end
        end
    
end
% tic
% Ctild = simplify(Ctild);
% toc
%toc
 %Gam_RF is tensor which makes up the components of the redfield tensor
 % The symbolic component from int_0^inf C_j(t) exp(iwt) is held seperately
% tic
 clear R_red Gam_RF
for e1 = 1:length(Hex)
    for e2 = 1:length(Hex)
        for e3 = 1:length(Hex)
            for e4 = 1:length(Hex)
                R_red{e1,e2,e3,e4} = 0; %preallocate while doing this loop     
                for lp = 1:length(Hex) %site loop, each site can have diff bath
                    Gam_RF{e1,e2,e3,e4}(lp)...
                    =ex_st(e1,lp)*ex_st(e2,lp)*ex_st(e1,lp)*ex_st(e2,lp);
               
                end
            end
        end
    end
end
%toc
%construct redfield tensor, won't be symbolic assuming all transition
%frequencies are given
%R_{abcd} = delta_{ab} Gam_{a,e,e,c}(om_ce) + delta_{ad} Gam^*_{b,e,e,d}(om_de)
%           - Gam_{d,b,a,c}(om_ca) - Gam^*_{c,a,b,d}(om_db)
%tic %this part is slow as subs is wank but tolerably so for small system
% For a large system I will have to sort this out
for e1 = 1:length(Hex)
    for e2 = 1:length(Hex)
        for e3 = 1:length(Hex)
            for e4 = 1:length(Hex)
 
                for lp = 1:length(Hex)
R_red{e1,e2,e3,e4}  = R_red{e1,e2,e3,e4} +double(-Gam_RF{e4,e2,e1,e3}(lp)*...
                        subs(Ctild(lp),deltaom,Hex(e3)-Hex(e1))...
                        -conj(Gam_RF{e3,e1,e2,e4}(lp)*...
                        subs(Ctild(lp),deltaom,Hex(e4)-Hex(e2))));   
               
                if e1==e2
                for e5 = 1:length(Hex) %site loop, each site can have diff bath
R_red{e1,e2,e3,e4}  = R_red{e1,e2,e3,e4}  + Gam_RF{e1,e5,e5,e3}(lp)*...
                        double(subs(Ctild(lp),deltaom,Hex(e3)-Hex(e5)));
                end
                elseif e1==e4
                for e5 = 1:length(Hex) %site loop, each site can have diff bath
R_red{e1,e2,e3,e4}  = R_red{e1,e2,e3,e4}  + conj(Gam_RF{e2,e5,e5,e4}(lp)*...
                         double(subs(Ctild(lp),deltaom,Hex(e4)-Hex(e5))));
                end
                end
                end
            end
        end
    end
end
%toc

%% express in Liouville space as super operator

% each density matrix element rho_{ab} decays with a term prop to
% -sum ( R_red(a,b,k1,k2) rho_{k1,k2})
% R_red(a,b,0,j) = R_red(a,b,k,0) = 0, i.e. terms which relate to ground
% state coherences etc
%
%Note reshape goes down columns then rows
if nargout == 2
R_red_op = zeros(length(Hex)^2); 
temp = zeros(length(Hex));
cnt1 = 0; cnt2 = 1;
for lp = 1:length(Hex)^2
    
    cnt1 = cnt1+1;
    if cnt1 > length(Hex)
        cnt1=1; cnt2 = cnt2+1;
    end
    temp = temp*0;
   for e1 = 1:length(Hex)
        for e2 = 1:length(Hex)
            temp(e1,e2) = R_red{cnt1,cnt2,e1,e2};
        end
   end
    R_red_op(lp,:) = reshape(temp,1,(length(Hex))^2);
    
end
end

%this can only be applied to the excited state part of the density matrix
%i.e. pckr = false(length(Hex)+1); pickr(2:end,2:end) = true;
% drho(pickr) /dt = R_red_op*rho(pickr)

%also include ground state

%test
% 
% rhotest = rand(length(Hex));
% rhodot1 = rhotest*0; rhodot2 = rhotest*0;
%    for e1 = 1:length(Hex)
%     for e2 = 1:length(Hex)
%         for e3 =1:length(Hex)
%            for e4 =1:length(Hex) 
%         rhodot1(e1,e2) = rhodot1(e1,e2) + R_red{e1,e2,e3,e4}*rhotest(e3,e4);
%            end
%         end
%     end
%    end
%         rhodot2 = R_red_op*reshape(rhotest,numel(rhotest),1);
% rhodot1-reshape(rhodot2.',length(Hex),length(Hex))
% rhodot1-reshape(rhodot2,length(Hex),length(Hex))
