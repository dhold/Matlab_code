function [R_red,R_red_op] = redfield_calc_backup(H_site,beta,gam_dru,lam_dru,gamma,lambda,om_0,use_markov)
%Calculate redfield tensor
%R_{abcd} = delta_{ab} Gam_{a,e,e,c}(om_ce) + delta_{ac} Gam^*_{b,e,e,d}
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
Ctilde = zeros(length(Hex),length(Hex)+1,length(Hex)+1);
%markovian calculation so int_0^inf C_j(t) exp(iwt) taken
syms deltaom
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
        fnname{lp} = strcat('Ctild_tmp',num2str(10000000*rand(1),7),'_fn');
        
        for pntlp = 1:10 %small finite chance to name same fn twice
            brkcnd = true;
            for lp2 = 1:lp-1
           if strcmp(fnname{lp},fnname{lp2})
              fnname{lp} = strcat('Ctild_tmp',num2str(10000000*rand(1),7),'_fn'); 
              brkcnd = false;
           end
            end
            if brkcnd
                break
            end
        end
sym_to_fn(strcat(fnname{lp},'.m'),vpa(Ctild(lp)),deltaom);

end

%preallocate with all values that it will take
for lp = 1:length(Hex)
    for a=1:length(Hex)
    Ctilde(lp,a+1,1) = eval(strcat(fnname{lp},'(Hex(a))'));
    for b = 1:length(Hex)
        Ctilde(lp,a+1,b+1) = eval(strcat(fnname{lp},'(Hex(a)-Hex(b))'));
        % array of all possible values of Ctilde
        if a==1
        Ctilde(lp,a+1,1) = eval(strcat(fnname{lp},'(-Hex(b))')) ;
        end
    end
    end
delete(fnname{lp}); %get rid of tmp func
end

 %Gam_RF is tensor which makes up the components of the redfield tensor
 % The symbolic component from int_0^inf C_j(t) exp(iwt) is held seperately
% tic
 clear R_red Gam_RF

 Gam_RF = zeros(length(Hex)+1,length(Hex)+1,length(Hex)+1,length(Hex)+1,length(Hex));
for a = 1:length(Hex)+1
    for b = 1:length(Hex)+1
        for c = 1:length(Hex)+1
            for d = 1:length(Hex)+1
               
                if a==1 || b==1 || c==1 || d==1
                    %Gam_RF(a,b,c,d)= zeros(length(Hex),1)  %preallocated anyway 
                else
                for lp = 1:length(Hex) %site loop, each site can have diff bath
                    Gam_RF(a,b,c,d,lp)...
                    =ex_st(a-1,lp)*ex_st(b-1,lp)*ex_st(a-1,lp)*ex_st(b-1,lp);
                end
                end
            end
        end
    end
end

%construct redfield tensor, won't be symbolic assuming all transition
%frequencies are given
%R_{abcd} = delta_{ab} Gam_{a,e,e,c}(om_ce) + delta_{ac} Gam^*_{b,e,e,d}(om_de)
%           - Gam_{d,b,a,c}(om_ca) - Gam^*_{c,a,b,d}(om_db)

% Note ma

R_red = zeros(length(Hex)+1,length(Hex)+1,length(Hex)+1,length(Hex)+1);

for a = 1:length(Hex)+1
    for b = 1:length(Hex)+1
        for c = 1:length(Hex)+1
            for d = 1:length(Hex)+1
 
                for lp = 1:length(Hex)+1
R_red(a,b,c,d)  = R_red(a,b,c,d) -(Gam_RF(d,b,a,c,lp)*Ctilde(lp,c,a) ...
                            +Gam_RF(c,a,b,d,lp)*Ctilde(lp,d,b));
               
                if a==b
                for e5 = 1:length(Hex)+1 %site loop, each site can have diff bath
R_red(a,b,c,d)  = R_red(a,b,c,d)  + Gam_RF(a,e5,e5,c,lp)*Ctilde(lp,c,e5);
                end
                elseif a==c
                for e5 = 1:length(Hex)+1 %site loop, each site can have diff bath
R_red(a,b,c,d)  = R_red(a,b,c,d)  + conj(Gam_RF(b,e5,e5,d,lp)*Ctilde(lp,d,e5));
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
R_red_op = zeros((length(Hex)+1)^2); 
temp = zeros(length(Hex)+1);
cnt1 = 0; cnt2 = 1;
for lp = 1:(length(Hex)+1)^2
    
    cnt1 = cnt1+1;
    if cnt1 > length(Hex)+1
        cnt1=1; cnt2 = cnt2+1;
    end
    temp = temp*0;
   for a = 1:length(Hex)+1
        for b = 1:length(Hex)+1
            temp(a,b) = R_red{cnt1,cnt2,a,b};
        end
   end
    R_red_op(lp,:) = reshape(temp,1,(length(Hex)+1)^2);
    
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
