function [R_red,R_red_op,Ctilde] = redfield_calc2(H_site,beta,gam_dru,...
                    lam_dru,gamma,lambda,om_0,use_markov)
%expand to include double excited manifold and vibrational states via just
%padding out the superoperator
                %Calculate redfield tensor
%R_{abcd} = delta_{ab} Gam_{a,e,e,c}(om_ce) + delta_{ac} Gam^*_{b,e,e,d}
%           - Gam_{d,b,a,c}(om_ca) - Gam^*_{c,a,b,d}(om_db)
%sum over e inplicit
%diagonalise the first excited state manifold into exciton states
% assumes no significant chance of transition gs
if nargin == 7
    use_markov = true;
end
ex_st = zeros(1+length(H_site));
[ex_st(2:end,2:end),Hex]=eig(H_site);

Hex = diag(Hex); % this assumes that the GS energy is set to zero
Hexfull = [0;Hex];
% assume on site baths so
% Gam_{abcd}(om) = sum_j Q_{j,ab} Q_{j,cd} * int_0^t dt' e^{i om t') C_j(t)
% no C_{ij}(t) cross terms etc, also take Markovian t->inf
% om is always om_dc = om_d - om_c
%Calculate Fourier Laplace transform of the site bath time correlations
if use_markov
Ctilde = zeros(length(Hex),length(Hexfull),length(Hexfull));
else
Ctilde = sym(zeros(length(Hex),length(Hexfull),length(Hexfull)));     
syms t real
end                    
%markovian calculation so int_0^inf C_j(t) exp(iwt) taken

%calculate C_{jj}(t) = sum_{k} b_k exp(-v_k t) (it will have this form with
%possibly complex v_k for the damped oscillator forms.
vn = 2*pi*(1:20)/beta;  %take 20 matsubara, probably more than enough
pm = [+1,-1];

   
    CC = zeros(length(Hex),max(cellfun(@length,gam_dru)) +...
                        2*max(cellfun(@length,gamma)));
    Cmatsu = zeros(length(Hex),length(vn));
    vv = CC;

for lp = 1:length(Hex)
    cnt =0;
    gam = gam_dru{lp}; lam = lam_dru{lp}; 
    for k =1:length(gam)  %drude mode
        cnt =cnt+1;
        CC(lp,cnt) =  lam(k)*gam(k)*(cot(beta*gam(k)/2)-1i);
        Cmatsu(lp,:) = Cmatsu(lp,:) + lam(k)*gam(k)*4*vn./(beta.*(vn.^2-gam(k)^2));
        vv(lp,cnt) = gam(k);
    end
    gam = gamma{lp}; lam = lambda{lp}; om0 = om_0{lp}; 
    for k =1:length(gam) %underdamped
        cnt =cnt+1;
            eta = sqrt(om0(k)^2-gam(k)^2/4); 
            phi = gam(k)/2 + pm*1i*eta;
            Phi = pm*(1i/2/eta) .* (cot(phi.*beta/2)-1i);
        CC(lp,cnt:cnt+1) = lam(k)*om0(k)^2*Phi;
        Cmatsu(lp,:) = Cmatsu(lp,:) -  4*lam(k)*om0(k)^2*gam(k)...
            .*vn./(beta.*((om0(k)^2+vn.^2).^2-gam(k)^2*vn.^2));
        vv(lp,cnt:cnt+1) = phi;
        cnt =cnt+1;
    end
end

for a=1:length(Hexfull)
    for b = 1:length(Hexfull)

            deltaom = Hexfull(a)-Hexfull(b);

            for lp = 1:length(Hex)
             if use_markov
             Ctilde(lp,a,b) = sum(-CC(lp,:)./(1i*deltaom-vv(lp,:)));
             Ctilde(lp,a,b) = Ctilde(lp,a,b)+ sum(-Cmatsu(lp,:)./(1i*deltaom-vn));
             else
             Ctilde(lp,a,b) = sum(CC(lp,:).*(exp((1i*deltaom-vv(lp,:))*t)-1)./(1i*deltaom-vv(lp,:)));
             Ctilde(lp,a,b) = Ctilde(lp,a,b)+ sum(Cmatsu(lp,:)...
                              .*(exp((1i*deltaom-vn)*t)-1)./(1i*deltaom-vn));                 
             end
            end
    end
end


 %Gam_RF is tensor which makes up the components of the redfield tensor
 % The symbolic component from int_0^inf C_j(t) exp(iwt) is held seperately
% tic

 Gam_RF = zeros(length(Hex)+1,length(Hex)+1,length(Hex)+1,length(Hex)+1,length(Hex));
for a = 1:length(Hex)+1
    for b = 1:length(Hex)+1
        for c = 1:length(Hex)+1
            for d = 1:length(Hex)+1
               
                if ~any([a,b,c,d]==1)
                for lp = 1:length(Hex) %site loop, each site can have diff bath
                    Gam_RF(a,b,c,d,lp)= Ctilde(lp,d,c).*...
                    ex_st(a,lp+1)*ex_st(b,lp+1)*ex_st(c,lp+1)*ex_st(d,lp+1);
                end
                end
            end
        end
    end
end

%Gam_RF =real(Gam_RF);
%construct redfield tensor, won't be symbolic assuming all transition
%frequencies are given
%R_{abcd} = delta_{db} Gam_{a,e,e,c}(om_ce) + delta_{ac} Gam^*_{b,e,e,d}(om_de)
%           - Gam_{d,b,a,c}(om_ca) - Gam^*_{c,a,b,d}(om_db)

% Note ma

R_red = zeros(length(Hex)+1,length(Hex)+1,length(Hex)+1,length(Hex)+1);

for a = 1:length(Hex)+1
    for b = 1:length(Hex)+1
        for c = 1:length(Hex)+1
            for d = 1:length(Hex)+1
 
                for lp = 1:length(Hex)
R_red(a,b,c,d)  = R_red(a,b,c,d) -(Gam_RF(d,b,a,c,lp) ...
                            +conj(Gam_RF(c,a,b,d,lp)));               
                if b==d
                for e5 = 2:length(Hex)+1 %site loop, each site can have diff bath
R_red(a,b,c,d)  = R_red(a,b,c,d)  + Gam_RF(a,e5,e5,c,lp);
                end
                end
                if a==c
                for e5 = 2:length(Hex)+1 %site loop, each site can have diff bath
                    %first one should be zero anyway
R_red(a,b,c,d)  = R_red(a,b,c,d)  + conj(Gam_RF(b,e5,e5,d,lp));
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
% -sum ( R_red(a,b,k1,k2) rho_{k1,k2}), in Liouville space, 
% rho(a,b) -> rho( a+ N*(b-1) ) so
% d rho( a+ N*(b-1) )/dt ~= -sum ( R_red(a,b,k1,k2) rho(k1+N*(k2-1)) )
% so the   a+ N*(b-1) th row of the Liouville space matrix is the reshaped
% R_red(a,b,:,:) 
% if N2 vibrational states are also present 
% d rho( N2*(a-1)+A+ N*N2*(N2*(b-1)+B) )/dt ~= -sum 
%( R_red(a,b,k1,k2) rho(N2*(k1-1)+A+N*N2*(N2*(k2-1)+B) )  )
%Note reshape goes down columns then rows
if nargout >= 2
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
            temp(a,b) = R_red(cnt1,cnt2,a,b);
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
% rhotest = rand(length(Hex)+1); rhotest  = rhotest + rhotest';
% rhodot1 = rhotest*0; 
%    for e1 = 1:length(Hex)+1
%     for e2 = 1:length(Hex)+1
%         for e3 =1:length(Hex)+1
%            for e4 =1:length(Hex)+1 
%         rhodot1(e1,e2) = rhodot1(e1,e2) + R_red(e1,e2,e3,e4)*rhotest(e3,e4);
%            end
%         end
%     end
%    end
%         rhodot2 = R_red_op*reshape(rhotest,numel(rhotest),1);
% rhodot1-reshape(rhodot2.',length(Hex)+1,length(Hex)+1)
% rhodot1-reshape(rhodot2,length(Hex)+1,length(Hex)+1)
