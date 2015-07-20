function [R_red,R_red_op,Ctilde, Gam_RF,Q_j] = redfield_calc_2(H_single,H_double...
            ,fock_space_rep,beta,gam_dru,lam_dru,gamma,lambda,om_0,use_markov)
%Calculate redfield tensor
%R_{abcd} = delta_{ab} Gam_{a,e,e,c}(om_ce) + delta_{ac} Gam^*_{b,e,e,d}
%           - Gam_{d,b,a,c}(om_ca) - Gam^*_{c,a,b,d}(om_db)
%sum over e inplicit
%diagonalise the first excited state manifold into exciton states
% assumes no significant chance of transition gs
%H_single is the first excited state basis Hamiltonian
%H_double should also include the ground state, which doublely excited
%states mix into via perm dipole moments...
if nargin == 9
    use_markov = true;
end
if ~isempty(gam_dru) && isempty(lam_dru) %alternative format passing
    %single N length cell array with n by 2 matrix in each bit
    tmp = gam_dru; N = length(tmp); clear gam_dru lam_dru
    for j = 1:N
    lam_dru{j} = tmp{j}(:,1);
    gam_dru{j} = tmp{j}(:,2);
    end
end
if ~isempty(gamma) && isempty(lambda) %alternative format passing
    %single N length cell array with n by 3 matrix in each bit
    tmp = gamma; N = length(tmp); clear gamma lambda om_0
    for j = 1:N
        if size(tmp{j},2)==4
        lg = tmp{j}(:,4) == 0; %only these modes
        else
         lg =  true(size(tmp{j},1),1);
        end
    lambda{j} = tmp{j}(lg,1);
    gamma{j} = tmp{j}(lg,2);
    om_0{j} = tmp{j}(lg,3);
    end
end


N = length(H_single);
[ex_st1,Hex] = eig(H_single);  %Hamiltonian is block diagonal in this way
if ~isempty(H_double)
[ex_st2,Hex_2]=eig(H_double); %ground state can mix into double

[ex_st,Hexfull] = eig(blkdiag(H_double(1,1),H_single,H_double(2:end,2:end)));
Hexfull = diag(Hexfull); 
else
 Hexfull = Hex;   ex_st = ex_st1;
end
LL = length(Hexfull);
%Hexfull = [tmp(1);diag(Hex);tmp(2:end)]; 
% assume on site baths so
% Gam_{abcd}(om) = sum_j Q_{j,ab} Q_{j,cd} * int_0^t dt' e^{i om t') C_j(t)
% no C_{ij}(t) cross terms etc, also take Markovian t->inf
% om is always om_dc = om_d - om_c
%Calculate Fourier Laplace transform of the site bath time correlations
if use_markov
Ctilde = zeros(length(Hex),LL,LL);
else
Ctilde = sym(zeros(length(Hex),LL,LL));     
syms t real
end                    
%markovian calculation so int_0^inf C_j(t) exp(iwt) taken

%calculate C_{jj}(t) = sum_{k} b_k exp(-v_k t) (it will have this form with
%possibly complex v_k for the damped oscillator forms.
%vn = 2*pi*(1:20)/beta;  %take 20 matsubara, probably more than enough
pm = [+1,-1];
vn = 2*pi*(1:5000)/beta; %this actually doesn't converge that fast, take more
    if ~isempty(gamma)
    CC = zeros(length(Hex),max(cellfun(@length,gam_dru)) +...
                        2*max(cellfun(@length,gamma)));
    else
       CC = zeros(length(Hex),max(cellfun(@length,gam_dru))); 
    end
    Cmatsu = zeros(length(Hex),length(vn));
    vv = CC;
    
  %  ww=sym('w','real'); ww2 = sym('W','real');
   % J_w =sym(zeros(length(H_single),1));  %total spec density
   % n_w =1/(exp(beta*ww)-1); %BEdist
    
for lp = 1:length(H_single)
    cnt =0;
    gam = gam_dru{lp}; lam = lam_dru{lp}; 
    for k =1:length(gam)  %drude mode
        cnt =cnt+1;
        CC(lp,cnt) =  lam(k)*gam(k)*(cot(beta*gam(k)/2)-1i);
        Cmatsu(lp,:) = Cmatsu(lp,:) + 4*lam(k)*gam(k)*vn./(beta.*(vn.^2-gam(k)^2));
        vv(lp,cnt) = gam(k);        
      %  J_w(lp) = J_w(lp) + 2*gam(k)*lam(k)*ww/(ww^2+gam(k)^2);
    end
    %pass empty values to not use underdamped modes in this
    if ~isempty(lambda) %don't include underdamped at all
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
       % J_w(lp) = J_w(lp) + 2*gam(k)*om0(k)^2*lam(k)*ww/((om0(k)^2-ww^2)^2+ww^2*gam(k)^2);
    end
    end
end

%real_part_fn = matlabFunction((n_w+1)*(J_w-subs(J_w,ww,-ww)));
%imag_fn = int(J_w/(ww2 - ww), ww, -inf, inf, 'PrincipalValue', true);
%CPV_fn = matlabFunction(imag_fn);

for a=1:LL
    for b = 1:LL

            deltaom = Hexfull(b)-Hexfull(a);
%note this assumes non correlated baths, I could include correlated baths
%with Ctilde(lp1,lp2,a,b) but I don't need to yet
            for lp = 1:N
             if use_markov %take time int to infinity
                 %first include contribution from poles in J(omega)
             Ctilde(lp,a,b) = sum(-CC(lp,:)./(1i*deltaom-vv(lp,:)));
             %now include matsubara freq terms
             Ctilde(lp,a,b) = Ctilde(lp,a,b)+ sum(-Cmatsu(lp,:)./(1i*deltaom-vn));
             %sum is analytic for drude only and deltaom = 0; could use
             %this
             
             
             
             else
             Ctilde(lp,a,b) = sum(CC(lp,:).*(exp((1i*deltaom-vv(lp,:))*t)-1)...
                                ./(1i*deltaom-vv(lp,:)));
             Ctilde(lp,a,b) = Ctilde(lp,a,b)+ sum(Cmatsu(lp,:)...
                              .*(exp((1i*deltaom-vn)*t)-1)./(1i*deltaom-vn));                 
             end            
            end
%             real_part = real_part_fn(deltaom);
%             if ~any(isnan(real_part)) 
%             Ctilde(:,a,b) = imag(Ctilde(:,a,b)) + real_part;  
%             end
    end
end
%toc
%%
 %Gam_RF is tensor which makes up the components of the redfield tensor
 % The symbolic component from int_0^inf C_j(t) exp(iwt) is held seperately
% tic

 Gam_RF = zeros(LL,LL,LL,LL,length(Hex));
        
 %find mixings <a | Q_j |b>   
 fock_space_rep  = double(fock_space_rep); %convert from logical
 %manifold_chk = round(sum(fock_space_rep,2));
 Q_j = zeros(LL,LL,N); %mixes only within the same excitation manifold
 %as I have defined it c_nk = ex_st(k,n)
 for lp = 1:N
     Q_j(:,:,lp) = ex_st'* diag(fock_space_rep(:,lp)) * ex_st;
     %{
     for lp2 = 2:LL
         occ_L = ex_st(:,lp2).*fock_space_rep(:,lp); %occupation of this state on bra
        for lp3 = 2:LL
            occ_R = ex_st(:,lp3).*fock_space_rep(:,lp); %occupation of this state on ket
            if manifold_chk(lp2) == manifold_chk(lp3) %else no mixing
                    Q_j(lp2,lp3,lp) =  occ_L'*occ_R;
            end
        end
     end
%}
 end

%Ctilde = permute(Ctilde,[3,2,1]);

%Need to calculate <a|Q_j|b>  <c|Q_j|d> C_jj (om_dc), note that last two
%indicies of C are the sumscripts of "om"

tmp = permute(Ctilde,[3,2,1]).*Q_j; %Ctilde_{b,a,lp}
%tic
for a = 1:LL
    for b = 1:LL
         tmp2 = tmp.*repmat(Q_j(a,b,:),LL,LL,1);
         Gam_RF(a,b,:,:,:)=tmp2;
    end
end
Gam_RF = sum(Gam_RF,5); %take the sum over sites

%construct redfield tensor, won't be symbolic assuming all transition
%frequencies are given
%R_{abcd} = delta_{db} Gam_{a,e,e,c}(om_ce) + delta_{ac} Gam^*_{b,e,e,d}(om_de)
%           - Gam_{d,b,a,c}(om_ca) - Gam^*_{c,a,b,d}(om_db)
    
%from the bottom line
R_red = -(permute(Gam_RF,[3,2,4,1])+conj(permute(Gam_RF,[2,3,1,4])));  

%next consider the other two expressions with delta functions
%contract the tensor along 2nd and 3rd indicies, in diff perms
tmp = diagsum(Gam_RF,2,3); %sum over middle two dimensions
 II = eye(LL); %rep of delta functions
R_red_tmp1 = zeros(size(R_red));  R_red_tmp2 = zeros(size(R_red));
for a=1:LL     %b=a; 
    for c = 1:LL         % d=c;
         R_red_tmp1(a,:,c,:) = tmp(a,c).*II; %could also do "diag"
         R_red_tmp2(:,a,:,c) = conj(tmp(a,c)).*II;          
    end
end
 %now add togther
R_red = R_red+R_red_tmp1+R_red_tmp2;

%% express in Liouville space as super operator if required

% each density matrix element rho_{ab} decays with a term prop to
% -sum ( R_red(a,b,k1,k2) rho_{k1,k2})
% R_red(a,b,0,j) = R_red(a,b,k,0) = 0, i.e. terms which relate to ground
% state coherences etc
%
%Note reshape goes down columns then rows
%{
if nargout >= 2
R_red_op = zeros((length(Hexfull))^2); 
temp = zeros(length(Hexfull));
cnt1 = 0; cnt2 = 1;
for lp = 1:(length(Hexfull))^2
    
    cnt1 = cnt1+1;
    if cnt1 > length(Hexfull)
        cnt1=1; cnt2 = cnt2+1;
    end
    temp = temp*0;
   for a = 1:length(Hexfull)
        for b = 1:length(Hexfull)
            temp(a,b) = R_red(cnt1,cnt2,a,b);
        end
   end
    R_red_op(lp,:) = reshape(temp,1,(length(Hexfull))^2);
    
end
end
%}
if nargout >= 2
R_red_op = zeros((length(Hexfull))^2); 

for a = 1:length(Hexfull)
    for b = 1:length(Hexfull)

    temp = squeeze(R_red(a,b,:,:));
    
    lp = (b-1)*length(Hexfull)+a;

    R_red_op(lp,:) = reshape(temp,1,(length(Hexfull))^2);
    end
end
end

%{
LL = sz1*Nv;
R_red_op = zeros(LL^2); 

for a = 1:sz1
    for b = 1:sz1

    temp = squeeze(R_red(a,b,:,:));
        for av = 1:Nv
            for bv = 1:Nv
                temp2 = sparse(av,bv,1,Nv,Nv);
                temp3 = kron(temp,temp2);
                
    lp = ((b-1)*Nv + bv-1)*sz1*Nv + (a-1)*Nv + av;

    R_red_op(lp,:) = reshape(temp3.',[LL^2,1]).';
            end
        end
    end
end
%}
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
