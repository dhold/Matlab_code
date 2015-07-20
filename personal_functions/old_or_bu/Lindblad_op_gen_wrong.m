function [Lindblad_op,Lindblad_op_L,Lindblad_op_indiv] = ...
            Lindblad_op_gen(B,om_vib,numvib,gamma,disp,sz1,M_prj);
% Generates lindblad form decoherence operator.  The Hamiltonian is assumed
% to be of the form H_a tensor H_vib with H_a of size "sz1" and H_vib is
% made up of vibrations (om_vib) each truncated at a level given in the
% vector "num_vib".  gamma is the damping and B = thermo beta.
% If M_prj is given, this will project into a new basis
% L(a) V = a V a^dag - 1/2* (a^dag a V + V a^dag a)
% V L(a) = a^dag V a - 1/2* (a^dag a V + V a^dag a) = L(a^dag) V
%Lind = gamma{k}*((1+nav)*L(a) + nav*L(a^dag))
% reshape(A * C, NN_1*NN_2,1) = kron(eye(NN_2),A)*reshape(C, NN_1*NN_2,1)
% reshape(C * B, NN_1*NN_2,1) = kron(A.',eye(NN_1))*reshape(C,NN_1*NN_2,1)
% reshape(A * C * B, NN_1*NN_2,1) = kron(B.',A)*reshape(C, NN_1*NN_2,1)

if ~isempty(om_vib) && isempty(gamma)
   %alternative format pass 
    tmp = om_vib; N = length(tmp); clear om_vib
    cnt  = 0 ;
    for j = 1:N
        if size(tmp{j},2)==4
        lg = tmp{j}(:,4) ~= 0; %only these modes
        else
         lg =  true(size(tmp{j},1),1);
        end
    cnt2 = cnt + sum(double(lg));    
    %lambda(cnt+1:cnt2) = tmp{j}(lg,1); %reorg energy
    gamma(cnt+1:cnt2) = tmp{j}(lg,2); %damping factor
    om_vib(cnt+1:cnt2) = tmp{j}(lg,3); %frequency of vibrations
    numvib(cnt+1:cnt2) = tmp{j}(lg,4); %truncation level
    disp(cnt+1:cnt2) =  sqrt(2*tmp{j}(lg,1)./tmp{j}(lg,3));  %oscillator displacement 
    %from equilibrium, this operator assumes the basis used is the same for
    %all the oscillators, which is around their GS equilibrium values
    cnt = cnt2;
    end    
end

if length(sz1)==1;   pdmat = speye(sz1);
else; pdmat = sparse(sz1); sz1 = length(sz1);
end

L_so = @(aa) (kron((aa').',aa) - 1/2 *  kron(speye(size(aa)),(aa')*aa)...
               - 1/2 *  kron(((aa')*aa).',speye(size(aa))))  ;
Lindblad_op =  sparse(sz1^2*prod(numvib.^2),sz1^2*prod(numvib.^2));
if nargout > 1
Lindblad_op_L =  sparse(sz1^2*prod(numvib.^2),sz1^2*prod(numvib.^2)); %left acting
end
%H_vib_op = zeros(length(H_ex_vib));
for k = 1:length(numvib)
nav = exp(-B*om_vib(k))/(1-exp(-B*om_vib(k))); %thermal occupancy

%express the annihilation op for this mode in the full hilbert space
aa = kron(pdmat,speye(prod(numvib(1:k-1))));   Nk=numvib(k);
tmp = sparse(1:Nk-1,2:Nk,sqrt(1:Nk-1),Nk,Nk); %unshifted hat{a}
%make new basis by diagonalisation
tmp = tmp - disp*speye(Nk);
aa = kron(aa, tmp);
aa = kron(aa,speye(prod(numvib(k+1:end))));

adag = aa'; %conjugate
if nargin == 7
%project into new basis, note for excited state relaxation this will still
%relax to the state of <adag a> = K_b T / omega_0 ignoring elec coupling,
%make sure you project to the shifted basis of the GROUND STATE HO
%eigenstates
aa = M_prj'*aa*M_prj;  adag = M_prj'*adag*M_prj;
end
Lindblad_indiv =  sparse(gamma(k)*((nav+1)*L_so(aa) + nav*L_so(adag))); 

%H_vib_op  = H_vib_op  + (adag*aa+eye(size(aa))/2)*om_0{1};
Lindblad_op = Lindblad_op  + Lindblad_indiv;

if nargout>1
Lindblad_op_L =  Lindblad_op_L + gamma(k)*((nav+1)*L_so(adag) + nav*L_so(aa));
    if nargout ==3
        Lindblad_op_indiv{k} = Lindblad_indiv;
    end
end
end    
Lindblad_op = sparse(Lindblad_op);
if nargout>1
Lindblad_op_L = sparse(Lindblad_op_L);
end