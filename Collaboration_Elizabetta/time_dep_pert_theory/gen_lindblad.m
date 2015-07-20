 function Lindblad_op = gen_lindblad(beta,om_vib,numvib,M_prj,gamma,sz1,sz2)
  %%   Calculate standard decoherence term on vibration   
%L(a) V = a V a^dag - 1/2* (a^dag a V + V a^dag a)
%Lind = gamma{k}*((1+nav)*L(a) + nav*L(a^dag))
% reshape(A * C, NN_1*NN_2,1) = kron(eye(NN_2),A)*reshape(C, NN_1*NN_2,1)
% reshape(C * B, NN_1*NN_2,1) = kron(A.',eye(NN_1))*reshape(C,NN_1*NN_2,1)
% reshape(A * C * B, NN_1*NN_2,1) = kron(B.',A)*reshape(C, NN_1*NN_2,1)
L_so = @(aa) (kron((aa').',aa) - 1/2 *  kron(eye(size(aa)),(aa')*aa)...
               - 1/2 *  kron(((aa')*aa).',eye(size(aa))))  ;

Lindblad_op =  zeros(sz1^2*sz2^2);
%H_vib_op = zeros(length(H_ex_vib));
for k = 1:length(numvib)
nav = exp(-beta*om_vib(k))/(1-exp(-beta*om_vib(k)));

%express the annihilation op in the full hilbert space
aa = kron(sparse(eye(sz1)),sparse(eye(prod(numvib(1:k-1)))));
aa = kron(aa, sparse(diag(sqrt(1:numvib(k)-1),1))); 
aa = kron(aa,sparse(eye(prod(numvib(k+1:end))))); 
adag = aa'; %conjugate
%project into exciton basis (not ex vib basis!) as Redfield op in ex basis
%This shouldn't actually be needed as it doesn't change anything
%aa = M_prj'*aa*M_prj;  adag = M_prj'*adag*M_prj;

Lindblad_op_indiv =  gamma{k}*((nav+1)*L_so(aa) + nav*L_so(adag));  
%H_vib_op  = H_vib_op  + (adag*aa+eye(size(aa))/2)*om_0{1};
Lindblad_op = Lindblad_op  + Lindblad_op_indiv;
end    
Lindblad_op = sparse(Lindblad_op);