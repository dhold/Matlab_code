function [total_prop_op]=H_prop_gen...
            (H_set,V_set,QQ,const_factor,coup_com,coup_acom)
   
H_g = H_set{1};  H_e = H_set{2};  H_f = H_set{3};
V_g = V_set{1};  V_e = V_set{2};  V_f = V_set{3};
       clear H_set V_set
        
%% Construct Ltot, the total Louville for self coupling density matricies

% For N X N A and B we have the following
% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho

 Lgg = -1i*(kron(sparse(eye(length(H_g))),sparse(H_g))...
            -kron(sparse(H_g).',sparse(eye(length(H_g)))) );
 Lee = -1i*(kron(sparse(eye(length(H_e))),sparse(H_e))...
            -kron(sparse(H_e).',sparse(eye(length(H_e)))) );
 Lff = -1i*(kron(sparse(eye(length(H_f))),sparse(H_f))...
            -kron(sparse(H_f).',sparse(eye(length(H_f)))) );        
 Lge = -1i*(kron(sparse(eye(length(H_e))),sparse(H_g))...
            -kron(sparse(H_e).',sparse(eye(length(H_g)))) );
 Lef = -1i*(kron(sparse(eye(length(H_f))),sparse(H_e))...
            -kron(sparse(H_f).',sparse(eye(length(H_e)))) );     
 Lgf = -1i*(kron(sparse(eye(length(H_f))),sparse(H_g))...
            -kron(sparse(H_f).',sparse(eye(length(H_g)))) );            
 
       
      %include decay factor first  
sz_gg = length(coup_com{1})*length(Lgg);  temp = kron(const_factor,ones(length(Lgg),1));  
prop_gg = - sparse(1:length(temp),1:length(temp),temp);
sz_ge = length(coup_com{1})*length(Lge);  temp = kron(const_factor,ones(length(Lge),1)); 
prop_ge = - sparse(1:length(temp),1:length(temp),temp);
sz_ee = length(coup_com{1})*length(Lee);  temp = kron(const_factor,ones(length(Lee),1)); 
prop_ee = - sparse(1:length(temp),1:length(temp),temp);
sz_ef = length(coup_com{1})*length(Lef);  temp = kron(const_factor,ones(length(Lef),1)); 
prop_ef = - sparse(1:length(temp),1:length(temp),temp);
sz_ff = length(coup_com{1})*length(Lff);  temp = kron(const_factor,ones(length(Lff),1)); 
prop_ff = - sparse(1:length(temp),1:length(temp),temp);
sz_gf = length(coup_com{1})*length(Lgf);  temp = kron(const_factor,ones(length(Lgf),1)); 
prop_gf = - sparse(1:length(temp),1:length(temp),temp); 
% QQ as input should be the constants in the correction 
if size(QQ,1) == 1 %assume same for every site
    QQ = repmat(QQ,N*viblvls,1);
end
 Q = LL*0;
% 
% %if include_truncation_correction

for j = 1:N %site loop
    
    
    %reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
    %[Q_i,[Q_i,rho_n]] = (1 X Q_i*Q_i + Q_i^T*Q_i^T X 1 - 2 Q_i^T X Q_i ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
    %[Q_uv ,rho_uv] = Quu rho_uv - rho_uv Q_vv 
    
    Qcom = -1i*(kron(eyeQ,Qjsec) - kron(Qjsec.',eyeQ)); %commtator
    %keep factor of -1i from -1i tilde(V)^X term
    Qacom =  (kron(eyeQ,Qjsec) + kron(Qjsec.',eyeQ)); %anti commutator
    %factor of -1i cancels with +1i in front of imag part
    Q = Q + QQ(j,1) .*Qcom*Qcom + QQ(j,2) .*Qacom*Qcom;
  
    prop_ge = prog_ge + kron(sparse(coup_com{j}),Qcom) + kron(sparse(coup_acom{j}),Qacom);
    prop_ee = prog_ee + kron(sparse(coup_com{j}),Qcom) + kron(sparse(coup_acom{j}),Qacom);
    prop_ff = prog_ff + kron(sparse(coup_com{j}),Qcom) + kron(sparse(coup_acom{j}),Qacom);
    prop_ge = prog_ge + kron(sparse(coup_com{j}),Qcom) + kron(sparse(coup_acom{j}),Qacom);
    prop_ge = prog_ge + kron(sparse(coup_com{j}),Qcom) + kron(sparse(coup_acom{j}),Qacom);
    
    
    total_prop_op = total_prop_op + kron(sparse(coup_com{j}),Qcom)...
                        + kron(sparse(coup_acom{j}),Qacom);
    
    
    % .*(kron(Vj.',eyeQ) + kron(eyeQ,Vj)); 
end

 
    
