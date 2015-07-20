

LL = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));

QQ = LL*0;

%Kappa = truncation parameter, assumed that for any e^(-v_m)t with m>Kappa
%will be ~delta(t)/v_m and treated markovianly
% treat_markov logical values length of cc2 of the values which will not be
% included in the heirarchy

lg = imag(vv2)~=0; %find complex frequencies

%if include_truncation_correction
for j = 1:N
    
    Vj = zeros(N);
    Vj(j,j) = 1;  Vj = kron(Vj,eye(viblvls)); eyeQ = eye(length(Vj));
      
    QQ = QQ + dot(cc1(j,Kappa+1:end),1./vv1(Kappa+1,2:end))...
        .*(kron(Vj.',eyeQ) + kron(eyeQ,Vj) - 2*kron(Vj.',Vj));
    %just take a really big sum to deal with the truncation, i.e. a lot of
    %cc1 and vv1 (converges pretty fast) that I then truncate
    
    %reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
    %[Q_i,[Q_i,rho_n]] = (1 X Q_i + Q_i^T X 1 -2 Q_i^T X Q_i ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
    
     QQ = QQ  + dot(cc2R(j,treat_markov),1./real(vv2(j,treat_markov)))...
        .*(kron(Vj.',eyeQ) + kron(eyeQ,Vj) - 2*kron(Vj.',Vj)); 
     QQ = QQ  - 1i*dot(cc2I(j,treat_markov),1./real(vv2(j,treat_markov)))...
                .*(kron(Vj.',eyeQ) + kron(eyeQ,Vj)); 
    
    
end
%end