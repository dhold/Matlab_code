NN=5;

temp = rand(NN);
temp2 = rand(NN);

ans1 = reshape(temp*temp2,NN^2,1);
ans11 = kron(eye(NN),temp)*reshape(temp2,NN^2,1);
[ans1-ans11]

ans2 = reshape(temp2*temp,NN^2,1);
ans22 = kron(temp.',eye(NN))*reshape(temp2,NN^2,1);
[ans2-ans22]

ans3 = reshape(temp*temp2*temp,NN^2,1);
ans33 =  kron(temp.',temp)*reshape(temp2,NN^2,1);
[ans3-ans33]

%Therefore for N X N A and B

% reshape(A * B, N^2,1) = kron(eye(N),A)*reshape(B, N^2,1)
% reshape(B * A, N^2,1) = kron(A.',eye(N))*reshape(B, N^2,1)
% reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
%%
ans1 = reshape((temp*temp2).',NN^2,1);
ans11 = kron(temp,eye(NN))*reshape(temp2.',NN^2,1);
[ans1-ans11]

ans2 = reshape((temp2*temp).',NN^2,1);
ans22 = kron(eye(NN),temp.')*reshape(temp2.',NN^2,1);
[ans2-ans22]

ans3 = reshape((temp*temp2*temp).',NN^2,1);
ans33 =  kron(temp,temp.')*reshape(temp2.',NN^2,1);
[ans3-ans33]

%Therefore for N X N A and B

% reshape((A * B).', N^2,1) = kron(A,eye(N))*reshape(B.', N^2,1)
% reshape((B * A).', N^2,1) = kron(eye(N),A')*reshape(B.', N^2,1)
% reshape((A * B * A).', N^2,1) = kron(A,A.')*reshape(B.', N^2,1)

%%
ans1 = reshape(temp*(temp*temp2),NN^2,1);
ans11 = kron(eye(NN),temp^2)*reshape(temp2,NN^2,1);
[ans1-ans11]

ans2 = reshape((temp2*temp)*temp,NN^2,1);
ans22= kron((temp^2).',eye(NN))*reshape(temp2,NN^2,1);
[ans2-ans22]