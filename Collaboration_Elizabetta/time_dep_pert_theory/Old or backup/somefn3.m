test_var = sym('A',[sz2*N,sz2]); test_var2 = sym('B',[sz2,sz2*N]);
test_var_full = sym(zeros(sz2*sz1)); test_var_full2 = sym(zeros(sz2*sz1));
test_var_full(sz2+1:(sz2*(N+1)),1:sz2) = test_var;
test_var_full2(1:sz2,sz2+1:(sz2*(N+1))) = test_var2;

test2a = supop*reshape(test_var_full,sz1*sz2*sz1*sz2,1);
test2b = supop*reshape(test_var_full2,sz1*sz2*sz1*sz2,1);

lga=test2a==0;  lgb=test2b==0;
test3a = test2a(~lga);  test3b = test2b(~lgb); 
%test_var_rs = reshape(test_var,sz2*N*sz2,1);
test4a = matlabFunction(test3a,'vars',{test_var});
test4b = matlabFunction(test3b,'vars',{test_var2});

test5a = mu_hilb{2}(sz2+1:(sz2*(N+1)),1:sz2); 
test6a = zeros(size(lg)); 
test6a(~lga) = test4a(test5a);
test5b = mu_hilb{2}(1:sz2,sz2+1:(sz2*(N+1))); 
test6b = zeros(size(lg)); 
test6b(~lgb) = test4b(test5b);

sum(sum(abs(reshape(test6a,sz1*sz2,sz1*sz2)-reshape(test6b,sz1*sz2,sz1*sz2)')))
