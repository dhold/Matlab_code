
test = rand(3); test = test+test'; test2 = zeros(3);
for k =1:3
    for j = 1:3
        
        test2 = test2 + R_red(:,:,k,j)*test(k,j);
        
    end
end
R_redd = permute(reshape(R_red,3,3,9),[2,3,1]);
test3 = mtimesx(R_redd,'G',reshape(test,9,1));
test3 = reshape(test3,3,3);
test2 - test3
%%

R_redd2 = permute(reshape(R_red2,4,4,16),[2,3,1]);
test = rand(4); test = test+test'; test2 = zeros(4);
for k =1:4
    for j = 1:4
        
        test2 = test2 + R_red2(:,:,k,j)*test(k,j);
        
    end
end

test3 = mtimesx(R_redd2,'G',reshape(test,16,1));
test2 - reshape(test3,4,4)