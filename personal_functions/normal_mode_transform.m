function proj_mat = normal_mode_transform(N)
%function returns matrix which projects to normal mode coordinates

proj_mat = zeros(N);

proj_mat(:,1) = 1/sqrt(N);

for k = 2:N
    
    proj_mat(k,k) = 1;
    proj_mat(1:k-1,k) = -1/(k-1);
    
    proj_mat(:,k) = proj_mat(:,k)*sqrt((k-1)/k);
    
end
end