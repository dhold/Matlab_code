function out = orientation_tensor_average(A,E)

szs = cellfun(@size,A,'UniformOutput',0);
szs2 = cellfun(@size,E,'UniformOutput',0);
rank = 0; rank2 = 0;

for k =1:length(szs)
    if any(szs{k}==0)
    else
        rank = rank +sum(double(szs{k}~=1));
    end
    if any(szs{k}==0)
    else
        rank2 = rank2 +sum(double(szs2{k}~=1));
    end    
    
end

if rank == 2
    
    if length(szs)==2 %two vectors
    out = dot(A{1},A{2})*dot(E{1},E{2})/3;
    else
    out = 'write something for single rank2 tensor average';   
    end

elseif rank == 3
    
    out = 'write something for rank 3'
elseif rank ==4 && length(szs)==4
    
    mat1 = [4,-1,-1;-1,4,-1;-1,-1,4]/30;
vec1 = [dot(A{1},A{2})*dot(A{3},A{4});dot(A{1},A{3})*dot(A{2},A{4});dot(A{1},A{4})*dot(A{2},A{3})];
vec2 = mat1*vec1;
out = 0;

k1rng = find(E{1}); k2rng = find(E{2}); k3rng = find(E{3}); k4rng = find(E{4}); 

for k1= k1rng
    for k2=k2rng
        for k3 = k3rng
            for k4=k4rng

                prefct = E{1}(k1)*E{2}(k2)*E{3}(k3)*E{4}(k4);
  kron_prefc = double([k1==k2 && k3==k4, k1==k3 && k2==k4, k1==k4 && k2==k3]);   
    out = out + prefct*kron_prefc * vec2;
    
            end
        end
    end
end

elseif rank == 5
       out = 'write something for rank 5'
else
    rank
    error('no routine for tensor of this rank')
    
end