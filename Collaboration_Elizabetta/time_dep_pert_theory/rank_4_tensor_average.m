function out = rank_4_tensor_average(a,b,c,d,e1,e2,e3,e4)
%calculates (inefficiently) rank 4 orientation averages
mat1 = [4,-1,-1;-1,4,-1;-1,-1,4]/30;
vec1 = [dot(a,b)*dot(c,d);dot(a,c)*dot(b,d);dot(a,d)*dot(b,c)];
vec2 = mat1*vec1;

if isempty(e1) %not given, give vector output based on the size of a
out = zeros(size(a));
k2rng = find(e2); k3rng = find(e3); k4rng = find(e4); 

for k1= 1:length(a)
    for k2=k2rng
        for k3 = k3rng
            for k4=k4rng

                prefct = e2(k2)*e3(k3)*e4(k4);
  kron_prefc = double([k1==k2 && k3==k4, k1==k3 && k2==k4, k1==k4 && k2==k3]);   
    out(k1) = out(k1) + prefct*kron_prefc * vec2;
    
            end
        end
    end
end
else
out = 0;
k1rng = find(e1); k2rng = find(e2); k3rng = find(e3); k4rng = find(e4); 

for k1= k1rng
    for k2=k2rng
        for k3 = k3rng
            for k4=k4rng

                prefct = e1(k1)*e2(k2)*e3(k3)*e4(k4);
  kron_prefc = double([k1==k2 && k3==k4, k1==k3 && k2==k4, k1==k4 && k2==k3]);   
    out = out + prefct*kron_prefc * vec2;
    
            end
        end
    end
end
end