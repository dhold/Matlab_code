function partlist = partgen(N,binnum)

partlist{1} = uint16(zeros(1,binnum));

for k = 1:N
    for j = 1:floor(k/2)
   temp = partlist{k-j+1};
   %add one to every element
   
   
   
    end
    
end


for k = kstart:N
%k  %uncomment to track output
    tempmat = uint16(zeros(table(k+1),min(ceil(k/min(allowednumrange(allowednumrange~=0))),binnum))); 
    if any(allowednumrange==k)
    tempmat(1,1)=k;
    count = 1;
    else
        count = 0;
    end
    
    for j = allowednumrange(allowednumrange<=floor(k/2) & allowednumrange>0)
       temp = partlist{k-j+1};
       lg = temp>j-1 | temp == 0; %all non zero elements of temp 
       lg2 = all(lg==1,2);

       lg2 = lg2 & (any(temp==0,2)|size(temp,2)<binnum); 
       %must have at least one space for another

       temp = temp(lg2,1:min(size(temp,2),binnum-1));    
       %[temp,ones(size(temp,1),1)]

        tempmat(count+1:count+size(temp,1),1:size(temp,2)+1) = [temp,j*ones(size(temp,1),1)];

        count = count + size(temp,1);
    end

    tempmat = sort(tempmat,2,'descend');
    partlist{k+1} = tempmat(1:count,:);

end

if ~any(allowednumrange==0)  %for technical reasons this is done at the end
    %it is really better to just included zero in the allowed range
    for k =0:N
        tmp = partlist{k+1};
       partlist{k+1} =  tmp(all(tmp~=0,2),:); %only keep ones with no zeros
    end
end
end