%this was used to test the matrix element generator for the coupling terms
%in the heirarchy
% This outputs

tic

if isempty(Kap1)
    Kap1 = inf; %no truncation based on Kap1
end
    
cc = [cc2R,cc1]; ccI = [cc2I,cc1*0];   vv = [vv2, vv1];
N2 = size(cc,2); N = size(cc,1);

cc =  reshape(cc.',1,numel(cc)); 
ccI =  reshape(ccI.',1,numel(ccI));
vv = reshape(vv.',1,numel(vv));

if Kap2>0  %else trivial

numwithn = zeros(1,Kap2+1); %size(cc2,2) is number of poles in J(omega)
for kk = 0:Kap2
numwithn(kk+1) = nchoosek((Kappa+size(cc2R,2))*N-1+kk,kk); 
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+size(cc2,2)))
nn{kk+1} = zeros(numwithn(kk+1),kk); %first is empty
end

 ncnt = cumsum(numwithn); 
 %hashfn = uint64(length(cc)+1).^uint64(1:Kap2);
 hashfn = ((length(cc)+1).^(1:Kap2)).';
 %these are used to remove repeated elements
 clear coup_p1 coup_m1 coup_m1_acom
 for j =1:N  %this method may be memory harsh
coup_p1{j} = zeros(ncnt(end));  %up coupling coefficient
coup_m1{j} =  zeros(ncnt(end)) ; %down coupling commutator coefficient
coup_m1_acom{j} =  zeros(ncnt(end));  %down coupling anti commutator coefficent
 end
totfreq = real(reshape(vv.',numel(vv),1));
count1 = 1; %last element used as a generator
%cn1 = count1; %same thing within tier
count2 = 1; %last element produced
%cn2 = count2; %same thing within tier

adding_vec = zeros(length(cc),1);
adding_vec(1:end) = 1:length(cc); current_tier =[];
for k =1:Kap2
        
        cn1 = 1; cn2 = 0; %these reset each time
        lower_tier  = current_tier;
        current_tier = nn{k+1};
        lgtest = repmat(adding_vec,1,k);
        
        if ~isempty(lower_tier)
        tmp = repmat(lower_tier(1,:),length(cc),1);
        else
        tmp = zeros(length(cc),0);
        end
        tmp2 = sort([tmp , adding_vec],2); %generate new element in tier above
        uniquetest = uint64(zeros(numwithn(k+1),1)); 
        %uniquetest(1:length(cc)) = uint64(tmp2) * hashfn(1:k);
        uniquetest(1:length(cc)) = uint64(tmp2 * hashfn(1:k));
        
        %must be unique as it is the first generated
        current_tier(1:size(tmp2,1),:) = sort(tmp2,2); %sort order
        %each of these elements has an assosiated upward and downward
        %coupling
        tmp3 = sqrt(sum(tmp2 == lgtest,2)); %sqrt(n_{k,j}+1)
        scale_fct = sqrt(abs(cc) + abs(ccI));
        map2 = count2+1:count2+length(cc);
        for j=1:N
            rng = (1+(j-1)*N2:j*N2);
        coup_p1{j}(count1,map2(rng)) = tmp3(rng).' .* scale_fct(rng) ;
        % scaling w.r.t. sqrt(abs(cc) + abs(ccI)) is somewhat unconvential 
        % down coupling at the transpose of these positions,
        % Note for matsubara terms these will be the same
        coup_m1{j}(map2(rng),count1) = tmp3(rng).' .* cc(rng) ./ scale_fct(rng); 
        coup_m1_acom{j}(map2(rng),count1) = tmp3(rng).' .* ccI(rng) ./ scale_fct(rng);   
        end                    
           count1=count1+1;  cn1 = cn1+1;
           count2 = count2 + length(cc);   cn2 = cn2+length(cc);  
           
        for kk = 2:size(lower_tier,1)
            
        tmp = repmat(lower_tier(kk,:),length(cc),1);
        tmp2 = sort([tmp , adding_vec],2); 
        tmp3 = sqrt(sum(tmp2 == lgtest,2)); 
        %test which of these elements create a unique element
        %uniquetmp = uint64(tmp2) * hashfn(1:k);
        uniquetmp = uint64(tmp2 * hashfn(1:k));
        %must be unique as it is the first generated
        uniquelg = true(size(uniquetmp)); 
        map2 = (1+count2):(count2+size(tmp2,1));
            for jj = 1:length(cc)
                tmplg = uniquetmp(jj) == uniquetest(1:cn2);
                uniquelg(jj) = ~any(tmplg);
                if ~uniquelg(jj) %will map to another element alread there
                    map2(jj) = ncnt(k) + find(tmplg); 
                    %uniquelg only measures from last tier
                end
            end

        tmp2 = tmp2(uniquelg,:);
        current_tier((cn2+1):(cn2 + size(tmp2,1)),:) = tmp2;
        uniquetest((cn2+1):(cn2+size(tmp2,1))) = uniquetmp(uniquelg);
        
        scale_fct = sqrt(abs(cc) + abs(ccI));
        for j=1:N
            rng = (1+(j-1)*N2:j*N2);
        coup_p1{j}(count1,map2(rng)) = tmp3(rng).' .* scale_fct(rng) ;
        coup_m1{j}(map2(rng),count1) = tmp3(rng).' .* cc(rng) ./ scale_fct(rng); 
        coup_m1_acom{j}(map2(rng),count1) = tmp3(rng).' .* ccI(rng) ./ scale_fct(rng);   
        end                   
            
            count1 = count1 + 1; count2 = count2+size(tmp2,1);
            cn1=cn1+1; cn2=cn2+size(tmp2,1);
        end
        
        nn{k+1} = current_tier;

end
else
    nn{1} =[]; count2=1;
end

const_factor = zeros(count2,1); count1 = 1;
for j= 1:Kap2
   
    tmp = nn{j+1};   tmp2 = zeros(size(nn{j+1},1),1);
    for jj = 1:j
        tmp2 = tmp2 + vv(tmp(:,j)).';
    end
    
    const_factor((count1+1):(count1+size(tmp2,1))) = tmp2;
    count1=count1 +size(tmp2,1) ;
end

toc