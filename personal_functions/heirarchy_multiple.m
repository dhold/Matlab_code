function out = heirarchy_multiple(A,B,nn,opt)
%this multiplies B by the matrix A, this assumes that B is density matrix
%vector for the HEOM which is flattened out such that 
% B = [rho_000;rho_100;rho_010;...] with any extra dimensions including
% time dependence or some shit
szB = size(B);

use_mtimesx = opt(1);  conj_mult = opt(2); 
keep_orig_shape = opt(3);

sz1 = szB(1)/size(nn,1); %length of hilbert space squared

if keep_orig_shape
    out = zeros(size(B));
else
    out = zeros([size(nn,1)*sqrt(sz1),sqrt(sz1),szB(2:end)]);
end

for lp = 1:size(nn,1) %normal matrix multiplication for each one
    
    Btmp = reshape(B(sz1*(lp-1)+1:sz1*lp,:),[sqrt(sz1),sqrt(sz1),szB(2:end)]);
    
    
        if use_mtimesx
            if conj_mult
            tmp = mtimesx(Btmp,A);
            else
            tmp = mtimesx(A,Btmp);    
            end
        else
            if conj_mult
            tmp = Btmp*A; 
            else
            tmp = A*Btmp;     
            end
        end
    if keep_orig_shape
    out((lp-1)*sz1+1:lp*sz1,:) = tmp;
    else
    out((lp-1)*sqrt(sz1)+1:lp*sqrt(sz1),:) = tmp;
    end    
end


