function coeff = coeffmat_gen(nn,jj,A0in)

persistent coeffmat A_0 

if nargin == 3
%recursively generate coeffmat
% if nn > length(coeffmat)+1 || jj > length(coeffmat)+1 || ~isequal(A_0,A0in)
% clear coeffmat
% end
A_0 = A0in; %past sum(cc,2) for this one; %in general length N vector

coeffmat = -inf*ones(max(nn+1,jj+1)); %set all elements to - inf to see which aren't generated
coeffmat(1,1) = 1; coeffmat(1,2:end) = 0;
%recursion relation c(n + 1,j + 1) = -c(n,j) + (n-1) A_0 c(n-1,j+1)
coeffmat(2,2) = -1;  coeffmat(3,3) = 1;
coeffmat(3,1) = A_0;
for nn2=nn:-1:0
    for jj2 = jj:-1:0
    coeffmat_gen(nn2,jj2);
    end
end
    coeff = coeffmat;
    return
end
if nn>=0 && jj>=0
   if isfinite(coeffmat(nn+1,jj+1))
       coeff = coeffmat(nn+1,jj+1); return
   end
else %at least matrix element is negative
    coeff = 0; %not is general true but probably sufficient
    return
end


if nn >=2 && jj >= 1
coeff = -coeffmat_gen(nn-1,jj-1) + (nn-1)*A_0*coeffmat_gen(nn-2,jj);
elseif nn >= 2 && jj<1
coeff = (nn-1)*A_0*coeffmat_gen(nn-2,jj);    
elseif nn ==1 && jj>=1
   coeff = -coeffmat_gen(nn-1,jj-1) ; 
else
    coeff = 0;
end
coeffmat(nn+1,jj+1) = coeff;

end
