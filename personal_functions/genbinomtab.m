function binomvtot = genbinomtab(nmax,~)

%this function generates all the possible binomial coefficients of the form
% nchoosek(n,k) with n <= nmax note that only k <= n are considered
% n and k both start from zero but matlab arrays start from one hence
% nchoosek(n,k) = binomvtot(n+1,k+1)
% if a second argument is given it will generate the log of this function
binomvtot = zeros(nmax+1); 

if nargin==1

 %use n choose k = [(n-1) choose k-1 ] + [n choose k -1 ]
binomvtot(:,1)=1;
 for n = 2:nmax+1  %note n=a but index is offset by 1
        binomvtot(n,2:n) = binomvtot(n-1,2:n)+binomvtot(n-1,1:n-1);
 end
 
else %take logs, note that biomial coeffs are given by exp of this
     %this uses the recurrsion relation valid for k \le n-1
 % nchoosek(n,k) = nchoosek(n-1,k) *n/(n-k)    n/(n-k) = 1/(1-k/n)
 % and nchoosek(a,a) == 1

 for n = 1:nmax  %note n=a but index is offset by 1
    binomvtot(n+1,1:n) = binomvtot(n,1:n)-log(1:-1/n:1/n);   
     %binomvtot(n+1,n+1:end) =-inf;  
 end        
  binomvtot=  binomvtot-tril(inf*ones(nmax+1),-1).'; %these elements are zero
end