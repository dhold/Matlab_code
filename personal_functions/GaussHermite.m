function [x, w] = GaussHermite(n,refine)

% This function determines the abscisas (x) and weights (w) for the
% Gauss-Hermite quadrature of order n>1, on the interval [-INF, +INF].
    % This function is valid for any degree n>=2, as the companion matrix
    % (of the n'th degree Hermite polynomial) is constructed as a
    % symmetrical matrix, guaranteeing that all the eigenvalues (roots)
    % will be real.
    
    
% © Geert Van Damme
% geert@vandamme-iliano.be
% February 21, 2010    

%Edit by me, refine is used to refine the eigenvalues to that accuracy

% Building the companion matrix CM
    % CM is such that det(xI-CM)=L_n(x), with L_n the Hermite polynomial
    % under consideration. Moreover, CM will be constructed in such a way
    % that it is symmetrical.
ii   = 1:n-1;
a   = sqrt(ii/2);
CM  = diag(a,1) + diag(a,-1);

% Determining the abscissas (x) and weights (w)
    % - since det(xI-CM)=L_n(x), the abscissas are the roots of the
    %   characteristic polynomial, i.d. the eigenvalues of CM;
    % - the weights can be derived from the corresponding eigenvectors.
[V, L]   = eig(CM);
[x, ind] = sort(diag(L));
V       = V(:,ind)';
w       = sqrt(pi) * V(:,1).^2;

if nargin == 2
    %refine should be the tolerance to refine to
if mod(n,2)==1 
x(ceil(n/2))=0; 
end

z=x(x>=0); 
success = false; 
for its=1:50 %maxit=50, usually will take 2 
p1 = pi^(-1/4); 
p2=0; 
for j=1:n %make hermite we need 
p3=p2; 
p2=p1; 
p1=z.*sqrt(2/j).*p2-sqrt((j-1)/j).*p3; 
end 
pp = sqrt(2*n).*p2; 
z1=z; 
z=z1-p1./pp; 
if all(abs(z-z1)< refine) %eps(20) is probably about correct
success = true; 
break 
end 
end 
if ~success 
warning('failed to converge to desired accuracy') 
end 
w(x>=0) = 2./pp.^2; 
w(x<=0) = flipud(w(x>=0)); 
x(x>=0) = z; 
x(x<=0) = -flipud(z);
end