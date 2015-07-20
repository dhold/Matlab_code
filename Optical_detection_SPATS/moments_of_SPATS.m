%moments of single photon added thermal state
% rho = (a^{dag})^m exp(-beta a^{dag} a ) a^m)
% with beta = const/T and <n> = 1/(exp(beta)-1) 
  num_mom = 3;  beta = 1;
  num_mom_tot = 2*(num_mom-1);
  n_av = (exp(beta)-1)^(-1); %thermal average photon number
  m = 1; %defined as part of state
  binom = genbinomtab(num_mom_tot);
  %value of <a^p a^{dag}^p > = ((<n>+1)^p*(m+p)!/m!)
  apadp =zeros(1,num_mom_tot+1); apadp(1) = 1;
  for p = 1:num_mom_tot
      
     apadp(p+1) =  (n_av+1)^p*factorial(m+p)/factorial(m);
      
  end

  m_k = zeros(1,num_mom_tot+1); m_k(1) = 1; %identity exp value
for p = 1:num_mom_tot
    for r= 0:p
        m_k(p+1) = m_k(p+1) + (-1)^(r+p)*factorial(p)^2*(apadp(r+1))/...
                                (factorial(r)^2*factorial(p-r));
    end
end

mm = zeros(num_mom);
for k = 1:num_mom
    mm(k,:) = m_k(k:k+num_mom-1);
end

% to get <(a^{dag} a)^(p+1) > = < a^{dag} sum_{k=0}^p nchoosek(p+1,p-k) *
%                                   (a^{dag} a)^(p)   a >

bb = zeros(num_mom_tot+1);  bb(1,1) = 1;
%  <(a^{dag} a)^n > = b_{0n} + b_{1n} <a^{dag}^1 a^1 > +...
%                   b_{jn} <a^{dag}^j a^j > + ... +b_{nn} <a^{dag}^n a^n >
% bb are the coefficients b, but 
for n = 1:num_mom_tot
   for j = 1:n
       %b_{j,n} = sum_{k=j-1}^n-1 (n chose n-1-k) b_{jk}
      bb(j+1,n+1) = dot(binom(n+1,n+1-(j:n)) ,bb(j,j:n));
   end    
end
mu_k = zeros(1,num_mom_tot+1); mu_k(1) = 1;
for p = 1:num_mom_tot
    mu_k(p+1) = dot(bb(:,p+1),m_k);
end

mumu = zeros(num_mom);
for k = 1:num_mom
    mumu(k,:) = mu_k(k:k+num_mom-1);
end

A_3 = det(mm(1:3,1:3))/(det(mumu(1:3,1:3))-det(mm(1:3,1:3)));