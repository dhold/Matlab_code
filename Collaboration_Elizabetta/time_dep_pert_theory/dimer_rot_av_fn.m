function out = dimer_rot_av_fn(order,kres,R,mu)
%order is either 1 or 3, this is the number of interactions, the order with
%respect to k.r is taken to be one, all assumed polarized along x
%kres is resultant k vector in exp(i R  dot (kres) ) 
% for a pump probe setup the only terms that survive the RWA 
%R is therefore R_12, the displacement between sites for a dimer.  Although
%in principe crazy (2*(R_1-R_0)*(k_r+k_u)) etc terms could crop up w/o RWA
% mu is a cell array filled with the dipole components
% Can take symbolic values I suppose
if isa(kres,'sym') || isa(R,'sym') || isa(mu,'sym') 
  out = sym(zeros(3,1));  
else
out = zeros(3,1);
end
% note this function is NOT symmetric with exchange of the last 
% two variables in mu for general R and kres

if order==1 
mu1 = mu{1}; mu2 = mu{2};
out(1) =  (mu1(1)*mu2(1))/3 + (mu1(2)*mu2(2))/3 + (mu1(3)*mu2(3))/3;
%if any(kres)
out(2) =  -(kres(3)*(R(1)*mu1(2)*mu2(3) - R(1)*mu1(3)*mu2(2) - R(2)*mu1(1)*mu2(3) + R(2)*mu1(3)*mu2(1) + R(3)*mu1(1)*mu2(2) - R(3)*mu1(2)*mu2(1))*1i)/6;
out(3) =   (kres(2)*(R(1)*mu1(2)*mu2(3) - R(1)*mu1(3)*mu2(2) - R(2)*mu1(1)*mu2(3) + R(2)*mu1(3)*mu2(1) + R(3)*mu1(1)*mu2(2) - R(3)*mu1(2)*mu2(1))*1i)/6;
%end
elseif order == 3
    mu1 = mu{1}; mu2 = mu{2}; mu3 = mu{3}; mu4 = mu{4};
out(1) = (mu1(1)*mu2(1)*mu3(1)*mu4(1))/5 + (mu1(1)*mu2(1)*mu3(2)*mu4(2))/15 + (mu1(1)*mu2(2)*mu3(1)*mu4(2))/15 + (mu1(1)*mu2(2)*mu3(2)*mu4(1))/15 + (mu1(2)*mu2(1)*mu3(1)*mu4(2))/15 + (mu1(2)*mu2(1)*mu3(2)*mu4(1))/15 + (mu1(2)*mu2(2)*mu3(1)*mu4(1))/15 + (mu1(1)*mu2(1)*mu3(3)*mu4(3))/15 + (mu1(1)*mu2(3)*mu3(1)*mu4(3))/15 + (mu1(1)*mu2(3)*mu3(3)*mu4(1))/15 + (mu1(2)*mu2(2)*mu3(2)*mu4(2))/5 + (mu1(3)*mu2(1)*mu3(1)*mu4(3))/15 + (mu1(3)*mu2(1)*mu3(3)*mu4(1))/15 + (mu1(3)*mu2(3)*mu3(1)*mu4(1))/15 + (mu1(2)*mu2(2)*mu3(3)*mu4(3))/15 + (mu1(2)*mu2(3)*mu3(2)*mu4(3))/15 + (mu1(2)*mu2(3)*mu3(3)*mu4(2))/15 + (mu1(3)*mu2(2)*mu3(2)*mu4(3))/15 + (mu1(3)*mu2(2)*mu3(3)*mu4(2))/15 + (mu1(3)*mu2(3)*mu3(2)*mu4(2))/15 + (mu1(3)*mu2(3)*mu3(3)*mu4(3))/5;
%if any(kres)
out(2) =  -(kres(3)*(R(1)*mu1(2)*mu2(1)*mu3(1)*mu4(3)*1i + R(1)*mu1(2)*mu2(1)*mu3(3)*mu4(1)*1i + R(1)*mu1(2)*mu2(3)*mu3(1)*mu4(1)*1i - R(1)*mu1(3)*mu2(1)*mu3(1)*mu4(2)*1i - R(1)*mu1(3)*mu2(1)*mu3(2)*mu4(1)*1i - R(1)*mu1(3)*mu2(2)*mu3(1)*mu4(1)*1i - R(2)*mu1(1)*mu2(1)*mu3(1)*mu4(3)*1i - R(2)*mu1(1)*mu2(1)*mu3(3)*mu4(1)*1i - R(2)*mu1(1)*mu2(3)*mu3(1)*mu4(1)*1i + R(2)*mu1(3)*mu2(1)*mu3(1)*mu4(1)*3*1i + R(3)*mu1(1)*mu2(1)*mu3(1)*mu4(2)*1i + R(3)*mu1(1)*mu2(1)*mu3(2)*mu4(1)*1i + R(3)*mu1(1)*mu2(2)*mu3(1)*mu4(1)*1i - R(3)*mu1(2)*mu2(1)*mu3(1)*mu4(1)*3*1i + R(1)*mu1(2)*mu2(2)*mu3(2)*mu4(3)*1i + R(1)*mu1(2)*mu2(2)*mu3(3)*mu4(2)*1i + R(1)*mu1(2)*mu2(3)*mu3(2)*mu4(2)*1i - R(1)*mu1(3)*mu2(2)*mu3(2)*mu4(2)*3*1i - R(2)*mu1(1)*mu2(2)*mu3(2)*mu4(3)*1i - R(2)*mu1(1)*mu2(2)*mu3(3)*mu4(2)*1i - R(2)*mu1(1)*mu2(3)*mu3(2)*mu4(2)*1i + R(2)*mu1(3)*mu2(1)*mu3(2)*mu4(2)*1i + R(2)*mu1(3)*mu2(2)*mu3(1)*mu4(2)*1i + R(2)*mu1(3)*mu2(2)*mu3(2)*mu4(1)*1i + R(3)*mu1(1)*mu2(2)*mu3(2)*mu4(2)*3*1i - R(3)*mu1(2)*mu2(1)*mu3(2)*mu4(2)*1i - R(3)*mu1(2)*mu2(2)*mu3(1)*mu4(2)*1i - R(3)*mu1(2)*mu2(2)*mu3(2)*mu4(1)*1i + R(1)*mu1(2)*mu2(3)*mu3(3)*mu4(3)*3*1i - R(1)*mu1(3)*mu2(2)*mu3(3)*mu4(3)*1i - R(1)*mu1(3)*mu2(3)*mu3(2)*mu4(3)*1i - R(1)*mu1(3)*mu2(3)*mu3(3)*mu4(2)*1i - R(2)*mu1(1)*mu2(3)*mu3(3)*mu4(3)*3*1i + R(2)*mu1(3)*mu2(1)*mu3(3)*mu4(3)*1i + R(2)*mu1(3)*mu2(3)*mu3(1)*mu4(3)*1i + R(2)*mu1(3)*mu2(3)*mu3(3)*mu4(1)*1i + R(3)*mu1(1)*mu2(2)*mu3(3)*mu4(3)*1i + R(3)*mu1(1)*mu2(3)*mu3(2)*mu4(3)*1i + R(3)*mu1(1)*mu2(3)*mu3(3)*mu4(2)*1i - R(3)*mu1(2)*mu2(1)*mu3(3)*mu4(3)*1i - R(3)*mu1(2)*mu2(3)*mu3(1)*mu4(3)*1i - R(3)*mu1(2)*mu2(3)*mu3(3)*mu4(1)*1i))/30;
out(3) =   (kres(2)*(R(1)*mu1(2)*mu2(1)*mu3(1)*mu4(3)*1i + R(1)*mu1(2)*mu2(1)*mu3(3)*mu4(1)*1i + R(1)*mu1(2)*mu2(3)*mu3(1)*mu4(1)*1i - R(1)*mu1(3)*mu2(1)*mu3(1)*mu4(2)*1i - R(1)*mu1(3)*mu2(1)*mu3(2)*mu4(1)*1i - R(1)*mu1(3)*mu2(2)*mu3(1)*mu4(1)*1i - R(2)*mu1(1)*mu2(1)*mu3(1)*mu4(3)*1i - R(2)*mu1(1)*mu2(1)*mu3(3)*mu4(1)*1i - R(2)*mu1(1)*mu2(3)*mu3(1)*mu4(1)*1i + R(2)*mu1(3)*mu2(1)*mu3(1)*mu4(1)*3*1i + R(3)*mu1(1)*mu2(1)*mu3(1)*mu4(2)*1i + R(3)*mu1(1)*mu2(1)*mu3(2)*mu4(1)*1i + R(3)*mu1(1)*mu2(2)*mu3(1)*mu4(1)*1i - R(3)*mu1(2)*mu2(1)*mu3(1)*mu4(1)*3*1i + R(1)*mu1(2)*mu2(2)*mu3(2)*mu4(3)*1i + R(1)*mu1(2)*mu2(2)*mu3(3)*mu4(2)*1i + R(1)*mu1(2)*mu2(3)*mu3(2)*mu4(2)*1i - R(1)*mu1(3)*mu2(2)*mu3(2)*mu4(2)*3*1i - R(2)*mu1(1)*mu2(2)*mu3(2)*mu4(3)*1i - R(2)*mu1(1)*mu2(2)*mu3(3)*mu4(2)*1i - R(2)*mu1(1)*mu2(3)*mu3(2)*mu4(2)*1i + R(2)*mu1(3)*mu2(1)*mu3(2)*mu4(2)*1i + R(2)*mu1(3)*mu2(2)*mu3(1)*mu4(2)*1i + R(2)*mu1(3)*mu2(2)*mu3(2)*mu4(1)*1i + R(3)*mu1(1)*mu2(2)*mu3(2)*mu4(2)*3*1i - R(3)*mu1(2)*mu2(1)*mu3(2)*mu4(2)*1i - R(3)*mu1(2)*mu2(2)*mu3(1)*mu4(2)*1i - R(3)*mu1(2)*mu2(2)*mu3(2)*mu4(1)*1i + R(1)*mu1(2)*mu2(3)*mu3(3)*mu4(3)*3*1i - R(1)*mu1(3)*mu2(2)*mu3(3)*mu4(3)*1i - R(1)*mu1(3)*mu2(3)*mu3(2)*mu4(3)*1i - R(1)*mu1(3)*mu2(3)*mu3(3)*mu4(2)*1i - R(2)*mu1(1)*mu2(3)*mu3(3)*mu4(3)*3*1i + R(2)*mu1(3)*mu2(1)*mu3(3)*mu4(3)*1i + R(2)*mu1(3)*mu2(3)*mu3(1)*mu4(3)*1i + R(2)*mu1(3)*mu2(3)*mu3(3)*mu4(1)*1i + R(3)*mu1(1)*mu2(2)*mu3(3)*mu4(3)*1i + R(3)*mu1(1)*mu2(3)*mu3(2)*mu4(3)*1i + R(3)*mu1(1)*mu2(3)*mu3(3)*mu4(2)*1i - R(3)*mu1(2)*mu2(1)*mu3(3)*mu4(3)*1i - R(3)*mu1(2)*mu2(3)*mu3(1)*mu4(3)*1i - R(3)*mu1(2)*mu2(3)*mu3(3)*mu4(1)*1i))/30;
%end
end
end