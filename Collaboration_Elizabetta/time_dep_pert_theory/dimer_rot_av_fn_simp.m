function out = dimer_rot_av_fn_simp(order,kres,R,mu)
%order is either 1 or 3, this is the number of interactions, the order with
%respect to k.r is taken to be one, all assumed polarized along x
%kres is resultant k vector in exp(i [0,0;R]T' * (kres) ) 
% R is the length of R_12
% Can take symbolic values I suppose
% I should probably also calculate the averages with the different
% interactions being along Y as well as just x
if isa(kres,'sym') || isa(R,'sym') || isa(mu,'sym') 
  out = sym(zeros(3,1));  
else
out = zeros(3,1);
end
if length(order)==1
    order(2) = 1;
end
% note this function is NOT exactly symmetric with exchange of the last 
% two variables in mu for general R and kres
if order(2) == 1;
    if order(1)==1 
mu1 = mu{1}; mu2 = mu{2};
out(1) =  dot(mu1,mu2)/3;
out(2) =   -(R*kres(3)*(mu1(1)*mu2(2) - mu1(2)*mu2(1))*1i)/6;
out(3) =    (R*kres(2)*(mu1(1)*mu2(2) - mu1(2)*mu2(1))*1i)/6;
%end
    elseif order == 3
    mu1 = mu{1}; mu2 = mu{2}; mu3 = mu{3}; mu4 = mu{4};
  out(1) =  (mu1(1)*mu2(1)*mu3(1)*mu4(1))/5 + (mu1(1)*mu2(1)*mu3(2)*mu4(2))/15 + (mu1(1)*mu2(2)*mu3(1)*mu4(2))/15 + (mu1(1)*mu2(2)*mu3(2)*mu4(1))/15 + (mu1(2)*mu2(1)*mu3(1)*mu4(2))/15 + (mu1(2)*mu2(1)*mu3(2)*mu4(1))/15 + (mu1(2)*mu2(2)*mu3(1)*mu4(1))/15 + (mu1(1)*mu2(1)*mu3(3)*mu4(3))/15 + (mu1(1)*mu2(3)*mu3(1)*mu4(3))/15 + (mu1(1)*mu2(3)*mu3(3)*mu4(1))/15 + (mu1(2)*mu2(2)*mu3(2)*mu4(2))/5 + (mu1(3)*mu2(1)*mu3(1)*mu4(3))/15 + (mu1(3)*mu2(1)*mu3(3)*mu4(1))/15 + (mu1(3)*mu2(3)*mu3(1)*mu4(1))/15 + (mu1(2)*mu2(2)*mu3(3)*mu4(3))/15 + (mu1(2)*mu2(3)*mu3(2)*mu4(3))/15 + (mu1(2)*mu2(3)*mu3(3)*mu4(2))/15 + (mu1(3)*mu2(2)*mu3(2)*mu4(3))/15 + (mu1(3)*mu2(2)*mu3(3)*mu4(2))/15 + (mu1(3)*mu2(3)*mu3(2)*mu4(2))/15 + (mu1(3)*mu2(3)*mu3(3)*mu4(3))/5;
  out(2) =  -(R*kres(3)*(mu1(1)*mu2(1)*mu3(1)*mu4(2) + mu1(1)*mu2(1)*mu3(2)*mu4(1) + mu1(1)*mu2(2)*mu3(1)*mu4(1) - 3*mu1(2)*mu2(1)*mu3(1)*mu4(1) + 3*mu1(1)*mu2(2)*mu3(2)*mu4(2) - mu1(2)*mu2(1)*mu3(2)*mu4(2) - mu1(2)*mu2(2)*mu3(1)*mu4(2) - mu1(2)*mu2(2)*mu3(2)*mu4(1) + mu1(1)*mu2(2)*mu3(3)*mu4(3) + mu1(1)*mu2(3)*mu3(2)*mu4(3) + mu1(1)*mu2(3)*mu3(3)*mu4(2) - mu1(2)*mu2(1)*mu3(3)*mu4(3) - mu1(2)*mu2(3)*mu3(1)*mu4(3) - mu1(2)*mu2(3)*mu3(3)*mu4(1))*1i)/30;
  out(3) = (R*kres(2)*(mu1(1)*mu2(1)*mu3(1)*mu4(2) + mu1(1)*mu2(1)*mu3(2)*mu4(1) + mu1(1)*mu2(2)*mu3(1)*mu4(1) - 3*mu1(2)*mu2(1)*mu3(1)*mu4(1) + 3*mu1(1)*mu2(2)*mu3(2)*mu4(2) - mu1(2)*mu2(1)*mu3(2)*mu4(2) - mu1(2)*mu2(2)*mu3(1)*mu4(2) - mu1(2)*mu2(2)*mu3(2)*mu4(1) + mu1(1)*mu2(2)*mu3(3)*mu4(3) + mu1(1)*mu2(3)*mu3(2)*mu4(3) + mu1(1)*mu2(3)*mu3(3)*mu4(2) - mu1(2)*mu2(1)*mu3(3)*mu4(3) - mu1(2)*mu2(3)*mu3(1)*mu4(3) - mu1(2)*mu2(3)*mu3(3)*mu4(1))*1i)/30;
    end
elseif order(2) == 2 %take to next order up
    if order(1)==1     
    mu1 = mu{1}; mu2 = mu{2};
out(1) = (mu1(1)*mu2(1))/3 + (mu1(2)*mu2(2))/3 + (mu1(3)*mu2(3))/3 - (R^2*(kres(1)^2*mu1(1)*mu2(1) + 2*kres(2)^2*mu1(1)*mu2(1) + kres(1)^2*mu1(2)*mu2(2) + 2*kres(3)^2*mu1(1)*mu2(1) + 2*kres(2)^2*mu1(2)*mu2(2) + 3*kres(1)^2*mu1(3)*mu2(3) + 2*kres(3)^2*mu1(2)*mu2(2) + kres(2)^2*mu1(3)*mu2(3) + kres(3)^2*mu1(3)*mu2(3)))/30;
out(2) =  (kres(1)*kres(2)*(mu1(1)*mu2(1) + mu1(2)*mu2(2) - 2*mu1(3)*mu2(3))*R^2)/30 + kres(3)*(mu1(1)*mu2(2) - mu1(2)*mu2(1))*R*(-1i/6);
out(3) =   (kres(1)*kres(3)*(mu1(1)*mu2(1) + mu1(2)*mu2(2) - 2*mu1(3)*mu2(3))*R^2)/30 + kres(2)*(mu1(1)*mu2(2) - mu1(2)*mu2(1))*R*(1i/6);
    elseif order(1) == 3
        mu1 = mu{1}; mu2 = mu{2}; mu3 = mu{3}; mu4 = mu{4};
out(1) =   (mu1(1)*mu2(1)*mu3(1)*mu4(1))/5 - (R^2*(3*kres(1)^2*mu1(1)*mu2(1)*mu3(1)*mu4(1) + 9*kres(2)^2*mu1(1)*mu2(1)*mu3(1)*mu4(1) + kres(1)^2*mu1(1)*mu2(1)*mu3(2)*mu4(2) + kres(1)^2*mu1(1)*mu2(2)*mu3(1)*mu4(2) + kres(1)^2*mu1(1)*mu2(2)*mu3(2)*mu4(1) + kres(1)^2*mu1(2)*mu2(1)*mu3(1)*mu4(2) + kres(1)^2*mu1(2)*mu2(1)*mu3(2)*mu4(1) + kres(1)^2*mu1(2)*mu2(2)*mu3(1)*mu4(1) + 9*kres(3)^2*mu1(1)*mu2(1)*mu3(1)*mu4(1) + 3*kres(2)^2*mu1(1)*mu2(1)*mu3(2)*mu4(2) + 3*kres(2)^2*mu1(1)*mu2(2)*mu3(1)*mu4(2) + 3*kres(2)^2*mu1(1)*mu2(2)*mu3(2)*mu4(1) + 3*kres(2)^2*mu1(2)*mu2(1)*mu3(1)*mu4(2) + 3*kres(2)^2*mu1(2)*mu2(1)*mu3(2)*mu4(1) + 3*kres(2)^2*mu1(2)*mu2(2)*mu3(1)*mu4(1) + 3*kres(1)^2*mu1(1)*mu2(1)*mu3(3)*mu4(3) + 3*kres(1)^2*mu1(1)*mu2(3)*mu3(1)*mu4(3) + 3*kres(1)^2*mu1(1)*mu2(3)*mu3(3)*mu4(1) + 3*kres(1)^2*mu1(2)*mu2(2)*mu3(2)*mu4(2) + 3*kres(1)^2*mu1(3)*mu2(1)*mu3(1)*mu4(3) + 3*kres(1)^2*mu1(3)*mu2(1)*mu3(3)*mu4(1) + 3*kres(1)^2*mu1(3)*mu2(3)*mu3(1)*mu4(1) + 3*kres(3)^2*mu1(1)*mu2(1)*mu3(2)*mu4(2) + 3*kres(3)^2*mu1(1)*mu2(2)*mu3(1)*mu4(2) + 3*kres(3)^2*mu1(1)*mu2(2)*mu3(2)*mu4(1) + 3*kres(3)^2*mu1(2)*mu2(1)*mu3(1)*mu4(2) + 3*kres(3)^2*mu1(2)*mu2(1)*mu3(2)*mu4(1) + 3*kres(3)^2*mu1(2)*mu2(2)*mu3(1)*mu4(1) + 2*kres(2)^2*mu1(1)*mu2(1)*mu3(3)*mu4(3) + 2*kres(2)^2*mu1(1)*mu2(3)*mu3(1)*mu4(3) + 2*kres(2)^2*mu1(1)*mu2(3)*mu3(3)*mu4(1) + 9*kres(2)^2*mu1(2)*mu2(2)*mu3(2)*mu4(2) + 2*kres(2)^2*mu1(3)*mu2(1)*mu3(1)*mu4(3) + 2*kres(2)^2*mu1(3)*mu2(1)*mu3(3)*mu4(1) + 2*kres(2)^2*mu1(3)*mu2(3)*mu3(1)*mu4(1) + 3*kres(1)^2*mu1(2)*mu2(2)*mu3(3)*mu4(3) + 3*kres(1)^2*mu1(2)*mu2(3)*mu3(2)*mu4(3) + 3*kres(1)^2*mu1(2)*mu2(3)*mu3(3)*mu4(2) + 3*kres(1)^2*mu1(3)*mu2(2)*mu3(2)*mu4(3) + 3*kres(1)^2*mu1(3)*mu2(2)*mu3(3)*mu4(2) + 3*kres(1)^2*mu1(3)*mu2(3)*mu3(2)*mu4(2) + 2*kres(3)^2*mu1(1)*mu2(1)*mu3(3)*mu4(3) + 2*kres(3)^2*mu1(1)*mu2(3)*mu3(1)*mu4(3) + 2*kres(3)^2*mu1(1)*mu2(3)*mu3(3)*mu4(1) + 9*kres(3)^2*mu1(2)*mu2(2)*mu3(2)*mu4(2) + 2*kres(3)^2*mu1(3)*mu2(1)*mu3(1)*mu4(3) + 2*kres(3)^2*mu1(3)*mu2(1)*mu3(3)*mu4(1) + 2*kres(3)^2*mu1(3)*mu2(3)*mu3(1)*mu4(1) + 2*kres(2)^2*mu1(2)*mu2(2)*mu3(3)*mu4(3) + 2*kres(2)^2*mu1(2)*mu2(3)*mu3(2)*mu4(3) + 2*kres(2)^2*mu1(2)*mu2(3)*mu3(3)*mu4(2) + 2*kres(2)^2*mu1(3)*mu2(2)*mu3(2)*mu4(3) + 2*kres(2)^2*mu1(3)*mu2(2)*mu3(3)*mu4(2) + 2*kres(2)^2*mu1(3)*mu2(3)*mu3(2)*mu4(2) + 15*kres(1)^2*mu1(3)*mu2(3)*mu3(3)*mu4(3) + 2*kres(3)^2*mu1(2)*mu2(2)*mu3(3)*mu4(3) + 2*kres(3)^2*mu1(2)*mu2(3)*mu3(2)*mu4(3) + 2*kres(3)^2*mu1(2)*mu2(3)*mu3(3)*mu4(2) + 2*kres(3)^2*mu1(3)*mu2(2)*mu3(2)*mu4(3) + 2*kres(3)^2*mu1(3)*mu2(2)*mu3(3)*mu4(2) + 2*kres(3)^2*mu1(3)*mu2(3)*mu3(2)*mu4(2) + 3*kres(2)^2*mu1(3)*mu2(3)*mu3(3)*mu4(3) + 3*kres(3)^2*mu1(3)*mu2(3)*mu3(3)*mu4(3)))/210 + (mu1(1)*mu2(1)*mu3(2)*mu4(2))/15 + (mu1(1)*mu2(2)*mu3(1)*mu4(2))/15 + (mu1(1)*mu2(2)*mu3(2)*mu4(1))/15 + (mu1(2)*mu2(1)*mu3(1)*mu4(2))/15 + (mu1(2)*mu2(1)*mu3(2)*mu4(1))/15 + (mu1(2)*mu2(2)*mu3(1)*mu4(1))/15 + (mu1(1)*mu2(1)*mu3(3)*mu4(3))/15 + (mu1(1)*mu2(3)*mu3(1)*mu4(3))/15 + (mu1(1)*mu2(3)*mu3(3)*mu4(1))/15 + (mu1(2)*mu2(2)*mu3(2)*mu4(2))/5 + (mu1(3)*mu2(1)*mu3(1)*mu4(3))/15 + (mu1(3)*mu2(1)*mu3(3)*mu4(1))/15 + (mu1(3)*mu2(3)*mu3(1)*mu4(1))/15 + (mu1(2)*mu2(2)*mu3(3)*mu4(3))/15 + (mu1(2)*mu2(3)*mu3(2)*mu4(3))/15 + (mu1(2)*mu2(3)*mu3(3)*mu4(2))/15 + (mu1(3)*mu2(2)*mu3(2)*mu4(3))/15 + (mu1(3)*mu2(2)*mu3(3)*mu4(2))/15 + (mu1(3)*mu2(3)*mu3(2)*mu4(2))/15 + (mu1(3)*mu2(3)*mu3(3)*mu4(3))/5;
out(2) =    (kres(1)*kres(2)*(3*mu1(1)*mu2(1)*mu3(1)*mu4(1) + mu1(1)*mu2(1)*mu3(2)*mu4(2) + mu1(1)*mu2(2)*mu3(1)*mu4(2) + mu1(1)*mu2(2)*mu3(2)*mu4(1) + mu1(2)*mu2(1)*mu3(1)*mu4(2) + mu1(2)*mu2(1)*mu3(2)*mu4(1) + mu1(2)*mu2(2)*mu3(1)*mu4(1) + 3*mu1(1)*mu2(1)*mu3(3)*mu4(3) + 3*mu1(1)*mu2(3)*mu3(1)*mu4(3) + 3*mu1(1)*mu2(3)*mu3(3)*mu4(1) + 3*mu1(2)*mu2(2)*mu3(2)*mu4(2) - 4*mu1(3)*mu2(1)*mu3(1)*mu4(3) - 4*mu1(3)*mu2(1)*mu3(3)*mu4(1) - 4*mu1(3)*mu2(3)*mu3(1)*mu4(1) + 3*mu1(2)*mu2(2)*mu3(3)*mu4(3) + 3*mu1(2)*mu2(3)*mu3(2)*mu4(3) + 3*mu1(2)*mu2(3)*mu3(3)*mu4(2) - 4*mu1(3)*mu2(2)*mu3(2)*mu4(3) - 4*mu1(3)*mu2(2)*mu3(3)*mu4(2) - 4*mu1(3)*mu2(3)*mu3(2)*mu4(2) - 6*mu1(3)*mu2(3)*mu3(3)*mu4(3))*R^2)/210 + kres(3)*(mu1(1)*mu2(1)*mu3(1)*mu4(2) + mu1(1)*mu2(1)*mu3(2)*mu4(1) + mu1(1)*mu2(2)*mu3(1)*mu4(1) - 3*mu1(2)*mu2(1)*mu3(1)*mu4(1) + 3*mu1(1)*mu2(2)*mu3(2)*mu4(2) - mu1(2)*mu2(1)*mu3(2)*mu4(2) - mu1(2)*mu2(2)*mu3(1)*mu4(2) - mu1(2)*mu2(2)*mu3(2)*mu4(1) + mu1(1)*mu2(2)*mu3(3)*mu4(3) + mu1(1)*mu2(3)*mu3(2)*mu4(3) + mu1(1)*mu2(3)*mu3(3)*mu4(2) - mu1(2)*mu2(1)*mu3(3)*mu4(3) - mu1(2)*mu2(3)*mu3(1)*mu4(3) - mu1(2)*mu2(3)*mu3(3)*mu4(1))*R*(-1i/30);
out(3) =      (kres(1)*kres(3)*(3*mu1(1)*mu2(1)*mu3(1)*mu4(1) + mu1(1)*mu2(1)*mu3(2)*mu4(2) + mu1(1)*mu2(2)*mu3(1)*mu4(2) + mu1(1)*mu2(2)*mu3(2)*mu4(1) + mu1(2)*mu2(1)*mu3(1)*mu4(2) + mu1(2)*mu2(1)*mu3(2)*mu4(1) + mu1(2)*mu2(2)*mu3(1)*mu4(1) + 3*mu1(1)*mu2(1)*mu3(3)*mu4(3) + 3*mu1(1)*mu2(3)*mu3(1)*mu4(3) + 3*mu1(1)*mu2(3)*mu3(3)*mu4(1) + 3*mu1(2)*mu2(2)*mu3(2)*mu4(2) - 4*mu1(3)*mu2(1)*mu3(1)*mu4(3) - 4*mu1(3)*mu2(1)*mu3(3)*mu4(1) - 4*mu1(3)*mu2(3)*mu3(1)*mu4(1) + 3*mu1(2)*mu2(2)*mu3(3)*mu4(3) + 3*mu1(2)*mu2(3)*mu3(2)*mu4(3) + 3*mu1(2)*mu2(3)*mu3(3)*mu4(2) - 4*mu1(3)*mu2(2)*mu3(2)*mu4(3) - 4*mu1(3)*mu2(2)*mu3(3)*mu4(2) - 4*mu1(3)*mu2(3)*mu3(2)*mu4(2) - 6*mu1(3)*mu2(3)*mu3(3)*mu4(3))*R^2)/210 + kres(2)*(mu1(1)*mu2(1)*mu3(1)*mu4(2) + mu1(1)*mu2(1)*mu3(2)*mu4(1) + mu1(1)*mu2(2)*mu3(1)*mu4(1) - 3*mu1(2)*mu2(1)*mu3(1)*mu4(1) + 3*mu1(1)*mu2(2)*mu3(2)*mu4(2) - mu1(2)*mu2(1)*mu3(2)*mu4(2) - mu1(2)*mu2(2)*mu3(1)*mu4(2) - mu1(2)*mu2(2)*mu3(2)*mu4(1) + mu1(1)*mu2(2)*mu3(3)*mu4(3) + mu1(1)*mu2(3)*mu3(2)*mu4(3) + mu1(1)*mu2(3)*mu3(3)*mu4(2) - mu1(2)*mu2(1)*mu3(3)*mu4(3) - mu1(2)*mu2(3)*mu3(1)*mu4(3) - mu1(2)*mu2(3)*mu3(3)*mu4(1))*R*(1i/30);
                                                                                                                                               
    end
end



