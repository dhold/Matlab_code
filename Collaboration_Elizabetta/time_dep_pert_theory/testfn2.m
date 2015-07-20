
t=t_int_rng;
z = exp(-2*pi*t/B); 
gam = gam_dru{1}(1); lambda = lam_dru{1}(1);   

[tmp1, tmp11]=hypergeometric2F1ODE(1,v1,1+v1,[0,max(z)]);
tmp11 = (interp1(tmp1,tmp11(:,1),z,'pchip'))/v1;
[tmp2, tmp22]=hypergeometric2F1ODE(1,v2,1+v2,[0,max(z)]);
tmp22 = (interp1(tmp2,tmp22(:,1),z,'pchip'))/v2;
tmp_b = lambda*(z/pi/gam).*(tmp11 +tmp22);
tmp_b = tmp_b + lambda*(2/pi/gam)*log(1-z);
v1 = 1-B*gam/2*pi; v2 = 1+B*gam/2*pi;

tmp = lambda*(z/pi/gam).*(hyp2f1(1,v1,1+v1,z)/v1 + hyp2f1(1,v2,1+v2,z)/v2);
tmp = tmp + lambda*(2/pi/gam)*log(1-z);

mat_dec_term = repmat(exp(-2*pi*t/B),500,1);
mat_dec_term(1,:) = z;
vn = 2*pi*(1:500)/B;

mat_dec_term = cumprod(mat_dec_term);  
Gmat_n  = 4*gam*lambda/B *((1./(vn.*(vn.^2-gam^2)))*mat_dec_term);