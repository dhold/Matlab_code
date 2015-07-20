function [g_t,g_mat_t]  = line_broad_fn_deriv(...
            B,gam_rng,om_0,lam_rng,gam_dru,lam_dru,t,tol,n_matsu)
if nargin == 6
    t = sym('t','real');
g_t = sym(0*t);
else
g_t = zeros(size(t));    
g_mat_t = g_t;

end
if nargin <8
   tol = 1e-8; %standard tolerence
end
% n=1:nmax; vn = 2*pi*n/B;
% mat_dec_term = repmat(exp(-2*pi*t/B),nmax,1);
% mat_dec_term = cumprod(mat_dec_term)-1;
% n=1:nmax; vn = 2*pi*n/B;
% mat_dec_term = zeros(nmax,length(t));
% for lp = 1:nmax
% mat_dec_term(lp,:) = expm1(-vn(lp)*t);
% end
eul_gam = 0.577215664901532860606512090082;
for k = 1:length(gam_rng)
    
    gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);

xi = sqrt(om0^2-gam^2/4);

g11 = -(1i*lambda.*(xi.*(gam.^2 + 4*xi.^2)+...
exp(-(gam*t/2)).*((-4*gam*xi^2*sin(t*xi) +xi*(gam.^2 - 4*xi.^2)*cos(t*xi))...
 -gam/2*(4*gam*xi*cos(t*xi) +(gam.^2 - 4*xi.^2)*sin(t*xi)))))/(xi*4*om0.^2);
        

g22 = (1/xi)*1i*lambda*om0.^2.*...
 (-(gam - 2*1i*xi).^(-2).*(2*(1i*xi-gam/2)*exp(t*(1i*xi-gam/2)) + (gam - 2i*xi))...
 .*cot(B*(gam - 2*1i*xi)/4)+... 
  (gam + 2*1i*xi).^(-2).*(-2*(gam/2 + 1i*xi)*exp(-t*(gam/2 + 1i*xi)) + (gam + 2i*xi))...
  .*cot(B*(gam + 2*1i*xi)/4));

g_t = g_t + g11+g22;

end
for k = 1:length(gam_dru)
    gam = gam_dru(k); lambda = lam_dru(k);
    
G11=-((1i.*lambda.*(-gam*exp(-gam.*t) + gam))/gam);

G22=lambda.*(-gam*exp(-gam*t) + gam)*cot(B*gam/2)/gam;

%Gmat_n = 4*B*gam*lambda *(1./(B^2*gam^2*(n.^0)-4*(n.^2)*pi^2))*mat_dec_term;
%lambda/(B*gam*pi).*(pi.*(2 - B*gam*cot((B*gam)/2)));

g_t = g_t + G11+G22;

end

non_exp_term_sum = zeros(size(t));
%analytic expressions for the terms without the exponential decay (whole
%function at long times)
for k = 1:length(gam_dru)
    %psi2 is the euler gamma function
gam = gam_dru(k); lambda = lam_dru(k);        
non_exp_term_sum =  non_exp_term_sum + (lambda/(pi))*(...
        + psi2(1+B*gam/2/pi)-psi2(1-B*gam/2/pi) ); 
end
for k = 1:length(gam_rng)   
gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);
xi = sqrt(om0^2-gam^2/4);

rooTs=[-B*gam - 2i*B*xi - 4*pi, B*gam - 2i*B*xi - 4*pi, ...
    -B*gam + 2i*B*xi- 4*pi, B*gam + 2i*B*xi - 4*pi]/(4*pi);

prefct = lambda*gam/2/pi/om0^2;
denom = B^2*om0^4+2*pi^2*(2*om0^2-gam^2)*(1+rooTs).^2;

non_exp_term_sum = non_exp_term_sum -...
    2*pi*prefct*B*om0^4*sum(psi2(-rooTs).*(1+rooTs)./denom);
end


%compute mastubara terms up until a tolerence
n_start = 1; change_q = inf; old_g_mat_t =  zeros(size(t));
n_step = 200; 
%split time as early and late times will converge at different rates
t_non_conv = true(size(t)); %these points are not yet converged
t_zero_lg = t ==0; 
t_non_conv(t_zero_lg) = false; %at zero g is zero
old_g_mat_t =  old_g_mat_t(t_non_conv);
%precompute terms decaying w/ matsubara frequencies


while change_q > tol
n_end = n_start + n_step-1;
n=n_start:n_end; vn = 2*pi*n/B;
tt = t(t_non_conv);
mat_dec_term = repmat(exp(-2*pi*tt/B),n_step,1);
mat_dec_term(1,:) = exp(-2*pi*n_start*tt/B);
mat_dec_term = cumprod(mat_dec_term);

for k = 1:length(gam_rng)   
gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);
gmat_n = 4*gam*lambda*om0^2*(1./(B.*(om0^4+(2*om0^2-gam^2)*vn.^2+vn.^4))*mat_dec_term);
g_mat_t(t_non_conv) = g_mat_t(t_non_conv)+ gmat_n;
end
for k = 1:length(gam_dru)
gam = gam_dru(k); lambda = lam_dru(k);    
Gmat_n  = -4*gam*lambda/B *(1./(vn.^2-gam^2))*mat_dec_term;
g_mat_t(t_non_conv) = g_mat_t(t_non_conv)+Gmat_n;
end %late time stuff will probably converge much faster
%g_mat_t(end)
if nargin == 9 %specific number to use given  
    n_start = n_start + n_step;
    if n_start >= n_matsu
        break
    end  
else
 n_start = n_start + n_step;
if n_start > 300000
   %warning('failed to converge with a sum over 300000 matsubara frequencies')
   %change_q
   break
end
change_qq =abs(old_g_mat_t-g_mat_t(t_non_conv))./(abs(g_mat_t(t_non_conv))+eps);

t_non_conv(change_qq<tol) = false; %these points have converged
%sum(t_non_conv)
change_q = max(change_qq); %biggest change

old_g_mat_t = g_mat_t(t_non_conv);
end
end

%Phi(z,1,v) = v^(-1) Hypergeometric2F1(1,v;1+v;z) for Abs(z)<1
% z = exp(-2*pi*t/B), v=1+/- B*gamma/2pi

%z = exp(-2*pi*t/B); 

% for k = 1:length(gam_dru) %in principle analytic but fucking useless
% gam = gam_dru(k); lambda = lam_dru(k);    
% v1 = 1-B*gam/2*pi; v2 = 1+B*gam/2*pi;
% %g_mat_t = g_mat_t +lambda*(z/pi).*(hyp2f1(1,v1,1+v1,z)/v1 - hyp2f1(1,v2,1+v2,z)/v2);
% 
% [tmp1, tmp11]=hypergeometric2F1ODE(1,v1,1+v1,[0,max(z)]);
% tmp11 = (interp1(tmp1,tmp11(:,2),z,'pchip'))/v1;
% [tmp2, tmp22]=hypergeometric2F1ODE(1,v2,1+v2,[0,max(z)]);
% tmp22 = (interp1(tmp2,tmp22(:,2),z,'pchip'))/v2;
% g_mat_t = g_mat_t +lambda*(z/pi).*(tmp11 - tmp22);
% end %late time stuff will probably converge much faster

g_mat_t = g_mat_t + non_exp_term_sum;
g_mat_t(t_zero_lg) = 0; %this term is exactly zero

if nargout == 1 %combine both into one
    g_t = g_t + g_mat_t;
end