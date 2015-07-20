function [g_t,g_mat_t,n_end,non_exp_term_sum]  = ...
    line_broad_fn_full(B,gam_rng,om_0,lam_rng,gam_dru,lam_dru,t,tol)

%this converges matsubara terms to a tolerancec "tol", outputing three args
%will give you the number of steps required to achieve this
if nargin==6
    t = sym('t','real'); g_t = t*0; g_mat_t = t*0;
else
g_t = zeros(size(t)); g_mat_t=g_t;
end
if nargin <8
   tol = 1e-8; %standard tolerence
end

eul_gam = 0.577215664901532860606512090082;

for k = 1:length(gam_rng)
    
    gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);

xi = sqrt(om0^2-gam^2/4);

g11 = -(1i*lambda.*(xi.*(-4*gam + gam.^2*t + 4*t*xi.^2)+...
exp(-(gam*t/2)).*(4*gam*xi*cos(t*xi) +(gam.^2 - 4*xi.^2)*sin(t*xi))))...
/(xi*4*om0.^2);

if ~isa(t,'sym') %apparently just to piss me off expm1 doesn't work w/sym
g22 = (1/xi)*1i*lambda*om0.^2.*...
 (-(gam - 2*1i*xi).^(-2).*(2*expm1(1i*t*xi-gam*t/2) + t*(gam - 2i*xi)).*cot(B*(gam - 2*1i*xi)/4)+... 
  (gam + 2*1i*xi).^(-2).*(2*expm1(-t*(gam/2 + 1i*xi)) + t*(gam + 2i*xi)).*cot(B*(gam + 2*1i*xi)/4));
else
g22 = (1/xi)*1i*lambda*om0.^2.*...
 (-(gam - 2*1i*xi).^(-2).*(2*(exp(t*(1i*xi-gam/2))-1) + t*(gam - 2i*xi)).*cot(B*(gam - 2*1i*xi)/4)+... 
  (gam + 2*1i*xi).^(-2).*(2*(exp(-t*(gam/2 + 1i*xi))-1) + t*(gam + 2i*xi)).*cot(B*(gam + 2*1i*xi)/4));    
end

g_t = g_t + g11+g22;

end
for k = 1:length(gam_dru)
    gam = gam_dru(k); lambda = lam_dru(k);
if ~isa(t,'sym')    
G11=-((1i.*lambda.*(expm1(-gam.*t) + gam.*t))/gam);
G22=lambda.*(expm1(-gam*t) + gam.*t)*cot(B*gam/2)/gam;
else
G11=-((1i.*lambda.*(exp(-gam.*t)-1 + gam.*t))/gam);
G22=lambda.*(exp(-gam*t)-1 + gam.*t)*cot(B*gam/2)/gam;   
end

g_t = g_t + G11+G22;

end
%compute mastubara terms up until a tolerence

non_exp_term_sum = zeros(size(t));
%analytic expressions for the terms without the exponential decay (whole
%function at long times)
for k = 1:length(gam_dru)
    %psi2 is the euler gamma function
gam = gam_dru(k); lambda = lam_dru(k);        
non_exp_term_sum =  non_exp_term_sum + (lambda/(pi*gam))*(2*eul_gam+...
        (1-gam*t)*psi2(1-B*gam/2/pi) + (1+gam*t)*psi2(1+B*gam/2/pi)); 
end
for k = 1:length(gam_rng)   
gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);
xi = sqrt(om0^2-gam^2/4);

rooTs=[-B*gam - 2i*B*xi - 4*pi, B*gam - 2i*B*xi - 4*pi, ...
    -B*gam + 2i*B*xi- 4*pi, B*gam + 2i*B*xi - 4*pi]/(4*pi);

prefct = lambda*gam/2/pi/om0^2;
denom = B^2*om0^4+2*pi^2*(2*om0^2-gam^2)*(1+rooTs).^2;

non_exp_term_sum = non_exp_term_sum +...
    prefct*(4*eul_gam + B^2*om0^4*sum(psi2(-rooTs)./denom));

non_exp_term_sum = non_exp_term_sum -...
    2*pi*t*prefct*B*om0^4*sum(psi2(-rooTs).*(1+rooTs)./denom);
end

n_start = 1; change_q = inf; old_g_mat_t =  zeros(size(t));
n_step = 200; %100 matsubara steps
%split time as early and late times will converge at different rates
t_non_conv = true(size(t)); 
t_zero_lg = t ==0; 
t_non_conv(t_zero_lg) = false; %at zero g is zero
old_g_mat_t =  old_g_mat_t(t_non_conv);

%precompute terms decaying w/ matsubara frequencies
while change_q > tol
    tt = t(t_non_conv);
n_end = n_start + n_step-1;
n=n_start:n_end; vn = 2*pi*n/B;

mat_dec_term = repmat(exp(-2*pi*tt/B),n_step,1);
mat_dec_term(1,:) = exp(-2*pi*n_start*tt/B);

mat_dec_term = cumprod(mat_dec_term);
for k = 1:length(gam_rng)   
gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);
gmat_n = -4*gam*lambda*om0^2*(1./(B*vn.*(om0^4+(2*om0^2-gam^2)*vn.^2+vn.^4))*mat_dec_term);
g_mat_t(t_non_conv) = g_mat_t(t_non_conv)+ gmat_n;
end
for k = 1:length(gam_dru)
gam = gam_dru(k); lambda = lam_dru(k);    
Gmat_n  = 4*gam*lambda/B *((1./(vn.*(vn.^2-gam^2)))*mat_dec_term);
g_mat_t(t_non_conv) = g_mat_t(t_non_conv) + Gmat_n;
end %late time stuff will probably converge much faster
%(g_mat_t(end-1)-g_mat_t(end))/(t(end-1)-t(end))
% norm_save = trapz(t,abs(old_g_mat_t).^2);
% change_q = trapz(t,abs(old_g_mat_t-g_mat_t).^2)/norm_save;
 n_start = n_start + n_step;
if n_start > 300000
   warning('failed to converge with a sum over 300000 matsubara frequencies')
   change_q
   break
end
change_qq =abs(old_g_mat_t-g_mat_t(t_non_conv))./(abs(g_mat_t(t_non_conv))+eps);

t_non_conv(change_qq<tol) = false; %these points have converged
%sum(t_non_conv)
change_q = max(change_qq); %biggest change

old_g_mat_t = g_mat_t(t_non_conv);
end

%Phi(z,1,v) = v^(-1) Hypergeometric2F1(1,v;1+v;z) for Abs(z)<1
% z = exp(-2*pi*t/B), v=1+/- B*gamma/2pi

% z = exp(-2*pi*t/B); 
% for k = 1:length(gam_dru) %can be performed ok analytically
% gam = gam_dru(k); lambda = lam_dru(k);    
% v1 = 1-B*gam/2*pi; v2 = 1+B*gam/2*pi;
% [tmp1, tmp11]=hypergeometric2F1ODE(1,v1,1+v1,[0,max(z)]);
% tmp11 = (interp1(tmp1,tmp11(:,1),z,'pchip'))/v1;
% [tmp2, tmp22]=hypergeometric2F1ODE(1,v2,1+v2,[0,max(z)]);
% tmp22 = (interp1(tmp2,tmp22(:,1),z,'pchip'))/v2;
% tmp = lambda*(z/pi/gam).*(tmp11 +tmp22);
% %tmp = lambda*(z/pi/gam).*(hyp2f1(1,v1,1+v1,z)/v1 + hyp2f1(1,v2,1+v2,z)/v2);
% tmp = tmp + lambda*(2/pi/gam)*log(1-z);
% g_mat_t = g_mat_t +tmp;
% end 

g_mat_t = g_mat_t + non_exp_term_sum;
g_mat_t(t_zero_lg) = 0; %this term is exactly zero

if nargout == 1 %add g_mat_t to g_t
    g_t = g_t + g_mat_t;
end