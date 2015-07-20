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
n_step = 500; 
%precompute terms decaying w/ matsubara frequencies
while change_q > tol
n_end = n_start + n_step-1;
n=n_start:n_end; vn = 2*pi*n/B;

mat_dec_term = repmat(exp(-2*pi*t*n_start/B),n_step,1);
mat_dec_term = cumprod(mat_dec_term)-1 + kron(t,vn');
for k = 1:length(gam_rng)   
gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);
gmat_n = -4*gam*lambda*om0^2*(1./(B*vn.*(om0^4+(2*om0^2-gam^2)*vn.^2+vn.^4))*mat_dec_term);
g_mat_t = g_mat_t+ gmat_n;
end
for k = 1:length(gam_dru)
gam = gam_dru(k); lambda = lam_dru(k);    
Gmat_n  = 4*gam*lambda/B *((1./(vn.*(vn.^2-gam^2)))*mat_dec_term);
g_mat_t = g_mat_t+Gmat_n;
end %late time stuff will probably converge much faster
%(g_mat_t(end-1)-g_mat_t(end))/(t(end-1)-t(end))
norm_save = trapz(t,abs(old_g_mat_t).^2);
change_q = trapz(t,abs(old_g_mat_t-g_mat_t).^2)/norm_save;
n_start = n_start + n_step;
if n_start > 100000
   warning('failed to converge with a sum over 100000 matsubara frequencies')
   change_q
   break
end
old_g_mat_t = g_mat_t;
end

if nargout == 1 %add g_mat_t to g_t
    g_t = g_t + g_mat_t;
end