function [g_t,g_mat_t]  = line_broad_fn_sec_der(B,gam_rng,om_0,lam_rng,...
                            gam_dru,lam_dru,t,tol,n_matsu)
%call with two outputs to see components seperately, otherwise Mastsu freq
%added to g_t
if nargin==6
    t = sym('t','real'); g_t = t*0; g_mat_t = t*0;
else
g_t = zeros(size(t)); g_mat_t=g_t;
end

%precompute terms decaying w/ matsubara frequencies
if nargin <8
   tol = 1e-8; %standard tolerence
end
% n=1:nmax;  vn = 2*pi*n/B;
% mat_dec_term = repmat(exp(-2*pi*t/B),nmax,1);
% mat_dec_term = cumprod(mat_dec_term);

for k = 1:length(gam_rng)
    
    gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);

xi = sqrt(om0^2-gam^2/4);

g11 = -1i*lambda/(xi*4*om0.^2).*(...
exp(-(gam*t/2)).*((-4*gam*xi^3*cos(t*xi) -xi^2*(gam.^2 - 4*xi.^2)*sin(t*xi))...
 -gam*(-4*gam*xi^2*sin(t*xi) +xi*(gam.^2 - 4*xi.^2)*cos(t*xi))...
 -gam^2/4*(4*gam*xi*cos(t*xi) +(gam.^2 - 4*xi.^2)*sin(t*xi))));
        
g22 = (1/xi)*1i*lambda*om0.^2.*...
 (-(gam - 2*1i*xi).^(-2).*(2*(1i*xi-gam/2)^2*exp(t*(1i*xi-gam/2)))...
 .*cot(B*(gam - 2*1i*xi)/4)+... 
  (gam + 2*1i*xi).^(-2).*(2*(gam/2 + 1i*xi)^2*exp(-t*(gam/2 + 1i*xi)))...
  .*cot(B*(gam + 2*1i*xi)/4));


%assumes matsubara frequencies are large enough that the exponentials decay
%rapidly and can be removed from the sum.  This causes some problems at
%very low time when this function should be zero

%rooTs=[-B*gam - 2i*B*xi - 4*pi, B*gam - 2i*B*xi - 4*pi, ...
%    -B*gam + 2i*B*xi- 4*pi, B*gam + 2i*B*xi - 4*pi]/(4*pi);

%gmat=-(1/(om0.^2*pi))*2*gam*lambda*( psi2(1) + 1/4*B*om0^4*...
%    sum(-B*psi2(-rooTs) + 2*pi*t*psi2(-rooTs)*(1+rooTs)...
%    ./(B.^2*om0^4 - 2*gam.^2*pi.^2 + 4*om0.^2*pi.^2 - 4*gam.^2*pi.^2*rooTs... 
% + 8*om0.^2*pi.^2*rooTs - 2*gam.^2*pi.^2*rooTs.^2 + 4*om0.^2*pi.^2*rooTs.^2)));



%gmat_n= -(8*om0.^2*gam*lambda*pi*n./...
%        (B^4*om0^4 - 4*(B^2*gam^2-2*om0^2)*pi^2+16*n.^4*pi^4))*mat_dec_term; 

%gmat=-(1/(om0.^2*pi))*2*gam*lambda*(1/4*B*om0^4*2*pi*t*...
%    sum(psi2(-rooTs)*(1+rooTs)...
%    ./(B.^2*om0^4 - 2*gam.^2*pi.^2 + 4*om0.^2*pi.^2 - 4*gam.^2*pi.^2*rooTs... 
% + 8*om0.^2*pi.^2*rooTs - 2*gam.^2*pi.^2*rooTs.^2 + 4*om0.^2*pi.^2*rooTs.^2)));

% root(1)=(-B*gam - sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);
% root(2)=(B*gam - sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);
% root(3)=(-B*gam + sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);
% root(4)=(B*gam + sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);

g_t = g_t + g11+g22;

end
for k = 1:length(gam_dru)
    gam = gam_dru(k); lambda = lam_dru(k);
    
G11=-((1i.*lambda.*(gam^2*exp(-gam.*t)))/gam);

G22=lambda.*(gam^2*exp(-gam*t))*cot(B*gam/2)/gam;
%Gmat_n= -(8*gam*lambda*pi*n./(B^2*gam^2-4*n.^2*pi^2))*mat_dec_term; 

g_t = g_t + G11+G22;

end

%compute mastubara terms up until a tolerence

n_start = 1; change_q = inf; old_g_mat_t = 0*t;
n_step = 500;
%precompute terms decaying w/ matsubara frequencies
while change_q > tol
n_end = n_start + n_step-1;
n=n_start:n_end; vn = 2*pi*n/B;

mat_dec_term = repmat(exp(-2*pi*t*n_start/B),n_step,1);
mat_dec_term = cumprod(mat_dec_term);
for k = 1:length(gam_rng)   
gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);
gmat_n = -4*gam*lambda*om0^2.*(vn./(B.*(om0^4+(2*om0^2-gam^2)*vn.^2+vn.^4))*mat_dec_term);
g_mat_t = g_mat_t+ gmat_n;
end
for k = 1:length(gam_dru)
gam = gam_dru(k); lambda = lam_dru(k);    
Gmat_n  = 4*gam*lambda/B *((vn./(vn.^2-gam^2))*mat_dec_term);
g_mat_t = g_mat_t+Gmat_n;
end %late time stuff will probably converge much faster
if nargin == 9 %specific number to use given  
    n_start = n_start + n_step;
    if n_start >= n_matsu
        break
    end  
else
norm_save = trapz(t,abs(old_g_mat_t).^2)+eps;
change_q = trapz(t,abs(old_g_mat_t-g_mat_t).^2)/norm_save;
n_start = n_start + n_step;
if n_start > 100000
   warning('failed to converge with a sum over 100000 matsubara frequencies')
   change_q
   break
end
old_g_mat_t = g_mat_t;
end
end
    
if nargout == 1
    g_t = g_t + g_mat_t;
end