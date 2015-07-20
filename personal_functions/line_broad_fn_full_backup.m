function [g_t,g_mat_t]  = line_broad_fn_full_backup(B,gam_rng,om_0,lam_rng,gam_dru,lam_dru,t)

if nargin==6
    t = sym('t','real'); g_t = t*0; g_mat_t = t*0;
else
g_t = zeros(size(t));
if nargout ==2
    g_mat_t = zeros(size(t));
end
end

%vn = 2*pi*n/B; vn*t is decay term
% n_max = round(10*B/2/pi./t); %max number to include in the sum for good convergence
% %if any of these are just too large use the zeno approx 
% %~t^2*(1i*C'(0)+C''(0))
% zeno_term = n_max>1e5 | t==0; 
% t_zeno = t(zeno_term);
% t_inter = t(n_max > 0);
% t = t(n_max==0);

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

%assumes matsubara frequencies are large enough that the exponentials decay
%rapidly and can be removed from the sum.  This causes some problems at
%very low time when this function should be zero

rooTs=[-B*gam - 2i*B*xi - 4*pi, B*gam - 2i*B*xi - 4*pi, ...
    -B*gam + 2i*B*xi- 4*pi, B*gam + 2i*B*xi - 4*pi]/(4*pi);

%gmat=-(1/(om0.^2*pi))*2*gam*lambda*( psi2(1) + 1/4*B*om0^4*...
%    sum(-B*psi2(-rooTs) + 2*pi*t*psi2(-rooTs)*(1+rooTs)...
%    ./(B.^2*om0^4 - 2*gam.^2*pi.^2 + 4*om0.^2*pi.^2 - 4*gam.^2*pi.^2*rooTs... 
% + 8*om0.^2*pi.^2*rooTs - 2*gam.^2*pi.^2*rooTs.^2 + 4*om0.^2*pi.^2*rooTs.^2)));

% linear order terms in t
gmat= t*B^3*gam*lambda*om0^2*sum(psi2(-rooTs)./(-B^2*gam^2+2*B^2*om0^2 + 8*pi^2 -...
    (B^2*gam^2+ 2*B^2*om0^2+24*pi^2)*rooTs + 24*pi^2*rooTs.^2 + 8*pi^2*rooTs.^2))/2/pi^2;

gmat_cnst = -pi*gam*lambda/(2)*(4*eul_gam/om0^2 + B^2*om0^2 *sum(psi2(-rooTs)./...
               (B^2*om0^4/pi^2 + (4*om0^2- 2*gam^2) +(8*om0^2-4*gam^2)*rooTs + ...
               (4*om0^2- 2*gam^2)*rooTs.^2)));

% n=1:nmax; vn = 2*pi*n/B;
% Gmat_n  = -2*B^4*gam*lambda*om0^2 *(expm1(vn*t))...
%             /(B^4*om0^2*n*pi-4*B^2*n^3*pi^3*(gam^2-2*om0^2)+16*n^5*pi^5);


%gmat=-(1/(om0.^2*pi))*2*gam*lambda*(1/4*B*om0^4*2*pi*t*...
%    sum(psi2(-rooTs)*(1+rooTs)...
%    ./(B.^2*om0^4 - 2*gam.^2*pi.^2 + 4*om0.^2*pi.^2 - 4*gam.^2*pi.^2*rooTs... 
% + 8*om0.^2*pi.^2*rooTs - 2*gam.^2*pi.^2*rooTs.^2 + 4*om0.^2*pi.^2*rooTs.^2)));

% root(1)=(-B*gam - sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);
% root(2)=(B*gam - sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);
% root(3)=(-B*gam + sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);
% root(4)=(B*gam + sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);

g_t = g_t + g11+g22;
g_mat_t = g_mat_t+gmat+gmat_cnst*ones(size(t));
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
Gmat = lambda/(B*gam*pi).*(pi.*t*(2 - B*gam*cot((B*gam)/2)));% + ...
%B*(HarmonicNumber(-((B*gam)/(2*pi))) + HarmonicNumber((B*gam)/(2*pi))));
%Gmat_n  = -2*B^2*gam*lambda *(expm1(vn*t)+vn*t)/(B^2*gam^2*n*pi-4*n^3*pi^3);
Gmat_cnst = lambda/pi/gam*(2*eul_gam+psi2(1+B*gam/2/pi)+psi2(1-B*gam/2/pi));
g_t = g_t + G11+G22;
g_mat_t = g_mat_t+Gmat+Gmat_cnst*ones(size(t));



% Gmat_zeno = (-((0.577215664901532860606512090082*gam*lambda)/pi)+...
%     (lambda*pi)/(3*B^2*gam)+...
% (2*lambda*pi*zeta(-1,1-(B*gam)/(2*pi)))/(B^2*gam)+...
% (2*lambda*pi*zeta(-1,1+(B*gam)/(2*pi)))/(B^2*gam)+...
% (2*lambda*zeta(0,1-(B*gam)/(2*pi)))/B-...
% (2*lambda*zeta(0,1+(B*gam)/(2*pi)))/B-(gam*lambda*Log((2*pi)/B))/pi-...
% (gam*lambda*Log(t))/pi-(gam*lambda*psi2(1-(B*gam)/(2*pi)))/(2*pi)-...
% (gam*lambda*psi2(1+(B*gam)/(2*pi)))/(2*pi))*t^2;


end