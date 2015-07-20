function [g_cnst,g_mat,G_mat] = line_broad_fn_markov(B,gam_rng,om_0,lam_rng,gam_dru,lam_dru)
%takes only terms of g(t) linear in t, this gives the constant g with g(t)
%= g * t
% B = beta value
g_cnst = 0; g_mat = zeros(length(lam_rng),1) ; G_mat =  zeros(length(lam_dru),1); 
%g_offset = zeros(size(k)); %offset for mastubara
for k = 1:length(gam_rng)
    
    gam = gam_rng(k); om0 = om_0(k);  lambda = lam_rng(k);

xi = sqrt(om0^2-gam^2/4);

g11 = -1i;%*om0.^2./(gam.^2 + 4*xi.^2); %this bit

g22 = (1/xi)*1i*om0.^2.* ...
     ( cot(B*(gam + 2*1i*xi)/4)./(gam + 2i*xi) - ... 
        cot(B*(gam - 2*1i*xi)/4)./(gam - 2i*xi));

rooTs=[-B*gam - 2i*B*xi - 4*pi, ...
        B*gam - 2i*B*xi - 4*pi, ...
        -B*gam + 2i*B*xi - 4*pi, ...
        B*gam + 2i*B*xi - 4*pi]/(4*pi);

%just take the linear order terms in t
gmat= B^3*gam*om0^2*sum(psi2(-rooTs)./(-B^2*gam^2+2*B^2*om0^2 + 8*pi^2 -...
    (B^2*gam^2+ 2*B^2*om0^2+24*pi^2)*rooTs + 24*pi^2*rooTs.^2 + 8*pi^2*rooTs.^2))/2/pi^2;

g_cnst = g_cnst + lambda*(g11+g22+gmat);

g_mat(k) = lambda*gmat;
end
for k = 1:length(gam_dru)
    gam = gam_dru(k); lambda = lam_dru(k);
    
G11= -1i;
G22= cot(B*gam/2);
Gmat = (2 - B*gam*cot(B*gam/2))./(B*gam);

g_cnst = g_cnst + lambda.*(G11+G22*Gmat);
G_mat(k) = Gmat;
end