g11 = -(4*1i*exp(-((gam*t)/2))+lambda*om0.^2*(exp((gam*t)/2)...
+xi*(-4*gam + gam.^2*t + 4*t*xi.^2)+4*gam*xi*cos(t*xi) +...
(gam.^2 - 4*xi.^2)*sin(t*xi)))/(xi*(gam.^2 + 4*xi.^2).^2);

g22 = (1/xi)*1i*lambda*om0.^2*(-1/(gam - 2*1i*xi).^2*(2*exp(1i*t*xi) +... 
(-2 + gam*t - 2*1i*t*xi))*cot(1/4*B*(gam - 2*1i*xi)) +... 
1/(gam + 2*1i*xi).^2*exp(-(1/2)*t*(gam + 2*1i*xi))*...
(2 +exp(1/2*t*(gam + 2*1i*xi))*(-2 + gam*t + 2*1i*t*xi))*cot(1/4*B*(gam + 2*1i*xi)));

%assumes matsubara frequencies are large enough that the exponentials decay
%rapidly and can be removed from the sum.  This causes some problems at
%very low time when this function should be zero

rooTs=[-B*gam - sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi, ...
        B*gam - sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi, ...
    -B*gam + sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi, ...
    B*gam + sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi]/(4*pi);

%gmat=-(1/(om0.^2*pi))*2*gam*lambda*( psi(1) + 1/4*B*om0^4*...
%    sum(-B*psi(-rooTs) + 2*pi*t*psi(-rooTs)*(1+rooTs)...
%    ./(B.^2*om0^4 - 2*gam.^2*pi.^2 + 4*om0.^2*pi.^2 - 4*gam.^2*pi.^2*rooTs... 
% + 8*om0.^2*pi.^2*rooTs - 2*gam.^2*pi.^2*rooTs.^2 + 4*om0.^2*pi.^2*rooTs.^2)));
gmat=-(1/(om0.^2*pi))*2*gam*lambda*(1/4*B*om0^4*2*pi*t*...
    sum(psi(-rooTs)*(1+rooTs)...
    ./(B.^2*om0^4 - 2*gam.^2*pi.^2 + 4*om0.^2*pi.^2 - 4*gam.^2*pi.^2*rooTs... 
 + 8*om0.^2*pi.^2*rooTs - 2*gam.^2*pi.^2*rooTs.^2 + 4*om0.^2*pi.^2*rooTs.^2)));

% root(1)=(-B*gam - sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);
% root(2)=(B*gam - sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);
% root(3)=(-B*gam + sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);
% root(4)=(B*gam + sqrt(B.^2*gam.^2 - 4*B.^2*om0.^2) - 4*pi)/(4*pi);

G11=-((1i*lambda*(-1 + exp(-gam*t)+gam*t))/gam);

G22=(exp(-gam*t)*lambda*(1 + exp(gam*t)*(-1 + gam*t))*cot((B*gam)/2))/gam;

Gmat = 1/(B*gam*pi)*lambda*(pi*t*(2 - B*gam*cot((B*gam)/2)));% + ...
%B*(HarmonicNumber(-((B*gam)/(2*pi))) + HarmonicNumber((B*gam)/(2*pi))));