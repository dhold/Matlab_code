function [g1,g2,gmatsu] = line_broadening_fn(lambda,om_0,gamma, beta)
if nargin == 0
syms lambda om_0 gamma real 
end
syms beta t1 t2 t v_n
if ~isempty(om_0)
xi = sqrt(om_0^2-gamma^2/4);
phi = gamma/2 + 1i*xi; phi_p = gamma/2 - 1i*xi;

C1 = -1/2/xi * exp(-gamma*t1/2)*sin(xi*t1);
C2 = 1/4/xi * (coth(1i*phi_p*beta/2)*exp(-phi_p*t1) - ...
            coth(1i*phi*beta/2)*exp(-phi8t1)); 
C_mat = - 2*gamma/beta *v_n*exp(-v_n*t1)/((om_0^2+v_n^2)^2-gamma^2*v_n^2);

g1 = int(C1,t1,0,t2);
g1 = 2i*om_0^2*lambda*int(g1,t2,0,t); %imaginary component
g2 = int(C2,t1,0,t2);
g2 = 2*om_0^2*lambda*int(g2,t2,0,t); %real
gmatsu = int(C_mat,t1,0,t2);
gmatsu = 2*om_0^2*lambda*int(gmatsu,t2,0,t); %summed to infinity
%v_n = 2*pi n /hbar/beta
else
C1 = -1i*gamma*lambda*exp(-gamma*t1);    
C2 = lambda*gamma*cot(beta*gamma/2)*exp(-gamma*t1);

C_mat = 4*gamma*lambda/beta*v_n*exp(-v_n*t1)/(v_n^2-gamma^2);
g1 = int(C1,t1,0,t2);
g1 = int(g1,t2,0,t); %imaginary component
g2 = int(C2,t1,0,t2);
g2 = int(g2,t2,0,t); %real
gmatsu = int(C_mat,t1,0,t2);
gmatsu = int(gmatsu,t2,0,t); %summed to infinity
end