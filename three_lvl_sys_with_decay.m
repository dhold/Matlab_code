function [rho_out,t_out] = three_lvl_sys_with_decay

syms E_1 E_2 rb_1 rb_2 g1 g2 t
%E_1 = 0; E_2 = 10;
%rb_1 = 0; rb_2 = 1; 

%three_lvl_LB = zeros(3*3);
Ham = [0,rb_1,0;rb_1,E_1,rb_2;0,rb_2,E_2];
%Ham = [0,0,0;0,0,rb_2;0,rb_2,E_2]; %simple form
Louiv = 1i*(-kron(eye(3),Ham) + kron(Ham,eye(3)));

%tmp = -1/2*[0,g1,g2;0,0,g1+g2;0,0,0];
%tmp = tmp + tmp'; %off diag

three_lvl_LB =diag(-1/2*[0,g1,g2,g1,0,g1+g2,g2,g1+g2,0]);
three_lvl_LB = three_lvl_LB  + diag([0,0,0,0,-g1,0,0,0,-g2]);
three_lvl_LB(1,5) = g1;
three_lvl_LB(5,9) = g2;
%this is just the coherence decay

%t_ev_op = expm(Louiv*t+three_lvl_LB*t);
time_prop({[E_1,E_2,rb_1,rb_2,g1,g2,t],[0,10,0,50,1,100,t]},Louiv+three_lvl_LB);

rho0 = [0,0,0,0,1,0,0,0,0];
options = odeset('RelTol',1e-5,'AbsTol',ones(9,1)*[1e-6]);
[t_out,rho_out] = ode45(@time_prop,[0,4],rho0,options);
end
function  drho = time_prop(tt,rho)
persistent propop
if iscell(tt)
propop = double(subs(rho,tt{1},tt{2}));
                            drho = pi; return
else
drho = propop*rho;
end
end