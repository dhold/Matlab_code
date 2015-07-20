syms R1 R2 R3 real
RR = [R1;R2;R3];

uu = [R2;-R1;0]/sqrt(R1^2+R2^2);

ucross = [0,-uu(3),uu(2);uu(3),0,-uu(1);-uu(2),uu(1),0];
ukron = kron(uu.',uu);
%cos(theta) = (R3/sqrt(R1^2+R2^2+R3^2));
%sin(theta) = (1 - R3^2/(R1^2 + R2^2 + R3^2))^(1/2);
Trot = (R3/(R1^2+R2^2+R3^2)^(1/2))*eye(3) + ...
         (1 - R3^2/(R1^2 + R2^2 + R3^2))^(1/2)*ucross + ...
         (1-R3/(R1^2+R2^2+R3^2)^(1/2))*ukron;
     
simplify(Trot*RR,'IgnoreAnalyticConstraints',true)     

syms u1 u2 u3 real
u = [u1,u2,u3];

uprime = u*(Trot^(-1));

Trot1 = subs(Trot,{R1,R2,R3},{uprime(1),uprime(2),uprime(3)});

%%
% 
% syms a b y real 
% %a,b,y are the euler angles (normally alpha beta gamma or phi theta xi)
% r1 = [cos(y),-sin(y),0;sin(y),cos(y),0;0,0,1];
% r2 = [cos(a),0,sin(a);0,1,0;-sin(a),0,cos(a)]; %a runs from 0 to pi
% r3 = [cos(b),-sin(b),0;sin(b),cos(b),0;0,0,1];
% TT = r3*r2*r1;
% 
% Rrand = rand(3,1)-1/2;
% test = solve(TT*Rrand==norm(Rrand)*[0;0;1],a,b,y)
% 
% syms R1 R2 R3 real
% RR = [R1;R2;R3];
% out = solve(TT*RR==norm(RR)*[0;0;1],a,b,y)