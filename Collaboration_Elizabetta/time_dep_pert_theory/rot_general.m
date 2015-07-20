mu = sym('mu',[4,3]); mu = mu.';  
k = sym('k',[3,1]);
%R = sym('R',[3,1]); IT WOULD be easier to take R along z, but I haven't

%these are general placeholders for T_0 mu_k where T_0 is the rotation
%which maps mu_k T_1^(-1) to the z axis, see those notes you made
syms a b y R real
%a,b,y are the euler angles (normally alpha beta gamma or phi theta xi)

%Z1Y2Z3 choice
r1 = [cos(y),-sin(y),0;sin(y),cos(y),0;0,0,1];
r2 = [cos(a),0,sin(a);0,1,0;-sin(a),0,cos(a)]; %a runs from 0 to pi
r3 = [cos(b),-sin(b),0;sin(b),cos(b),0;0,0,1];
%r1(y)*r1(-y) = I
%(r3(b)*r2(a)*r1(y))  * ( r1(-y)*r2(-a)*r3(-b) )= I

TT = r3*r2*r1;
TTinv = subs(r1*r2*r3,{a,b,y},{-a,-b,-y} );

mutran = TT*mu;
%first order terms contain averages of mu(j,A)*mu(k,B)
%where j=1,2,3 etc is the output polarization direction, k is the direction
%of the polarization of the light, if this is not along one of the chosen
%axis there will be multiple non zero elements


expfct = 1i*(TT*[0;0;R]).' * k; %take R along z
%expfct_ap = 1 + expfct + expfct^2/2 + expfct^3/6; %3rd order
facto = [1,1,1/2,1/6];
%%  1st order terms
AV1 = sym(zeros(3,1)); 
for n = 0:2 %1st order and zeroth only
    tic
for j = 1:3

    tmp = int(expfct^n*facto(n+1)*mutran(j,1)*mutran(1,2),y,0,2*pi);
    tmp = int(tmp/pi,b,0,2*pi);
    tmp = int(tmp*sin(a)/pi,a,0,pi);

        AV1(j) = AV1(j) + simplify(tmp/8,'IgnoreAnalyticConstraints',true);

end
    toc
end
%%  3rd order terms
AV3 = sym(zeros(3,1));
main_fct = mutran(1,2)*mutran(1,3)*mutran(1,4);  %assume all fields along x
for n = 0:2 %first order is enough this is fucking slow
    tic
for j = 1:3
    
    tmp = int(expfct^n*facto(n+1)*mutran(j,1)*main_fct,y,0,2*pi);
    tmp = int(tmp/pi,b,0,2*pi);
    tmp = int(tmp*sin(a)/pi,a,0,pi);
    
        AV3(j) =AV3(j)+ simplify(tmp/8,'IgnoreAnalyticConstraints',true);%/pi^(3/2);
end
    toc
end

%% Save because of the long time taken
notee = strcat('this assumed all incoming beams were polarized in the x direction'...
        ,' but otherwise makes no assumptions about the resultant wavevec k '...
        ,'except that k.R is small and so can be taken to first order');
save('first_and_third_order_dimer_averages.mat','AV1','AV3','notee')

