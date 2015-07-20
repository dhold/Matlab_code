N = 2; %number of chromophores

kr = sym('kr','real');  ku = sym('ku','real'); %probe and pump
k1 = sym([0;0;kr]); %take k1 along z w.l.o.g
e1 = [1;0;0]; %linear polarized along x;

theta_pu = sym('theta_pu','real');%pump linear polarizated at an angle theta_pu to z
k2 = sym(ku*[0;sin(theta_pu);cos(theta_pu)]); 
e2 = [1;0;0]; %also linear polarized along x;

mu = sym('mu',[N,3]);   R = sym('R',[N,3]);  mu = mu.';  R=R.';
%mu is an 3 X N matrix containing all the dipole elements vector elements
%stacked on top of one another
% R is another 3 X N matrix consisting of vectors connecting the sites to
% the centre of mass
syms a b y real 
%a,b,y are the euler angles (normally alpha beta gamma or phi theta xi)

% TT = [cos(b)*cos(a)*cos(y)-sin(b)*sin(y),...
%       sin(b)*cos(a)*cos(y)+cos(b)*sin(y),...
%     -sin(a)*cos(y);-cos(b)*cos(a)*sin(y)-sin(b)*cos(y),...
%     -sin(b)*cos(a)*sin(y)-cos(b)*cos(y),sin(a)*sin(y);...
%     cos(b)*sin(a),sin(b)*sin(a),cos(a)];

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

for A = 1:N
    for B = 1:N
        temp = sym(zeros(3));   DR = TT*(R(:,A)-R(:,B)) ; 
        exfct = exp(1i*dot(k1,DR))*sin(a);
        for j=1:3
            for k = 1:3
                %CAN integrate over b and y without worrying about
                %sin(a)*exp(icos(a)*(k_pr dot (R_A-R_B))) prefact
                if e1(k) ~=0
                      
                tmp = int(exfct*mutran(j,A)*mutran(k,B)*e1(k)/8/pi^2,y,0,2*pi);
                tmp = int(tmp,b,0,2*pi);
                tmp = int(tmp,a,0,pi);
                
                temp(j,k) = tmp;
                end
            end
        end
        AV{A,B} = temp;
    end
end


