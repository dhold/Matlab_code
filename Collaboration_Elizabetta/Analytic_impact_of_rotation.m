%phase impact of rotation

syms a b c d k
%k vec along z direction
ry = [cos(a),0,sin(a);0,1,0;-sin(a),0,cos(a)];
rz = [cos(b),-sin(b),0;sin(b),cos(b),0;0,0,1];

R = [0;0;1]; %start along z, average over orientations

rot_op = rz*ry;
%rot_op = rz*ry*rz; %neither

mu = [sin(c)*cos(d);sin(c)*sin(d);cos(c)];
mu2 = rot_op*mu;

expfct = exp(1i*k*cos(a)); 
%factor due to the local phase exp(i k dot (R+R_cm)) for this orientation
%at the position R_cm

%field polarized along z... unusual (massive focus of light w angular mom will do this)
tmpzx = int(sin(a)*mu2(3)*mu2(1),b,0,2*pi); %integrate over b dependence
tmpzy = int(sin(a)*mu2(3)*mu2(2),b,0,2*pi);
tmpzz = int(sin(a)*mu2(3)^2,b,0,2*pi);

%note that int_0^2pi  da exp(1i*k*cos(a)) sin^(3-n)(a) cos^n(a) 
% is non zero only if n = 1, also
% int_0^2pi  da exp(1i*k*cos(a)) sin^(2)(a)  is non zero but sin(a)cos(a)
% vanishes

%field polarized along y... 
tmpyx = int(sin(a)*mu2(2)*mu2(1),b,0,2*pi); %integrate over b dependence
tmpyy = int(sin(a)*mu2(2)*mu2(2),b,0,2*pi);
tmpyz = int(sin(a)*mu2(2)*mu2(3),b,0,2*pi);

%field polarized along x... 
tmpxx = int(sin(a)*mu2(1)*mu2(1),b,0,2*pi); %integrate over b dependence
tmpxy = int(sin(a)*mu2(1)*mu2(2),b,0,2*pi);
tmpxz = int(sin(a)*mu2(1)*mu2(3),b,0,2*pi);

int(tmpxx*exp(1i*k*cos(a)),a,0,pi)
int(tmpyy*exp(1i*k*cos(a)),a,0,pi)
int(tmpzz*exp(1i*k*cos(a)),a,0,pi)
