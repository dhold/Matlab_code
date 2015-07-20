%this code was written to be used for calculate the state of a density
%matrix after a pump pulse needs parameters
% tstart, polvec, mu,R, L0, 
% L0 = -1i*(kron(eye(length(Htot)),Htot)-kron(Htot.',eye(length(Htot))));

kk = sym('k','real');  %pump mean wavevector
kvec = kk*[sin(pu_angle),0,cos(pu_angle)];
t=sym('t','real'); 
om = sym('om','real'); %pump resonant freq
phi = sym('phi','real');  %phase, occurs from k dot r =phi terms taken

syms a b c %Euler angles
ry = [cos(a),0,sin(a);0,1,0;-sin(a),0,cos(a)];
rz = [cos(b),-sin(b),0;sin(b),cos(b),0;0,0,1];
ry2 = [cos(c),0,sin(c);0,1,0;-sin(c),0,cos(c)];
rot_op = ry2*rz*ry; 

mu = Tmat*mu;  R = Tmat*R; %transform with specified rotation

   muop(1,k+1) = dot(mu(:,k),polvec);

%to select the correct phase component of the output polarization

sigma = FWHM/(2*sqrt(2*ln(2)));
E_env = exp(-t^2/2/sigma^2)/sigma/sqrt(2*pi); %normalised envelope

LE_fwd = sym(zeros(N+1)); LE_bak = LE_fwd;
for k = 1 : N
LE_fwd(1,k+1) = dot(mu(:,k),polvec)*exp( 1i*(dot(R(:,k),kvec))); 
LE_bak(1,k+1) = dot(mu(:,k),polvec)*exp(-1i*(dot(R(:,k),kvec)));
end
LE_fwd = LE_fwd+LE_fwd.'; LE_bak=LE_bak+LE_bak.';
%note these are not actually hermitian but the sum is
LE_fwd = kron(eye(length(LE_fwd)),LE_fwd)-kron(LE_fwd.',eye(length(LE_fwd)));
LE_bak = kron(eye(length(LE_fwd)),LE_bak)-kron(LE_bak.',eye(length(LE_fwd)));
%commutator super operators in Liouville space

TT = diag(exp(-diag(L0)*t)); %diagonal elements of system Liouvillian

barescalefct = double(TT^(-1)*diff(TT,t)); %should be just a number!
L0_test = L0 - barescalefct; %should have no diagonal elements
L0 = L0 - diag(diag(L0));
%note this is not the actual interaction picture liouvillian for
%computational reasons with the HEOM

temp = TT*(exp( 1i*(phi -om*t))*ones(size(L0)))*TT^(-1);
temp2 = TT*(exp(-1i*(phi -om*t))*ones(size(L0)))*TT^(-1);
LE_tmp1 = sym(zeros(size(temp))); LE_tmp2 = LE_tmp1;
%remove rapidly oscillating terms, i.e. oscillating fast than RWA_cutoff
for lp1 = 1:N+1  %must be a less fucking stupid way to do this....
    for lp2 = 1:N+1
        tmp = express_exp_series(temp(lp1,lp2));
        tmpp1 = sym(0);
        for lp3 = 1:size(tmp,2)
            tmpp = subs(diff(tmp{2,lp3},t),t,0); %picks linear comp
            if abs(imag(tmpp)) < RWA_cutoff
                tmpp1 = tmp{1,lp3}*exp(tmp{2,lp3});
            end
        end
        LE_tmp1(lp1,lp2) = tmpp1;        
        tmp2 = express_exp_series(temp2(lp1,lp2));
        tmpp1 = sym(0);
        for lp3 = 1:size(tmp,2)
            tmpp = subs(diff(tmp2{2,lp3},t),t,0); %picks linear comp
            if abs(imag(tmpp)) < RWA_cutoff
                tmpp1 = tmp2{1,lp3}*exp(tmp2{2,lp3});
            end
        end
        LE_tmp2(lp1,lp2) = tmpp1;
    end
end
LE_fwd = LE_fwd.*LE_tmp1;%this actually IS in the interaction picture as 
%LE is symbolic and time dependent 
LE_bak = LE_bak.*LE_tmp2;


