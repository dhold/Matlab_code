%ALL THIS DOES IS PICK OUT Y COMPONENT AND THEN SHIFT IT A BIT IN ORDER TO
%SEE VARIATION
epsil = sym('ep');
X = sym('X');
L = sym('L'); alpha = sym('a'); eta = sym('n'); delta = sym('d'); 

MBS = [cos(X),1i*sin(X);1i*sin(X),cos(X)]; %45 degrees + retardation 2X
Manal = [sin(epsil)^2,sin(epsil)*cos(epsil);sin(epsil)*cos(epsil),cos(epsil)^2];

Msample = exp(-alpha*L/2).*[cosh(eta/4+1i/2*delta),-1i*sinh(eta/4+1i/2*delta);...
                            1i*sinh(eta/4+1i/2*delta),cosh(eta/4+1i/2*delta)];
                        
Mtot = Manal*MBS*Msample;
Eout = Mtot*[1;0];

syms z real
Msample2 = exp(-alpha*z/2).*[cosh(eta*z/4+1i/2*delta*z),-1i*sinh(eta*z/4+1i/2*delta*z);...
                            1i*sinh(eta*z/4+1i/2*delta*z),cosh(eta*z/4+1i/2*delta*z)];

test = taylor(Eout, [eta, delta,X,epsil],[0,0,0,0], 'Order',2)

test2 = taylor(Eout.*conj(Eout), [eta, delta,X,epsil],[0,0,0,0], 'Order',4)