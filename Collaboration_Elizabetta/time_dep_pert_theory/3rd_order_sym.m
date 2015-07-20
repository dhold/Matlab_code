syms f1 f2 f3 f4 g1 g2 g3 g4 h1 h2 h3 h4 i1 i2 i3 i4
V1 = [0,f1,f2,0;f1',0,0,f3;f2',0,0,f4;0,f3',f4',0];
V2 = subs(V1,{f1 f2 f3 f4},{g1 g2 g3 g4});
V3 = subs(V1,{f1 f2 f3 f4},{h1 h2 h3 h4});
V4 = subs(V1,{f1 f2 f3 f4},{i1 i2 i3 i4});

tmp = 2*imag(trace(V2*V3*V4*V1*rho_eq + V1*V3*V4*V2*rho_eq+...
                        V1*V2*V4*V3*rho_eq+V4*V3*V2*V1*rho_eq));
                    
%f1,f2,f3,f4,.... etc are all functions of 

% R1_t_tmp(tlp1,tlp2,tlp3) = trace(V2*V3*V4*V1*rho_eq);
% R2_t_tmp(tlp1,tlp2,tlp3) = trace(V1*V3*V4*V2*rho_eq); 
% R3_t_tmp(tlp1,tlp2,tlp3) = trace(V1*V2*V4*V3*rho_eq); 
% R4_t_tmp(tlp1,tlp2,tlp3) = trace(V4*V3*V2*V1*rho_eq); 

%% finds analytic solution for sum of pathways in terms of (time dep) sym vars
N=2; rho_eq = zeros(1+N+N*(N-1)/2); rho_eq(1,1) = 1;
A1 = sym('Fge', [1 N]); B1= sym('Fef', [N,N*(N-1)/2 ]) ;
A2 = sym('Gge', [1 N]); B2= sym('Gef', [N,N*(N-1)/2 ]) ;
A3 = sym('Hge', [1 N]); B3= sym('Hef', [N,N*(N-1)/2 ]) ;
A4 = sym('Ige', [1 N]); B4= sym('Ief', [N,N*(N-1)/2 ]) ;
V1 = sym(zeros(N*(N+1)/2+1)); 
V1(1,2:N+1) = A1; V1(2:N+1,N+2:end) = B1; V1=V1+V1';
V2 = sym(zeros(N*(N+1)/2+1)); 
V2(1,2:N+1) = A2; V2(2:N+1,N+2:end) = B2; V2=V2+V2';
V3 = sym(zeros(N*(N+1)/2+1)); 
V3(1,2:N+1) = A3; V3(2:N+1,N+2:end) = B3;  V3=V3+V3';
V4 = sym(zeros(N*(N+1)/2+1));
V4(1,2:N+1) = A4; V4(2:N+1,N+2:end) = B4;  V4=V4+V4';

tmp= 2*imag(trace(V2*V3*V4*V1*rho_eq + V1*V3*V4*V2*rho_eq+...
                        V1*V2*V4*V3*rho_eq+V4*V3*V2*V1*rho_eq));
              
% Fit functions to time dependent mu operators in order to be able 
% to evaluate these analytically                   