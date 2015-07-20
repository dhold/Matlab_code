mu = sym('mu',[4,3]); mu = mu.';  
%these are general placeholders for T_0 mu_k where T_0 is the rotation
%which maps mu_k T_1^(-1) to the z axis, see those notes you made
syms a b y x real 
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

AV1 = sym(zeros(3,1)); AV3 = sym(zeros(3,1));
syms holdvr holdvr2
%%  1st order terms

for j = 1:3

    tmp = int(mutran(j,1)*mutran(1,2),y,0,2*pi);
    tmp = int(tmp,b,0,2*pi);
    tmp = tmp*sin(a);
    %those two integrals can be performed quite easily analytically... the
    %final integral requires a little bit more thought
    %collect terms with sin^n(theta) cos^m(theta) prefactor
    if tmp~=0         
    tmp2 = 0;
    difvar =  subs(tmp,sin(a),holdvr);
    
            %find strings with powers of sin(a)^n for logic
    for n = 1:5

        difvar = diff(difvar,holdvr)/n; % pick terms of order sin(a)*n
        %divide by n to counter factor from differentiation
            if difvar == 0
                break %no higher terms
            end
        difvar2 = subs(difvar,holdvr,0); 
        difvar3 = subs(difvar2,cos(a),holdvr2);
        
        %deal with zero order terms in cos(a) first
        tmpp=subs(difvar3,holdvr2,0);
        prfct = (2/x)^(n/2)* gamma(1/2+n/2);
        if tmpp~=0
            
        %tmp2 = tmp2 + tmpp*beta((1+n)/2,(1+m)/2)*(1+x^2*(m+1)/2/(m+n+2) + ...
        %                            x^4*(m+1)*(m+3)/24/(m+n+2)/(m+n+4));
        tmp2 = tmp2 + tmpp*prfct*besselj(n/2,x);
        end                       

        
        for m =1:5
            %use same technique to pick out variables with cos(a)
            
            difvar3 = diff(difvar3,holdvr2)/m; %divide by m to counter factor
            if difvar3 == 0
                break %no higher terms
            end
            difvar4 = subs(difvar3,holdvr2,0); % prefactor to term with order cos^m 
            if difvar4~=0
                
if m==1
tmp2 = tmp2 + 1i*difvar4*prfct*besselj(1+n/2,x);  
elseif m==2
tmp2 = tmp2 + difvar4*prfct*(besselj(1+n/2,x)/x-besselj(2+n/2,x));      
elseif m==3
tmp2 = tmp2 + 1i*difvar4*prfct*(3*besselj(2+n/2,x)/x-besselj(3+n/2,x));      
elseif m==4
tmp2 = tmp2 + difvar4*prfct*(3*((6+n)/x^2-2)*besselj(3+n/2,x)/x...
            +(1-3/x^2)*besselj(4+n/2,x));      
elseif m==5
tmp2 = tmp2 + 1i*difvar4*prfct*(5*(3*(8+n)/x^2-2)*besselj(4+n/2,x)/x...
            +(1-15/x^2)*besselj(5+n/2,x));   
end
                
                
%                 if 2*floor(m/2)==m %m even
%             tmp2 = tmp2 + beta((1+n)/2,(1+m)/2)*difvar4*(1+x^2*(m+1)/2/(m+n+2) + ...
%                                     x^4*(m+1)*(m+3)/24/(m+n+2)/(m+n+4));
%                                 %take to 4th order in x which is k dot r
%                                 % This should be ample as x~1/1000
%                  else
%              tmp2 = tmp2 + beta((1+n)/2,(2+m)/2)*1i*x*difvar4*(1+...
%                  x^2*(m+2)/6/(m+n+3) + x^4*(m+2)*(m+4)/120/(m+n+3)/(m+n+5));
%                  end           
            end        
        end
    end
        AV1(j) = simplify(tmp2/8,'IgnoreAnalyticConstraints',true);%/pi^(3/2);
        %include factor of pi^(-3/2) later
    end
end

%%  3rd order terms

main_fct = mutran(1,2)*mutran(1,3)*mutran(1,4);  %assume all fields along x

for j = 1:3

    tmp = int(mutran(j,1)*main_fct,y,0,2*pi);
    tmp = int(tmp,b,0,2*pi);
    tmp = tmp*sin(a);
    %those two integrals can be performed quite easily analytically... the
    %final integral requires a little bit more thought
    %collect terms with sin^n(theta) cos^m(theta) prefactor
    if tmp~=0         
    tmp2 = 0;
    difvar =  subs(tmp,sin(a),holdvr);
  
            %find strings with powers of sin(a)^n for logic
            
    for n = 1:5

        difvar = diff(difvar,holdvr)/n; % pick terms of order sin(a)*n
        %divide by n to counter factor from differentiation
            if difvar == 0
                break %no higher terms
            end
        difvar2 = subs(difvar,holdvr,0); 
        difvar3 = subs(difvar2,cos(a),holdvr2);
        
        %deal with zero order terms in cos(a) first
        tmpp=subs(difvar3,holdvr2,0);
        prfct = (2/x)^(n/2)* gamma(1/2+n/2); %leave sqrt pi factor for a bit
        if tmpp~=0
            
        %tmp2 = tmp2 + tmpp*beta((1+n)/2,(1+m)/2)*(1+x^2*(m+1)/2/(m+n+2) + ...
        %                            x^4*(m+1)*(m+3)/24/(m+n+2)/(m+n+4));
        tmp2 = tmp2 + tmpp*prfct*besselj(n/2,x);
        end                       

        
        for m =1:5
            %use same technique to pick out variables with cos(a)
            
            difvar3 = diff(difvar3,holdvr2)/m; %divide by m to counter factor
            if difvar3 == 0
                break %no higher terms
            end
            difvar4 = subs(difvar3,holdvr2,0); % prefactor to term with order cos^m 
            if difvar4~=0
                
if m==1
tmp2 = tmp2 + 1i*difvar4*prfct*besselj(1+n/2,x);  
elseif m==2
tmp2 = tmp2 + difvar4*prfct*(besselj(1+n/2,x)/x-besselj(2+n/2,x));      
elseif m==3
tmp2 = tmp2 + 1i*difvar4*prfct*(3*besselj(2+n/2,x)/x-besselj(3+n/2,x));      
elseif m==4
tmp2 = tmp2 + difvar4*prfct*(3*((6+n)/x^2-2)*besselj(3+n/2,x)/x...
            +(1-3/x^2)*besselj(4+n/2,x));      
elseif m==5
tmp2 = tmp2 + 1i*difvar4*prfct*(5*(3*(8+n)/x^2-2)*besselj(4+n/2,x)/x...
            +(1-15/x^2)*besselj(5+n/2,x));   
end
            end
        end
    end
            
%     for n = 1:5
% 
%         difvar = diff(difvar,holdvr)/n; % pick terms of order sin(a)*n
%             if difvar == 0
%                 break %no higher terms
%             end  
%         difvar2 = subs(difvar,holdvr,0); 
%         difvar3 = subs(difvar2,cos(a),holdvr2);
%         
%         %deal with zero order terms in cos(a) first
%         tmpp=subs(difvar3,holdvr2,0);
%         if tmpp~=0
%             
%         tmp2 = tmp2 + tmpp*beta((1+n)/2,(1+m)/2)*(1+x^2*(m+1)/2/(m+n+2) + ...
%                                     x^4*(m+1)*(m+3)/24/(m+n+2)/(m+n+4));
%         end                       
%                     %   (2/x)^(n/2)*sqrt(pi)* gamma(1/2+n/2)*besselj(n/2,x);
%         
%         for m =1:5
%             %use same technique to pick out variables with cos(a)
%             
%             difvar3 = diff(difvar3,holdvr2)/m; 
%             if difvar3 == 0
%                 break %no higher terms
%             end  
%             difvar4 = subs(difvar3,holdvr2,0); % prefactor to term with order cos^m 
%             if difvar4~=0
%                 if 2*floor(m/2)==m %m even
%             tmp2 = tmp2 + beta((1+n)/2,(1+m)/2)*difvar4*(1+x^2*(m+1)/2/(m+n+2) + ...
%                                     x^4*(m+1)*(m+3)/24/(m+n+2)/(m+n+4));
%                                 %take to 4th order in x which is k dot r
%                                 % This should be ample as x~1/1000
%                  else
%              tmp2 = tmp2 + beta((1+n)/2,(2+m)/2)*1i*x*difvar4*(1+...
%                  x^2*(m+2)/6/(m+n+3) + x^4*(m+2)*(m+4)/120/(m+n+3)/(m+n+5));
%                  end           
%             end
%         end
%     end
        AV3(j) = simplify(tmp2/8,'IgnoreAnalyticConstraints',true);%/pi^(3/2);
        %include factor of pi^(-3/2) later
    end
end

%% expand to fifth order
AV1app = AV1*0;  AV3app = AV3*0;  
for j = 1:3
    
    AV1app(j) = simplify(taylor(AV1(j),x,0,'order',5));
    AV3app(j) = simplify(taylor(AV3(j),x,0,'order',5));
end

%% this sections gets the projection matrix for a general vector to the z axis

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
