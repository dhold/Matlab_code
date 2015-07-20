%This scripts finds the poles in a J(omega) which is a sum of underdamped
%modes
% 2 lambda * gamma*om_0*om / ( ( om_0^2 - om^2)^2 + gamma^2*om^2)
Temp =300; %units Kelvin
    light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); 
Beta = (2 * pi * hbar * light_speed * length_unit)/ ( Temp * boltz_const);
lambda = [1,5,6]; %amplitude
 gamma = [10,3,20]; %damping
 om_0 = [3,10,50]; %resonant value
 
 if 1==0
 om_rng = linspace(0,2*max(om_0+gamma));
 J = om_rng*0;
 for k = 1:length(lambda)
 J = J + (lambda(k).*gamma(k).*om_0(k).^2 ).*(om_rng./( (om_0(k).^2 ...
        - om_rng.^2).^2 + gamma(k).^2.*om_rng.^2));
 end
 J = J*2/pi;
 end
 % pole locations
 %poles = zeros(4,length(lambda)); 
 poles2 = zeros(4,length(lambda));
 lhpoles = zeros(2,length(lambda)); uhpoles = lhpoles;
 pmmat = [1,1,-1,-1;1,-1,1,-1];
for k = 1:length(lambda)
    %tmp = roots([1,0,gamma(k)^2 - 2*om_0(k)^2,0,om_0(k)^4]);
     %lower half plane important
    % lhplane(:,k) = (imag(tmp)<0);
    %poles(:,k) = tmp;
    tmp2 = pmmat(1,:).*sqrt(-gamma(k)^2+2*om_0(k)^2+ pmmat(2,:).*...
        sqrt(gamma(k)^2*(gamma(k)^2-4*om_0(k)^2)))/sqrt(2) ;
    [~,tmp3] = sort(real(tmp2)); tmp2= tmp2(tmp3);
    poles2(:,k) = tmp2; %numerically better
    lhpoles(:,k) = tmp2(imag(tmp2)<0);
    uhpoles(:,k) = tmp2(imag(tmp2)>0);
end

%calculate residues at the poles via 
% 2 \pi i F(p) = int_{cont} dw F(w)/(w-p) 

lg = abs(lhpoles(1,:)-lhpoles(2,:))<=eps(abs(lhpoles(1,:))); %repeated roots occur at critical damping
if any(lg)
   warning('I havent made this code work for critically damped modes yet, this will bug') 
end
res = lhpoles*0;
 for k = 1:length(lambda)

     tmp = (lambda(k).*gamma(k).*om_0(k).^2 );
 res(1,k) =  tmp.*lhpoles(1,k)/((1-exp(-beta*lhpoles(1,k))).*(...
        (lhpoles(1,k)-lhpoles(2,k)).*(lhpoles(1,k)-uhpoles(2,k))...
        .*(lhpoles(1,k)-uhpoles(1,k))));
 res(2,k)  =   tmp.*lhpoles(2,k)/((1-exp(-beta*lhpoles(2,k))).*(...
        (lhpoles(2,k)-lhpoles(1,k)).*(lhpoles(2,k)-uhpoles(2,k))...
        .*(lhpoles(2,k)-uhpoles(1,k))));
 end
res = res*(4i); %2 pi i *2/pi  

%note sums of pole resides should be real for underdamped modes.  Potentially this
%should also work for overdamped modes
xi = zeros(size(lambda));
phipm = zeros(2,length(lambda)); Phi_pm = phipm;
 for k = 1:length(lambda)
     
     xi(k) = sqrt(om_0(k)^2-gamma(k)^2/4);
     phipm(:,k) = gamma(k)/2 + 1i*[1;-1]*xi(k);
    Phi_pm(:,k) = ([1;-1].*1i*lambda(k)*om_0(k)^2/(2*xi(k))) .*...
                 (cot(phipm(:,k)*beta/2)-1i);
    %Phi_k(k) = - (4*gamma(k)*lambda(k)*om_0(k)^2/beta).*
                    
 end
 
 %% Construct operator with coefficients from poles

 op1 = zeros(numel(lhpoles),1); % first coeff is commutator component second is
 op2 = op1;                                           %anti commutator
 
 lg = abs(real(lhpoles))<10*eps(abs(real(lhpoles))); %find pure imaginary poles

  lg2 = abs(real(1i*phipm))<10*eps(abs(real(1i*lhpoles)));
  
  op1(reshape(lg,numel(op1),1)) = real(res(lg));
  op2(reshape(lg,numel(op2),1)) = imag(res(lg));
 
  %now the complex poles
  nlg = abs(real(lhpoles))>=10*eps(abs(real(lhpoles)));
  op1(reshape(nlg,numel(op1),1)) = [real(res(1,nlg(1,:))+res(2,nlg(1,:)))...
        +1i*imag(res(1,nlg(1,:))-res(2,nlg(1,:))), real(res(1,nlg(1,:))+res(2,nlg(1,:)))...
        -1i*imag(res(1,nlg(1,:))-res(2,nlg(1,:)))];
  op2(reshape(nlg,numel(op2),1)) = [imag(res(1,nlg(1,:))+res(2,nlg(1,:)))...
        +1i*real(res(1,nlg(1,:))-res(2,nlg(1,:))),imag(res(1,nlg(1,:))+res(2,nlg(1,:)))...
        -1i*real(res(1,nlg(1,:))-res(2,nlg(1,:)))];
  
  
 %% Poles from the 1/(1-exp(beta*omega)) terms
 Kappa = 10; %number of matsubara frequencies included
 matsu_freq = 2*pi*(1:Kappa)/beta;
  om_rng = -1i*matsu_freq; %pole locations
 J = om_rng*0;
 for k = 1:length(lambda)
 J = J + (lambda(k).*gamma(k).*om_0(k).^2 ).*(om_rng./( (om_0(k).^2 ...
        - om_rng.^2).^2 + gamma(k).^2.*om_rng.^2));
 end
 J = J*4i; %pole residues