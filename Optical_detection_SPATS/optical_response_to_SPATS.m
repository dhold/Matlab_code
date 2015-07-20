function  [Pph,Kav,Aav,KAav,Aavsq,Anosat,Asqnosat,KAavnosat] = ...
    optical_response_to_SPATS(madd,ntherm,eta,Asat,A0,sigA,sigD,etaD,nmax,sqs)
%optical response for SPATS
%%
% IN variables, all currents in  peco Amps
%nmax = 300; %max n considered
%Arange = range of photocurrents to include
%madd = 1; % mode added
%Beta = 1; %inverse temperate
%eta; %overall quantum efficiency of detector
%Asat; %saturation current,
%A0;  %is the average, 
%sigA; %is the standard deviation of the photocurrent amplitude in 
%response to a single isomerization,
%sigD; %is the standard deviation of the photocurrent in darkness.
% etaD; %efficiency of detector
% in PRL 109, 113601 (2012) 
% Asat in range 18,26 pA,  0.5<A0<0.85 pA,  0.25<sigA<0.6, 0.1<sigD<0.2 
% 0.007 < eta < 0.01 (really really low) 
%
% OUT variables
%  Pph is the original input density matrix on diagonal
% Kav = average count on detector,Aav = Average photocurrent from rod
% KAav = coincidence average,  Aavsq =  Averaged squared of photocurrent
% 
%  rho = (a^dag)^m exp(-beta*a^dag*a)*a^m *(1-e^(-beta))^(m+1)/m! 
biomtab = genbinomtab(nmax);
%ntherm = exp(-Beta)/(1-exp(-Beta));

n = 0:nmax;  %now make n the whole range
%Pph =zeros(nmax+1,1); 
%Pph(n(n>=madd)+1) = (biomtab(n(n>=madd)+1,madd+1).').*exp(-Beta*(n(n>=madd)-madd));
Pph =zeros(nmax+1,1);

if nargin == 10 %squeezed state input if given a 10th argument
   r = sqs; u = tanh(r); z = sqrt(ntherm)*(1+tanh(r))/sqrt(2*tanh(r));   
   %tabulate hermite polynomials with prefactor (u/2)^(n/2) / sqrt(n!)
   Hp = zeros(nmax+1); Hp(1) = 1; Hp(2) = sqrt(2*u)*z;
   for nlp = 2:nmax
       Hp(nlp+1) = sqrt(2*u/(nlp))*Hp(nlp) - u*sqrt((nlp-1)/nlp)*Hp(nlp-1);
   end
   Pph = exp(-ntherm(1+u)).*Hp.^2/cosh(r);
    
else
Pph(n(n>=madd)+1) = ((biomtab(n(n>=madd)+1,madd+1)).')...
.*[1,cumprod( ones(1,nmax-madd).*ntherm/(1+ntherm))];
Pph = Pph./(1+ntherm)^(madd+1);
end

norm = sum(Pph) %print out norm so I can see how much I'm losing

Pph = Pph/norm; %renormalise to account for truncation   

%%Calculate expectation values

Kav = etaD*dot(Pph,n)/2; %factor of 2 from BS

%integral over this times sign(A)*min(|A_s|,|A|) gives

    s = 2*(sigD.^2+n*sigA.^2);
    
 expAcomp = (n*A0).*(erf((Asat+n.*A0)./sqrt(s)) + erf((Asat-n.*A0)./sqrt(s)) ) /2 ...
                + (Asat/2) .*(erfc((Asat+n*A0)./sqrt(s)) - erfc((Asat-n*A0)./sqrt(s))) ...
               - sqrt(s./(4*pi)).*(exp(-(Asat-n.*A0).^2 ./s) -exp(-(Asat+n.*A0).^2 ./s));

            
 expAsatqcomp = (s + 2*n.^2*A0^2).*(erf( (Asat + n*A0)./sqrt(s)) + erf( (Asat - n*A0)./sqrt(s)) ) /4 ... 
                    + (Asat^2/2) .*(erfc((Asat+n*A0)./sqrt(s)) - erfc((Asat-n*A0)./sqrt(s))) ...
        + sqrt(s./(4*pi)).*((n*A0-Asat).*exp(-(Asat+n.*A0).^2 ./s) - (Asat + n*A0).*exp(-(Asat-n.*A0).^2 ./s));
    
    %without saturation
 expAnosat = n*A0;                    
 expAsqnosat = s/2 + A0^2*n.^2;                   

%     expAcomp = (n*A0/2).*(1+erf((Asat-n.*A0)./(2*(sigD.^2+n*sigA.^2)))) - ...
%     sqrt((sigD.^2+n.*sigA.^2)./(2*pi)).*exp(-(Asat-n.*A0).^2./(2*(sigD.^2+n*sigA.^2))) ...
%      + (Asat/2) .* erfc((Asat-n*A0)./sqrt(2*(sigD.^2+n*sigA.^2)));    
%  
%  expAsatqcomp = ((sigD.^2+n*sigA.^2)/2+n.^2*A0^2/2).*...
%                 (1+erf((Asat-n.*A0)./(2*(sigD.^2+n*sigA.^2)))) + ...
%                sqrt((sigD.^2+n*sigA.^2)/2).*(n*A0-Asat).*exp((Asat-n*A0).^2/(2*(sigD.^2+n*sigA.^2)) - ...                
%     A0*n.*sqrt((sigD.^2+n.*sigA.^2)./pi).*exp(-(Asat-n.*A0).^2./(2*(sigD.^2+n*sigA.^2))) ...
%      + (Asat^2/2) .* erfc((Asat-n*A0)./sqrt(2*(sigD.^2+n*sigA.^2)));    
 
 Aav = 0; Aavsq = 0; 
 KAav = 0;   %combined average
 if nargout >= 6
 Anosat = 0; Asqnosat = 0; KAavnosat = 0;
 end
 
 etap = eta.^(0:nmax);   onem_etap = (1-eta).^(0:nmax);
 
for n = 0:nmax
    prefc = Pph(n+1)*2^(-n);
    for m=0:n
        
        mm = 0:m;

        %temp = dot(expAcomp(mm+1),...
        %            biomtab(m+1,mm+1).*eta.^(mm) .*(1-eta).^(m-mm));  

        temp = expAcomp(mm+1)*(biomtab(m+1,mm+1).*etap(mm+1) .*onem_etap(m-mm+1)).'; 
        temp2 =expAsatqcomp(mm+1)*(biomtab(m+1,mm+1).*etap(mm+1) .*onem_etap(m-mm+1))'; 
                
        Aav = Aav + prefc*biomtab(n+1,m+1)*temp;
        Aavsq = Aavsq + prefc*biomtab(n+1,m+1)*temp2;
        KAav =  KAav + prefc*biomtab(n+1,m+1)*etaD*(n-m)*temp; 
        
        if nargout >= 6
        tmp1 = expAcompnosat(mm+1)*(biomtab(m+1,mm+1).*etap(mm+1) .*onem_etap(m-mm+1)).';
        tmp2 = expAsatqcompnosat(mm+1)*(biomtab(m+1,mm+1).*etap(mm+1) .*onem_etap(m-mm+1))';
        Anosat = Anosat + prefc*biomtab(n+1,m+1)*tmp1;        
        Asqnosat = Anosat + prefc*biomtab(n+1,m+1)*tmp2; 
        KAavnosat = KAavnosat + prefc*biomtab(n+1,m+1)*etaD*(n-m)*tmp1; 
        end       

    end            
end

