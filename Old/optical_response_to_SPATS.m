function  [PA,PI,Pph,PK,stn,g2] = optical_response_to_SPATS(m,Beta,eta,Asat,A0,sigA,sigD,Arange,etaD,nmax)
%optical response for SPATS
%%
% IN variables, all currents in  peco Amps
%nmax = 300; %max n considered
%Arange = range of photocurrents to include
%m = 1; % mode added
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
%  PA is probability of current A
%  PI is probability of isomerization of n molecules
%  Pph is density matrix on diagonal
%  rho = (a^dag)^m exp(-beta*a^dag*a)*a^m *(1-e^(-beta))^(m+1)/m! 
biomtab = genbinomtab(nmax);
%ntherm = exp(-Beta)/(1-exp(-Beta));

%n! = factlist(n+2)
n = 0:nmax;  %now make n the whole range
Pph =zeros(nmax+1,1); 
Pph(n(n>=m)+1) = (biomtab(n(n>=m)+1,m+1).').*exp(-Beta*(n(n>=m)-m));
%Pph2 =zeros(nmax+1,1);
%Pph2(n(n>=m)+1) = (biomtab(n(n>=m)+1,m+1)).'.*(ntherm/(1+ntherm)).^(n(n>=m)-m);
%Pph-Pph2
Pph(1:m-1) = 0;

Pph = Pph.*(1-exp(-Beta))^(m+1);
%sum(Pph)
Pph = Pph/sum(Pph); %renormalise to account for truncation

 %Pisom(mm+1,nn+1) is prob of nn detections given mm photons
[nn,mm] = meshgrid(0:nmax);
Pisom = biomtab((0:nmax)+1,n+1).*eta.^(nn).*(1-eta).^(mm-nn); %P(n|m)

%PI(n) = sum_mm P(n,mm) Pph(mm)
PI = Pisom'*Pph; 
% for kk = 1:length(Pisom)
% PI2(kk) = dot(Pisom(:,kk),Pph);
% end

%PI = PI/(sum(PI));
if isempty(Arange)
Arange = linspace(0,Asat,10000);
elseif length(Arange) == 1 && floor(Arange)==Arange %integer test
 Arange = linspace(0,Asat,Arange);     
end

[AA,nn] = meshgrid(Arange,n);
%size(AA)
%size(nn)
%PA_given_n = is prob of current A from n isomerisations
% PA_given_n (n+1, k) = P_s(A(k) | n) 

%C2= C1.*(erf( A0*n./(sqrt(2)*sqrt(n*sigA^2 + sigD^2))) + ...
%        erf( (Asat - A0*n)./(sqrt(2)*sqrt(n*sigA^2 + sigD^2))) ) ;
    
PA_given_n = exp(-(AA-nn*A0).^2 ./ (2*(sigD.^2+nn*sigA.^2))) ./ ...
            (1 + erf( A0*nn./(sqrt(2)*sqrt(nn*sigA^2 + sigD^2)))) .* ...
            sqrt(2)./sqrt(pi*(sigD.^2+nn*sigA.^2));
PA_given_n(:,end) = 1-trapz(Arange(1:end-1),PA_given_n(:,1:end-1),2); %saturated               
%Note that I am not explicitly making this so trapz will give the correct 
% integral over all the range, the last value MUST be treated explicitly

%normcheck = trapz(Arange(1:end-1),PA_given_n(:,1:end-1),2)+PA_given_n(:,end);
%normcheck should be ~ 1 for all n
%normcheck
%sum(normcheck-1)

PA = PA_given_n'*PI;
%PA(end) = 1 - trapz(Arange(1:end-1),PA(1:end-1)); 
%saturation condition now done earlier

%calculate statisatical moments

Abar = trapz(Arange,PA.'.*Arange);
Abarsq = trapz(Arange,PA.'.*Arange.^2);
varA = Abarsq - Abar.^2;
stn = Abar ./ sqrt(varA); %signal to noise

%calculate statistical distribution of attenuated light
[K,mm] = meshgrid(0:nmax);
%PK_given_n(m+1,K+1) = P(K|m)
PK_given_n = (biomtab((0:nmax)+1,(0:nmax)+1).*etaD.^(K).*(1-etaD).^(mm-K));
PK = PK_given_n.' * Pph; 

for lp = 
PKtest = PK*0;
end
%PK is slightly different from a thermal distribution, the lower the
%component of the thermal distribution compared to single photon mode the
%more significant.
%average photon number is etaD * Pph average
%calculate g^(2), correlation function

g2 = 0*PA_given_n(1,:); 

for m = 0:nmax
    for K = 0:nmax
        for n = 0:nmax

            g2 = g2 + K.*Pph(m+1).*PK_given_n(m+1,K+1).*Pisom(m+1,n+1).*PA_given_n(n+1,:);
            
        end
    end
end

g2= (trapz(Arange(1:end-1),Arange(1:end-1).*g2(1:end-1)) + Asat*g2(end));
g2 = g2 /dot(0:nmax,PK);
g2 = g2 /(trapz(Arange(1:end-1),Arange(1:end-1).*(PA(1:end-1).'))+Asat*PA(end));


% size(repmat(0:nmax,nmax+1,1).*PK_given_n)
% P_K_A_m = repmat(0:nmax,nmax+1,1).*PK_given_n .* (Pisom.')*trapz(Arange,repmat(Arange,nmax+1,1).*PA_given_n,2);
% size(P_K_A_m)
% P_KA = (P_K_A_m).' * Pph;

%Possible message to the boss
%Hope everything is going well over in Columbia. Just had a couple of questions I wanted to ask.

%I am trying to follow through the method of the PRL by Nigel Sim et al (109, 113601 (2012)) with single photon added thermal states. I can follow everything through up to the point where they give an equation for g^(2), which is eq. (6). It seems to me that K and A are independent and P_s (A) and P(K) are functions are A and K only, so the whole thing should just separate.

%Clearly this is not the case as both of them are dependent on the input photon distribution, in our case P(K) is not simply a rescaling of P_{ph}(m) because an attenuated SPATS is not another SPATS but otherwise the situation is very much the same. I see how A and K function of "m" from the input photon distribution "P_{ph}(m)", but it's not clear that evaluating int_0^infinity dA(m) sum_{k=1}^infinity K(m) A(m) P(A(m)) P(K(m)) should change 
 