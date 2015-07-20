function  [PA,PI,Pph,PK,stn,g2] = optical_response_to_SPATS(madd,Beta,eta,Asat,A0,sigA,sigD,Arange,etaD,nmax)
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
%  PA is probability of current A
%  PI is probability of isomerization of n molecules
%  Pph is density matrix on diagonal
%  rho = (a^dag)^m exp(-beta*a^dag*a)*a^m *(1-e^(-beta))^(m+1)/m! 
biomtab = genbinomtab(nmax);
%ntherm = exp(-Beta)/(1-exp(-Beta));

n = 0:nmax;  %now make n the whole range
Pph =zeros(nmax+1,1); 
Pph(n(n>=madd)+1) = (biomtab(n(n>=madd)+1,madd+1).').*exp(-Beta*(n(n>=madd)-madd));
%Pph2 =zeros(nmax+1,1);
%Pph2(n(n>=m)+1) = (biomtab(n(n>=m)+1,m+1)).'.*(ntherm/(1+ntherm)).^(n(n>=m)-m);
%cumprod( ones(0:nmax-m).*ntherm/(1+ntherm))
%Pph-Pph2
Pph(1:madd-1) = 0;

Pph = Pph.*(-expm1(-Beta))^(madd+1);
%sum(Pph)
Pph = Pph/sum(Pph); %renormalise to account for truncation

 %Pisom(mm+1,nn+1) = P(nn|mm) is prob of nn detections given mm photons
[nn,mm] = meshgrid(0:nmax,0:nmax);
Pisom = biomtab((0:nmax)+1,n+1).*eta.^(nn).*(1-eta).^(mm-nn); %P(n|m)

%PI(n) = sum_m P(n|m) P_{ph}(m) = sum_m  Pisom(m+1,n)*Pph(m+1)
PI = (Pisom')*Pph; 
% for nnn = 0:nmax
% PI2(nnn+1) = dot(Pisom(:,nnn+1),Pph(:));
% end
% PI-PI2'
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

Abar = trapz(Arange(1:end-1),(PA(1:end-1).').*Arange(1:end-1)) + Asat*PA(end);
Abarsq = trapz(Arange(1:end-1),(PA(1:end-1).').*Arange(1:end-1).^2) + Asat^2*PA(end);
varA = Abarsq - Abar.^2;
stn = Abar ./ sqrt(varA); %signal to noise

%calculate statistical distribution of attenuated light
[K,mm] = meshgrid(0:nmax);
%PK_given_n(m+1,K+1) = P(K|m)
PK_given_n = (biomtab((0:nmax)+1,(0:nmax)+1).*etaD.^(K).*(1-etaD).^(mm-K));
PK = (PK_given_n.') * Pph; 
%dot(0:nmax,PK)
%etaD*dot(0:nmax,Pph)

%PK is slightly different from a thermal distribution, the lower the
%component of the thermal distribution compared to single photon mode the
%more significant.
%average photon number is etaD * Pph average
%calculate g^(2), correlation function

g2 = 0*PA_given_n(1,:);  %g22 = g2;  
PA_given_m = Pisom * PA_given_n;

for m = 0:nmax
    for K = 0:nmax
        tic
       % for n = 0:nmax

       %     g2 = g2 + K.*Pph(m+1).*PK_given_n(m+1,K+1).*Pisom(m+1,n+1).*PA_given_n(n+1,:);
            
       % end
       % toc
        %tic
        g2 = g2 + K.*Pph(m+1).*PK_given_n(m+1,K+1).*PA_given_m(m+1,:);
       % toc
    end
end

g2= (trapz(Arange(1:end-1),Arange(1:end-1).*g2(1:end-1)) + Asat*g2(end));
g2 = g2 /dot(0:nmax,PK);
g2 = g2 /(trapz(Arange(1:end-1),Arange(1:end-1).*(PA(1:end-1).'))+Asat*PA(end));


% size(repmat(0:nmax,nmax+1,1).*PK_given_n)
% P_K_A_m = repmat(0:nmax,nmax+1,1).*PK_given_n .* (Pisom.')*trapz(Arange,repmat(Arange,nmax+1,1).*PA_given_n,2);
% size(P_K_A_m)
% P_KA = (P_K_A_m).' * Pph;
