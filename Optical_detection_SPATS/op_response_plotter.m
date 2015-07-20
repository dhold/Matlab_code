
navrange = 10.^(linspace(-2,2,15));
%navrange = [0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50]; %thermal average photon number
maddrange = [0,1,2,4,8,16,25,50,100];
g2save = zeros(length(navrange),length(maddrange)); stnsave = g2save; %Asatrange = [18,26];

for lp = 1:length(navrange)
    for lp2 = 1:length(maddrange)
nav = navrange(lp); % thermal average, i.e. average with NO added photons
madd = maddrange(lp2); %added single photon mode, 0 is standard thermal state
%note actual photon average is (exp(-beta) + m)/(1-exp(-beta))
%                              =  n_therm + m(1+n_therm)
%actualnav = nav + madd*(1+nav);
Beta = log((1+nav)/nav);
Asat = 18;
A0 = 0.5; sigA = 0.25; sigD = 0.1;
nmax = round(madd*(1+nav)+nav)*5+10;
%nmax = round(2*nav)+20; 
etaD =0.5; etaI = 0.01;
% ranges 18<Asat<26 pA,  0.5<A0<0.85 pA,  0.25<sigA<0.6, 0.1<sigD<0.2 
% 0.007 < eta < 0.01 (really really low) 
[Pph,Kav,Aav,KAav,Aavsq] = optical_response_to_SPATS(madd,nav,etaI,Asat,0.5,0.25,0.1,etaD,nmax);

g2 = KAav/Kav/Aav; %g2nosat = KAavnosat/Kav/Anosat;
stn = Aav/sqrt(Aavsq-Aav.^2); %stnnosat = Anosat/sqrt(Asqnosat-Anosat.^2);

g2save(lp,lp2) = g2;  %g2savenosat(lp) = g2nosat;  
stnsave(lp,lp2) = stn;  %stnsavenosat(lp) = stnnosat;  
    end
end
%%
figure0 = figure;
axes0 = axes('Parent',figure0,'XScale','log','XMinorTick','on',...
    'Position',[0.13 0.132231404958678 0.775 0.792768595041322],...
    'FontSize',12);

%plot(navrange,stnsave,'r')%,'XScale','log','XMinorTick','on')
%hold on 
[navr,maddr] =meshgrid(navrange,maddrange);
%actualnav = (navrange + madd*(1+navrange));
actualnav = (navr+ maddr.*(1+navr));
%expnsq = 3*madd*navrange.*(navrange+1) + madd^2*(1+navrange).^2 + navrange.*(2*navrange+1);
%expnsq_minus_nsq = (1+madd)*(navrange+1).*navrange;
expnsq_minus_nsq = (1+maddr).*(navr+1).*navr;
%plot(navrange,A0*etaI*actualnav/2./sqrt(sigD^2 + etaI.*actualnav/2.*...
%        (sigA^2 + A0^2*(1-etaI/2)) +etaI^2*A0^2.*expnsq_minus_nsq/2) ,'g')
 nosatSTN = A0*etaI*actualnav/2./sqrt(sigD^2 + etaI.*actualnav/2.*...
        (sigA^2 + A0^2*(1-etaI/2)) +etaI^2*A0^2.*expnsq_minus_nsq/2);
semilogx0 = semilogx(navrange,[stnsave.';nosatSTN]...
    ,'Parent',axes0,'LineWidth',2);
%set(semilogx1(1),'DisplayName','A_s = 18 pA');
%set(semilogx1(2),'LineStyle','--','Color',[0 1 0],...
%    'DisplayName','No saturation');    
    
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on',...
    'Position',[0.13 0.132231404958678 0.775 0.792768595041322],...
    'FontSize',12);

box(axes1,'on');
hold(axes1,'all');

nosatg2 = (2*navr.^2+maddr.^2 .*(1+navr).^2+maddr.*...
                (-1+2.*navr+3*navr.^2))./(maddr.*(1+navr)+navr).^2;

% Create multiple lines using matrix input to semilogx
semilogx1 = semilogx(navrange,[g2save.';nosatg2]...
    ,'Parent',axes1,'LineWidth',2);
set(semilogx1(1),'DisplayName','A_s = 18 pA');
set(semilogx1(2),'LineStyle','--','Color',[0 1 0],...
    'DisplayName','No saturation');

% Create xlabel
xlabel('Mean thermal photon number','FontSize',14);

% Create ylabel
ylabel('$g^{(2)}$','Interpreter','latex','FontSize',16);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.482911392405064 0.186294765840221 0.350632911392405 0.153236914600551]);

%end
if 1==0
plot(Arange,PA)
end
if 1==0 %checking difference between 

biomtab = genbinomtab(nmax);
ntherm = etaD*exp(-Beta)/(1-exp(-Beta));
n = 0:nmax;  %now make n the whole range

Pph2 =zeros(nmax+1,1);
Pph2(n(n>=madd)+1) = (biomtab(n(n>=madd)+1,madd+1)).'.*(ntherm/(1+ntherm)).^(n(n>=madd)-madd);
Pph2 = Pph2/sum(Pph2);

figure
plot([PK,Pph2])
hold on 

end