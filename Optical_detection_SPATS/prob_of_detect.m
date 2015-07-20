
%parameters used in Teich et all 1982 
% eta = 0.2
% M = 0.5
% <d> ~ 30 

Drange = [0,0.5,10,15,20];
Nrange = 0:70;
etarange = [0.1,0.2,0.3,0.4]; %take larger value of eta
Trange = [1,2,3,5,8,10,12,14,17,20];
Mrange = [0,0.5,1,1.5,2];
maxn = 250;
binomtab = genbinomtab(maxn);
Psee_poisson = zeros(length(Drange),length(Nrange),length(Trange),length(Mrange),length(etarange));
Psee_fock = Psee_poisson;
%av_pois = zeros(length(Drange),length(Nrange));
%av_fock = zeros(length(Drange),length(Nrange));

for lp3 = 1:length(Trange) %threshold
    T = Trange(lp3);
    for lp4=1:length(Mrange) %amplification coefficient
        M = Mrange(lp4);
        if M ~=0
        inc_gam_save = gammainc(T,M*(0:maxn),'upper');
        else %assume NO multiplication noise
        inc_gam_save = zeros(maxn+1,1); inc_gam_save(T+1:end)=1;
        M=1;
        end
        %this value is equal to sum( exp(-Mn)*(Mn)^s/s!,{s,T,inf})
        
        for lp2 = 1:length(Nrange)  %average photon number
            N = Nrange(lp2);
            for lp1 = 1:length(Drange)
                D =Drange(lp1); %noise average
                for lp5 = 1:length(etarange)
                    eta = etarange(lp5);
           % poisson_input = exp(-(N+D)*eta)*[1,cumprod((N+D)*eta./(1:maxn))];
           %don't scale dark light with eta just makes it confusing
           poisson_input = exp(-(N*eta+D))*[1,cumprod((N*eta+D)./(1:maxn))];
            %this is convolved with the noise already
            atten_fock_input = zeros(maxn+1,1);
            atten_fock_input(1:N+1) = (binomtab(N+1,1:N+1)).'.*(1-eta)^N...
                                        .*[1;cumprod(ones(N,1)*eta./(1-eta))];
           % noise_input = exp(-(D)*eta)*[1,cumprod((D)*eta./(1:maxn))];
           noise_input = exp(-D)*[1,cumprod(D./(1:maxn))];
            %convolve to get distn
            
            noisefock_input = zeros(maxn+1,1);
            
            for nn = 0:maxn
                ll=max(0,nn-N):nn;
                noisefock_input(nn+1) = noise_input(ll+1)*atten_fock_input(nn+1-ll);
            end
            
            %av_pois(lp1,lp2) =  dot(poisson_input,0:maxn);
            %av_fock(lp1,lp2) =  dot(noisefock_input,0:maxn);
            
            Psee_poisson(lp1,lp2,lp3,lp4,lp5) = dot(poisson_input,inc_gam_save);
            Psee_fock(lp1,lp2,lp3,lp4,lp5) = dot(noisefock_input,inc_gam_save);
            
            %now calculate the probability of seeing
                end
            end
        end
    end
end
%% line graphs
           
figure
plot(Nrange*etarange(3),100*[Psee_poisson(2,:,10,4,3);Psee_fock(2,:,10,4,3)])  
xlabel('Average number of photons absorbed','FontSize',16);
ylabel('Percent chance of seeing','FontSize',16);
hold on
plot(Nrange*etarange(3),100*[Psee_poisson(3,:,10,4,3);Psee_fock(3,:,10,4,3)],'--')  
plot(Nrange*etarange(3),100*[Psee_poisson(4,:,10,4,3);Psee_fock(4,:,10,4,3)],'.-') 
%% difference between vacuum state input, no noise
% Create figure
figure1 = figure('PaperSize',[20.98404194812 29.67743169791]);
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.118357487922705 0.841285892634207 0.822680247926351],...
    'FontSize',14);
 xlim(axes1,[-0.169466329716072 10.3144046380259]);
 ylim(axes1,[-0.0159443427929058 1]);
box(axes1,'on');
hold(axes1,'all');

plot1 = plot(Nrange,[reshape(Psee_fock(1,:,1:4,1,end),size(Psee_fock,2),4),...
    reshape(Psee_poisson(1,:,1:4,1,end),size(Psee_fock,2),4)],'Parent',axes1);
set(plot1(1),'Marker','square','LineStyle','--','DisplayName','T=1');
set(plot1(2),'Marker','square','LineStyle','--','DisplayName','T=2');
set(plot1(3),'Marker','square','LineStyle','--','DisplayName','T=3');
set(plot1(4),'Marker','square','LineStyle','--','Color',[1 0 1],...
    'DisplayName','T=5');
set(plot1(5),'DisplayName','T=1');
set(plot1(6),'DisplayName','T=2');
set(plot1(7),'DisplayName','T=3');
set(plot1(8),'Color',[1 0 1],'DisplayName','T=5');

xlabel('Mean/exact number of photons in input pulse','FontSize',14);
ylabel('Probability of seeing','FontSize',14);
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.140498002171765 0.525927255187916 0.0998751560549313 0.380110062893085],...
    'FontSize',10);
%% small noise
figure2 = figure('PaperSize',[20.98404194812 29.67743169791]);
axes2 = axes('Parent',figure2,...
    'Position',[0.13 0.118357487922705 0.841285892634207 0.822680247926351],...
    'FontSize',14);
 xlim(axes2,[-0.169466329716072 10.3144046380259]);
 ylim(axes2,[-0.0159443427929058 1]);
box(axes2,'on');
hold(axes2,'all');
plot2 = plot(Nrange,[reshape(Psee_fock(2,:,1:4,1,end),size(Psee_fock,2),4),...
        reshape(Psee_poisson(2,:,1:4,1,end),size(Psee_fock,2),4)]);

set(plot2(1),'Marker','square','LineStyle','--','DisplayName','T=1');
set(plot2(2),'Marker','square','LineStyle','--','DisplayName','T=2');
set(plot2(3),'Marker','square','LineStyle','--','DisplayName','T=3');
set(plot2(4),'Marker','square','LineStyle','--','Color',[1 0 1],...
    'DisplayName','T=5');
set(plot2(5),'DisplayName','T=1');
set(plot2(6),'DisplayName','T=2');
set(plot2(7),'DisplayName','T=3');
set(plot2(8),'Color',[1 0 1],'DisplayName','T=5');
xlabel('Mean/exact number of photons in input pulse','FontSize',14);
ylabel('Probability of seeing','FontSize',14);
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.140498002171765 0.525927255187916 0.0998751560549313 0.380110062893085],...
    'FontSize',10);
%%

figure
plot(Nrange,reshape(Psee_fock(1,:,1:6,1,end),size(Psee_fock,2),6))
hold on
plot(Nrange,reshape(Psee_poisson(1,:,1:6,1,end),size(Psee_fock,2),6),'--')
%amplification noise only
plot(Nrange,reshape(Psee_fock(1,:,1:6,3,end),size(Psee_fock,2),6),'-')



%% realistic noise and thresh
figure
plot(Nrange,Psee_fock([3,4,5],:,end,1,end),'--')
hold on
plot(Nrange,Psee_poisson([3,4,5],:,end,1,end))
%%
figure
plot(Nrange,max(1-reshape(Psee_fock(3,:,:,1,end),length(Nrange),length(Trange)),eps(10)),'--')
hold on
plot(Nrange,max(1-reshape(Psee_poisson(3,:,:,1,end),length(Nrange),length(Trange)),eps(10)))


figure
plot(Trange,reshape(Psee_fock(3,2,:,1,end)-Psee_fock(3,1,:,1,end),length(Trange),1),'--')
hold on
plot(Trange,reshape(Psee_poisson(3,2,:,1,end)-Psee_poisson(3,1,:,1,end),length(Trange),1))

%%
figure
plot(Nrange,reshape(Psee_fock(1,:,1:6,1,1),size(Psee_fock,2),6),'Marker','square','LineStyle','--')
hold on
plot(Nrange,reshape(Psee_poisson(1,:,1:6,1,1),size(Psee_fock,2),6))
%%
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'FontSize',14);

box(axes1,'on');
hold(axes1,'all');

% Create surface
surface(Nrange*etarange(1),Drange,100*(Psee_poisson(:,:,4,4,1)-Psee_fock(:,:,4,4,1)))

% Create xlabel
xlabel('Average number of photons absorbed','FontSize',16);
% Create ylabel
ylabel('Noise input','FontSize',16);
% Create colorbar
%colorbar('peer',axes1,'FontSize',14);  

figure
surface(Nrange*etarange(2),Drange,100*(Psee_poisson(:,:,4,4,2)-Psee_fock(:,:,4,4,2)))
xlabel('Average number of photons absorbed','FontSize',16);
ylabel('Noise input','FontSize',16);
figure
surface(Nrange*etarange(3),Drange,100*(Psee_poisson(:,:,4,4,3)-Psee_fock(:,:,4,4,3)))  
xlabel('Average number of photons absorbed','FontSize',16);
ylabel('Noise input','FontSize',16);

%%
figure
surface(Nrange*etarange(1),Drange(4:end),200*(Psee_poisson(4:end,:,4,4,1)-Psee_fock(4:end,:,4,4,1))...
    ./(Psee_poisson(4:end,:,4,4,1)+Psee_fock(4:end,:,4,4,1)))
xlabel('Average number of photons absorbed','FontSize',16);
ylabel('Noise input','FontSize',16);
figure
surface(Nrange*etarange(2),Drange(4:end),200*(Psee_poisson(4:end,:,4,4,2)-Psee_fock(4:end,:,4,4,2))...
    ./(Psee_poisson(4:end,:,4,4,2)+Psee_fock(4:end,:,4,4,2)))
xlabel('Average number of photons absorbed','FontSize',16);
ylabel('Noise input','FontSize',16);
figure
surface(Nrange*etarange(3),Drange(4:end),200*(Psee_poisson(4:end,:,4,4,3)-Psee_fock(4:end,:,4,4,3))...
    ./(Psee_poisson(4:end,:,4,4,3)+Psee_fock(4:end,:,4,4,3)))
xlabel('Average number of photons absorbed','FontSize',16);
ylabel('Noise input','FontSize',16);