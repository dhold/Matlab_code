figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'FontSize',14);

 xlim(axes1,[700 900]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1,'LineWidth',3);

% Create xlabel
xlabel('Wavelength (nm)','FontSize',16);

% Create ylabel
ylabel('CD signal (relative to absorption)','FontSize',14);


%%
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'FontSize',14);
xlim(axes1,[700 900]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1,'LineWidth',3);

% Create xlabel
xlabel('Wavelength (nm)','FontSize',14);

% Create ylabel
ylabel('Linear absorption (AU)','FontSize',14);


%% cut throughs and fft
j = 1
pnts = [peaks,peaks+8,troughs,troughs+8]; 
to_plot2=zeros(length(pnts),floor(length(t_delay_range_fs)/2));
to_plot2y = to_plot2;
for k = 1:length(pnts)
  %  tmp = real(Spp_alpha(pnts(k),:,j))- mean(real(Spp_alpha(pnts(k),:,j)));
 %   tmp2 = real(Spp_CD(pnts(k),:,j))- mean(real(Spp_CD(pnts(k),:,j)));    
    tmp = real(Spp_f(pnts(k),:,j))- mean(real(Spp_f(pnts(k),:,j)));
    tmp2 = real(Spp_CDf(pnts(k),:,j))- mean(real(Spp_CDf(pnts(k),:,j)));    
 [to_plot1,to_plot2(k,:)] = ezfft(t_sep_rng,tmp);
[to_plot1y,to_plot2y(k,:)] = ezfft(t_sep_rng,tmp2);
end

rng = to_plot1<3000;
to_plot1 = to_plot1(rng); to_plot1y = to_plot1y(rng); 
to_plot2 = sqrt(to_plot2(:,rng)); to_plot2y = sqrt(to_plot2y(:,rng)); 

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XTick',[0 400 800 1200 1600 2000],...
    'XMinorTick','on',...
    'FontSize',14);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
%plot1 = plot(t_delay_range_fs,Spp_alpha(pnts,:,j),'Parent',axes1,'LineWidth',2);
plot1 = plot(to_plot1,to_plot2,'Parent',axes1,'LineWidth',2);
for el = 1:length(pnts)
set(plot1(el),'DisplayName',strcat('\lambda = ',num2str(lam_nm(pnts(el))),'nm'));
end

% Create xlabel
xlabel('angular frequency, cm^{-1}','FontSize',14);

% Create ylabel
ylabel('Fourier component amplitude','FontSize',14);

% Create title
title('Fourier components in the absorption signal');

% Create legend
legend(axes1,'show');


% Create figure
figure2 = figure;

% Create axes
axes2 = axes('Parent',figure2,'XTick',[0 400 800 1200 1600 2000],...
    'XMinorTick','on',...
    'FontSize',14);
 xlim(axes2,[0 2000]);

box(axes2,'on');
hold(axes2,'all');

% Create multiple lines using matrix input to plot
plot2 = plot(to_plot1y,to_plot2y,'Parent',axes2,'LineWidth',2);
for el = 1:length(pnts)
set(plot2(el),'DisplayName',strcat('\lambda = ',num2str(lam_nm(pnts(el))),'nm'));
end

% Create xlabel
xlabel('angular frequency, cm^{-1}','FontSize',14);

% Create ylabel
ylabel('Fourier component amplitude','FontSize',14);

title('Fourier components in the CD signal');

% Create legend
legend(axes2,'show');

