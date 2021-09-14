function fig_S1

% fig_S1.m
%
% Generates Figure S1 in "Does differential mortality after parental
% investment affect sex ratio evolution? No."

% Generate or load data
if(~exist('sexRatioSweep.mat','file'))
    sexRatioSweep;
end
load('sexRatioSweep.mat')

% Set figure dimensions
figure(3)
clf
set(gcf,'color','w')
set(gcf,'PaperUnits','centimeters')
xSize = 8; ySize = 8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 100 xSize*50 ySize*50])

% Plot data
subplot(2,2,1)
imagesc(1-mean(FREQ1(:,:,2,:),4),[0,0.55]);
set(gca,'ydir','normal')
set(gca,'xtick',0.5+linspace(0,res2,5),'xticklabel',UMJ(1:5:end)/uFJ)
set(gca,'ytick',0.5+linspace(0,res1,6),'yticklabel',linspace(0,1,6))
title('$q>0$','interpreter','latex','fontsize',14)
x1=xlabel('Relative juvenile male mortality, $\mu_{MJ}/\mu_{FJ}$','interpreter','latex','fontsize',14);

subplot(2,2,2)
imagesc(1-mean(FREQ1(:,:,1,:),4),[0,0.55]);
set(gca,'ydir','normal')
set(gca,'xtick',0.5+linspace(0,res2,5),'xticklabel',UMJ(1:5:end)/uFJ)
set(gca,'ytick',0.5+linspace(0,res1,6),'yticklabel',linspace(0,1,6))
title('$q=0$','interpreter','latex','fontsize',14)

subplot(2,2,3)
imagesc(1-mean(FREQ2(:,:,2,:),4),[0,0.55]);
set(gca,'ydir','normal')
set(gca,'xtick',0.5+linspace(0,res2,5),'xticklabel',UMJ(1:5:end)/uFJ)
set(gca,'ytick',0.5+linspace(0,res1,6),'yticklabel',linspace(0,1,6))
title('$q>0$','interpreter','latex','fontsize',14)
x2=xlabel('Relative adult male mortality, $\mu_{MA}/\mu_{FA}$','interpreter','latex','fontsize',14);
y1=ylabel('sex ratio at birth (proportion male), $s_2$','interpreter','latex','fontsize',14);

subplot(2,2,4)
imagesc(1-mean(FREQ1(:,:,1,:),4),[0,0.55]);
set(gca,'ydir','normal')
set(gca,'xtick',0.5+linspace(0,res2,5),'xticklabel',UMJ(1:5:end)/uFJ)
set(gca,'ytick',0.5+linspace(0,res1,6),'yticklabel',linspace(0,1,6))
title('$q=0$','interpreter','latex','fontsize',14)

% Add labels and adjust figures for presentation
labs = {'A','B','C','D'};
for i=1:4
    subplot(2,2,i)
    temp = get(gca,'position');
    
    temp(1) = temp(1) - 0.03*(2-mod(i,2));
    set(gca,'position',temp)
    text(-0.05,1.1*res1,strcat('\bf{',labs{i},'}'),'interpreter','latex','fontsize',14);
end

C = colorbar('EastOutside');

temp = get(x1,'position');
temp(1) = temp(1) + 0.3*res1;
set(x1,'position',temp)
temp = get(x2,'position');
temp(1) = temp(1) + 0.3*res1;
set(x2,'position',temp)
temp = get(y1,'position');
temp(2) = temp(2) + 1.32*res2;
set(y1,'position',temp)

temp = get(C,'position');
temp(1) = temp(1) + 0.12;
temp(4) = 0.815;
set(C,'position',temp)
set(C,'ytick',0:0.1:0.5)
ylabel(C,'$s_2$ frequency','interpreter','latex','fontsize',14);

% Update colormap
if(exist('cvidis.mat','file'))
    load('cvidis.mat')
    colormap(cvidis)
end

% Save to PDF
if(exist('save2pdf','file'))
    save2pdf('fig_S1.pdf')
end