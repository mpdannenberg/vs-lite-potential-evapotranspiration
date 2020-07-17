% VS-Lite conceptual figure
T = -5:30;
M = 0:0.01:0.5;

T1 = 5;
T2 = 20;
M1 = 0.03;
M2 = 0.3;

gT = (T-T1) / (T2-T1); gT(T<T1) = 0; gT(T>T2) = 1;
gM = (M-M1) / (M2-M1); gM(M<M1) = 0; gM(M>M2) = 1;

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 2];

subplot(1,2,1)
plot(T, gT, 'k-', 'LineWidth',1)
box off;
set(gca, 'TickDir','out', 'TickLength',[0.02 0])
xlabel(['Mean temperature (',char(176),'C)'])
ylabel('\itg_{T}')
hold on;
plot([T1 T1], [0 1], '--', 'Color',[0.4 0.4 0.4], 'LineWidth',0.5)
plot([T2 T2], [0 1], '--', 'Color',[0.4 0.4 0.4], 'LineWidth',0.5)
text(T1-0.2, 0.9, 'T1', 'HorizontalAlignment','right')
text(T2+0.2, 0.1, 'T2')
text(-13, 1, 'a', 'FontSize',12)

subplot(1,2,2)
plot(M, gM, 'k-', 'LineWidth',1)
box off;
set(gca, 'TickDir','out', 'TickLength',[0.02 0], 'XLim',[min(M) max(M)])
xlabel('Soil moisture (vol/vol)')
ylabel('\itg_{M}')
hold on;
plot([M1 M1], [0 1], '--', 'Color',[0.4 0.4 0.4], 'LineWidth',0.5)
plot([M2 M2], [0 1], '--', 'Color',[0.4 0.4 0.4], 'LineWidth',0.5)
text(M1+0.01, 0.9, 'M1')
text(M2+0.01, 0.1, 'M2')
text(-0.12, 1, 'b', 'FontSize',12)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/VS-Lite-conceptual.tif')
close all;
