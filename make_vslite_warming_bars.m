%% Calculate mean/SE of PET, soil moisture, and gM in each ecoregion
load ./data/ITRDB_simulations.mat;
alphabet = 'abcdefghijklmnopqrstuvwxyz';

clr = wesanderson('fantasticfox1');
idx = ~cellfun(@isempty, {ITRDB.EcoL1_Code});

ITRDB = ITRDB(idx);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

pet_T2_m = NaN(length(ecos), 4);
pet_T2_s = NaN(length(ecos), 4);
pet_T2_c = NaN(length(ecos), 2, 4);
M_T2_m = NaN(length(ecos), 4);
M_T2_s = NaN(length(ecos), 4);
M_T2_c = NaN(length(ecos), 2, 4);
gM_T2_m = NaN(length(ecos), 4);
gM_T2_s = NaN(length(ecos), 4);
gM_T2_c = NaN(length(ecos), 2, 4);
n = NaN(length(ecos), 1);

for i = 1:length(ecos)
    
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    n(i) = length(ITRDB_sub);
    
    % Thornthwaite
    model = [ITRDB_sub.Th];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet_T2_m(i, 1) = median(12*[T2.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T2_s(i, 1) = std(12*[T2.PET] - 12*[T0.PET]) / sqrt(length(model));
    pet_T2_c(i, :, 1) = bootci(1000, @median, (12*[T2.PET] - 12*[T0.PET]) );
    M_T2_m(i, 1) = median([T2.M] - [T0.M]);
    M_T2_s(i, 1) = std([T2.M] - [T0.M]) / sqrt(length(model));
    M_T2_c(i, :, 1) = bootci(1000, @median, ([T2.M] - [T0.M]) );
    gM_T2_m(i, 1) = median([T2.gM] - [T0.gM]);
    gM_T2_s(i, 1) = std([T2.gM] - [T0.gM]) / sqrt(length(model));
    gM_T2_c(i, :, 1) = bootci(1000, @median, ([T2.gM] - [T0.gM]) );
    
    % Hargreaves
    model = [ITRDB_sub.Hg];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet_T2_m(i, 2) = median(12*[T2.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T2_s(i, 2) = std(12*[T2.PET] - 12*[T0.PET]) / sqrt(length(model));
    pet_T2_c(i, :, 2) = bootci(1000, @median, (12*[T2.PET] - 12*[T0.PET]) );
    M_T2_m(i, 2) = median([T2.M] - [T0.M]);
    M_T2_s(i, 2) = std([T2.M] - [T0.M]) / sqrt(length(model));
    M_T2_c(i, :, 2) = bootci(1000, @median, ([T2.M] - [T0.M]) );
    gM_T2_m(i, 2) = median([T2.gM] - [T0.gM]);
    gM_T2_s(i, 2) = std([T2.gM] - [T0.gM]) / sqrt(length(model));
    gM_T2_c(i, :, 2) = bootci(1000, @median, ([T2.gM] - [T0.gM]) );
    
    % Priestly-Taylor
    model = [ITRDB_sub.PT];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet_T2_m(i, 3) = median(12*[T2.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T2_s(i, 3) = std(12*[T2.PET] - 12*[T0.PET]) / sqrt(length(model));
    pet_T2_c(i, :, 3) = bootci(1000, @median, (12*[T2.PET] - 12*[T0.PET]) );
    M_T2_m(i, 3) = median([T2.M] - [T0.M]);
    M_T2_s(i, 3) = std([T2.M] - [T0.M]) / sqrt(length(model));
    M_T2_c(i, :, 3) = bootci(1000, @median, ([T2.M] - [T0.M]) );
    gM_T2_m(i, 3) = median([T2.gM] - [T0.gM]);
    gM_T2_s(i, 3) = std([T2.gM] - [T0.gM]) / sqrt(length(model));
    gM_T2_c(i, :, 3) = bootci(1000, @median, ([T2.gM] - [T0.gM]) );
    
    % Penman-Monteith
    model = [ITRDB_sub.PM];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet_T2_m(i, 4) = median(12*[T2.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T2_s(i, 4) = std(12*[T2.PET] - 12*[T0.PET]) / sqrt(length(model));
    pet_T2_c(i, :, 4) = bootci(1000, @median, (12*[T2.PET] - 12*[T0.PET]) );
    M_T2_m(i, 4) = median([T2.M] - [T0.M]);
    M_T2_s(i, 4) = std([T2.M] - [T0.M]) / sqrt(length(model));
    M_T2_c(i, :, 4) = bootci(1000, @median, ([T2.M] - [T0.M]) );
    gM_T2_m(i, 4) = median([T2.gM] - [T0.gM]);
    gM_T2_s(i, 4) = std([T2.gM] - [T0.gM]) / sqrt(length(model));
    gM_T2_c(i, :, 4) = bootci(1000, @median, ([T2.gM] - [T0.gM]) );
    
end

%% Bargraph of +2C effect
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 6];

subplot(3,1,1)
b = bar(ecos,pet_T2_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
% for k1 = 1:size(pet_T2_m,2)
%     ctr(k1,:) = bsxfun(@plus, b(k1).XData, b(k1).XOffset');   % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
%     ydt(k1,:) = b(k1).YData;                                     % Individual Bar Heights
% end
% hold on
% errorbar(ctr', ydt', ydt'-squeeze(pet_T2_c(:,1,:)), squeeze(pet_T2_c(:,2,:))-ydt', '.k')
set(gca, 'TickDir','out','TickLength',[0.02 0.05],'XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'})
ylabel('\DeltaPET (mm)', 'FontSize',11);
ylim = get(gca, 'Ylim');
text(5,ylim(2),'a', 'FontSize',12, 'HorizontalAlignment','right', 'VerticalAlignment','top');
lgd = legend('Thornthwaite','Hargreaves','Priestly-Taylor','Penman-Monteith');
lgd.Position = [0.5965 0.89 0.2854 0.1068];
lgd.FontSize = 7;
legend('boxoff');

subplot(3,1,2)
b = bar(ecos,M_T2_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','YLim',[-0.015 0],'TickLength',[0.02 0.05],'XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'}, 'XAxisLocation','bottom', 'YDir','reverse')
ylabel('\DeltaSoil moisture (vol/vol)', 'FontSize',11);
ylim = get(gca, 'Ylim');
text(5,ylim(1),'b', 'FontSize',12, 'HorizontalAlignment','right', 'VerticalAlignment','top');

subplot(3,1,3)
b = bar(ecos,gM_T2_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05], 'XAxisLocation','bottom','XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'}, 'YDir','reverse')
ylabel('\Deltag_{M}', 'FontSize',11);
ylim = get(gca, 'Ylim');
text(5,ylim(1),'c', 'FontSize',12, 'HorizontalAlignment','right', 'VerticalAlignment','top');
xlabel('Ecoregion', 'FontSize',11);


set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/vslite-t+2-bars.tif')
close all;

