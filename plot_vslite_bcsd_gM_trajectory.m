% Make Novick-like figure showing projected change in a) average number of
% water-limited months and b) average magnitude of water stress during
% those months

load ./output/ITRDB_projections.mat;
ITRDB = ITRDB(cellfun(@ischar, {ITRDB.EcoL1_Code}));

alphabet = 'abcdefghijklmnopqrstuvwxyz';
clr = wesanderson('moonrise4');
clr=clr([5 2 4 3],:);

ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));
models = ITRDB(1).RCP85.Models;

syear = 1950;
eyear = 2100;
years = syear:eyear;

% RCP8.5
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 5];

ax = tight_subplot(3,3,0.07,[0.1 0.02], [0.1 0.04]);

for i = 1:length(ecos)
    
    axes(ax(i))
    
    temp = ITRDB(ecol1 == ecos(i));
    n = length(temp);
    
    gM_time_hist = NaN(n, length(models), 4);
    gM_magn_hist = NaN(n, length(models), 4);
    gM_time_proj = NaN(n, length(models), 4);
    gM_magn_proj = NaN(n, length(models), 4);
    
    for k=1:n
        
        % Thornthwaite
        gT = repmat(temp(k).Th.Tplus0.gTseas', 1, length(models));
        gM_hist = squeeze(mean(temp(k).RCP85.Th.gM(:,years>=1981 & years<=2010,:), 2));
        gM_proj = squeeze(mean(temp(k).RCP85.Th.gM(:,years>=2071 & years<=2099,:), 2));
        gM_time_hist(k,:,1) = sum(gM_hist < gT);
        gM_time_proj(k,:,1) = sum(gM_proj < gT);
        for j=1:length(models)
            gM_magn_hist(k,j,1) = mean(gM_hist(gM_hist(:,j) < gT(:,j), j));
            gM_magn_proj(k,j,1) = mean(gM_proj(gM_proj(:,j) < gT(:,j), j));
        end
        
        % Hargreaves
        gT = repmat(temp(k).Hg.Tplus0.gTseas', 1, length(models));
        gM_hist = squeeze(mean(temp(k).RCP85.Hg.gM(:,years>=1981 & years<=2010,:), 2));
        gM_proj = squeeze(mean(temp(k).RCP85.Hg.gM(:,years>=2071 & years<=2099,:), 2));
        gM_time_hist(k,:,2) = sum(gM_hist < gT);
        gM_time_proj(k,:,2) = sum(gM_proj < gT);
        for j=1:length(models)
            gM_magn_hist(k,j,2) = mean(gM_hist(gM_hist(:,j) < gT(:,j), j));
            gM_magn_proj(k,j,2) = mean(gM_proj(gM_proj(:,j) < gT(:,j), j));
        end
        
        % Priestly-Taylor
        gT = repmat(temp(k).PT.Tplus0.gTseas', 1, length(models));
        gM_hist = squeeze(mean(temp(k).RCP85.PT.gM(:,years>=1981 & years<=2010,:), 2));
        gM_proj = squeeze(mean(temp(k).RCP85.PT.gM(:,years>=2071 & years<=2099,:), 2));
        gM_time_hist(k,:,3) = sum(gM_hist < gT);
        gM_time_proj(k,:,3) = sum(gM_proj < gT);
        for j=1:length(models)
            gM_magn_hist(k,j,3) = mean(gM_hist(gM_hist(:,j) < gT(:,j), j));
            gM_magn_proj(k,j,3) = mean(gM_proj(gM_proj(:,j) < gT(:,j), j));
        end
        
        % Penman-Monteith
        gT = repmat(temp(k).PM.Tplus0.gTseas', 1, length(models));
        gM_hist = squeeze(mean(temp(k).RCP85.PM.gM(:,years>=1981 & years<=2010,:), 2));
        gM_proj = squeeze(mean(temp(k).RCP85.PM.gM(:,years>=2071 & years<=2099,:), 2));
        gM_time_hist(k,:,4) = sum(gM_hist < gT);
        gM_time_proj(k,:,4) = sum(gM_proj < gT);
        for j=1:length(models)
            gM_magn_hist(k,j,4) = mean(gM_hist(gM_hist(:,j) < gT(:,j), j));
            gM_magn_proj(k,j,4) = mean(gM_proj(gM_proj(:,j) < gT(:,j), j));
        end
        
    end
    clear gT gM_hist gM_proj k j;
    
    gM_magn_hist = reshape(gM_magn_hist, [], 4);
    gM_magn_proj = reshape(gM_magn_proj, [], 4);
    gM_time_hist = reshape(gM_time_hist, [], 4);
    gM_time_proj = reshape(gM_time_proj, [], 4);
    
    %% Plot historical 
    t = linspace(0, 2*pi);
    a = nanstd(gM_magn_hist);
    b = nanstd(gM_time_hist);
    c = nanmean(gM_magn_hist);
    d = nanmean(gM_time_hist);
    
    scatter(d(1),c(1),25, clr(1,:), 'filled')
    hold on;
    scatter(d(2),c(2),25, clr(2,:), 'filled')
    scatter(d(3),c(3),25, clr(3,:), 'filled')
    scatter(d(4),c(4),25, clr(4,:), 'filled')
    
    %% Plot projected
    f = nanmean(gM_magn_proj);
    g = nanmean(gM_time_proj);
    
    q1 = arrow('Start',[d(1) c(1)], 'Stop',[g(1) f(1)], 'Length',15,...
        'FaceColor',clr(1,:), 'EdgeColor',clr(1,:));
    q2 = arrow('Start',[d(2) c(2)], 'Stop',[g(2) f(2)], 'Length',15,...
        'FaceColor',clr(2,:), 'EdgeColor',clr(2,:));
    q3 = arrow('Start',[d(3) c(3)], 'Stop',[g(3) f(3)], 'Length',15,...
        'FaceColor',clr(3,:), 'EdgeColor',clr(3,:));
    q4 = arrow('Start',[d(4) c(4)], 'Stop',[g(4) f(4)], 'Length',15,...
        'FaceColor',clr(4,:), 'EdgeColor',clr(4,:));
    
    %% Labels and legend
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    text(xlim(1)+0.05*diff(xlim), ylim(1)+0.95*diff(ylim), alphabet(i), 'FontSize',12)
    set(gca, 'TickDir','out','TickLength',[0.025 0])
    if i==8; xlabel('Mean number of water-limited months per year', 'FontSize',11); end
    if i==4; ylabel('Mean {\itg_{M}} during water-limited months', 'FontSize',11); end
    
    if i==3 
        legend([q1 q2 q3 q4], 'Th','Hg','PT','PM', 'Location','southwest', 'FontSize',8);
        legend('boxoff')
    end
    
end

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/bcsd-gM-trajectory-rcp85.tif')
close all;

%% RCP4.5
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 5];

ax = tight_subplot(3,3,0.07,[0.1 0.02], [0.1 0.04]);

for i = 1:length(ecos)
    
    axes(ax(i))
    
    temp = ITRDB(ecol1 == ecos(i));
    n = length(temp);
    
    gM_time_hist = NaN(n, length(models), 4);
    gM_magn_hist = NaN(n, length(models), 4);
    gM_time_proj = NaN(n, length(models), 4);
    gM_magn_proj = NaN(n, length(models), 4);
    
    for k=1:n
        
        % Thornthwaite
        gT = repmat(temp(k).Th.Tplus0.gTseas', 1, length(models));
        gM_hist = squeeze(mean(temp(k).RCP45.Th.gM(:,years>=1981 & years<=2010,:), 2));
        gM_proj = squeeze(mean(temp(k).RCP45.Th.gM(:,years>=2071 & years<=2099,:), 2));
        gM_time_hist(k,:,1) = sum(gM_hist < gT);
        gM_time_proj(k,:,1) = sum(gM_proj < gT);
        for j=1:length(models)
            gM_magn_hist(k,j,1) = mean(gM_hist(gM_hist(:,j) < gT(:,j), j));
            gM_magn_proj(k,j,1) = mean(gM_proj(gM_proj(:,j) < gT(:,j), j));
        end
        
        % Hargreaves
        gT = repmat(temp(k).Hg.Tplus0.gTseas', 1, length(models));
        gM_hist = squeeze(mean(temp(k).RCP45.Hg.gM(:,years>=1981 & years<=2010,:), 2));
        gM_proj = squeeze(mean(temp(k).RCP45.Hg.gM(:,years>=2071 & years<=2099,:), 2));
        gM_time_hist(k,:,2) = sum(gM_hist < gT);
        gM_time_proj(k,:,2) = sum(gM_proj < gT);
        for j=1:length(models)
            gM_magn_hist(k,j,2) = mean(gM_hist(gM_hist(:,j) < gT(:,j), j));
            gM_magn_proj(k,j,2) = mean(gM_proj(gM_proj(:,j) < gT(:,j), j));
        end
        
        % Priestly-Taylor
        gT = repmat(temp(k).PT.Tplus0.gTseas', 1, length(models));
        gM_hist = squeeze(mean(temp(k).RCP45.PT.gM(:,years>=1981 & years<=2010,:), 2));
        gM_proj = squeeze(mean(temp(k).RCP45.PT.gM(:,years>=2071 & years<=2099,:), 2));
        gM_time_hist(k,:,3) = sum(gM_hist < gT);
        gM_time_proj(k,:,3) = sum(gM_proj < gT);
        for j=1:length(models)
            gM_magn_hist(k,j,3) = mean(gM_hist(gM_hist(:,j) < gT(:,j), j));
            gM_magn_proj(k,j,3) = mean(gM_proj(gM_proj(:,j) < gT(:,j), j));
        end
        
        % Penman-Monteith
        gT = repmat(temp(k).PM.Tplus0.gTseas', 1, length(models));
        gM_hist = squeeze(mean(temp(k).RCP45.PM.gM(:,years>=1981 & years<=2010,:), 2));
        gM_proj = squeeze(mean(temp(k).RCP45.PM.gM(:,years>=2071 & years<=2099,:), 2));
        gM_time_hist(k,:,4) = sum(gM_hist < gT);
        gM_time_proj(k,:,4) = sum(gM_proj < gT);
        for j=1:length(models)
            gM_magn_hist(k,j,4) = mean(gM_hist(gM_hist(:,j) < gT(:,j), j));
            gM_magn_proj(k,j,4) = mean(gM_proj(gM_proj(:,j) < gT(:,j), j));
        end
        
    end
    clear gT gM_hist gM_proj k j;
    
    gM_magn_hist = reshape(gM_magn_hist, [], 4);
    gM_magn_proj = reshape(gM_magn_proj, [], 4);
    gM_time_hist = reshape(gM_time_hist, [], 4);
    gM_time_proj = reshape(gM_time_proj, [], 4);
    
    %% Plot historical 
    t = linspace(0, 2*pi);
    a = nanstd(gM_magn_hist);
    b = nanstd(gM_time_hist);
    c = nanmean(gM_magn_hist);
    d = nanmean(gM_time_hist);
    
    scatter(d(1),c(1),25, clr(1,:), 'filled')
    hold on;
    scatter(d(2),c(2),25, clr(2,:), 'filled')
    scatter(d(3),c(3),25, clr(3,:), 'filled')
    scatter(d(4),c(4),25, clr(4,:), 'filled')
    
    %% Plot projected
    f = nanmean(gM_magn_proj);
    g = nanmean(gM_time_proj);
    
    q1 = arrow('Start',[d(1) c(1)], 'Stop',[g(1) f(1)], 'Length',15,...
        'FaceColor',clr(1,:), 'EdgeColor',clr(1,:));
    q2 = arrow('Start',[d(2) c(2)], 'Stop',[g(2) f(2)], 'Length',15,...
        'FaceColor',clr(2,:), 'EdgeColor',clr(2,:));
    q3 = arrow('Start',[d(3) c(3)], 'Stop',[g(3) f(3)], 'Length',15,...
        'FaceColor',clr(3,:), 'EdgeColor',clr(3,:));
    q4 = arrow('Start',[d(4) c(4)], 'Stop',[g(4) f(4)], 'Length',15,...
        'FaceColor',clr(4,:), 'EdgeColor',clr(4,:));
    
    %% Labels and legend
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    text(xlim(1)+0.05*diff(xlim), ylim(1)+0.95*diff(ylim), alphabet(i), 'FontSize',12)
    set(gca, 'TickDir','out','TickLength',[0.025 0])
    if i==8; xlabel('Mean number of water-limited months per year', 'FontSize',11); end
    if i==4; ylabel('Mean {\itg_{M}} during water-limited months', 'FontSize',11); end
    
    if i==3 
        legend([q1 q2 q3 q4], 'Th','Hg','PT','PM', 'Location','southwest', 'FontSize',8);
        legend('boxoff')
    end
    
end

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/bcsd-gM-trajectory-rcp45.tif')
close all;

