% Plot projected changes across PET models

load ./output/ITRDB_projections.mat;
ITRDB = ITRDB(cellfun(@ischar, {ITRDB.EcoL1_Code}));

alphabet = 'abcdefghijklmnopqrstuvwxyz';
clr = wesanderson('moonrise4');
clr=clr([5 2 4 3],:);

ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

syear = 1950;
eyear = 2100;
years = syear:eyear;

%% Trend in tree-ring width
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 4];

ax = tight_subplot(3,3,0.03,[0.08 0.02], [0.1 0.04]);

for i = 1:length(ecos)
    
    axes(ax(i))
    
    temp = ITRDB(ecol1 == ecos(i));
    temp = [temp.RCP85];
    
    % Thornthwaite
    temp2 = [temp.Th];
    temp2 = [temp2.TRW]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(1,:), 'LineWidth',1.5)
    hold on;
    
    % Hargreaves
    temp2 = [temp.Hg];
    temp2 = [temp2.TRW]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(2,:), 'LineWidth',1.5)
    
    % Priestly-Taylor 
    temp2 = [temp.PT];
    temp2 = [temp2.TRW]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(3,:), 'LineWidth',1.5)
    
    % Penman-Monteith
    temp2 = [temp.PM];
    temp2 = [temp2.TRW]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(4,:), 'LineWidth',1.5)
    
    box off;
    set(gca, 'TickDir','out', 'TickLength',[0.02 0], 'YLim',[0.8 1.05],...
        'FontSize',8)
    
    if rem(i, 3) ~= 1
        set(gca, 'YTickLabel','')
    end
    if i < 7
        set(gca, 'XTickLabel','')
    end
    xtickangle(-30)
    if i==4
        ylabel('Relative ring width', 'FontSize',11)
    end
    if i==3
        legend('Th','Hg','PT','PM', ...
            'Location','southwest','FontSize',8);
        legend('boxoff')
    end
    
    ylim = get(gca, 'YLim');
    text(1955, ylim(2), alphabet(i), 'FontSize',12)
    
    
end
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/bcsd-trend-trw.tif')
close all;
    
%% Trend in gM
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 4];

ax = tight_subplot(3,3,0.03,[0.08 0.02], [0.1 0.04]);

for i = 1:length(ecos)
    
    axes(ax(i))
    
    temp = ITRDB(ecol1 == ecos(i));
    temp = [temp.RCP85];
    
    % Thornthwaite
    temp2 = [temp.Th];
    temp2 = [temp2.gM_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(1,:), 'LineWidth',1.5)
    hold on;
    
    % Hargreaves
    temp2 = [temp.Hg];
    temp2 = [temp2.gM_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(2,:), 'LineWidth',1.5)
    
    % Priestly-Taylor - interesting... some of the climate models didn't
    % work here?
    temp2 = [temp.PT];
    temp2 = [temp2.gM_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(3,:), 'LineWidth',1.5)
    
    % Penman-Monteith
    temp2 = [temp.PM];
    temp2 = [temp2.gM_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(4,:), 'LineWidth',1.5)
    
    box off;
    set(gca, 'TickDir','out', 'TickLength',[0.02 0], 'YLim',[-15 5],...
        'FontSize',8)
    
    if rem(i, 3) ~= 1
        set(gca, 'YTickLabel','')
    end
    if i < 7
        set(gca, 'XTickLabel','')
    end
    xtickangle(-30)
    if i==4
        ylabel('\Deltag_{M} (%)', 'FontSize',11)
    end
    if i==3
        legend('Th','Hg','PT','PM', ...
            'Location','southwest','FontSize',8);
        legend('boxoff')
    end
    
    ylim = get(gca, 'YLim');
    text(1955, ylim(2), alphabet(i), 'FontSize',12)
    
    
end
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/bcsd-trend-gM.tif')
close all;

%% Trend in soil moisture
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 4];

ax = tight_subplot(3,3,0.03,[0.08 0.02], [0.1 0.04]);

for i = 1:length(ecos)
    
    axes(ax(i))
    
    temp = ITRDB(ecol1 == ecos(i));
    temp = [temp.RCP85];
    
    % Thornthwaite
    temp2 = [temp.Th];
    temp2 = [temp2.M_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(1,:), 'LineWidth',1.5)
    hold on;
    
    % Hargreaves
    temp2 = [temp.Hg];
    temp2 = [temp2.M_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(2,:), 'LineWidth',1.5)
    
    % Priestly-Taylor - interesting... some of the climate models didn't
    % work here?
    temp2 = [temp.PT];
    temp2 = [temp2.M_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(3,:), 'LineWidth',1.5)
    
    % Penman-Monteith
    temp2 = [temp.PM];
    temp2 = [temp2.M_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(4,:), 'LineWidth',1.5)
    
    box off;
    set(gca, 'TickDir','out', 'TickLength',[0.02 0], 'YLim',[-15 5],...
        'FontSize',8)
    
    if rem(i, 3) ~= 1
        set(gca, 'YTickLabel','')
    end
    if i < 7
        set(gca, 'XTickLabel','')
    end
    xtickangle(-30)
    if i==4
        ylabel('\DeltaM (%)', 'FontSize',11)
    end
    if i==3
        legend('Th','Hg','PT','PM', ...
            'Location','southwest','FontSize',8);
        legend('boxoff')
    end
    
    ylim = get(gca, 'YLim');
    text(1955, ylim(2), alphabet(i), 'FontSize',12)
    
    
end
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/bcsd-trend-soilmoisture.tif')
close all;

%% Trend in PET
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 4];

ax = tight_subplot(3,3,0.03,[0.08 0.02], [0.1 0.04]);

for i = 1:length(ecos)
    
    axes(ax(i))
    
    temp = ITRDB(ecol1 == ecos(i));
    temp = [temp.RCP85];
    
    % Thornthwaite
    temp2 = [temp.Th];
    temp2 = [temp2.PET_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(1,:), 'LineWidth',1.5)
    hold on;
    
    % Hargreaves
    temp2 = [temp.Hg];
    temp2 = [temp2.PET_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(2,:), 'LineWidth',1.5)
    
    % Priestly-Taylor - interesting... some of the climate models didn't
    % work here?
    temp2 = [temp.PT];
    temp2 = [temp2.PET_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(3,:), 'LineWidth',1.5)
    
    % Penman-Monteith
    temp2 = [temp.PM];
    temp2 = [temp2.PET_trend]; temp2(end, :) = NaN; % Remove 2100 since some models have it and some don't
    xbar = repmat(mean(temp2(years>=1981 & years<=2010, :)), length(years), 1);
    temp2 = 100 * (temp2 - xbar) ./ xbar;
    trw = nanmedian(temp2, 2);
    plot(years, trw, '-', 'Color',clr(4,:), 'LineWidth',1.5)
    
    box off;
    set(gca, 'TickDir','out', 'TickLength',[0.02 0], 'YLim',[-10 40],...
        'FontSize',8)
    
    if rem(i, 3) ~= 1
        set(gca, 'YTickLabel','')
    end
    if i < 7
        set(gca, 'XTickLabel','')
    end
    xtickangle(-30)
    if i==4
        ylabel('\DeltaPET (%)', 'FontSize',11)
    end
    if i==3
        lgd = legend('Th','Hg','PT','PM', ...
            'Location','north','FontSize',8);
        legend('boxoff')
        lgd.Position(2) = 0.84;
    end
    
    ylim = get(gca, 'YLim');
    text(1955, ylim(2), alphabet(i), 'FontSize',12)
    
    
end
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/bcsd-trend-pet.tif')
close all;


