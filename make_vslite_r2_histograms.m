% Plot distributions of validation metrics

load ./data/ITRDB_vslite.mat;
alphabet = 'abcdefghijklmnopqrstuvwxyz';

%% Add EcoRegions 
ecoL3 = shaperead('D:\Data_Analysis\EcoRegions\NA_CEC_Eco_Level3_GEO.shp', 'UseGeoCoords',true);

for i = 1:length(ecoL3)
    [IN, ON] = inpolygon([ITRDB.LAT], [ITRDB.LON], ecoL3(i).Lat, ecoL3(i).Lon);
    
    if sum(IN)>0 | sum(ON)>0 
        [ITRDB(IN==1 | ON==1).EcoL3_Code] = deal(ecoL3(i).NA_L3CODE);
        [ITRDB(IN==1 | ON==1).EcoL3_Name] = deal(ecoL3(i).NA_L3NAME);
        [ITRDB(IN==1 | ON==1).EcoL2_Code] = deal(ecoL3(i).NA_L2CODE);
        [ITRDB(IN==1 | ON==1).EcoL2_Name] = deal(ecoL3(i).NA_L2NAME);
        [ITRDB(IN==1 | ON==1).EcoL1_Code] = deal(ecoL3(i).NA_L1CODE);
        [ITRDB(IN==1 | ON==1).EcoL1_Name] = deal(ecoL3(i).NA_L1NAME);
    end
    
end
clear IN ON ecoL3;

idx = ~cellfun(@isempty, {ITRDB.EcoL1_Code});

ITRDB = ITRDB(idx);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

%% Make figure
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 4];

clr = wesanderson('fantasticfox1');

T = table({'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'}',...
    {'Northern Forests','Northwestern Forested Mountains',...
    'Marine West Coast Forest','Eastern Temperate Forests','Great Plains',...
    'North American Deserts','Mediterranean California',...
    'Southern Semi-Arid Highlands','Temperate Sierras'}', 'VariableNames',{'Code','Name'});
n = NaN(9,1);
Th = cell(9,1);
Hg = cell(9,1);
PT = cell(9,1);
PM = cell(9,1);

for i = 1:length(ecos)
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    subplot(3,3,i)
    x = 0.05:0.1:0.8;
    n(i) = length(ITRDB_sub);
    
    % Thornthwaite
    r2 = [ITRDB_sub.Th];
    r2 = [r2.r2_val];
    [f] = histcounts(r2, 0:0.1:0.8);
    p1 = plot(x, f, '-', 'Color',clr(1,:), 'LineWidth',2);
    hold on;
    scatter(x, f, 20, clr(1,:), 'filled')
    set(gca, 'XLim',[0 0.8], 'XTick',0:0.2:0.8, 'TickDir','out','TickLength',[0.02 0.02])
    box off;
    ci = bootci(1000,@median,r2);
    Th{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % Hargreaves
    r2 = [ITRDB_sub.Hg];
    r2 = [r2.r2_val];
    [f] = histcounts(r2, 0:0.1:0.8);
    p2 = plot(x, f, '-', 'Color',clr(2,:), 'LineWidth',2);
    scatter(x, f, 20, clr(2,:), 'filled')
    ci = bootci(1000,@median,r2);
    Hg{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % Priestly-Taylor
    r2 = [ITRDB_sub.PT];
    r2 = [r2.r2_val];
    [f] = histcounts(r2, 0:0.1:0.8);
    p3 = plot(x, f, '-', 'Color',clr(3,:), 'LineWidth',2);
    scatter(x, f, 20, clr(3,:), 'filled')
    ci = bootci(1000,@median,r2);
    PT{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % Penman-Monteith
    r2 = [ITRDB_sub.PM];
    r2 = [r2.r2_val];
    [f] = histcounts(r2, 0:0.1:0.8);
    p4 = plot(x, f, '-', 'Color',clr(4,:), 'LineWidth',2);
    scatter(x, f, 20, clr(4,:), 'filled')
    ci = bootci(1000,@median,r2);
    PM{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    if i == 3
        lgd = legend([p1 p2 p3 p4],'Th','Hg','PT','PM');
        legend('boxoff')
        lgd.FontSize = 8;
        lgd.Position = [0.76 0.8 0.2035 0.1497];
    end
    
    if i == 8
        xlabel('R^{2}', 'FontSize',12)
    end
    
    if i == 4
        ylabel('Count', 'FontSize',12)
    end
    
    if i<7
        set(gca, 'XTickLabels','');
    end
    
    ax = gca;
    set(ax, 'YLim',[ax.YLim(1) ax.YLim(2)*1.2])
    text(0.05, ax.YLim(2), alphabet(i), 'FontSize',12)

end

T.n = n;
T.Th = Th;
T.Hg = Hg;
T.PT = PT;
T.PM = PM;

%% All
T.Name(10) = {'All'};
T.n(10) = length(ITRDB);

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.r2_val];
ci = bootci(1000,@median,r2);
T.Th{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.r2_val];
ci = bootci(1000,@median,r2);
T.Hg{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Priestly-Taylor
r2 = [ITRDB.PT];
r2 = [r2.r2_val];
ci = bootci(1000,@median,r2);
T.PT{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.r2_val];
ci = bootci(1000,@median,r2);
T.PM{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

%% Save figure and table
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/vslite-r2-histogram-byEcoregion.tif')
close all;

writetable(T, './output/vslite-r2-median-ci.csv');

