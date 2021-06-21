% Show maps of skill during calibration/validation periods

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
ITRDB = ITRDB(ecol1 > 0);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

%% Make maps
states = shaperead('cb_2017_us_nation_5m','UseGeoCoords', true);
latlim = [25.1 49.5];
lonlim = [-125 -66];

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 7];

lat = [ITRDB.LAT];
lon = [ITRDB.LON];

clr = wesanderson('fantasticfox1');
clr = make_cmap([1 1 1; clr(3,:).^1.5; clr(3,:).^4; clr(3,:).^8], 7);

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.r2_val];
subplot(3,18,1:9)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 10, r2, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
caxis([0 0.7]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'a', 'FontSize',12)
title('Thornthwaite', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.05;

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.r2_val];
subplot(3,18,10:18)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 10, r2, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
caxis([0 0.7]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'b', 'FontSize',12)
title('Hargreaves', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.55;

% Priestley-Taylor
r2 = [ITRDB.PT];
r2 = [r2.r2_val];
subplot(3,18,19:27)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 10, r2, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
caxis([0 0.7]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'c', 'FontSize',12)
title('Priestley-Taylor', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.05;

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.r2_val];
subplot(3,18,28:36)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 10, r2, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
caxis([0 0.7]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'d', 'FontSize',12)
title('Penman-Monteith', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.55;

cb = colorbar('southoutside');
cb.Position = [0.1 0.37 0.8 0.02];
cb.TickLength = 0.027;
xlabel(cb, 'R^{2}', 'FontSize',10);
cbarrow('right');

%% Boxplots by ecoregion
clr = wesanderson('fantasticfox1');
ttl = {'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'};
for i = 1:length(ecos)
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    
    ind = [i*2+35 i*2+36];
    
    n = length(ITRDB_sub);
    
    r2 = [ITRDB_sub.Th];
    dat(:, 1) = [r2.r2_val];
    r2 = [ITRDB_sub.Hg];
    dat(:, 2) = [r2.r2_val];
    r2 = [ITRDB_sub.PT];
    dat(:, 3) = [r2.r2_val];
    r2 = [ITRDB_sub.PM];
    dat(:, 4) = [r2.r2_val];
    
    if i==1
        axes('Position',[0.02,0.03,0.12,0.22])
    else
        axes('Position',[0.1*(i-1)+0.07,0.03,0.07,0.22])
    end
    bxp = boxplot(dat, 'PlotStyle','compact', 'OutlierSize',3);
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    t = get(a,'tag');   % List the names of all the objects 
    set(a(17), 'Color', clr(4,:));   % Set the color of the first box to green
    set(a(18), 'Color', clr(3,:));   % Set the color of the first box to green
    set(a(19), 'Color', clr(2,:));   % Set the color of the first box to green
    set(a(20), 'Color', clr(1,:));   % Set the color of the first box to green
    set(a(21:24), 'Color', 'k');   % Set the color of the first box to green
    box off;
    ax = gca;
    set(ax, 'TickDir','out', 'TickLength',[0.02 0.02], 'YTick',0:0.2:1, 'YLim',[0 1])
    if i>1; set(ax, 'YTickLabel',''); end
    if i==1; ylabel('R^{2}'); end
    text(1,0.98,alphabet(i+4), 'FontSize',12);
    text(2.5, -0.06, ttl{i}, 'FontSize',10, 'HorizontalAlignment','center');
    
    clear dat;

end

h = findobj(gcf,'tag','Outliers');
for iH = 1:length(h)
    h(iH).MarkerEdgeColor = 'k';
end    
h = findobj(gcf,'tag','MedianOuter');
for iH = 1:length(h)
    h(iH).MarkerEdgeColor = 'k';
end    

hLegend = legend(findall(gca,'Tag','Box'), {'Penman-Monteith','Priestley-Taylor','Hargreaves','Thornthwaite'}, 'Orientation','horizontal');
hLegend.Position = [0.14    0.27    0.7179    0.0253];
hLegend.FontSize = 10;
legend('boxoff')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/vslite-r2-maps.tif')
close all;

