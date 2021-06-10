% Make maps of VS-Lite M1 and M2 parameters

%% load tree-ring data
load ./data/ITRDB_simulations;
ITRDB = ITRDB(cellfun(@ischar, {ITRDB.EcoL1_Code}));
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ITRDB = ITRDB(ecol1 > 0);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

n = length(ITRDB);

%% Variables for figure
states = shaperead('cb_2017_us_nation_5m','UseGeoCoords', true);
latlim = [25.1 49.5];
lonlim = [-125 -66];
lat = [ITRDB.LAT];
lon = [ITRDB.LON];

%% Make maps of M1 parameter
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 7];

clr = wesanderson('fantasticfox1');
clr = make_cmap([1 1 1; clr(3,:).^3; clr(3,:).^8], 5);

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.M1];
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
caxis([0 0.1]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'a', 'FontSize',12)
title('Thornthwaite', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.05;

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.M1];
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
caxis([0 0.1]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'b', 'FontSize',12)
title('Hargreaves', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.55;

% Priestly-Taylor
r2 = [ITRDB.PT];
r2 = [r2.M1];
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
caxis([0 0.1]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'c', 'FontSize',12)
title('Priestly-Taylor', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.05;

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.M1];
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
caxis([0 0.1]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'d', 'FontSize',12)
title('Penman-Monteith', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.55;

cb = colorbar('southoutside');
cb.Position = [0.1 0.37 0.8 0.02];
cb.TickLength = 0.027;
cb.Ticks = 0:0.02:0.1;
xlabel(cb, 'M1 (vol/vol)', 'FontSize',10);
cbarrow('right');

% Add Boxplots by ecoregion
clr = wesanderson('fantasticfox1');
ttl = {'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'};
alphabet = 'abcdefghijklmnopqrstuvwxyz';
for i = 1:length(ecos)
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    
    ind = [i*2+35 i*2+36];
    
    n = length(ITRDB_sub);
    
    r2 = [ITRDB_sub.Th];
    dat(:, 1) = [r2.M1];
    r2 = [ITRDB_sub.Hg];
    dat(:, 2) = [r2.M1];
    r2 = [ITRDB_sub.PT];
    dat(:, 3) = [r2.M1];
    r2 = [ITRDB_sub.PM];
    dat(:, 4) = [r2.M1];
    
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
    set(ax, 'TickDir','out', 'TickLength',[0.02 0.02], 'YTick',0:0.02:0.12, 'YLim',[0 0.12])
    if i>1; set(ax, 'YTickLabel',''); end
    if i==1; ylabel('M1 (vol/vol)'); end
    text(1,ax.YLim(2),alphabet(i+4), 'FontSize',12);
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

hLegend = legend(findall(gca,'Tag','Box'), {'Penman-Monteith','Priestly-Taylor','Hargreaves','Thornthwaite'}, 'Orientation','horizontal');
hLegend.Position = [0.14    0.27    0.7179    0.0253];
hLegend.FontSize = 10;
legend('boxoff')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-vslite-M1-maps.tif')
close all;

%% Make maps of M2 parameter
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 7];

clr = wesanderson('fantasticfox1');
clr = make_cmap([1 1 1; clr(3,:).^1.5; clr(3,:).^4; clr(3,:).^8], 7);

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.M2];
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
caxis([0.15 0.5]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'a', 'FontSize',12)
title('Thornthwaite', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.05;

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.M2];
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
caxis([0.15 0.5]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'b', 'FontSize',12)
title('Hargreaves', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.55;

% Priestly-Taylor
r2 = [ITRDB.PT];
r2 = [r2.M2];
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
caxis([0.15 0.5]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'c', 'FontSize',12)
title('Priestly-Taylor', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.05;

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.M2];
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
caxis([0.15 0.5]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'d', 'FontSize',12)
title('Penman-Monteith', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.55;

cb = colorbar('southoutside');
cb.Position = [0.1 0.37 0.8 0.02];
cb.TickLength = 0.027;
cb.Ticks = 0.15:0.05:0.5;
xlabel(cb, 'M2 (vol/vol)', 'FontSize',10);
cbarrow;

% Add Boxplots by ecoregion
clr = wesanderson('fantasticfox1');
ttl = {'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'};
alphabet = 'abcdefghijklmnopqrstuvwxyz';
for i = 1:length(ecos)
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    
    ind = [i*2+35 i*2+36];
    
    n = length(ITRDB_sub);
    
    r2 = [ITRDB_sub.Th];
    dat(:, 1) = [r2.M2];
    r2 = [ITRDB_sub.Hg];
    dat(:, 2) = [r2.M2];
    r2 = [ITRDB_sub.PT];
    dat(:, 3) = [r2.M2];
    r2 = [ITRDB_sub.PM];
    dat(:, 4) = [r2.M2];
    
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
    set(ax, 'TickDir','out', 'TickLength',[0.02 0.02], 'YTick',0.05:0.05:0.55, 'YLim',[0.05 0.55])
    if i>1; set(ax, 'YTickLabel',''); end
    if i==1; ylabel('M2 (vol/vol)'); end
    text(1,ax.YLim(2),alphabet(i+4), 'FontSize',12);
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

hLegend = legend(findall(gca,'Tag','Box'), {'Penman-Monteith','Priestly-Taylor','Hargreaves','Thornthwaite'}, 'Orientation','horizontal');
hLegend.Position = [0.14    0.27    0.7179    0.0253];
hLegend.FontSize = 10;
legend('boxoff')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-vslite-M2-maps.tif')
close all;

