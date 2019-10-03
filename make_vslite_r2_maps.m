% Show maps of skill during calibration/validation periods

load ./data/ITRDB_vslite.mat;

states = shaperead('cb_2017_us_nation_5m','UseGeoCoords', true);
latlim = [25.1 49.5];
lonlim = [-125 -66];

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 7];

lat = [ITRDB.LAT];
lon = [ITRDB.LON];

clr = cbrewer('seq','PuRd',10);

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.r2_cal];
subplot(4,2,1)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 8, r2, 'filled');
caxis([0 1]);
colormap(clr);
subplotsqueeze(gca, 1.3)
text(-0.38,0.86,'A', 'FontSize',12, 'FontWeight','bold')
text(-0.46,0.68,'Thornthwaite','FontSize',11,'Rotation',90,...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');
title('Calibration', 'FontSize',12);

r2 = [ITRDB.Th];
r2 = [r2.r2_val];
subplot(4,2,2)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 8, r2, 'filled');
caxis([0 1]);
colormap(clr);
subplotsqueeze(gca, 1.3)
text(-0.38,0.86,'B', 'FontSize',12, 'FontWeight','bold')
title('Validation', 'FontSize',12);

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.r2_cal];
subplot(4,2,3)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 8, r2, 'filled');
caxis([0 1]);
colormap(clr);
subplotsqueeze(gca, 1.3)
text(-0.38,0.86,'C', 'FontSize',12, 'FontWeight','bold')
text(-0.46,0.68,'Hargreaves','FontSize',11,'Rotation',90,...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');

r2 = [ITRDB.Hg];
r2 = [r2.r2_val];
subplot(4,2,4)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 8, r2, 'filled');
caxis([0 1]);
colormap(clr);
subplotsqueeze(gca, 1.3)
text(-0.38,0.86,'D', 'FontSize',12, 'FontWeight','bold')

% Priestly-Taylor
r2 = [ITRDB.PT];
r2 = [r2.r2_cal];
subplot(4,2,5)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 8, r2, 'filled');
caxis([0 1]);
colormap(clr);
subplotsqueeze(gca, 1.3)
text(-0.38,0.86,'E', 'FontSize',12, 'FontWeight','bold')
text(-0.46,0.68,'Priestly-Taylor','FontSize',11,'Rotation',90,...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');

r2 = [ITRDB.PT];
r2 = [r2.r2_val];
subplot(4,2,6)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 8, r2, 'filled');
caxis([0 1]);
colormap(clr);
subplotsqueeze(gca, 1.3)
text(-0.38,0.86,'F', 'FontSize',12, 'FontWeight','bold')

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.r2_cal];
subplot(4,2,7)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 8, r2, 'filled');
caxis([0 1]);
colormap(clr);
subplotsqueeze(gca, 1.3)
text(-0.38,0.86,'G', 'FontSize',12, 'FontWeight','bold')
text(-0.46,0.68,'Penman-Monteith','FontSize',11,'Rotation',90,...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');

r2 = [ITRDB.PM];
r2 = [r2.r2_val];
subplot(4,2,8)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 8, r2, 'filled');
caxis([0 1]);
colormap(clr);
subplotsqueeze(gca, 1.3)
text(-0.38,0.86,'H', 'FontSize',12, 'FontWeight','bold')

cb = colorbar('southoutside');
cb.Position = [0.1 0.06 0.8 0.0124];
cb.TickLength = 0.02;
xlabel(cb, 'R^{2}', 'FontSize',10);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/vslite-r2-maps.tif')
close all;

