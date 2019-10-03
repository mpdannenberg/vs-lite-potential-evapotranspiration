% Map L1 ecoregions w/ tree-ring sites and AmeriFlux sites

load ./data/ITRDB_vslite.mat;

load ./data/conus_mask;
ecol1 = double(geotiffread('./data/us_EcoL1_4km.tif'));
ecol1(ecol1==255 | ecol1==0 | ecol1==150) = NaN;
ecol1 = ecol1.*CONUS;

clear CONUS;

ecos = unique( ecol1(isfinite(ecol1))' );

ecol1_reclass = NaN(size(ecol1));
for i=1:length(ecos)
    
    ecol1_reclass(ecol1 == ecos(i)) = i;
    
end

ecolabs = {'5.0 - Northern Forests','6.0 - Northwestern Forested Mountains',...
    '7.0 - Marine West Coast Forest','8.0 - Eastern Temperate Forests','9.0 - Great Plains',...
    '10.0 - North American Deserts','11.0 - Mediterranean California',...
    '12.0 - Southern Semi-Arid Highlands','13.0 - Temperate Sierras'};

clr = wesanderson('fantasticfox1');
cbrew = [clr(4,:)
    (clr(4,:).^(1/3))
    clr(3,:).^3
    clr(3,:)
    clr(1,:).^2
    clr(1,:)
    clr(2,:)
    sqrt(clr(2,:))
    clr(5,:)];

tlat = [ITRDB.LAT];
tlon = [ITRDB.LON];

% Time to make a map...
load ./data/prism_latlon;

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [24 49];
lonlim = [-125 -67];

h=figure('Color','w');
h.Units = 'inches';
h.Position = [1 4 6.5 3];

subplot(1,3,[1 2])
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',[28,34,40,46],'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica', 'GColor',[0.8 0.8 0.8],...
        'FLineWidth',1, 'FontColor',[0.4 0.4 0.4], 'MLabelLocation',20,...
        'MLabelParallel',24.01, 'FontSize',7);
surfm(PRISMlat, PRISMlon, ecol1_reclass);
geoshow(states,'FaceColor','none','LineWidth',0.4)
plotm(tlat, tlon, 'k.','MarkerSize',10);
caxis([0.5 9.5])
colormap(cbrew)
axis off;
axis image;
set(ax, 'Position',[-0.08 0.08 0.82 0.87]);

h = colorbar('eastoutside');
pos = get(h, 'Position');
set(h, 'Position',[0.64 0.15 0.03 0.7],'FontName','Times New Roman',...
    'XTick',1:1:9, 'XTickLabel',ecolabs, 'FontSize',9, 'TickLength',0);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/site-map.tif')
close all;


