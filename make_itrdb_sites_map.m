% Map L1 ecoregions w/ tree-ring sites and AmeriFlux sites

load ./data/ITRDB_vslite.mat;
% Add EcoRegions 
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
clear IN ON ecoL3 i;
ITRDB = ITRDB(~cellfun(@isempty,{ITRDB.EcoL1_Code}));

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

spc = {ITRDB.SPECIES};
unSpc = sort(unique(spc));
nSpc = NaN(size(unSpc));
for i = 1:length(unSpc)
    nSpc(i) = sum(strcmp(unSpc{i}, spc));
end


% Time to make a map...
load ./data/prism_latlon;

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [24 49];
lonlim = [-125 -67];

h=figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 4.5];

% Bar graph of species
nt = 5;
subplot(3,3,7:9)
[B,I] = sort(nSpc, 'descend');
bar(B(B>nt), 'FaceColor',[0.8 0.8 0.8]);
box off;
set(gca, 'YLim',[0 200], 'XTick',1:length(B(B>nt)), 'TickDir','out',...
    'XTickLabels',unSpc(I(B>nt)), 'XTickLabelRotation',90, 'FontSize',9,...
    'XLim',[0 length(B(B>nt))+1]);
ylabel('\itN', 'Rotation',0, 'HorizontalAlignment','right', 'VerticalAlignment','middle')
ax = gca;
ax.Position(1) = 0.1;
ax.Position(2) = 0.12;
ax.Position(3) = 0.85;
ax.Position(4) = 0.28;

text(-3, 180, 'b', 'FontSize',12);

% Map of ecoregions and site locations
subplot(3, 3, [1 2 4 5])
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
set(ax, 'Position',[-0.08 0.2 0.82 0.9]);

h = colorbar('eastoutside');
pos = get(h, 'Position');
set(h, 'Position',[0.64 0.42 0.03 0.46],'FontName','Times New Roman',...
    'XTick',1:1:9, 'XTickLabel',ecolabs, 'FontSize',9, 'TickLength',0);

text(-0.314,0.85,'a', 'FontSize',12)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/site-map.tif')
close all;

T = table(unSpc', nSpc', 'VariableNames',{'Code','N'});
writetable(T, './output/TableS1_Species.xlsx');

