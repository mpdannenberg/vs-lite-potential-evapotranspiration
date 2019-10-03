% Plot distributions of validation metrics

load ./data/ITRDB_vslite.mat;

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

%% Make figure
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 4];

clr = wesanderson('fantasticfox1');

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.r2_val];
[f,x] = histcounts(r2, 0:0.05:1);
plot(0.025:0.05:1, f, '-', 'Color',clr(1,:), 'LineWidth',2)
hold on;

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.r2_val];
[f,x] = histcounts(r2, 0:0.05:1);
plot(0.025:0.05:1, f, '-', 'Color',clr(2,:), 'LineWidth',2)

% Priestly-Taylor
r2 = [ITRDB.PT];
r2 = [r2.r2_val];
[f,x] = histcounts(r2, 0:0.05:1);
plot(0.025:0.05:1, f, '-', 'Color',clr(3,:), 'LineWidth',2)

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.r2_val];
[f,x] = histcounts(r2, 0:0.05:1);
plot(0.025:0.05:1, f, '-', 'Color',clr(4,:), 'LineWidth',2)

