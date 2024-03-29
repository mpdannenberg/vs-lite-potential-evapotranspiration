% Simulate asymmetric diurnal warming

load ./data/ITRDB_vslite;

%% Read climate data
ppt = matfile('D:\Data_Analysis\PRISM\PRISM_PPT');
tmax = matfile('D:\Data_Analysis\PRISM\PRISM_TMAX');
tmin = matfile('D:\Data_Analysis\PRISM\PRISM_TMIN');
tdmean = matfile('D:\Data_Analysis\PRISM\PRISM_TDMEAN.mat');
lat = ppt.lat;
lon = ppt.lon;
year = ppt.year;
nx = length(lon);
ny = length(lat);
prismLatLon = [reshape(repmat(lat', 1, nx), [], 1) reshape(repmat(lon, ny, 1), [], 1)];
syear = min(year);
eyear = max(year);

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

%% Get monthly climatology and simulate +2 and +4 C

n = length(ITRDB);
COAST = geotiffread('./data/us_CoastalBoundary_4km.tif');
COAST(COAST<-1000) = NaN;

for i = 1:n
    
    phi = ITRDB(i).LAT;
    elev = ITRDB(i).ELEV;
    
    xy = [ITRDB(i).LAT ITRDB(i).LON];
    DistDeg = distance(xy(1), xy(2), prismLatLon(:,1), prismLatLon(:,2));
    DistKM = distdim(DistDeg, 'deg', 'km');
    xy = prismLatLon(DistKM == min(DistKM), :);
    xind = find(lon == xy(1,2));
    yind = find(lat == xy(1,1));
    
    coast = COAST(yind, xind);
    P = squeeze(ppt.PPT(yind, xind, :, :))';
    Tmin = squeeze(tmin.tmin(yind, xind, :, :))';
    Tmax = squeeze(tmax.tmax(yind, xind, :, :))';
    Tdmean = squeeze(tdmean.tdmean(yind, xind, :, :))';
    
    % Thornthwaite
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','Th');
    ITRDB(i).Th.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin+1,Tmax+3,Tdmean+1,P,coast,elev, 'pet_model','Th');
    ITRDB(i).Th.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus2.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin+3,Tmax+1,Tdmean+3,P,coast,elev, 'pet_model','Th');
    ITRDB(i).Th.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus4.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    % Hargreaves
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','Hg');
    ITRDB(i).Hg.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin+1,Tmax+3,Tdmean+1,P,coast,elev, 'pet_model','Hg');
    ITRDB(i).Hg.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus2.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin+3,Tmax+1,Tdmean+3,P,coast,elev, 'pet_model','Hg');
    ITRDB(i).Hg.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus4.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    % Priestley-Taylor
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','PT');
    ITRDB(i).PT.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin+1,Tmax+3,Tdmean+1,P,coast,elev, 'pet_model','PT');
    ITRDB(i).PT.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus2.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin+3,Tmax+1,Tdmean+3,P,coast,elev, 'pet_model','PT');
    ITRDB(i).PT.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus4.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    % Penman-Monteith
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','PM');
    ITRDB(i).PM.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin+1,Tmax+3,Tdmean+1,P,coast,elev, 'pet_model','PM');
    ITRDB(i).PM.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus2.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,~] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin+3,Tmax+1,Tdmean+3,P,coast,elev, 'pet_model','PM');
    ITRDB(i).PM.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus4.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
end
clear coast DistDeg DistKM elev gM i ia ib m M ms mw n P PET phi rwi rwi_sim s T Tdmean Tmax Tmin w xind xy yind yr;

%% Calculate mean/SE of PET, soil moisture, and gM in each ecoregion
clr = wesanderson('fantasticfox1');
idx = ~cellfun(@isempty, {ITRDB.EcoL1_Code});

ITRDB = ITRDB(idx);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

pet_T2_m = NaN(length(ecos), 4);
pet_T2_s = NaN(length(ecos), 4);
M_T2_m = NaN(length(ecos), 4);
M_T2_s = NaN(length(ecos), 4);
gM_T2_m = NaN(length(ecos), 4);
gM_T2_s = NaN(length(ecos), 4);
pet_T4_m = NaN(length(ecos), 4);
pet_T4_s = NaN(length(ecos), 4);
M_T4_m = NaN(length(ecos), 4);
M_T4_s = NaN(length(ecos), 4);
gM_T4_m = NaN(length(ecos), 4);
gM_T4_s = NaN(length(ecos), 4);
n = NaN(length(ecos), 1);

for i = 1:length(ecos)
    
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    n(i) = length(ITRDB_sub);
    
    % Thornthwaite
    model = [ITRDB_sub.Th];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    T4 = [model.Tplus4];
    pet_T2_m(i, 1) = mean(12*[T2.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T2_s(i, 1) = std(12*[T2.PET] - 12*[T0.PET]) / sqrt(length(model));
    M_T2_m(i, 1) = mean([T2.M] - [T0.M]);
    M_T2_s(i, 1) = std([T2.M] - [T0.M]) / sqrt(length(model));
    gM_T2_m(i, 1) = mean([T2.gM] - [T0.gM]);
    gM_T2_s(i, 1) = std([T2.gM] - [T0.gM]) / sqrt(length(model));
    pet_T4_m(i, 1) = mean(12*[T4.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T4_s(i, 1) = std(12*[T4.PET] - 12*[T0.PET]) / sqrt(length(model));
    M_T4_m(i, 1) = mean([T4.M] - [T0.M]);
    M_T4_s(i, 1) = std([T4.M] - [T0.M]) / sqrt(length(model));
    gM_T4_m(i, 1) = mean([T4.gM] - [T0.gM]);
    gM_T4_s(i, 1) = std([T4.gM] - [T0.gM]) / sqrt(length(model));
    
    % Hargreaves
    model = [ITRDB_sub.Hg];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    T4 = [model.Tplus4];
    pet_T2_m(i, 2) = mean(12*[T2.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T2_s(i, 2) = std(12*[T2.PET] - 12*[T0.PET]) / sqrt(length(model));
    M_T2_m(i, 2) = mean([T2.M] - [T0.M]);
    M_T2_s(i, 2) = std([T2.M] - [T0.M]) / sqrt(length(model));
    gM_T2_m(i, 2) = mean([T2.gM] - [T0.gM]);
    gM_T2_s(i, 2) = std([T2.gM] - [T0.gM]) / sqrt(length(model));
    pet_T4_m(i, 2) = mean(12*[T4.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T4_s(i, 2) = std(12*[T4.PET] - 12*[T0.PET]) / sqrt(length(model));
    M_T4_m(i, 2) = mean([T4.M] - [T0.M]);
    M_T4_s(i, 2) = std([T4.M] - [T0.M]) / sqrt(length(model));
    gM_T4_m(i, 2) = mean([T4.gM] - [T0.gM]);
    gM_T4_s(i, 2) = std([T4.gM] - [T0.gM]) / sqrt(length(model));
    
    % Priestly-Taylor
    model = [ITRDB_sub.PT];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    T4 = [model.Tplus4];
    pet_T2_m(i, 3) = mean(12*[T2.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T2_s(i, 3) = std(12*[T2.PET] - 12*[T0.PET]) / sqrt(length(model));
    M_T2_m(i, 3) = mean([T2.M] - [T0.M]);
    M_T2_s(i, 3) = std([T2.M] - [T0.M]) / sqrt(length(model));
    gM_T2_m(i, 3) = mean([T2.gM] - [T0.gM]);
    gM_T2_s(i, 3) = std([T2.gM] - [T0.gM]) / sqrt(length(model));
    pet_T4_m(i, 3) = mean(12*[T4.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T4_s(i, 3) = std(12*[T4.PET] - 12*[T0.PET]) / sqrt(length(model));
    M_T4_m(i, 3) = mean([T4.M] - [T0.M]);
    M_T4_s(i, 3) = std([T4.M] - [T0.M]) / sqrt(length(model));
    gM_T4_m(i, 3) = mean([T4.gM] - [T0.gM]);
    gM_T4_s(i, 3) = std([T4.gM] - [T0.gM]) / sqrt(length(model));
    
    % Penman-Monteith
    model = [ITRDB_sub.PM];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    T4 = [model.Tplus4];
    pet_T2_m(i, 4) = mean(12*[T2.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T2_s(i, 4) = std(12*[T2.PET] - 12*[T0.PET]) / sqrt(length(model));
    M_T2_m(i, 4) = mean([T2.M] - [T0.M]);
    M_T2_s(i, 4) = std([T2.M] - [T0.M]) / sqrt(length(model));
    gM_T2_m(i, 4) = mean([T2.gM] - [T0.gM]);
    gM_T2_s(i, 4) = std([T2.gM] - [T0.gM]) / sqrt(length(model));
    pet_T4_m(i, 4) = mean(12*[T4.PET] - 12*[T0.PET]); % Mean monthly --> mean annual
    pet_T4_s(i, 4) = std(12*[T4.PET] - 12*[T0.PET]) / sqrt(length(model));
    M_T4_m(i, 4) = mean([T4.M] - [T0.M]);
    M_T4_s(i, 4) = std([T4.M] - [T0.M]) / sqrt(length(model));
    gM_T4_m(i, 4) = mean([T4.gM] - [T0.gM]);
    gM_T4_s(i, 4) = std([T4.gM] - [T0.gM]) / sqrt(length(model));
    
end

%% Bargraph of higher daytime warming effect
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 5.5];

axes('Position',[0.1,0.7,0.4,0.25])
b = bar(ecos,pet_T2_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05],'XTickLabels',...
    '','YLim',[-50 200])
text(2.75, 75, '\DeltaPET (mm)', 'FontSize',11, 'Rotation',90, 'HorizontalAlignment','center');
ylim = get(gca, 'Ylim');
text(5,ylim(2),'a', 'FontSize',12, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
text(9,ylim(2),'T_{max} > T_{min}', 'FontSize',12, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');

axes('Position',[0.1,0.4,0.4,0.25])
b = bar(ecos,M_T2_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05],'XTickLabels',...
    '', 'XAxisLocation','bottom', 'YDir','reverse', 'YLim',[-0.025 0.005])
text(2.75, -0.01, '\DeltaSoil moisture (vol/vol)', 'FontSize',11, 'Rotation',90, 'HorizontalAlignment','center');
ylim = get(gca, 'Ylim');
text(5,ylim(1),'b', 'FontSize',12, 'HorizontalAlignment','center', 'VerticalAlignment','middle');

axes('Position',[0.1,0.1,0.4,0.25])
b = bar(ecos,gM_T2_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05], 'XAxisLocation','bottom',...
    'XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'},...
    'YDir','reverse','YLim',[-0.06 0.01])
text(2.75, -0.025, '\Deltag_{M}', 'FontSize',11, 'Rotation',90, 'HorizontalAlignment','center');
ylim = get(gca, 'Ylim');
text(5,ylim(1),'c', 'FontSize',12, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
text(9,0.03,'Ecoregion', 'FontSize',11, 'HorizontalAlignment','center', 'VerticalAlignment','middle');

%% Bargraph of higher nighttime warming effect
axes('Position',[0.57,0.7,0.4,0.25])
b = bar(ecos,pet_T4_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05],'XTickLabels',...
    '','YLim',[-50 200])
ylim = get(gca, 'Ylim');
text(5,ylim(2),'d', 'FontSize',12, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
text(9,ylim(2),'T_{min} > T_{max}', 'FontSize',12, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');

axes('Position',[0.57,0.4,0.4,0.25])
b = bar(ecos,M_T4_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05],'XTickLabels',...
    '', 'XAxisLocation','bottom', 'YDir','reverse', 'YLim',[-0.025 0.005])
ylim = get(gca, 'Ylim');
text(5,ylim(1),'e', 'FontSize',12, 'HorizontalAlignment','center', 'VerticalAlignment','middle');

axes('Position',[0.57,0.1,0.4,0.25])
b = bar(ecos,gM_T4_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05], 'XAxisLocation','bottom',...
    'XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'},...
    'YDir','reverse','YLim',[-0.06 0.01])
ylim = get(gca, 'Ylim');
text(5,ylim(1),'f', 'FontSize',12, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
lgd = legend('Thornthwaite','Hargreaves','Priestly-Taylor','Penman-Monteith');
lgd.Position = [0.73 0.26 0.2854 0.1068];
lgd.FontSize = 7;
legend('boxoff');
text(9,0.03,'Ecoregion', 'FontSize',11, 'HorizontalAlignment','center', 'VerticalAlignment','middle');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/vslite-asymmetric-diurnal-bars.tif')
close all;

