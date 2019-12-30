% Simulate +2 and +4 C

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

cal_yrs = syear:1960;

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
for i = 1:n
    
    phi = ITRDB(i).LAT;
    elev = ITRDB(i).ELEV;
    coast = 0;
    yr = ITRDB(i).YEAR;
    rwi = ITRDB(i).STD;
    
    m = mean(rwi(yr>=syear & yr<=1960));
    s = std(rwi(yr>=syear & yr<=1960));
    
    xy = [ITRDB(i).LAT ITRDB(i).LON];
    DistDeg = distance(xy(1), xy(2), prismLatLon(:,1), prismLatLon(:,2));
    DistKM = distdim(DistDeg, 'deg', 'km');
    xy = prismLatLon(DistKM == min(DistKM), :);
    xind = find(lon == xy(1,2));
    yind = find(lat == xy(1,1));
    
    P = squeeze(ppt.PPT(yind, xind, :, :))';
    Tmin = squeeze(tmin.tmin(yind, xind, :, :))';
    Tmax = squeeze(tmax.tmax(yind, xind, :, :))';
    Tdmean = squeeze(tdmean.tdmean(yind, xind, :, :))';
    T = (Tmin+Tmax)/2;
    
    [~,ia,ib] = intersect(cal_yrs,year);
    
    % Thornthwaite
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin,Tmax,Tdmean,P,0,elev, 'pet_model','Th');
    w(1) = NaN;
    mw = nanmean(w(ia));
    ms = nanstd(w(ia));
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Th.Tplus0.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Th.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin+2,Tmax+2,Tdmean+2,P,0,elev, 'pet_model','Th');
    w(1) = NaN;
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Th.Tplus2.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Th.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus2.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin+4,Tmax+4,Tdmean+4,P,0,elev, 'pet_model','Th');
    w(1) = NaN;
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Th.Tplus4.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Th.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus4.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    % Hargreaves
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin,Tmax,Tdmean,P,0,elev, 'pet_model','Hg');
    w(1) = NaN;
    mw = nanmean(w(ia));
    ms = nanstd(w(ia));
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Hg.Tplus0.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Hg.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin+2,Tmax+2,Tdmean+2,P,0,elev, 'pet_model','Hg');
    w(1) = NaN;
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Hg.Tplus2.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Hg.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus2.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin+4,Tmax+4,Tdmean+4,P,0,elev, 'pet_model','Hg');
    w(1) = NaN;
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Hg.Tplus4.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Hg.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus4.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    % Priestley-Taylor
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin,Tmax,Tdmean,P,0,elev, 'pet_model','PT');
    w(1) = NaN;
    mw = nanmean(w(ia));
    ms = nanstd(w(ia));
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PT.Tplus0.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PT.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin+2,Tmax+2,Tdmean+2,P,0,elev, 'pet_model','PT');
    w(1) = NaN;
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PT.Tplus2.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PT.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus2.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin+4,Tmax+4,Tdmean+4,P,0,elev, 'pet_model','PT');
    w(1) = NaN;
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PT.Tplus4.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PT.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus4.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    % Penman-Monteith
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin,Tmax,Tdmean,P,0,elev, 'pet_model','PM');
    w(1) = NaN;
    mw = nanmean(w(ia));
    ms = nanstd(w(ia));
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PM.Tplus0.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PM.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin+2,Tmax+2,Tdmean+2,P,0,elev, 'pet_model','PM');
    w(1) = NaN;
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PM.Tplus2.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PM.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus2.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    
    [~,~,gM,~,~,M,PET,w] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin+4,Tmax+4,Tdmean+4,P,0,elev, 'pet_model','PM');
    w(1) = NaN;
    rwi_sim = (w-mw)/ms; % Scale series relative to ambient calibration period
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PM.Tplus4.rwi = mean(rwi_sim(year>=1981 & year<=2010));
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

%% Bargraph of +2C effect
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 6];

subplot(3,1,1)
b = bar(ecos,pet_T2_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05],'XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'})
ylabel('\DeltaPET (mm)', 'FontSize',11);
ylim = get(gca, 'Ylim');
text(5,ylim(2),'a', 'FontSize',12, 'HorizontalAlignment','right', 'VerticalAlignment','top');
lgd = legend('Thornthwaite','Hargreaves','Priestly-Taylor','Penman-Monteith');
lgd.Position = [0.5965 0.89 0.2854 0.1068];
lgd.FontSize = 7;
legend('boxoff');

subplot(3,1,2)
b = bar(ecos,M_T2_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05],'XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'}, 'XAxisLocation','bottom', 'YDir','reverse')
ylabel('\DeltaSoil moisture (vol/vol)', 'FontSize',11);
ylim = get(gca, 'Ylim');
text(5,ylim(1),'b', 'FontSize',12, 'HorizontalAlignment','right', 'VerticalAlignment','top');

subplot(3,1,3)
b = bar(ecos,gM_T2_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05], 'XAxisLocation','bottom','XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'}, 'YDir','reverse')
ylabel('\Deltag_{M}', 'FontSize',11);
ylim = get(gca, 'Ylim');
text(5,ylim(1),'c', 'FontSize',12, 'HorizontalAlignment','right', 'VerticalAlignment','top');
xlabel('Ecoregion', 'FontSize',11);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/vslite-t+2-bars.tif')
close all;

%% Bargraph of +4C effect
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 6];

subplot(3,1,1)
b = bar(ecos,pet_T4_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05],'XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'})
ylabel('\DeltaPET (mm)', 'FontSize',11);
ylim = get(gca, 'Ylim');
text(5,ylim(2),'a', 'FontSize',12, 'HorizontalAlignment','right', 'VerticalAlignment','top');
lgd = legend('Thornthwaite','Hargreaves','Priestly-Taylor','Penman-Monteith');
lgd.Position = [0.5965 0.89 0.2854 0.1068];
lgd.FontSize = 7;
legend('boxoff');

subplot(3,1,2)
b = bar(ecos,M_T4_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05],'XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'}, 'XAxisLocation','bottom', 'YDir','reverse')
ylabel('\DeltaSoil moisture (vol/vol)', 'FontSize',11);
ylim = get(gca, 'Ylim');
text(5,ylim(1),'b', 'FontSize',12, 'HorizontalAlignment','right', 'VerticalAlignment','top');

subplot(3,1,3)
b = bar(ecos,gM_T4_m);
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
b(4).FaceColor = clr(4,:);
box off;
set(gca, 'TickDir','out','TickLength',[0.02 0.05], 'XAxisLocation','bottom','XTickLabels',{'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'}, 'YDir','reverse')
ylabel('\Deltag_{M}', 'FontSize',11);
ylim = get(gca, 'Ylim');
text(5,ylim(1),'c', 'FontSize',12, 'HorizontalAlignment','right', 'VerticalAlignment','top');
xlabel('Ecoregion', 'FontSize',11);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/vslite-t+4-bars.tif')
close all;

