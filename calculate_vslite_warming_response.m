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

%% Get monthly climatology and simulate +2 and +4 C

n = length(ITRDB);
for i = 1:n
    
    phi = ITRDB(i).LAT;
    elev = ITRDB(i).ELEV;
    coast = 0;
    yr = ITRDB(i).YEAR;
    rwi = ITRDB(i).STD;
    
    m = mean(rwi(year>=syear & year<=1960));
    s = std(rwi(year>=syear & year<=1960));
    
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
    
    % Thornthwaite
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin,Tmax,Tdmean,P,0,elev, 'pet_model','Th');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Th.Tplus0.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Th.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin+2,Tmax+2,Tdmean+2,P,0,elev, 'pet_model','Th');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Th.Tplus2.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Th.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin+4,Tmax+4,Tdmean+4,P,0,elev, 'pet_model','Th');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Th.Tplus4.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Th.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    % Hargreaves
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin,Tmax,Tdmean,P,0,elev, 'pet_model','Hg');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Hg.Tplus0.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Hg.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin+2,Tmax+2,Tdmean+2,P,0,elev, 'pet_model','Hg');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Hg.Tplus2.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Hg.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin+4,Tmax+4,Tdmean+4,P,0,elev, 'pet_model','Hg');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).Hg.Tplus4.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).Hg.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    % Priestley-Taylor
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin,Tmax,Tdmean,P,0,elev, 'pet_model','PT');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PT.Tplus0.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PT.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin+2,Tmax+2,Tdmean+2,P,0,elev, 'pet_model','PT');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PT.Tplus2.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PT.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin+4,Tmax+4,Tdmean+4,P,0,elev, 'pet_model','PT');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PT.Tplus4.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PT.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    % Penman-Monteith
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin,Tmax,Tdmean,P,0,elev, 'pet_model','PM');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PM.Tplus0.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PM.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin+2,Tmax+2,Tdmean+2,P,0,elev, 'pet_model','PM');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PM.Tplus2.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PM.Tplus2.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus2.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
    [rwi_sim,~,gM,~,~,M] = VSLite_v2_3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin+4,Tmax+4,Tdmean+4,P,0,elev, 'pet_model','PM');
    rwi_sim = rwi_sim*s + m; % Back to original scale during calibration period
    ITRDB(i).PM.Tplus4.rwi = mean(rwi_sim(year>=1981 & year<=2010));
    ITRDB(i).PM.Tplus4.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus4.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));;
    
end

