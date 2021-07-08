% Run calibrated model with observed climate

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

ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ITRDB = ITRDB(ecol1 > 0);

%% Get monthly climate, run model, and calculate climatological means

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
    [~,gT,gM,~,~,M,PET,~] = VSLite_v3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','Th');
    ITRDB(i).Th.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Th.Tplus0.gMseas = mean(gM(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).Th.Tplus0.Mseas = mean(M(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).Th.Tplus0.PETseas = mean(PET(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).Th.Tplus0.gTseas = mean(gT(:,year>=1981 & year<=2010), 2)';
    
    % Hargreaves
    [~,gT,gM,~,~,M,PET,~] = VSLite_v3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','Hg');
    ITRDB(i).Hg.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).Hg.Tplus0.gMseas = mean(gM(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).Hg.Tplus0.Mseas = mean(M(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).Hg.Tplus0.PETseas = mean(PET(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).Hg.Tplus0.gTseas = mean(gT(:,year>=1981 & year<=2010), 2)';
    
    % Priestley-Taylor
    [~,gT,gM,~,~,M,PET,~] = VSLite_v3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','PT');
    ITRDB(i).PT.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PT.Tplus0.gMseas = mean(gM(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).PT.Tplus0.Mseas = mean(M(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).PT.Tplus0.PETseas = mean(PET(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).PT.Tplus0.gTseas = mean(gT(:,year>=1981 & year<=2010), 2)';
        
    % Penman-Monteith
    [~,gT,gM,~,~,M,PET,~] = VSLite_v3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','PM');
    ITRDB(i).PM.Tplus0.gM = mean(reshape(gM(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus0.M = mean(reshape(M(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus0.PET = mean(reshape(PET(:,year>=1981 & year<=2010), 1,[]));
    ITRDB(i).PM.Tplus0.gMseas = mean(gM(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).PM.Tplus0.Mseas = mean(M(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).PM.Tplus0.PETseas = mean(PET(:,year>=1981 & year<=2010), 2)';
    ITRDB(i).PM.Tplus0.gTseas = mean(gT(:,year>=1981 & year<=2010), 2)';
    
end
clear coast DistDeg DistKM elev gM i m M ms mw n P PET phi rwi rwi_sim s T Tdmean Tmax Tmin w xind xy yind yr;

save('./data/ITRDB_simulations.mat', 'ITRDB');

