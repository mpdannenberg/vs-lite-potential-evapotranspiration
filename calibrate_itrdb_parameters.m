% Calibrate VS-Lite model for multiple PET models

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

%% Load and process ITRDB data
load ./data/ITRDB;
ITRDB = ITRDB_TW;
clear ITRDB_EW ITRDB_LWadj ITRDB_TW;
ITRDB = ITRDB([ITRDB.LAT] > 28 & [ITRDB.LAT] <= 49 & [ITRDB.LON] > -125  & [ITRDB.LON] < -66);
idx = cellfun(@(x) length(regexpi(x, '[a-zA-Z]')), {ITRDB.SITE});
ITRDB = ITRDB(idx <= 3); clear idx;
ITRDB([ITRDB.END]>2019) = [];
ITRDB([ITRDB.END]<1980) = [];
ITRDB([ITRDB.START]>1895) = [];
ITRDB = ITRDB(~isnan([ITRDB.ELEV]));
ITRDB = ITRDB(~isnan([ITRDB.LAT]));
ITRDB = ITRDB(~isnan([ITRDB.LON]));

%% Calibrate parameters - takes a LONG time (even with parallel processing)

n = length(ITRDB);
COAST = geotiffread('./data/us_CoastalBoundary_4km.tif');
COAST(COAST<-1000) = NaN;

for i = 1:n
    
    phi = ITRDB(i).LAT;
    elev = ITRDB(i).ELEV;
    yr = ITRDB(i).YEAR;
    rwi = ITRDB(i).STD;
    
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
    T = (Tmin+Tmax)/2;
    
    if sum(sum(~isnan(T))) > 0
        [~,~,ia] = intersect(cal_yrs, year);
        [~,~,ib] = intersect(cal_yrs, yr);

        % Thornthwaite
        [T1,T2,M1,M2,T1dist,T2dist,M1dist,M2dist] = ...
            estimate_vslite_params_v2_3(Tmin(:,ia),Tmax(:,ia),Tdmean(:,ia),P(:,ia),...
            phi,coast,elev,rwi(ib)','pet_model','Th', 'pt_ests','mle', 'verbose',0,...
            'nsamp',2500, 'nbi',500);

        rwi_sim = VSLite_v2_3(syear, eyear, phi, T1,T2,M1,M2,0,0,...
            Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','Th');
        rwi_sim(1) = NaN;

        ITRDB(i).Th.VSLite = rwi_sim;
        ITRDB(i).Th.T1 = T1;
        ITRDB(i).Th.T2 = T2;
        ITRDB(i).Th.M1 = M1;
        ITRDB(i).Th.M2 = M2;
        ITRDB(i).Th.T1dist = T1dist;
        ITRDB(i).Th.T2dist = T2dist;
        ITRDB(i).Th.M1dist = M1dist;
        ITRDB(i).Th.M2dist = M2dist;

        % Hargreaves
        [T1,T2,M1,M2,T1dist,T2dist,M1dist,M2dist] = ...
            estimate_vslite_params_v2_3(Tmin(:,ia),Tmax(:,ia),Tdmean(:,ia),P(:,ia),...
            phi,coast,elev,rwi(ib)','pet_model','Hg', 'pt_ests','mle', 'verbose',0,...
            'nsamp',2500, 'nbi',500);

        rwi_sim = VSLite_v2_3(syear, eyear, phi, T1,T2,M1,M2,0,0,...
            Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','Hg');
        rwi_sim(1) = NaN;

        ITRDB(i).Hg.VSLite = rwi_sim;
        ITRDB(i).Hg.T1 = T1;
        ITRDB(i).Hg.T2 = T2;
        ITRDB(i).Hg.M1 = M1;
        ITRDB(i).Hg.M2 = M2;
        ITRDB(i).Hg.T1dist = T1dist;
        ITRDB(i).Hg.T2dist = T2dist;
        ITRDB(i).Hg.M1dist = M1dist;
        ITRDB(i).Hg.M2dist = M2dist;

        % Priestley-Taylor
        [T1,T2,M1,M2,T1dist,T2dist,M1dist,M2dist] = ...
            estimate_vslite_params_v2_3(Tmin(:,ia),Tmax(:,ia),Tdmean(:,ia),P(:,ia),...
            phi,coast,elev,rwi(ib)','pet_model','PT', 'pt_ests','mle', 'verbose',0,...
            'nsamp',2500, 'nbi',500);

        rwi_sim = VSLite_v2_3(syear, eyear, phi, T1,T2,M1,M2,0,0,...
            Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','PT');
        rwi_sim(1) = NaN;

        ITRDB(i).PT.VSLite = rwi_sim;
        ITRDB(i).PT.T1 = T1;
        ITRDB(i).PT.T2 = T2;
        ITRDB(i).PT.M1 = M1;
        ITRDB(i).PT.M2 = M2;
        ITRDB(i).PT.T1dist = T1dist;
        ITRDB(i).PT.T2dist = T2dist;
        ITRDB(i).PT.M1dist = M1dist;
        ITRDB(i).PT.M2dist = M2dist;

        % Penman-Monteith
        [T1,T2,M1,M2,T1dist,T2dist,M1dist,M2dist] = ...
            estimate_vslite_params_v2_3(Tmin(:,ia),Tmax(:,ia),Tdmean(:,ia),P(:,ia),...
            phi,coast,elev,rwi(ib)','pet_model','PM', 'pt_ests','mle', 'verbose',0,...
            'nsamp',2500, 'nbi',500);

        rwi_sim = VSLite_v2_3(syear, eyear, phi, T1,T2,M1,M2,0,0,...
            Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','PM');
        rwi_sim(1) = NaN;

        ITRDB(i).PM.VSLite = rwi_sim;
        ITRDB(i).PM.T1 = T1;
        ITRDB(i).PM.T2 = T2;
        ITRDB(i).PM.M1 = M1;
        ITRDB(i).PM.M2 = M2;
        ITRDB(i).PM.T1dist = T1dist;
        ITRDB(i).PM.T2dist = T2dist;
        ITRDB(i).PM.M1dist = M1dist;
        ITRDB(i).PM.M2dist = M2dist;
    end
end

idx = ~cellfun(@isempty, {ITRDB.Th});
ITRDB = ITRDB(idx);

save('./data/ITRDB_vslite.mat','ITRDB','syear','eyear','cal_yrs');


