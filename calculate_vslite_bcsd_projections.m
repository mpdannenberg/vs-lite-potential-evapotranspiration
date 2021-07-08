% Use down-scaled CMIP5 model ensemble to drive calibrated VS-Lite model
% with different PET formulations to see differences in projected PET, soil
% moisture, soil moisture stress, and growth under warming

load ./output/ITRDB_bcsd;
ITRDB = ITRDB(cellfun(@ischar, {ITRDB.EcoL1_Code}));
syear = 1950;
eyear = 2100;
bcsd_year = syear:eyear;

% Exclude sites with no BCSD data (likely on coasts)
temp = [ITRDB.RCP45];
idx = arrayfun(@(x) sum(~isnan(x.P(:))), temp);
ITRDB = ITRDB(idx > 0);
clear temp idx;

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

%% Get monthly climatology and simulate +2 and +4 C
n = length(ITRDB);
COAST = geotiffread('./data/us_CoastalBoundary_4km.tif');
COAST(COAST<-1000) = NaN;

for j = 1:n
    
    phi = ITRDB(j).LAT;
    elev = ITRDB(j).ELEV;
    
    xy = [ITRDB(j).LAT ITRDB(j).LON];
    DistDeg = distance(xy(1), xy(2), prismLatLon(:,1), prismLatLon(:,2));
    DistKM = distdim(DistDeg, 'deg', 'km');
    xy = prismLatLon(DistKM == min(DistKM), :);
    xind = find(lon == xy(1,2));
    yind = find(lat == xy(1,1));
    
    coast = COAST(yind, xind);
    
    % Calculate adjustment factor for estimating Tdmean from Tmin based on
    % PRISM data.
    Tmin = reshape(squeeze(tmin.tmin(yind, xind, :, :))', [], 1);
    Tdmean = reshape(squeeze(tdmean.tdmean(yind, xind, :, :))', [], 1);
    a = nanmean(Tmin-Tdmean);
    
    % Get 1981-2010 mean precip from observations (this is 
    % controlled in model simulations to isolate just the effect of
    % temperature on PET, M, gM, and growth - gT is also controlled at
    % 1981-2010 levels (in the for loop below because it differs slightly
    % for each PET model)
    P = squeeze(ppt.PPT(yind, xind, :, :))';
    P = repmat(nanmean(P(:, year>=1981 & year<=2010), 2), 1, length(syear:eyear));
    
    % Get RCP4.5 projections
    models = ITRDB(j).RCP45.Models;
    Tmin = ITRDB(j).RCP45.Tmin;
    Tmax = ITRDB(j).RCP45.Tmax;
    if a > 0
        Tdmean = Tmin - a;
    else
        Tdmean = Tmin;
    end
    
    % Loop through each model
    for k = 1:length(models)
        % Thornthwaite
        gT = repmat(ITRDB(j).Th.Tplus0.gTseas', 1, length(syear:eyear));
        [~,~,gM,~,~,M,PET,width] = VSLite_v3(syear, eyear, phi,...
            ITRDB(j).Th.T1, ITRDB(j).Th.T2,...
            ITRDB(j).Th.M1, ITRDB(j).Th.M2,0,0,...
            Tmin(:,:,k), Tmax(:,:,k), Tdmean(:,:,k), P, coast, elev, 'pet_model','Th', 'gT_0',gT);
        ITRDB(j).RCP45.Th.PET_trend(:,k) = sum(PET)';
        ITRDB(j).RCP45.Th.M_trend(:,k) = mean(M)'; ITRDB(j).RCP45.Th.M_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.Th.gM_trend(:,k) = mean(gM)'; ITRDB(j).RCP45.Th.gM_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.Th.TRW(:,k) = width / mean(width(bcsd_year>=1981 & bcsd_year<=2010)); 
        ITRDB(j).RCP45.Th.TRW(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.Th.gM(:,:,k) = gM;
        
        % Hargreaves
        gT = repmat(ITRDB(j).Hg.Tplus0.gTseas', 1, length(syear:eyear));
        [~,~,gM,~,~,M,PET,width] = VSLite_v3(syear, eyear, phi,...
            ITRDB(j).Hg.T1,ITRDB(j).Hg.T2,...
            ITRDB(j).Hg.M1,ITRDB(j).Hg.M2,0,0,...
            Tmin(:,:,k), Tmax(:,:,k), Tdmean(:,:,k), P, coast, elev, 'pet_model','Hg', 'gT_0',gT);
        ITRDB(j).RCP45.Hg.PET_trend(:,k) = sum(PET)';
        ITRDB(j).RCP45.Hg.M_trend(:,k) = mean(M)'; ITRDB(j).RCP45.Hg.M_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.Hg.gM_trend(:,k) = mean(gM)'; ITRDB(j).RCP45.Hg.gM_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.Hg.TRW(:,k) = width / mean(width(bcsd_year>=1981 & bcsd_year<=2010)); 
        ITRDB(j).RCP45.Hg.TRW(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.Hg.gM(:,:,k) = gM;
        
        % Priestley-Taylor
        gT = repmat(ITRDB(j).PT.Tplus0.gTseas', 1, length(syear:eyear));
        [~,~,gM,~,~,M,PET,width] = VSLite_v3(syear, eyear, phi,...
            ITRDB(j).PT.T1,ITRDB(j).PT.T2,...
            ITRDB(j).PT.M1,ITRDB(j).PT.M2,0,0,...
            Tmin(:,:,k), Tmax(:,:,k), Tdmean(:,:,k), P, coast, elev, 'pet_model','PT', 'gT_0',gT);
        ITRDB(j).RCP45.PT.PET_trend(:,k) = sum(PET)';
        ITRDB(j).RCP45.PT.M_trend(:,k) = mean(M)'; ITRDB(j).RCP45.PT.M_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.PT.gM_trend(:,k) = mean(gM)'; ITRDB(j).RCP45.PT.gM_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.PT.TRW(:,k) = width / mean(width(bcsd_year>=1981 & bcsd_year<=2010)); 
        ITRDB(j).RCP45.PT.TRW(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.PT.gM(:,:,k) = gM;
        
        % Penman-Monteith
        gT = repmat(ITRDB(j).PM.Tplus0.gTseas', 1, length(syear:eyear));
        [~,~,gM,~,~,M,PET,width] = VSLite_v3(syear, eyear, phi,...
            ITRDB(j).PM.T1,ITRDB(j).PM.T2,...
            ITRDB(j).PM.M1,ITRDB(j).PM.M2,0,0,...
            Tmin(:,:,k), Tmax(:,:,k), Tdmean(:,:,k), P, coast, elev, 'pet_model','PM', 'gT_0',gT);
        ITRDB(j).RCP45.PM.PET_trend(:,k) = sum(PET)';
        ITRDB(j).RCP45.PM.M_trend(:,k) = mean(M)'; ITRDB(j).RCP45.PM.M_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.PM.gM_trend(:,k) = mean(gM)'; ITRDB(j).RCP45.PM.gM_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.PM.TRW(:,k) = width / mean(width(bcsd_year>=1981 & bcsd_year<=2010));
        ITRDB(j).RCP45.PM.TRW(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP45.PM.gM(:,:,k) = gM;
        
    end
    
    % Get RCP8.5 projections
    models = ITRDB(j).RCP85.Models;
    Tmin = ITRDB(j).RCP85.Tmin;
    Tmax = ITRDB(j).RCP85.Tmax;
    if a > 0
        Tdmean = Tmin - a;
    else
        Tdmean = Tmin;
    end
    
    % Loop through each model
    for k = 1:length(models)
        % Thornthwaite
        gT = repmat(ITRDB(j).Th.Tplus0.gTseas', 1, length(syear:eyear));
        [~,~,gM,~,~,M,PET,width] = VSLite_v3(syear, eyear, phi,...
            ITRDB(j).Th.T1, ITRDB(j).Th.T2,...
            ITRDB(j).Th.M1, ITRDB(j).Th.M2,0,0,...
            Tmin(:,:,k), Tmax(:,:,k), Tdmean(:,:,k), P, coast, elev, 'pet_model','Th', 'gT_0',gT);
        ITRDB(j).RCP85.Th.PET_trend(:,k) = sum(PET)';
        ITRDB(j).RCP85.Th.M_trend(:,k) = mean(M)'; ITRDB(j).RCP85.Th.M_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.Th.gM_trend(:,k) = mean(gM)'; ITRDB(j).RCP85.Th.gM_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.Th.TRW(:,k) = width / mean(width(bcsd_year>=1981 & bcsd_year<=2010)); 
        ITRDB(j).RCP85.Th.TRW(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.Th.gM(:,:,k) = gM;
        test = ~isreal(ITRDB(j).RCP85.Th.TRW(:,k));
        if ~isreal(ITRDB(j).RCP85.Th.TRW(:,k))
            disp([ITRDB(j).SITE, ' - ',models{k}, ' - ', 'Thornthwaite']);
        end
        
        % Hargreaves
        gT = repmat(ITRDB(j).Hg.Tplus0.gTseas', 1, length(syear:eyear));
        [~,~,gM,~,~,M,PET,width] = VSLite_v3(syear, eyear, phi,...
            ITRDB(j).Hg.T1,ITRDB(j).Hg.T2,...
            ITRDB(j).Hg.M1,ITRDB(j).Hg.M2,0,0,...
            Tmin(:,:,k), Tmax(:,:,k), Tdmean(:,:,k), P, coast, elev, 'pet_model','Hg', 'gT_0',gT);
        ITRDB(j).RCP85.Hg.PET_trend(:,k) = sum(PET)';
        ITRDB(j).RCP85.Hg.M_trend(:,k) = mean(M)'; ITRDB(j).RCP85.Hg.M_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.Hg.gM_trend(:,k) = mean(gM)'; ITRDB(j).RCP85.Hg.gM_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.Hg.TRW(:,k) = width / mean(width(bcsd_year>=1981 & bcsd_year<=2010)); 
        ITRDB(j).RCP85.Hg.TRW(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.Hg.gM(:,:,k) = gM;
        if ~isreal(ITRDB(j).RCP85.Hg.TRW(:,k))
            disp([ITRDB(j).SITE, ' - ',models{k}, ' - ', 'Hargreaves']);
        end
        
        % Priestley-Taylor
        gT = repmat(ITRDB(j).PT.Tplus0.gTseas', 1, length(syear:eyear));
        [~,~,gM,~,~,M,PET,width] = VSLite_v3(syear, eyear, phi,...
            ITRDB(j).PT.T1,ITRDB(j).PT.T2,...
            ITRDB(j).PT.M1,ITRDB(j).PT.M2,0,0,...
            Tmin(:,:,k), Tmax(:,:,k), Tdmean(:,:,k), P, coast, elev, 'pet_model','PT', 'gT_0',gT);
        ITRDB(j).RCP85.PT.PET_trend(:,k) = sum(PET)';
        ITRDB(j).RCP85.PT.M_trend(:,k) = mean(M)'; ITRDB(j).RCP85.PT.M_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.PT.gM_trend(:,k) = mean(gM)'; ITRDB(j).RCP85.PT.gM_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.PT.TRW(:,k) = width / mean(width(bcsd_year>=1981 & bcsd_year<=2010)); 
        ITRDB(j).RCP85.PT.TRW(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.PT.gM(:,:,k) = gM;
        if ~isreal(ITRDB(j).RCP85.PT.TRW(:,k))
            disp([ITRDB(j).SITE, ' - ',models{k}, ' - ', 'Priestley-Taylor']);
        end

        % Penman-Monteith
        gT = repmat(ITRDB(j).PM.Tplus0.gTseas', 1, length(syear:eyear));
        [~,~,gM,~,~,M,PET,width] = VSLite_v3(syear, eyear, phi,...
            ITRDB(j).PM.T1,ITRDB(j).PM.T2,...
            ITRDB(j).PM.M1,ITRDB(j).PM.M2,0,0,...
            Tmin(:,:,k), Tmax(:,:,k), Tdmean(:,:,k), P, coast, elev, 'pet_model','PM', 'gT_0',gT);
        ITRDB(j).RCP85.PM.PET_trend(:,k) = sum(PET)';
        ITRDB(j).RCP85.PM.M_trend(:,k) = mean(M)'; ITRDB(j).RCP85.PM.M_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.PM.gM_trend(:,k) = mean(gM)'; ITRDB(j).RCP85.PM.gM_trend(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.PM.TRW(:,k) = width / mean(width(bcsd_year>=1981 & bcsd_year<=2010));
        ITRDB(j).RCP85.PM.TRW(1,k) = NaN; % Exclude first value since it hasn't equilibrated
        ITRDB(j).RCP85.PM.gM(:,:,k) = gM;
        if ~isreal(ITRDB(j).RCP85.PM.TRW(:,k))
            disp([ITRDB(j).SITE, ' - ',models{k}, ' - ', 'Penman-Monteith']);
        end

    end
    
end
clear a coast DistDeg DistKM elev gM j m M ms mw n P PET phi rwi rwi_sim s T Tdmean Tmax Tmin w xind xy yind yr;

save('./output/ITRDB_projections.mat', 'ITRDB', '-v7.3');

