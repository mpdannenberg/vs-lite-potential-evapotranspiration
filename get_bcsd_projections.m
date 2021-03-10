% Get downscaled RCP4.5 and RCP8.5 projections

load ./data/ITRDB_simulations;
n = length(ITRDB);

cd('D:\Data_Analysis\BCSD_DownscaledCMIP5')
models = {'ACCESS1-0','bcc-csm1-1','bcc-csm1-1-m','CanESM2',...
    'CCSM4','CESM1-BGC','CESM1-CAM5','CMCC-CM','CNRM-CM5',...
    'CSIRO-Mk3-6-0','FGOALS-g2','FIO-ESM','GFDL-CM3','GFDL-ESM2G',...
    'GFDL-ESM2M','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES',...
    'inmcm4','MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM-LR',...
    'MPI-ESM-MR','MRI-CGCM3','NorESM1-M'}; % Subject to change
% ACCESS1-3 does not have tasmin or tasmax
% IPSL models project Tmin > Tmax (nonsense)
lat = double(ncread(char(glob(['*_pr_Amon_',models{1},'_historical*'])), 'latitude'));
lon = double(ncread(char(glob(['*_pr_Amon_',models{1},'_historical*'])), 'longitude'))-360;
days_in_month = [31 28 31 30 31 30 31 31 30 31 30 31];
t0 = datenum(1950, 1, 1);
nx = length(lon);
ny = length(lat);
LatLon = [reshape(repmat(lat', 1, nx), [], 1) reshape(repmat(lon, ny, 1), [], 1)];

for i = 1:length(models)
    
    % Time vector
    t1 = ncread(char(glob(['*_pr_Amon_',models{i},'_historical*'])), 'time');
    t2 = ncread(char(glob(['*_pr_Amon_',models{i},'_rcp45*'])), 'time');
    [yr45, mo45] = datevec(t0 + [t1;t2]);
    t2 = ncread(char(glob(['*_pr_Amon_',models{i},'_rcp85*'])), 'time');
    [yr85, mo85] = datevec(t0 + [t1;t2]);
    
    % Precip
    hstrcl = ncread(char(glob(['*_pr_Amon_',models{i},'_historical*'])), 'pr');
    rcp45 = ncread(char(glob(['*_pr_Amon_',models{i},'_rcp45*'])), 'pr');
    rcp85 = ncread(char(glob(['*_pr_Amon_',models{i},'_rcp85*'])), 'pr');
    
    pr_rcp45 = cat(3, hstrcl, rcp45);
    pr_rcp85 = cat(3, hstrcl, rcp85);
    clear hstrcl rcp45 rcp85;
    
    pr_rcp45(pr_rcp45 == 1e+20) = NaN;
    pr_rcp45 = permute(pr_rcp45, [2 1 3]); % Switch lat & lon
    pr_rcp45 = pr_rcp45 .* permute(repmat(days_in_month', length(yr45)/12, size(pr_rcp45, 1), size(pr_rcp45, 2)), [2 3 1]);
    
    pr_rcp85(pr_rcp85 == 1e+20) = NaN;
    pr_rcp85 = permute(pr_rcp85, [2 1 3]); % Switch lat & lon
    pr_rcp85 = pr_rcp85 .* permute(repmat(days_in_month', length(yr85)/12, size(pr_rcp85, 1), size(pr_rcp85, 2)), [2 3 1]);
    
    % Tmin
    hstrcl = ncread(char(glob(['*_tasmin_Amon_',models{i},'_historical*'])), 'tasmin');
    rcp45 = ncread(char(glob(['*_tasmin_Amon_',models{i},'_rcp45*'])), 'tasmin');
    rcp85 = ncread(char(glob(['*_tasmin_Amon_',models{i},'_rcp85*'])), 'tasmin');
    
    tasmin_rcp45 = cat(3, hstrcl, rcp45);
    tasmin_rcp85 = cat(3, hstrcl, rcp85);
    clear hstrcl rcp45 rcp85;
    
    tasmin_rcp45(tasmin_rcp45 == 1e+20) = NaN;
    tasmin_rcp45 = permute(tasmin_rcp45, [2 1 3]); % Switch lat & lon
    
    tasmin_rcp85(tasmin_rcp85 == 1e+20) = NaN;
    tasmin_rcp85 = permute(tasmin_rcp85, [2 1 3]); % Switch lat & lon
    
    % Tmax
    hstrcl = ncread(char(glob(['*_tasmax_Amon_',models{i},'_historical*'])), 'tasmax');
    rcp45 = ncread(char(glob(['*_tasmax_Amon_',models{i},'_rcp45*'])), 'tasmax');
    rcp85 = ncread(char(glob(['*_tasmax_Amon_',models{i},'_rcp85*'])), 'tasmax');
    
    tasmax_rcp45 = cat(3, hstrcl, rcp45);
    tasmax_rcp85 = cat(3, hstrcl, rcp85);
    clear hstrcl rcp45 rcp85;
    
    tasmax_rcp45(tasmax_rcp45 == 1e+20) = NaN;
    tasmax_rcp45 = permute(tasmax_rcp45, [2 1 3]); % Switch lat & lon
    
    tasmax_rcp85(tasmax_rcp85 == 1e+20) = NaN;
    tasmax_rcp85 = permute(tasmax_rcp85, [2 1 3]); % Switch lat & lon
    
    
    for j = 1:n
        
        % Initialize projection matrix
        if i==1
            ITRDB(j).RCP45.Year = 1950:2100;
            ITRDB(j).RCP45.Models = models;
            ITRDB(j).RCP45.P = NaN(12, length(1950:2100), length(models));
            ITRDB(j).RCP45.Tmin = NaN(12, length(1950:2100), length(models));
            ITRDB(j).RCP45.Tmax = NaN(12, length(1950:2100), length(models));
            ITRDB(j).RCP85.Year = 1950:2100;
            ITRDB(j).RCP85.Models = models;
            ITRDB(j).RCP85.P = NaN(12, length(1950:2100), length(models));
            ITRDB(j).RCP85.Tmin = NaN(12, length(1950:2100), length(models));
            ITRDB(j).RCP85.Tmax = NaN(12, length(1950:2100), length(models));
        
        end
        
        xy = [ITRDB(j).LAT ITRDB(j).LON];
        DistDeg = distance(xy(1), xy(2), LatLon(:,1), LatLon(:,2));
        DistKM = distdim(DistDeg, 'deg', 'km');
        xy = LatLon(DistKM == min(DistKM), :);
        xind = find(lon == xy(1,2));
        yind = find(lat == xy(1,1));
        
        % RCP4.5
        tmin = reshape(squeeze(tasmin_rcp45(yind, xind, :)), 12, []);
        tmax = reshape(squeeze(tasmax_rcp45(yind, xind, :)), 12, []);
        p = reshape(squeeze(pr_rcp45(yind, xind, :)), 12, []);
        [~,ia,ib] = intersect(ITRDB(j).RCP45.Year, sort(unique(yr45)));
        ITRDB(j).RCP45.Tmin(:, ia, i) = tmin(:, ib);
        ITRDB(j).RCP45.Tmax(:, ia, i) = tmax(:, ib);
        ITRDB(j).RCP45.P(:, ia, i) = p(:, ib);
        
        % RCP8.5
        tmin = reshape(squeeze(tasmin_rcp85(yind, xind, :)), 12, []);
        tmax = reshape(squeeze(tasmax_rcp85(yind, xind, :)), 12, []);
        p = reshape(squeeze(pr_rcp85(yind, xind, :)), 12, []);
        [~,ia,ib] = intersect(ITRDB(j).RCP85.Year, sort(unique(yr85)));
        ITRDB(j).RCP85.Tmin(:, ia, i) = tmin(:, ib);
        ITRDB(j).RCP85.Tmax(:, ia, i) = tmax(:, ib);
        ITRDB(j).RCP85.P(:, ia, i) = p(:, ib);
        
    end
    
    
end

cd('D:\Publications\Dannenberg_VS-Lite_PET');
save('./output/ITRDB_bcsd.mat','ITRDB', '-v7.3')
