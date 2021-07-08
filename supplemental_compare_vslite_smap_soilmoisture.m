% Compare simulate soil moisture to SMAP

%% load tree-ring data
load ./data/ITRDB_simulations;
ITRDB = ITRDB(cellfun(@ischar, {ITRDB.EcoL1_Code}));
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ITRDB = ITRDB(ecol1 > 0);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

n = length(ITRDB);

%% load SMAP data
load ./data/SMAP_L4_SM_monthly.mat;
clear Tsoil_monthly SurfSM_monthly;
slat = lat'; slon = lon'; syr = yr; smo = mo; % add "s" prefix to distinguish from PRISM
smapLatLon = [reshape(repmat(slat', 1, length(slon)), [], 1) reshape(repmat(slon, length(slat), 1), [], 1)];
clear lat lon yr mo;

%% load PRISM data
ppt = matfile('D:\Data_Analysis\PRISM\PRISM_PPT');
tmax = matfile('D:\Data_Analysis\PRISM\PRISM_TMAX');
tmin = matfile('D:\Data_Analysis\PRISM\PRISM_TMIN');
tdmean = matfile('D:\Data_Analysis\PRISM\PRISM_TDMEAN.mat');
plat = ppt.lat;
plon = ppt.lon;
pyr = reshape(repmat(ppt.year, 12, 1),[],1);
pmo = reshape(repmat([1:12]',1,length(ppt.year)),[],1);
prismLatLon = [reshape(repmat(plat', 1, length(plon)), [], 1) reshape(repmat(plon, length(plat), 1), [], 1)];
syear = min(pyr);
eyear = max(pyr);

%% loop through each site, simulate soil moisture with leaky bucket model, and compare to SMAP soil moisture
COAST = geotiffread('./data/us_CoastalBoundary_4km.tif');
COAST(COAST<-1000) = NaN;

for i=1:n
    
    % Get info
    phi = ITRDB(i).LAT;
    elev = ITRDB(i).ELEV;
    
    % Find nearest PRISM cell
    xy = [ITRDB(i).LAT ITRDB(i).LON];
    DistDeg = distance(xy(1), xy(2), prismLatLon(:,1), prismLatLon(:,2));
    DistKM = distdim(DistDeg, 'deg', 'km');
    xy = prismLatLon(DistKM == min(DistKM), :);
    xind = find(plon == xy(1,2));
    yind = find(plat == xy(1,1));
    
    % Get PRISM data
    coast = COAST(yind, xind);
    P = squeeze(ppt.PPT(yind, xind, :, :))';
    Tmin = squeeze(tmin.Tmin(yind, xind, :, :))';
    Tmax = squeeze(tmax.Tmax(yind, xind, :, :))';
    Tdmean = squeeze(tdmean.Tdmean(yind, xind, :, :))';
    
    % Find nearest SMAP cell
    xy = [ITRDB(i).LAT ITRDB(i).LON];
    DistDeg = distance(xy(1), xy(2), smapLatLon(:,1), smapLatLon(:,2));
    DistKM = distdim(DistDeg, 'deg', 'km');
    xy = smapLatLon(DistKM == min(DistKM), :);
    xind = find(slon == xy(1,2));
    yind = find(slat == xy(1,1));
    smap = squeeze(RootSM_monthly(yind, xind, :));
    
    % Simulate soil moisture
    % Thornthwaite
    [~,~,~,~,~,M,PET,~] = VSLite_v3(syear, eyear, phi,...
        ITRDB(i).Th.T1,ITRDB(i).Th.T2,...
        ITRDB(i).Th.M1,ITRDB(i).Th.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','Th');
    ITRDB(i).Th.M = M;
    ITRDB(i).Th.PET = PET;
    M = reshape(M,[],1);
    [r,p] = corr(M(find(pyr==syr(1) & pmo==smo(1)):find(pyr==syr(end) & pmo==smo(end))), smap, 'rows','pairwise');
    ITRDB(i).Th.R_SMAP = r;
    ITRDB(i).Th.p_SMAP = p;
    
    % Hargreaves
    [~,~,~,~,~,M,PET,~] = VSLite_v3(syear, eyear, phi,...
        ITRDB(i).Hg.T1,ITRDB(i).Hg.T2,...
        ITRDB(i).Hg.M1,ITRDB(i).Hg.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','Hg');
    ITRDB(i).Hg.M = M;
    ITRDB(i).Hg.PET = PET;
    M = reshape(M,[],1);
    [r,p] = corr(M(find(pyr==syr(1) & pmo==smo(1)):find(pyr==syr(end) & pmo==smo(end))), smap, 'rows','pairwise');
    ITRDB(i).Hg.R_SMAP = r;
    ITRDB(i).Hg.p_SMAP = p;
    
    % Priestley-Taylor
    [~,~,~,~,~,M,PET,~] = VSLite_v3(syear, eyear, phi,...
        ITRDB(i).PT.T1,ITRDB(i).PT.T2,...
        ITRDB(i).PT.M1,ITRDB(i).PT.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','PT');
    ITRDB(i).PT.M = M;
    ITRDB(i).PT.PET = PET;
    M = reshape(M,[],1);
    [r,p] = corr(M(find(pyr==syr(1) & pmo==smo(1)):find(pyr==syr(end) & pmo==smo(end))), smap, 'rows','pairwise');
    ITRDB(i).PT.R_SMAP = r;
    ITRDB(i).PT.p_SMAP = p;
        
    % Penman-Monteith
    [~,~,~,~,~,M,PET,~] = VSLite_v3(syear, eyear, phi,...
        ITRDB(i).PM.T1,ITRDB(i).PM.T2,...
        ITRDB(i).PM.M1,ITRDB(i).PM.M2,0,0,...
        Tmin,Tmax,Tdmean,P,coast,elev, 'pet_model','PM');
    ITRDB(i).PM.M = M;
    ITRDB(i).PM.PET = PET;
    M = reshape(M,[],1);
    [r,p] = corr(M(find(pyr==syr(1) & pmo==smo(1)):find(pyr==syr(end) & pmo==smo(end))), smap, 'rows','pairwise');
    ITRDB(i).PM.R_SMAP = r;
    ITRDB(i).PM.p_SMAP = p;
    
    clear r p M PET xy DistDeg DistKM xind yind elev phi coast P Tmin Tmax Tdmean smap;
    
end

clear COAST syear eyear i plat plon pmo prismLatLon pyr RootSM_monthly slat slon smapLatLon smo syr tdmean tmax tmin; 

%% Make some maps
states = shaperead('cb_2017_us_nation_5m','UseGeoCoords', true);
latlim = [25.1 49.5];
lonlim = [-125 -66];

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 7];

lat = [ITRDB.LAT];
lon = [ITRDB.LON];

clr = wesanderson('fantasticfox1');
clr = make_cmap([1 1 1; clr(3,:).^1.5; clr(3,:).^4; clr(3,:).^8], 7);

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.R_SMAP];
subplot(3,18,1:9)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 10, r2, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
caxis([0.3 1.0]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'a', 'FontSize',12)
title('Thornthwaite', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.05;

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.R_SMAP];
subplot(3,18,10:18)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 10, r2, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
caxis([0.3 1.0]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'b', 'FontSize',12)
title('Hargreaves', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.55;

% Priestley-Taylor
r2 = [ITRDB.PT];
r2 = [r2.R_SMAP];
subplot(3,18,19:27)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 10, r2, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
caxis([0.3 1.0]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'c', 'FontSize',12)
title('Priestley-Taylor', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.05;

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.R_SMAP];
subplot(3,18,28:36)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'EdgeColor',[0.6 0.6 0.6])
axis off;
axis image;
p1 = scatterm(lat, lon, 10, r2, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
caxis([0.3 1.0]);
colormap(clr);
subplotsqueeze(gca, 1.1)
text(-0.38,0.86,'d', 'FontSize',12)
title('Penman-Monteith', 'FontSize',12);
ax = gca;
ax.Position(1) = 0.55;

cb = colorbar('southoutside');
cb.Position = [0.1 0.37 0.8 0.02];
cb.TickLength = 0.027;
xlabel(cb, 'Pearson''s R', 'FontSize',10);
cbarrow('left');

% Add Boxplots by ecoregion
clr = wesanderson('fantasticfox1');
ttl = {'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'};
alphabet = 'abcdefghijklmnopqrstuvwxyz';
for i = 1:length(ecos)
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    
    ind = [i*2+35 i*2+36];
    
    n = length(ITRDB_sub);
    
    r2 = [ITRDB_sub.Th];
    dat(:, 1) = [r2.R_SMAP];
    r2 = [ITRDB_sub.Hg];
    dat(:, 2) = [r2.R_SMAP];
    r2 = [ITRDB_sub.PT];
    dat(:, 3) = [r2.R_SMAP];
    r2 = [ITRDB_sub.PM];
    dat(:, 4) = [r2.R_SMAP];
    
    if i==1
        axes('Position',[0.02,0.03,0.12,0.22])
    else
        axes('Position',[0.1*(i-1)+0.07,0.03,0.07,0.22])
    end
    bxp = boxplot(dat, 'PlotStyle','compact', 'OutlierSize',3);
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    t = get(a,'tag');   % List the names of all the objects 
    set(a(17), 'Color', clr(4,:));   % Set the color of the first box to green
    set(a(18), 'Color', clr(3,:));   % Set the color of the first box to green
    set(a(19), 'Color', clr(2,:));   % Set the color of the first box to green
    set(a(20), 'Color', clr(1,:));   % Set the color of the first box to green
    set(a(21:24), 'Color', 'k');   % Set the color of the first box to green
    box off;
    ax = gca;
    set(ax, 'TickDir','out', 'TickLength',[0.02 0.02], 'YTick',0:0.2:1, 'YLim',[0 1])
    if i>1; set(ax, 'YTickLabel',''); end
    if i==1; ylabel('Pearson''s R'); end
    text(1,1.03,alphabet(i+4), 'FontSize',12);
    text(2.5, -0.06, ttl{i}, 'FontSize',10, 'HorizontalAlignment','center');
    
    clear dat;

end

h = findobj(gcf,'tag','Outliers');
for iH = 1:length(h)
    h(iH).MarkerEdgeColor = 'k';
end    
h = findobj(gcf,'tag','MedianOuter');
for iH = 1:length(h)
    h(iH).MarkerEdgeColor = 'k';
end    

hLegend = legend(findall(gca,'Tag','Box'), {'Penman-Monteith','Priestley-Taylor','Hargreaves','Thornthwaite'}, 'Orientation','horizontal');
hLegend.Position = [0.14    0.27    0.7179    0.0253];
hLegend.FontSize = 10;
legend('boxoff')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-vslite-smap-maps.tif')
close all;

clear hLegend h ax cb n iH i bxp a t r2 p1 ind;

%% Make table of VS-Lite/SMAP correlations
T = table({'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'}',...
    {'Northern Forests','Northwestern Forested Mountains',...
    'Marine West Coast Forest','Eastern Temperate Forests','Great Plains',...
    'North American Deserts','Mediterranean California',...
    'Southern Semi-Arid Highlands','Temperate Sierras'}', 'VariableNames',{'Code','Name'});
n = NaN(9,1);
Th = cell(9,1);
Hg = cell(9,1);
PT = cell(9,1);
PM = cell(9,1);

for i = 1:length(ecos)
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    n(i) = length(ITRDB_sub);
    
    % Thornthwaite
    r2 = [ITRDB_sub.Th];
    r2 = [r2.R_SMAP];
    ci = bootci(1000,@median,r2);
    Th{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % Hargreaves
    r2 = [ITRDB_sub.Hg];
    r2 = [r2.R_SMAP];
    ci = bootci(1000,@median,r2);
    Hg{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % Priestly-Taylor
    r2 = [ITRDB_sub.PT];
    r2 = [r2.R_SMAP];
    ci = bootci(1000,@median,r2);
    PT{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % Penman-Monteith
    r2 = [ITRDB_sub.PM];
    r2 = [r2.R_SMAP];
    ci = bootci(1000,@median,r2);
    PM{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

end

T.n = n;
T.Th = Th;
T.Hg = Hg;
T.PT = PT;
T.PM = PM;

% Add row for "All" sites
T.Name(10) = {'All'};
T.n(10) = length(ITRDB);

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.R_SMAP];
ci = bootci(1000,@median,r2);
T.Th{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.R_SMAP];
ci = bootci(1000,@median,r2);
T.Hg{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Priestly-Taylor
r2 = [ITRDB.PT];
r2 = [r2.R_SMAP];
ci = bootci(1000,@median,r2);
T.PT{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.R_SMAP];
ci = bootci(1000,@median,r2);
T.PM{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Save table
writetable(T, './output/supplemental-vslite-smap-median-ci.csv');
clear T n Th Hg PT PM r2 ci;

%% Add figures for mean seasonal soil moisture cycle and annual soil moisture averaged by ecoregion
clr = wesanderson('moonrise4');
clr=clr([5 2 4 3],:);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 9];

yr = ppt.year;

for i = length(ecos):-1:1
    
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    n = length(ITRDB_sub);
    
    x = 0.025 + (length(ecos)-i)*(0.95 / length(ecos));
    
    % Seasonal cycle
    subplot(9, 4, i*4-3)
    
    % Thornthwaite
    model = [ITRDB_sub.Th];
    T0 = [model.Tplus0];
    Mseas = reshape([T0.Mseas], 12, [])';
    plot(1:12, mean(Mseas), '-', 'Color',clr(1,:), 'LineWidth',2)
    hold on;
    
    % Hargreaves
    model = [ITRDB_sub.Hg];
    T0 = [model.Tplus0];
    Mseas = reshape([T0.Mseas], 12, [])';
    plot(1:12, mean(Mseas), '-', 'Color',clr(2,:), 'LineWidth',2)
    
    % Priestley-Taylor
    model = [ITRDB_sub.PT];
    T0 = [model.Tplus0];
    Mseas = reshape([T0.Mseas], 12, [])';
    plot(1:12, mean(Mseas), '-', 'Color',clr(3,:), 'LineWidth',2)
    
    % Penman-Monteith
    model = [ITRDB_sub.PM];
    T0 = [model.Tplus0];
    Mseas = reshape([T0.Mseas], 12, [])';
    plot(1:12, mean(Mseas), '-', 'Color',clr(4,:), 'LineWidth',2)
    hold off;
    
    box off;
    ax1 = gca;
    ax1.Position(1) = 0.08;
    ax1.Position(2) = x;
    ax1.Position(3) = 0.22;
    ax1.Position(4) = 0.08;
    ax1.YLim = [floor(ax1.YLim(1)*20)/20 ceil(ax1.YLim(2)*20)/20]; % Round up/down to nearest 0.05
    ax1.TickLength = [0.03 0];
    ax1.TickDir = 'out';
    ax1.XTick = 1:12;
    ax1.XLim = [0.5 12.5];
    if i<length(ecos)
        ax1.XTickLabels = '';
    else
        ax1.XTickLabels = {'J','F','M','A','M','J','J','A','S','O','N','D'};
    end
    ax1.FontSize = 7;
    
    text(1, ax1.YLim(2), alphabet(i*2-1), 'FontSize',11)
    if i == ceil(length(ecos)/2); ylabel('Mean soil moisture (vol/vol)', 'FontSize',11); end

    
    % Annual time series
    subplot(9, 4, [i*4-2 i*4])
    
    % Thornthwaite
    model = [ITRDB_sub.Th];
    M = mean(reshape(mean([model.M])', length(yr), [])'); M(1) = NaN;
    plot(yr, M, '-', 'Color',clr(1,:), 'LineWidth',2)
    hold on;
    
    % Hargreaves
    model = [ITRDB_sub.Hg];
    M = mean(reshape(mean([model.M])', length(yr), [])'); M(1) = NaN;
    plot(yr, M, '-', 'Color',clr(2,:), 'LineWidth',2)
    
    % Priestley-Taylor
    model = [ITRDB_sub.PT];
    M = mean(reshape(mean([model.M])', length(yr), [])'); M(1) = NaN;
    plot(yr, M, '-', 'Color',clr(3,:), 'LineWidth',2)
    
    % Penman-Monteith
    model = [ITRDB_sub.PM];
    M = mean(reshape(mean([model.M])', length(yr), [])'); M(1) = NaN;
    plot(yr, M, '-', 'Color',clr(4,:), 'LineWidth',2)
    
    box off;
    ax2 = gca;
    ax2.YLim = ax1.YLim;
    ax2.XLim = [min(yr) max(yr)];
    ax2.Position(2) = ax1.Position(2);
    ax2.Position(4) = ax1.Position(4);
    ax2.TickLength = [0.01 0];
    ax2.TickDir = 'out';
    ax2.FontSize = 7;
    ax2.YTickLabels = '';
    if i<length(ecos); ax2.XTickLabels = ''; end
    if i == 1
        lgd = legend('Th','Hg','PT','PM', 'Location','northeast');
        legend('boxoff')
        lgd.Position(2) = 0.94;
    end
    
    text(1897, ax2.YLim(2), alphabet(i*2), 'FontSize',11)
    
    
    
end

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-vslite-monthly-annual-soilmoisture.tif')
close all;


