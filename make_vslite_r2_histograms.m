% Make table of ecoregion-level validation metrics

load ./data/ITRDB_vslite.mat;
alphabet = 'abcdefghijklmnopqrstuvwxyz';

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

idx = ~cellfun(@isempty, {ITRDB.EcoL1_Code});

ITRDB = ITRDB(idx);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ITRDB = ITRDB(ecol1 > 0);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

%% Make table
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
    r2 = [r2.r2_val];
    ci = bootci(1000,@median,r2);
    Th{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % Hargreaves
    r2 = [ITRDB_sub.Hg];
    r2 = [r2.r2_val];
    ci = bootci(1000,@median,r2);
    Hg{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % Priestly-Taylor
    r2 = [ITRDB_sub.PT];
    r2 = [r2.r2_val];
    ci = bootci(1000,@median,r2);
    PT{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % Penman-Monteith
    r2 = [ITRDB_sub.PM];
    r2 = [r2.r2_val];
    ci = bootci(1000,@median,r2);
    PM{i} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

end

T.n = n;
T.Th = Th;
T.Hg = Hg;
T.PT = PT;
T.PM = PM;

%% Add row for "All" sites
T.Name(10) = {'All'};
T.n(10) = length(ITRDB);

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.r2_val];
ci = bootci(1000,@median,r2);
T.Th{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.r2_val];
ci = bootci(1000,@median,r2);
T.Hg{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Priestly-Taylor
r2 = [ITRDB.PT];
r2 = [r2.r2_val];
ci = bootci(1000,@median,r2);
T.PT{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.r2_val];
ci = bootci(1000,@median,r2);
T.PM{10} = [num2str(round(median(r2),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

%% Save table
writetable(T, './output/vslite-r2-median-ci.csv');

