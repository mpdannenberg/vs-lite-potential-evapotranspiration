%% Make tables with median responses and CIs
load ./data/ITRDB_simulations.mat;
idx = ~cellfun(@isempty, {ITRDB.EcoL1_Code});
ITRDB = ITRDB(idx);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

T = table({'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'}',...
    {'Northern Forests','Northwestern Forested Mountains',...
    'Marine West Coast Forest','Eastern Temperate Forests','Great Plains',...
    'North American Deserts','Mediterranean California',...
    'Southern Semi-Arid Highlands','Temperate Sierras'}', 'VariableNames',{'Code','Name'});
Th = cell(9,1);
Hg = cell(9,1);
PT = cell(9,1);
PM = cell(9,1);

for i = 1:length(ecos)
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    
    % Thornthwaite
    model = [ITRDB_sub.Th];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    dat = 12*[T2.PET] - 12*[T0.PET];
    ci = bootci(1000,@median,dat);
    Th{i} = [num2str(round(median(dat),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];

    % hargreaves
    model = [ITRDB_sub.Hg];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    dat = 12*[T2.PET] - 12*[T0.PET];
    ci = bootci(1000,@median,dat);
    Hg{i} = [num2str(round(median(dat),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];
    
    % Priestly-Taylor
    model = [ITRDB_sub.PT];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    dat = 12*[T2.PET] - 12*[T0.PET];
    ci = bootci(1000,@median,dat);
    PT{i} = [num2str(round(median(dat),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];
    
    % Penman-Monteith
    model = [ITRDB_sub.PM];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    dat = 12*[T2.PET] - 12*[T0.PET];
    ci = bootci(1000,@median,dat);
    PM{i} = [num2str(round(median(dat),2)),' [',num2str(round(ci(1),2)),', ',num2str(round(ci(2),2)),']'];
    
    clear dat;

end

T.Th = Th;
T.Hg = Hg;
T.PT = PT;
T.PM = PM;

% continute with other responses (soil moisture and gM)