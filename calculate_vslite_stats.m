% Calculate stats

load ./data/ITRDB_vslite.mat;
n = length(ITRDB);
cal_yrs = 1901:1960;
sim_yrs = syear:eyear;

for i = 1:n
    
    yr = ITRDB(i).YEAR;
    rwi = ITRDB(i).STD;
    [~,~,ib] = intersect(cal_yrs,yr);
    [~,~,ia] = intersect(cal_yrs,syear:eyear);
    [~,~,iv] = intersect(yr((max(ib)+1):end), syear:eyear);
    
    % Thornthwaite
    rwi_sim = ITRDB(i).Th.VSLite';
    z = (rwi_sim - nanmean(rwi_sim(ia))) / nanstd(rwi_sim(ia));
    rwi_sc = z*nanstd(rwi(ib)) + nanmean(rwi(ib));
    
    ITRDB(i).Th.VSLite_sc = rwi_sc';
    ITRDB(i).Th.year = sim_yrs;
    ITRDB(i).Th.r2_cal = corr(rwi(ib), rwi_sc(ia), 'rows','pairwise')^2;
    ITRDB(i).Th.r2_val = corr(rwi((max(ib)+1):end), rwi_sc(iv), 'rows','pairwise')^2;
    ITRDB(i).Th.rmse_cal = rmse(rwi(ib), rwi_sc(ia));
    ITRDB(i).Th.rmse_val = rmse(rwi((max(ib)+1):end), rwi_sc(iv));
    ITRDB(i).Th.bias_val = mean(rwi((max(ib)+1):end)-rwi_sc(iv));
    
    % Hargreaves
    rwi_sim = ITRDB(i).Hg.VSLite';
    z = (rwi_sim - nanmean(rwi_sim(ia))) / nanstd(rwi_sim(ia));
    rwi_sc = z*nanstd(rwi(ib)) + nanmean(rwi(ib));
    
    ITRDB(i).Hg.VSLite_sc = rwi_sc';
    ITRDB(i).Hg.year = sim_yrs;
    ITRDB(i).Hg.r2_cal = corr(rwi(ib), rwi_sc(ia), 'rows','pairwise')^2;
    ITRDB(i).Hg.r2_val = corr(rwi((max(ib)+1):end), rwi_sc(iv), 'rows','pairwise')^2;
    ITRDB(i).Hg.rmse_cal = rmse(rwi(ib), rwi_sc(ia));
    ITRDB(i).Hg.rmse_val = rmse(rwi((max(ib)+1):end), rwi_sc(iv));
    ITRDB(i).Hg.bias_val = mean(rwi((max(ib)+1):end)-rwi_sc(iv));
    
    % Priestly-Taylor
    rwi_sim = ITRDB(i).PT.VSLite';
    z = (rwi_sim - nanmean(rwi_sim(ia))) / nanstd(rwi_sim(ia));
    rwi_sc = z*nanstd(rwi(ib)) + nanmean(rwi(ib));
    
    ITRDB(i).PT.VSLite_sc = rwi_sc';
    ITRDB(i).PT.year = sim_yrs;
    ITRDB(i).PT.r2_cal = corr(rwi(ib), rwi_sc(ia), 'rows','pairwise')^2;
    ITRDB(i).PT.r2_val = corr(rwi((max(ib)+1):end), rwi_sc(iv), 'rows','pairwise')^2;
    ITRDB(i).PT.rmse_cal = rmse(rwi(ib), rwi_sc(ia));
    ITRDB(i).PT.rmse_val = rmse(rwi((max(ib)+1):end), rwi_sc(iv));
    ITRDB(i).PT.bias_val = mean(rwi((max(ib)+1):end)-rwi_sc(iv));
    
    % Penman-Monteith
    rwi_sim = ITRDB(i).PM.VSLite';
    z = (rwi_sim - nanmean(rwi_sim(ia))) / nanstd(rwi_sim(ia));
    rwi_sc = z*nanstd(rwi(ib)) + nanmean(rwi(ib));
    
    ITRDB(i).PM.VSLite_sc = rwi_sc';
    ITRDB(i).PM.year = sim_yrs;
    ITRDB(i).PM.r2_cal = corr(rwi(ib), rwi_sc(ia), 'rows','pairwise')^2;
    ITRDB(i).PM.r2_val = corr(rwi((max(ib)+1):end), rwi_sc(iv), 'rows','pairwise')^2;
    ITRDB(i).PM.rmse_cal = rmse(rwi(ib), rwi_sc(ia));
    ITRDB(i).PM.rmse_val = rmse(rwi((max(ib)+1):end), rwi_sc(iv));
    ITRDB(i).PM.bias_val = mean(rwi((max(ib)+1):end)-rwi_sc(iv));
    
    
end

save('./data/ITRDB_vslite.mat','ITRDB');

