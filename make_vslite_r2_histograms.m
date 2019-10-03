% Plot distributions of validation metrics

load ./data/ITRDB_vslite.mat;

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 4];

clr = wesanderson('cavalcanti');

% Thornthwaite
r2 = [ITRDB.Th];
r2 = [r2.rmse_val];
[f,x] = histcounts(r2, 0:0.05:1);
plot(0.025:0.05:1, f, '-', 'Color',clr(1,:), 'LineWidth',2)
hold on;

% Hargreaves
r2 = [ITRDB.Hg];
r2 = [r2.rmse_val];
[f,x] = histcounts(r2, 0:0.05:1);
plot(0.025:0.05:1, f, '-', 'Color',clr(2,:), 'LineWidth',2)

% Priestly-Taylor
r2 = [ITRDB.PT];
r2 = [r2.rmse_val];
[f,x] = histcounts(r2, 0:0.05:1);
plot(0.025:0.05:1, f, '-', 'Color',clr(4,:), 'LineWidth',2)

% Penman-Monteith
r2 = [ITRDB.PM];
r2 = [r2.rmse_val];
[f,x] = histcounts(r2, 0:0.05:1);
plot(0.025:0.05:1, f, '-', 'Color',clr(5,:), 'LineWidth',2)

