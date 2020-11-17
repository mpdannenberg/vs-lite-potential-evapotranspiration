%% Make monthly difference plots
load ./data/ITRDB_simulations.mat;
alphabet = 'abcdefghijklmnopqrstuvwxyz';

clr = wesanderson('fantasticfox1');
idx = ~cellfun(@isempty, {ITRDB.EcoL1_Code});

ITRDB = ITRDB(idx);
ecol1 = cellfun(@str2num, {ITRDB.EcoL1_Code});
ecos = sort(unique(ecol1));

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 9];

ha = tight_subplot(9, 3, [0.02 0.1], [0.05 0.04], [0.1 0.05]);

for i = 1:length(ecos)
    
    ITRDB_sub = ITRDB(ecol1 == ecos(i));
    n = length(ITRDB_sub);
    
    %% PET differences
    axes(ha(3*(i-1)+1))
    
    % Thornthwaite
    model = [ITRDB_sub.Th];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.PETseas], 12, []) - reshape([T0.PETseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(1,:), 'LineWidth',2)
    hold on;
    
    % Hargreaves
    model = [ITRDB_sub.Hg];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.PETseas], 12, []) - reshape([T0.PETseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(2,:), 'LineWidth',2)
    
    % Priestly-Taylor
    model = [ITRDB_sub.PT];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.PETseas], 12, []) - reshape([T0.PETseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(3,:), 'LineWidth',2)
    
    % Penman-Monteith
    model = [ITRDB_sub.PM];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.PETseas], 12, []) - reshape([T0.PETseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(4,:), 'LineWidth',2)
    hold off;
    
    set(gca, 'XLim',[1 12], 'TickDir','out', 'TickLength',[0.02 0],...
        'YLim',[0 12],'YTick',0:5:10)
    if i<length(ecos); set(gca, 'XTickLabel',''); end
    if i==length(ecos); xlabel('Month'); end
    ylim = get(gca, 'Ylim');
    if i == 1 
        ttl = title('\DeltaPET', 'FontSize',11); 
        ttl.Position(2) = ylim(2) + (ylim(2) - ylim(1))*0.1;
    end
    box off;
    ylabel('\DeltaPET (mm)')
    ylim = get(gca, 'Ylim');
    text(1.5, ylim(2), [alphabet(i),') ',sprintf('%1.1f',ecos(i))])

    %% M differences
    axes(ha(3*(i-1)+2))
    
    % Thornthwaite
    model = [ITRDB_sub.Th];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.Mseas], 12, []) - reshape([T0.Mseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(1,:), 'LineWidth',2)
    hold on;
    
    % Hargreaves
    model = [ITRDB_sub.Hg];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.Mseas], 12, []) - reshape([T0.Mseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(2,:), 'LineWidth',2)
    
    % Priestly-Taylor
    model = [ITRDB_sub.PT];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.Mseas], 12, []) - reshape([T0.Mseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(3,:), 'LineWidth',2)
    
    % Penman-Monteith
    model = [ITRDB_sub.PM];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.Mseas], 12, []) - reshape([T0.Mseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(4,:), 'LineWidth',2)
    hold off;
    
    set(gca, 'XLim',[1 12], 'TickDir','out', 'TickLength',[0.02 0],...
        'YLim',[-0.015 0], 'YTick',-0.015:0.005:0, 'YTickLabel',{'','-0.01','','0'})
    if i<length(ecos); set(gca, 'XTickLabel',''); end
    if i==length(ecos); xlabel('Month'); end
    ylim = get(gca, 'Ylim');
    if i == 1 
        ttl = title('\DeltaSoil moisture', 'FontSize',11); 
        ttl.Position(2) = ylim(2) + (ylim(2) - ylim(1))*0.1;
    end
    box off;
    ylabel('\DeltaM (vol/vol)')
    text(1.5, ylim(2), [alphabet(i+length(ecos)),') ',sprintf('%1.1f',ecos(i))])

    %% gM differences
    axes(ha(3*(i-1)+3))
    
    % Thornthwaite
    model = [ITRDB_sub.Th];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.gMseas], 12, []) - reshape([T0.gMseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(1,:), 'LineWidth',2)
    hold on;
    
    % Hargreaves
    model = [ITRDB_sub.Hg];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.gMseas], 12, []) - reshape([T0.gMseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(2,:), 'LineWidth',2)
    
    % Priestly-Taylor
    model = [ITRDB_sub.PT];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.gMseas], 12, []) - reshape([T0.gMseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(3,:), 'LineWidth',2)
    
    % Penman-Monteith
    model = [ITRDB_sub.PM];
    T0 = [model.Tplus0];
    T2 = [model.Tplus2];
    pet = reshape([T2.gMseas], 12, []) - reshape([T0.gMseas], 12, []);
    plot(1:12, median(pet, 2)', '-', 'Color',clr(4,:), 'LineWidth',2)
    hold off;
    
    set(gca, 'XLim',[1 12], 'TickDir','out', 'TickLength',[0.02 0],...
        'YLim',[-0.055 0], 'YTick',-0.05:0.01:0, 'YTickLabel',{'','-0.04','','-0.02','','0'})
    if i<length(ecos); set(gca, 'XTickLabel',''); end
    if i==length(ecos); xlabel('Month'); end
    ylim = get(gca, 'Ylim');
    if i == 1 
        ttl = title('\Deltag_{M}', 'FontSize',11); 
        ttl.Position(2) = ylim(2) + (ylim(2) - ylim(1))*0.1;
    end
    box off;
    ylabel('\Deltag_{M}')
    if i<length(ecos)
        text(1.5, ylim(2), [alphabet(i+length(ecos)*2),') ',sprintf('%1.1f',ecos(i))])
    else
        text(1.5, ylim(2), ['aa) ',sprintf('%1.1f',ecos(i))])
    end
    if i == 1
        lgd = legend('Th', 'Hg', 'PT','PM', 'Location','southeast');
        legend('boxoff')
        lgd.FontSize = 6;
        lgd.Position(1) = 0.735;
        lgd.Position(2) = 0.885;
    end

    
end

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/vslite-t+2-monthly.tif')
close all;


