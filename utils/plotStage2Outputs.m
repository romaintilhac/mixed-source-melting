% Plot the outputs of the stage 2 model

font_size = 10;

%% Major-element compositions

figure('Name','Major-element compositions')

ME_plotlist=split("SiO2 TiO2 Al2O3 Cr2O3 FeO MgO CaO Na2O");
for n=1:length(ME_plotlist)
    ME=string(ME_plotlist(n));
    subplot(4,2,n); title(ME);
    plot(ME_meltPdt.(ME),P,'b');hold on
    if PropPyr>0; plot(ME_meltPyr.(ME),P,'r'); end
    if PropEcl2>0; plot(ME_meltEcl.(ME),P,'g'); end
    plot(ME_residuePdt.(ME),P, '--b');
    if PropPyr>0; plot(ME_residuePyr.(ME),P,'--r'); end
    if PropEcl2>0; plot(ME_residueEcl.(ME),P,'--g'); end
    title(regexprep(ME,'(\d+)','_{$1}'),'Interpreter','tex')
    set(gca,'LineWidth',1,'TickLength',[0.01 0.01]); set(gca, 'YDir','reverse');xlabel('wt %');ylabel('P (GPa)');
end

% edit plot n (e.g. 2 = Ti)
subplot(4,2,2);
if PropEcl2 > 0 && PropPyr > 0
    legend('Peridotite melt', 'Pyroxenite melt', 'Eclogite melt', ...
           'Peridotite residue', 'Pyroxenite residue', 'Eclogite residue')
elseif PropEcl2 > 0 && PropPyr == 0
    legend('Peridotite melt', 'Eclogite melt', ...
           'Peridotite residue', 'Eclogite residue')
elseif PropEcl2 == 0 && PropPyr > 0
    legend('Peridotite melt', 'Pyroxenite melt', ...
           'Peridotite residue', 'Pyroxenite residue')
else
    legend('Peridotite melt', 'Peridotite residue')
end

set(findall(gcf,'type','axes'),'FontSize',font_size)

%% Sulfur systematics

figure('Name','Sulfur systematics')

subplot(1,3,1);

plot(SCSS_Pdt,P,'b');hold on;
if PropPyr>0; plot(SCSS_Pyr,P,'r'); end
if PropEcl2>0; plot(SCSS_Ecl,P,'g'); end
title('SCSS');
set(gca,'LineWidth',1,'TickLength',[0.01 0.01]); set(gca,'LineWidth',1,'TickLength',[0.01 0.01]);set(gca, 'YDir','reverse');xlabel('µg·g⁻¹');ylabel('P (GPa)');

subplot(1,3,2);

plot(S_residuePdt(:,1),P,'b');hold on;
if PropPyr>0; plot(S_residuePyr(:,1),P,'r'); end
if PropEcl2>0; plot(S_residueEcl(:,1),P,'g'); end
title('Sulfur in residue');
set(gca,'LineWidth',1,'TickLength',[0.01 0.01]); set(gca,'LineWidth',1,'TickLength',[0.01 0.01]);set(gca, 'YDir','reverse');xlabel('µg·g⁻¹');ylabel('P (GPa)');

subplot(1,3,3); 

plot(100.*nPdtSul,P,'b');hold on;
if PropPyr>0; plot(100.*nPyrSul,P,'r'); end
if PropEcl2>0; plot(100.*nEclSul,P,'g'); end
plot(100.*nPdtMSS,P,'--b');
if PropPyr>0; plot(100.*nPyrMSS,P,'--r'); end
if PropEcl2>0; plot(100.*nEclMSS,P,'--g'); end
title('Sulfide mode');
if PropEcl2 > 0 && PropPyr > 0
    legend('Sulf (Pdt)', 'Sulf (Pyr)', 'Sulf (Ecl)', ...
           'MSS (Pdt)', 'MSS (Pyr)', 'MSS (Ecl)')
elseif PropEcl2 > 0 && PropPyr == 0
    legend('Sulf (Pdt)', 'Sulf (Ecl)', ...
           'MSS (Pdt)', 'MSS (Ecl)')
elseif PropEcl2 == 0 && PropPyr > 0
    legend('Sulf (Pdt)', 'Sulf (Pyr)', ...
           'MSS (Pdt)', 'MSS (Pyr)')
else
    legend('Sulf (Pdt)', 'MSS (Pdt)')
end
set(gca,'LineWidth',1,'TickLength',[0.01 0.01]); set(gca,'LineWidth',1,'TickLength',[0.01 0.01]);set(gca, 'YDir','reverse');xlabel('wt %');ylabel('P (GPa)');

set(findall(gcf,'type','axes'),'FontSize',font_size)

%% Chalcophile trace-element compositions

figure('Name','Chalcophile trace-element compositions')

nCEmax = 1; % Se only

nplot=0;
for n=1:nCEmax
     
    CE=string(lst.CE(n));

    subplot(nCEmax,3,1); 
    plot(CE_residuePdt(:,1),P,'b');hold on;
    if PropPyr>0; plot(CE_residuePyr(:,1),P,'r'); end
    if PropEcl2>0; plot(CE_residueEcl(:,1),P,'g'); end
    title(append(CE, ' in residue'));
    set(gca,'LineWidth',1,'TickLength',[0.01 0.01]);set(gca, 'YDir','reverse');xlabel('ng·g⁻¹');ylabel('P (GPa)');

    subplot(nCEmax,3,2); 
    plot(CE_meltPdt(:,1),P,'b');hold on;
    if PropPyr>0; plot(CE_meltPyr(:,1),P,'r'); end
    if PropEcl2>0; plot(CE_meltEcl(:,1),P,'g'); end
    title(append(CE, ' in melt'));
    set(gca,'LineWidth',1,'TickLength',[0.01 0.01]);set(gca, 'YDir','reverse');xlabel('ng·g⁻¹');ylabel('P (GPa)');
    set(gca, 'XScale', 'log')
    
    subplot(nCEmax,3,3); 
    plot(DPdtCE(:,1),P,'b');hold on;
    if PropPyr>0; plot(DPyrCE(:,1),P,'r'); end
    if PropEcl2>0; plot(DEclCE(:,1),P,'g'); end
    title(append('Bulk Kd ', CE));
    set(gca,'LineWidth',1,'TickLength',[0.01 0.01]);set(gca, 'YDir','reverse');ylabel('P (GPa)');

end

set(findall(gcf,'type','axes'),'FontSize',font_size)