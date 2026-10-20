% Plot the inputs of the stage 2 model

font_size = 10;

%% Melting function parameters

figure('Name','Melting function parameters')

subplot(1,2,1)

plot(Tsol_Pdt, P,'b','LineWidth',2); hold on
plot(Tsol_Pyr, P,'r','LineWidth',2);
if PropEcl2>0; plot(Tsol_Ecl, P,'g','LineWidth',2); end
plot(Tadiabat,P,'--k');
plot(Temperature,P,'LineWidth',2,'Color','k');
plot(realT_Pdt(node_list_Pdt),realP_Pdt(node_list_Pdt)./constants.Pbar, '-ob');
if PropPyr>0; plot(realT_Pyr(node_list_Pyr),realP_Pyr(node_list_Pyr)./constants.Pbar, '-or'); end
if PropEcl2>0; plot(realT_Ecl(node_list_Ecl),realP_Ecl(node_list_Ecl)./constants.Pbar, '-og'); end

set(gca, 'YDir','reverse')
xlim([1100, 1500])
xlabel('T [ºC]')
ylabel('P [GPa]')

title('Thermal regime')
if PropPyr>0 && PropEcl2>0
    legend('Peridotite solidus','Pyroxenite solidus','Eclogite solidus', 'Adiabat','Actual P-T path', 'Peridotite (pMELTS)','Pyroxenite (pMELTS)','Eclogite (pMELTS)')
elseif PropPyr>0
    legend('Peridotite solidus','Pyroxenite solidus', 'Adiabat','Actual P-T path', 'Peridotite (pMELTS)','Pyroxenite (pMELTS)')
elseif PropEcl2>0
    legend('Peridotite solidus','Eclogite solidus', 'Adiabat','Actual P-T path', 'Peridotite (pMELTS)','Eclogite (pMELTS)')
else
    legend('Peridotite solidus', 'Adiabat','Actual P-T path', 'Peridotite (pMELTS)')
end

subplot(1,2,2)

plot(100.*MeltingD_Pdt,Pressure,'b'); hold on
scatter(100.*realF_Pdt(node_list_Pdt),realP_Pdt(node_list_Pdt)./constants.Pbar,'b');
if PropPyr>0
    plot(100.*MeltingD_Pyr,Pressure,'r'); 
    scatter(100.*realF_Pyr(node_list_Pyr),realP_Pyr(node_list_Pyr)./constants.Pbar,'r');
end
if PropEcl2>0
    plot(100.*MeltingD_Ecl,Pressure,'g');
    scatter(100.*realF_Ecl(node_list_Ecl),realP_Ecl(node_list_Ecl)./constants.Pbar,'g');
end

set(gca, 'YDir','reverse')
xlim([0,100])
xlabel('F [wt %]')
ylabel('P [GPa]')

title('Melting degree')
if PropPyr>0 && PropEcl2>0
    legend('Peridotite (interpolated)','Peridotite (pMELTS)', ...
           'Pyroxenite (interpolated)','Pyroxenite (pMELTS)', ...
           'Eclogite (interpolated)','Eclogite (pMELTS)')
elseif PropPyr>0
    legend('Peridotite (interpolated)','Peridotite (pMELTS)', ...
           'Pyroxenite (interpolated)','Pyroxenite (pMELTS)')
elseif PropEcl2>0
    legend('Peridotite (interpolated)','Peridotite (pMELTS)', ...
           'Eclogite (interpolated)','Eclogite (pMELTS)')
else
    legend('Peridotite (interpolated)','Peridotite (pMELTS)')
end

set(findall(gcf,'type','axes'),'FontSize',font_size)

%% Mixing function

figure('Name','Mixing function')

subplot(1,3,1)

plot(100.*fPdt,P,'b'); hold on
if PropPyr>0; plot(100.*fPyr,P,'r'); end
if PropEcl2>0; plot(100.*fEcl,P,'g'); end

set(gca, 'YDir','reverse')
xlabel('wt %')
ylabel('P [GPa]')

title('Melt increment')
if PropPyr>0 && PropEcl2>0
    legend('Peridotite','Pyroxenite','Eclogite')
elseif PropPyr>0
    legend('Peridotite','Pyroxenite')
elseif PropEcl2>0
    legend('Peridotite','Eclogite')
else
    legend('Peridotite')
end

subplot(1,3,2)

plot(100.*FracPdt,P','b'); hold on
if PropPyr>0; plot(100.*FracPyr,P,'r'); end
if PropEcl2>0; plot(100.*FracEcl,P,'g'); end

set(gca, 'YDir','reverse')
xlabel('wt %')
ylabel('P [GPa]')

title('Weighing functions')
if PropPyr>0 && PropEcl2>0
    legend('Peridotite','Pyroxenite','Eclogite')
elseif PropPyr>0
    legend('Peridotite','Pyroxenite')
elseif PropEcl2>0
    legend('Peridotite','Eclogite')
else
    legend('Peridotite')
end

subplot(1,3,3)

plot(v,P, 'c'); hold on
plot(Urate./(sum(Urate)),P, 'm');

set(gca, 'YDir','reverse')
ylabel('P [GPa]')

title('Relative weight in MORB')
legend('Increment volume','Upwelling rate')  

set(findall(gcf,'type','axes'),'FontSize',font_size)
