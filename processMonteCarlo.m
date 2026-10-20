% Process Monte Carlo simulation outputs

close all;
clear all;
addpath(fullfile(pwd,'utils'));
pth = userFolderPaths; % load folder paths
output_path = pth.outputs;
folder_name = pth.folder_name;

save_name = 'MonteCarlo_compile';       % .mat file where compiled results are saved or loaded from (if existing)
tab_name = 'tabulated data';                     % .xlsx file where tabulated results will be saved (overwritting if existing)

%% COMPILE

force_recompile = 0; % 1 to recompile even if save_name exists, 0 to compile only if save_name is missing

load(fullfile(output_path,folder_name,strcat('MonteCarlo_parameters.mat')));
MonteCarlo_outputs = load(fullfile(output_path,folder_name,strcat('MonteCarlo_outputs.mat')));

if  fopen(fullfile(output_path,folder_name,strcat(save_name,'.mat')))==-1 || force_recompile == 1

    iterations=settings.iterations;
    iterations_preadapt=settings.iterations_preadapt;
    
    MonteCarlo_compile=[];
    passID_compile=[];
    runID_compile=[];
    runID_preadapt=(1:iterations)';
    runID=(1:iterations_preadapt)';

    for i=1:settings.maxcycle + 1
        pass_string=string(['pass' num2str(i-1)]);
        if i==1
            runID_compile=[runID_compile; (1:iterations_preadapt)'];
            passID_compile=[passID_compile; repmat(i-1, iterations_preadapt,1)];
        else
            runID_compile=[runID_compile; (1:iterations)'];
            passID_compile=[passID_compile; repmat(i-1, iterations,1)];
        end
        MonteCarlo_compile=[MonteCarlo_compile; MonteCarlo_outputs.(pass_string)];
    end
    MonteCarlo_compile=[passID_compile runID_compile MonteCarlo_compile];
    
    % convert to table
    header=['Pass' 'Run' MonteCarlo_outputs.rand_inputs MonteCarlo_outputs.target_outputs MonteCarlo_outputs.passive_outputs 'Goodness'];
    header=replace(header,'.','_');
    MonteCarlo_compile=array2table(MonteCarlo_compile, VariableNames=header);

    % handle failed scenarios
    if option.stage1
        MonteCarlo_compile(MonteCarlo_compile.T1 == 0, :) = []; 
    end
    
    % check proportions
    if ismember('PropEcl2', MonteCarlo_compile.Properties.VariableNames)
        PropEcl2 = min(MonteCarlo_compile.PropEcl2);
    else
        PropEcl2 = 0;
    end
    if ismember('PropPyr', MonteCarlo_compile.Properties.VariableNames)
        PropPyr = min(MonteCarlo_compile.PropPyr);
    else
        PropPyr = 0;
    end

    % stage-2 contributions
    if option.stage2
        Maf_contrib = 0;
        PropMaf = 0;
        if PropPyr>0
            Pyr_contribFactor=MonteCarlo_compile.Pyr_contrib./MonteCarlo_compile.PropPyr;
            Maf_contrib=MonteCarlo_compile.Pyr_contrib;
            PropMaf=MonteCarlo_compile.PropPyr;
        end
        if PropEcl2>0
            Ecl_contribFactor=MonteCarlo_compile.Ecl_contrib./MonteCarlo_compile.PropEcl2;
            Maf_contrib= Maf_contrib + MonteCarlo_compile.Ecl_contrib;
            PropMaf=PropMaf +MonteCarlo_compile.PropEcl2; 
        end
        Maf_contribFactor=Maf_contrib./PropMaf;
        if PropEcl2>0
            isoCE_finalmelt_Se82_rel=abs(MonteCarlo_compile.isoCE_finalmelt_Se82-compo_init.isoCE_Pdt.Se82)./abs(compo_init.isoCE_Ecl.Se82-compo_init.isoCE_Pdt.Se82);
        else
            if PropPyr>0; isoCE_finalmelt_Se82_rel=abs(MonteCarlo_compile.isoCE_finalmelt_Se82-compo_init.isoCE_Pdt.Se82)./abs(compo_init.isoCE_Pyr.Se82-compo_init.isoCE_Pdt.Se82); end
        end
        
        MonteCarlo_compile=addvars(MonteCarlo_compile,PropMaf,'NewVariableNames',{'PropMaf'}, Before='Ur');
        MonteCarlo_compile=addvars(MonteCarlo_compile,Maf_contrib, Maf_contribFactor,isoCE_finalmelt_Se82_rel,'NewVariableNames',{'Maf_contrib', 'Maf_contribFactor', 'isoCE_finalmelt_Se82_rel'}, Before='Goodness');
    end
    
    % temporary handle complex number outputs
    varNames = MonteCarlo_compile.Properties.VariableNames;
    MonteCarlo_compile = varfun(@real, MonteCarlo_compile);
    MonteCarlo_compile.Properties.VariableNames = varNames;

    % save
    save(fullfile(output_path,folder_name,strcat(save_name,'.mat')));
    writetable(MonteCarlo_compile,fullfile(output_path,folder_name,strcat(tab_name,'.xlsx')))

else

    load(fullfile(output_path,folder_name,strcat(save_name,'.mat')));

end

%% FILTER

MonteCarlo_compile_raw = MonteCarlo_compile;
MonteCarlo_compile = sortrows(MonteCarlo_compile,'Goodness','ascend');

plot_filter = 0;        % 0 for all, 1 for pre-adaptative, 2 for adaptive, 3 for adaptive range
stat_filter = 2;         % 0 for all, 1 for pre-adaptative, 2 for adaptive, 3 for adaptive range
plot_split = 1;         % 0 none, 1 for adaptive/preadaptive, 1 for GOF
report_bestfit = 0; % 0 for mean+sdt, 1 for bestfit 

% adaptive range (only used if filter == 3)
pass_min = 5;
pass_max = 5;

% goodness filter
GOF_threshold_exclude = 1e-4;
GOF_threshold_filled = 1e-3;

% build full/preadapt/adapt/range tables
MonteCarlo_compile = MonteCarlo_compile(MonteCarlo_compile.Goodness >= GOF_threshold_exclude, :);
index = MonteCarlo_compile.Pass < 1;
MonteCarlo_compile_preadapt = MonteCarlo_compile(index,:);
index = MonteCarlo_compile.Pass >= 1;
MonteCarlo_compile_adapt = MonteCarlo_compile(index,:);
index = MonteCarlo_compile.Pass >= pass_min & MonteCarlo_compile.Pass <= pass_max;
MonteCarlo_compile_range = MonteCarlo_compile(index,:);

if plot_filter==0
    plot_str = ' (all outputs)';
    MonteCarlo_plot=MonteCarlo_compile;
elseif plot_filter==1
    plot_str = ' (pre-adaptive outputs)';
    MonteCarlo_plot=MonteCarlo_compile_preadapt;
elseif plot_filter==2
    plot_str = ' (adaptive outputs)';
    MonteCarlo_plot=MonteCarlo_compile_adapt;
elseif plot_filter==3
    plot_str = ' (filtered outputs)';
    MonteCarlo_plot=MonteCarlo_compile_range;
end

if stat_filter==0
    stat_str = ' (all outputs)';
    MonteCarlo_stat=MonteCarlo_compile;
elseif stat_filter==1
    stat_str = ' (pre-adaptive outputs)';
    MonteCarlo_stat=MonteCarlo_compile_preadapt;
elseif stat_filter==2
    stat_str = ' (adaptive outputs)';
    MonteCarlo_stat=MonteCarlo_compile_adapt;
elseif stat_filter==3
    stat_str = ' (filtered outputs)';
    MonteCarlo_stat=MonteCarlo_compile_range;
end

if plot_split==1
    MonteCarlo_plot1 =MonteCarlo_compile_adapt;
    MonteCarlo_plot2 =MonteCarlo_compile_preadapt;
elseif plot_split==2
    MonteCarlo_plot1 =MonteCarlo_compile(MonteCarlo_compile.Goodness >= GOF_threshold_filled, :);
    MonteCarlo_plot2 =MonteCarlo_compile(MonteCarlo_compile.Goodness < GOF_threshold_filled, :);
end

if report_bestfit==1
    MonteCarlo_report=MonteCarlo_stat(end,:);
else
    MonteCarlo_stat_array=table2array(MonteCarlo_stat);
    MonteCarlo_report=mean(MonteCarlo_stat_array, 1, 'omitnan');
    MonteCarlo_report_std=std(MonteCarlo_stat_array, 1, 'omitnan');
    MonteCarlo_report_min = min(MonteCarlo_stat_array, [], 1, 'omitnan');
    MonteCarlo_report_max = max(MonteCarlo_stat_array, [], 1, 'omitnan');
    MonteCarlo_report=array2table(MonteCarlo_report,VariableNames=MonteCarlo_stat. Properties. VariableNames);
    MonteCarlo_report_std=array2table(MonteCarlo_report_std,VariableNames=MonteCarlo_stat. Properties. VariableNames);
    MonteCarlo_report_min=array2table(MonteCarlo_report_min,VariableNames=MonteCarlo_stat. Properties. VariableNames);
    MonteCarlo_report_max=array2table(MonteCarlo_report_max,VariableNames=MonteCarlo_stat. Properties. VariableNames);
    MonteCarlo_report_full=[MonteCarlo_report;MonteCarlo_report_std; MonteCarlo_report_min; MonteCarlo_report_max];
end

%% PLOT

loadAnalyticalData

figure('Name', ['Model vs data'])
sgtitle(['Monte Carlo simulation: ', folder_name, plot_str, newline], 'Interpreter', 'none', 'FontWeight', 'bold');

Nplots=9; Nrows=3;
Ncols=Nplots/Nrows;
spot1=1; %boxplot
spot5=2; %normative
spot6=3; %report_filter
spot2=[4, 7]; %histogram
spot3=[5, 8]; %scatter1
spot4=[6, 9]; %scatter2

caxis='Goodness';
if plot_split == 0
    color=MonteCarlo_plot.(caxis);
else
    color1=MonteCarlo_plot1.(caxis);
    color2=MonteCarlo_plot2.(caxis);
end
transparency = 0.6;
transparency1 = 0.6;
transparency2 = 0.1;
size = 50;
size1 = 50;
size2 = 50;
size_target = 15;

Semin=0; Semax=700;
Se82min=-0.5; Se82max=0.1;

%% Box plots

title('Target analytical data');
subplot(Nrows,Ncols,spot1)

boxplot(data.Se82_MORB, 'Positions', 0.4,'Orientation','horizontal'); hold on
boxplot(data.Se82_Pyr, 'Positions', 0.1,'Orientation','horizontal')

h = findobj(gca,'Tag','Box');
patch(get(h(2),'XData'),get(h(2),'YData'),'g','FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),'r','FaceAlpha',.5);
xlim([Se82min Se82max])
xlabel(['δ^{82}Se (' char(8240) ')']);
set(gca,'ytick',[])
set(gca,'YTickLabel',[]);

legend('MORB (Yierpan {\itet al.} 2019)', 'Pyroxenites (This study)')

%% Histograms

title('Monte Carlo MSM3 model simulation');
subplot(Nrows,Ncols,spot2); hold on

bin_width = 0.005;
if option.stage2; histogram(MonteCarlo_plot.isoCE_finalmelt_Se82, 'BinWidth', bin_width, 'FaceColor','g'); end
if option.stage1; histogram(MonteCarlo_plot.isoCE_mixPyr1_Se82, 'BinWidth', bin_width, 'FaceColor','r'); end

xlim([Se82min Se82max])
xlabel(['δ^{82}Se (' char(8240) ')']);
ylabel('N');

labels = {};
if option.stage2; labels{end+1} = 'Stage-2 outputs (MORB)'; end
if option.stage1; labels{end+1} = 'Stage-1 outputs (pyroxenites)'; end
legend(labels{:});

%% Pyroxenites covariation

if option.stage1

    subplot(Nrows,Ncols,spot3); hold on
    title('Stage-1 outputs (pyroxenites)');

    % axes
    xaxis='CE_mixPyr1_Se'; xlabeltxt='Se (ppm)';
    yaxis='isoCE_mixPyr1_Se82'; ylabeltxt=['δ^{82}Se (' char(8240) ')'];
    xindex=find(MonteCarlo_outputs.target_outputs==xaxis);
    yindex=find(MonteCarlo_outputs.target_outputs==yaxis);
    
    % model
    if plot_split == 0
        scatter(MonteCarlo_plot.(xaxis), MonteCarlo_plot.(yaxis), size, color, 'filled', 'MarkerFaceAlpha', transparency);
    else
        scatter(MonteCarlo_plot2.(xaxis), MonteCarlo_plot2.(yaxis), size2, color2, 'MarkerEdgeAlpha', transparency2);
        scatter(MonteCarlo_plot1.(xaxis), MonteCarlo_plot1.(yaxis), size1, color1, 'filled', 'MarkerFaceAlpha', transparency1);
    end
    
    % target
    x=MonteCarlo_outputs.target(xindex);
    y=MonteCarlo_outputs.target(yindex);
    xneg=MonteCarlo_outputs.error(xindex);
    xpos=MonteCarlo_outputs.error(xindex);
    yneg=MonteCarlo_outputs.error(yindex);
    ypos=MonteCarlo_outputs.error(yindex);
    errorbar(x,y,yneg,ypos,xneg,xpos,'LineWidth',0.5,'Color','k','MarkerFaceColor','r','Marker','square','MarkerSize',size_target);

    % caption
    xlim([Semin Semax])
    ylim([Se82min Se82max])
    xlabel(xlabeltxt);
    ylabel(ylabeltxt);
    box on

    % legend
    if plot_split==0
        legend('Model outputs', 'Target analytical data')
    elseif plot_split==1
       legend('Model outputs (pre-adaptive)', 'Model outputs (adaptive)', 'Target analytical data')
    else
       legend({'Model outputs', sprintf('Model outputs (GOF > %.0e)', GOF_threshold_filled), 'Target analytical data'})
    end
    
end

%% MORB covariation

if option.stage2

    subplot(Nrows,Ncols,spot4); hold on
    title('Stage-2 outputs (MORB)');

    % axes    
    xaxis='CE_finalmelt_Se'; xlabeltxt='Se (ppm)';
    yaxis='isoCE_finalmelt_Se82'; ylabeltxt=['δ^{82}Se (' char(8240) ')'];
    xindex=find(MonteCarlo_outputs.target_outputs==xaxis);
    yindex=find(MonteCarlo_outputs.target_outputs==yaxis);

    % model
    if plot_split == 0
        scatter(MonteCarlo_plot.(xaxis), MonteCarlo_plot.(yaxis), size, color, 'filled', 'MarkerFaceAlpha', transparency);
    else
        scatter(MonteCarlo_plot2.(xaxis), MonteCarlo_plot2.(yaxis), size2, color2, 'MarkerEdgeAlpha', transparency2);
        scatter(MonteCarlo_plot1.(xaxis), MonteCarlo_plot1.(yaxis), size1, color1, 'filled', 'MarkerFaceAlpha', transparency1);
    end

    % target
    x=MonteCarlo_outputs.target(xindex);
    y=MonteCarlo_outputs.target(yindex);
    xneg=MonteCarlo_outputs.error(xindex);
    xpos=MonteCarlo_outputs.error(xindex);
    yneg=MonteCarlo_outputs.error(yindex);
    ypos=MonteCarlo_outputs.error(yindex);
    errorbar(x,y,yneg,ypos,xneg,xpos,'LineWidth',0.5,'Color','k','MarkerFaceColor','r','Marker','square','MarkerSize',size_target);
    
    % caption
    xlim([Semin Semax])
    ylim([Se82min Se82max])
    xlabel(xlabeltxt);
    ylabel(ylabeltxt);
    box on

    % legend
    if plot_split==0
        legend('Model outputs', 'Target analytical data')
    elseif plot_split==1
       legend('Model outputs (pre-adaptive)', 'Model outputs (adaptive)', 'Target analytical data')
    else
       legend({'Model outputs', sprintf('Model outputs (GOF > %.0e)', GOF_threshold_filled), 'Target analytical data'})
    end

end

%covariation color bar
colormap(jet)
cbar=colorbar('Ticks',[]);
cbar.Label.String = "Goodness of fit";
cbar.Label.Rotation = 270;
cbar.Label.VerticalAlignment = "bottom";
set(cbar,'XTick',[0:round(max(MonteCarlo_plot.(caxis))/5,2):round(max(MonteCarlo_plot.(caxis)),2)])

%% Normative ternary

if option.stage1

    subplot(Nrows,Ncols,spot5); hold on
    title('Fo-CaTs-Qz projection');
    
    % axes
    xaxis='normX_mixPyr1';
    yaxis='normY_mixPyr1';
    xmin=0; xmax=100;
    ymin=0; ymax=100;
    xindex=find(MonteCarlo_outputs.target_outputs==xaxis);
    yindex=find(MonteCarlo_outputs.target_outputs==yaxis);

    % model
    if plot_split == 0
        scatter(MonteCarlo_plot.(xaxis), MonteCarlo_plot.(yaxis), size, color, 'filled', 'MarkerFaceAlpha', transparency);
    else
        scatter(MonteCarlo_plot2.(xaxis), MonteCarlo_plot2.(yaxis), size2, color2, 'MarkerEdgeAlpha', transparency2);
        scatter(MonteCarlo_plot1.(xaxis), MonteCarlo_plot1.(yaxis), size1, color1, 'filled', 'MarkerFaceAlpha', transparency1);
    end
    
    % target
    x=MonteCarlo_outputs.target(3);
    y=MonteCarlo_outputs.target(4);
    xneg=MonteCarlo_outputs.error(3);
    xpos=MonteCarlo_outputs.error(3);
    yneg=MonteCarlo_outputs.error(4);
    ypos=MonteCarlo_outputs.error(4);
    errorbar(x,y,yneg,ypos,xneg,xpos,'LineWidth',0.5,'Color','k','MarkerFaceColor','r','Marker','square','MarkerSize',size_target);

    % caption
    xlim([xmin xmax]);
    ylim([ymin ymax])
    template=imread('resources/normTemplate.png');
    template2=image(xlim,[xmax xmin],template);
    uistack(template2,'bottom')
    axis tight;
    axis equal

end

%% REPORT

subplot(Nrows,Ncols,spot6);
xpos=0;
yshift=0;
axis off;

if option.stage1
    
    % title
    if report_bestfit == 1
        txt_title = ['Stage-1 pyroxenite formation - best-fit scenario', stat_str];
    else
        txt_title = ['Stage-1 pyroxenite formation', stat_str];
    end

    % values
    if report_bestfit==1
        txt_stage1=sprintf('P = %.1f GPa, T = %.0f ºC (%.0f ºC above Ecl. solidus)\n', MonteCarlo_report.P1, MonteCarlo_report.T1, MonteCarlo_report.T1_excess);
        txt_stage1F=sprintf('Melting degree: Eclogite = %.0f %%, Peridotite = %.0f %%\n\n', 100*MonteCarlo_report.F_Ecl1, 100*MonteCarlo_report.F_Pdt1);
        txt_stage1Ecl=sprintf('Ecl. composition: S = %.0f ppm, Se = %.0f ppm, δ^{82}Se = %.2f\n',MonteCarlo_report.S_Ecl, MonteCarlo_report.CE_Ecl_Se, MonteCarlo_report.isoCE_Ecl_Se82);
        txt_stage1Pyr=sprintf('Pyr. composition: S = %.0f ppm, Se = %.0f ppm, δ^{82}Se = %.2f\n',MonteCarlo_report.S_mixPyr1, MonteCarlo_report.CE_mixPyr1_Se, MonteCarlo_report.isoCE_mixPyr1_Se82);
        txt_stage1contrib=sprintf('Contributions: Ecl. melt = %.0f %%, Pdt. melt = %.0f %%, Pdt. residue = %.0f %%\n', 100*MonteCarlo_report.PropEcl_melt, 100*MonteCarlo_report.PropPdt_melt, 100*MonteCarlo_report.PropPdt_residue);
    else
        txt_stage1=sprintf('P = %.1f ± %.1f GPa, T = %.0f ± %.0f ºC (%.0f ± %.0f ºC above Ecl. solidus)\n', MonteCarlo_report.P1, MonteCarlo_report_std.P1, MonteCarlo_report.T1, MonteCarlo_report_std.T1, MonteCarlo_report.T1_excess, MonteCarlo_report_std.T1_excess);
        txt_stage1F=sprintf('Melting degree: Eclogite = %.0f ± %.0f %%, Peridotite = %.0f ± %.0f %%\n', 100*MonteCarlo_report.F_Ecl1, 100*MonteCarlo_report_std.F_Ecl1, 100*MonteCarlo_report.F_Pdt1, 100*MonteCarlo_report_std.F_Pdt1);
        txt_stage1Ecl=sprintf('Ecl. composition: S = %.0f ppm, Se = %.0f ppm, δ^{82}Se = %.2f\n',MonteCarlo_report.S_Ecl, MonteCarlo_report.CE_Ecl_Se, MonteCarlo_report.isoCE_Ecl_Se82);
        txt_stage1Pyr=sprintf('Pyr. composition: S = %.0f ppm, Se = %.0f ppm, δ^{82}Se = %.2f\n',MonteCarlo_report.S_mixPyr1, MonteCarlo_report.CE_mixPyr1_Se, MonteCarlo_report.isoCE_mixPyr1_Se82);
        txt_stage1contrib=sprintf('Contributions: Ecl. melt = %.0f ± %.0f %%, Pdt. melt = %.0f ± %.0f %%, Pdt. residue = %.0f ± %.0f %%\n\n', 100*MonteCarlo_report.PropEcl_melt, 100*MonteCarlo_report_std.PropEcl_melt, 100*MonteCarlo_report.PropPdt_melt, 100*MonteCarlo_report_std.PropPdt_melt, 100*MonteCarlo_report.PropPdt_residue, 100*MonteCarlo_report_std.PropPdt_residue);
    end

    % print
    text(xpos, yshift + 1.12, txt_title, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'BackgroundColor', 'none');
    text(xpos, yshift+0.7, [txt_stage1, txt_stage1F, txt_stage1contrib txt_stage1Ecl, txt_stage1Pyr], 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'BackgroundColor', 'none');

end

if option.stage2

    % title
    if report_bestfit == 1
        txt_title = ['Stage-2 MORB formation - best-fit scenario', stat_str];
    else
        txt_title = ['Stage-2 MORB formation', stat_str];
    end

    % values
    if report_bestfit == 1

        txt_stage2 = sprintf('Potential temperature = %.0f ºC\n', MonteCarlo_report.Tp);
    
        % Source composition
        src = {};
        if PropPyr > 0; src{end+1} = sprintf('%.0f %% pyroxenite',100*MonteCarlo_report.PropPyr); end
        if PropEcl2 > 0; src{end+1} = sprintf('%.0f %% eclogite',100*MonteCarlo_report.PropEcl2); end
        txt_stage2source = sprintf('Source composition: %s\n\n', strjoin(src, ', '));
    
        % Contributions
        contrib = {};
        if PropPyr > 0; contrib{end+1} = sprintf('%.0f %% Pyr-derived',100*MonteCarlo_report.Pyr_contrib); end
        if PropEcl2 > 0; contrib{end+1} = sprintf('%.0f %% Ecl-derived',100*MonteCarlo_report.Ecl_contrib); end
        txt_stage2contrib = sprintf('Contributions: %s\n', strjoin(contrib, ', '));
        txt_stage2MORB = sprintf('MORB composition: Se = %.0f ppm, δ^{82}Se = %.2f\n', MonteCarlo_report.CE_finalmelt_Se, MonteCarlo_report.isoCE_finalmelt_Se82);
    
    else

        txt_stage2 = sprintf('Potential temperature = %.0f ± %.0f ºC\n', MonteCarlo_report.Tp, MonteCarlo_report_std.Tp);
    
        % Source composition
        src = {};
        if PropPyr > 0; src{end+1} = sprintf('%.0f ± %.0f %% pyroxenite', 100*MonteCarlo_report.PropPyr,100*MonteCarlo_report_std.PropPyr); end
        if PropEcl2 > 0; src{end+1} = sprintf('%.0f ± %.0f %% eclogite', 100*MonteCarlo_report.PropEcl2,100*MonteCarlo_report_std.PropEcl2); end
        txt_stage2source = sprintf('Source composition: %s\n', strjoin(src, ', '));
    
        % Contributions
        contrib = {};
        if PropPyr > 0; contrib{end+1} = sprintf('%.0f ± %.0f %% Pyr-derived', 100*MonteCarlo_report.Pyr_contrib,100*MonteCarlo_report_std.Pyr_contrib); end
        if PropEcl2 > 0; contrib{end+1} = sprintf('%.0f ± %.0f %% Ecl-derived', 100*MonteCarlo_report.Ecl_contrib,100*MonteCarlo_report_std.Ecl_contrib); end
        txt_stage2contrib = sprintf('Contributions: %s\n\n', strjoin(contrib, ', '));
        txt_stage2MORB = sprintf('MORB composition: Se = %.0f ppm, δ^{82}Se = %.2f\n', MonteCarlo_report.CE_finalmelt_Se, MonteCarlo_report.isoCE_finalmelt_Se82);
    
    end

    % print
    text(xpos, yshift + 0.35, txt_title,'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'BackgroundColor', 'none');
    text(xpos, yshift, [txt_stage2,txt_stage2source, txt_stage2contrib, txt_stage2MORB], 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'BackgroundColor', 'none');

end