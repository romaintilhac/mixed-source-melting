% Plot outputs of Monte Carlo simulations

%offline=1;

rand_inputs=MonteCarlo_save.rand_inputs;
target_outputs=MonteCarlo_save.target_outputs;

if exist('offline','var')

    close all; %close figures
    pth=userFolderPaths; % load folder paths
    
    folder_name='Se-Se82-normX-Y_Stage-1inversion';
    setup_name='';
    parameterID= 'MonteCarlo_parameters.mat';
    pass_list=string(fieldnames(MonteCarlo_save));
    pass_list=pass_list(contains(pass_list,'pass'));
    %pass_list=split("pass0 pass1 pass2 pass20");
    stage1only=1;
    if stage1only==1
        output=fullfile(pth.outputs1,folder_name,setup_name);
    else
        output=fullfile(pth.outputs,folder_name,setup_name);
    end
    pth.outputsID=fullfile(output,parameterID);
    load(pth.outputsID);
    Nrand_inputs=length(rand_inputs);
    Ntarget_outputs=length(target_outputs);

else

    pass_list=string(str_pass);
    
end

title_list(1)="Normative Fo-CaTs-Qz stage-1 pyroxenite compositions";
title_list(2)="Stage-1 pyroxenites";
title_list(3)="Final aggregated MORB";

xaxis_list(1)="normX_mixPyr1";
yaxis_list(1)="normY_mixPyr1";
xaxis_list(2)="CE_mixPyr1_Se";
yaxis_list(2)="isoCE_mixPyr1_Se82";
xaxis_list(3)="CE_finalmelt_Se";
yaxis_list(3)="isoCE_finalmelt_Se82";

xlabel_list(1)="";
ylabel_list(1)="";
xlabel_list(2)="Se [ppm]";
ylabel_list(2)="dSe82";
xlabel_list(3)="Se [ppm]";
ylabel_list(3)="dSe82";

Nplots=length(xaxis_list);

if pass==0 || exist('offline','var');fig_fit=figure('Name','Monte Carlo targets vs outputs');end

for i=1:length(pass_list)

    str_pass=pass_list(i);
    pass=double(regexp(str_pass,'\d*','Match'));
    
        for n=1:Nplots
            
            xindexvec=find(strcmp(target_outputs, xaxis_list(n)));
            yindexvec=find(strcmp(target_outputs, yaxis_list(n)));
            xindex=xindexvec+MonteCarlo_save.Nrand_inputs;
            yindex=yindexvec+MonteCarlo_save.Nrand_inputs;
    
            if pass==0 || exist('offline','var')
                plots.(strcat('plot',num2str(n)))=subplot(1,Nplots,n);
                x=MonteCarlo_save.target(xindexvec);
                y=MonteCarlo_save.target(yindexvec);
                xneg=MonteCarlo_save.error(xindexvec);
                xpos=MonteCarlo_save.error(xindexvec);
                yneg=MonteCarlo_save.error(yindexvec);
                ypos=MonteCarlo_save.error(yindexvec);
                errorbar(x,y,yneg,ypos,xneg,xpos,'o'); hold on
                title(title_list(n));
                if n>1
                    xlabel(xlabel_list(n));
                    ylabel(ylabel_list(n));
                end
            else
                figure(fig_fit);
            end
    
            x=MonteCarlo_save.(str_pass)(:,xindex);
            y=MonteCarlo_save.(str_pass)(:,yindex);
            scatter(plots.(strcat('plot',num2str(n))), x, y);
         
            % format normative plot
            if exist('offline','var') && n==1 || pass==0 && n==1
                xlim([0 100]); ylim([0 100])
                template=imread('resources/normTemplate.png');
                template2=image(xlim,[100 0],template);
                uistack(template2,'bottom')
                axis tight; axis equal
            end
    
        end
        if exist('offline','var'); keyboard; end
end