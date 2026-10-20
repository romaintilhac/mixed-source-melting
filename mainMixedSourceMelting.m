%% MAIN MIXED-SOURCE MELTING
% Adaptative MonteCarlo simulations based on the mixed-source melting model MSM3
% By Romain Tilhac, October 2026
% Contact: romain.tilhac@umontpellier.fr

%% Load environment

% clear
clear; % workspace
clc; % command
close all; % figures
fclose('all'); % files
delete(findall(0,'type','figure','tag','TMWWaitbar')); % waitbars

% load paths
restoredefaultpath
addpath(fullfile(pwd,'utils'));
pth=userFolderPaths; %load folder paths

% build tables
forcebuildTable=0;
cd(pth.resources)
if ~isfile('TableIsotope.mat') || forcebuildTable==1; buildTableIsotope(pth); clear getTableIso; end
if ~isfile('TableMineral.mat') || forcebuildTable==1; buildTableMineral(pth); clear getTableMineral; end
if ~isfile('TableOxide.mat') || forcebuildTable==1; buildTableOxide(pth); clear getTableOxide; end
load TableIsotope.mat  TableIso
load TableMineral.mat  TableMin
load TableOxide.mat  TableOxi
cd(pth.default)

% load model
lst = userModelComponents; %

[lst, option, settings, parameters_init, parameters_references, compo_init, compo_references, Kd, MonteCarlo, targets] = initializeSetup(pth.setup_name, lst, pth); %initialize setup

closevar; %close all opened variables

%% Initialize simulation

% handle single run vs MC
if option.MonteCarlo==1
    
    % build Monte Carlo save
    [MonteCarlo_save, Ntarget_outputs] = MonteCarloTarget(targets, lst);

    % waitbar
    h1 = waitbar(0,'Running Monte Carlo pass', 'Position', [588 500 360 75], 'windowstyle', 'modal');
    frames = java.awt.Frame.getFrames();frames(end).setAlwaysOnTop(1); 

else
    
    % prepare for single run
    Niterations=1;
    settings.maxcycle=0;
    runID=1;
    run_output=pth.run;

end

% create/clear pth.folder
if exist(pth.folder, 'dir')
    warning("The folder " + pth.folder + " already exists and will be deleted.")
    keyboard
    rmdir(pth.folder,"s");
end
mkdir(pth.folder);

% set init parameters & compo
compo=compo_init;
parameters=parameters_init;
compo.ME_PdtRef=compo_references.ME_PdtRef;
compo.ME_PyrRef=compo_references.ME_PyrRef;
compo.ME_EclRef=compo_references.ME_EclRef;

%% Run simulation

t0 = datetime('now');
pass=0; stagnation_count=0;
while pass <= settings.maxcycle
    
    %% Prepare pass
    if option.MonteCarlo==1
        
        % pass#
        passID = num2str(pass);
        str_pass=strcat('pass', passID);
        str_output=strcat('output',passID);
        str_array=strcat('array', passID);        

        %adaptive/pre-adaptive phase
        if pass==0
            Niterations=settings.iterations_preadapt;
        else
            Niterations=settings.iterations;
        end

        % store & plot randomized inputs
        [MonteCarlo_save] = MonteCarloRand(MonteCarlo,MonteCarlo_save,lst,str_pass,str_array);
        plotMCinputs

        % waitbar
        waitbar(pass/settings.maxcycle,h1,['Running Monte Carlo passes: ' num2str(round(100*pass/settings.maxcycle)), '%']);
        disp(['///////////////////////// MONTE CARLO PASS ' num2str(pass) ' /////////////////////////']);
        if pass==0
            h2 = waitbar(0,'Running Monte Carlo simulations', 'Position', [588 400 360 75], 'windowstyle', 'modal');
            frames = java.awt.Frame.getFrames();frames(end).setAlwaysOnTop(1); 
        end

    end

    %% Run iteration(s)
    for runID=1:Niterations
        tic

        % initialize MonteCarlo
        if option.MonteCarlo==1

            % print run header
            disp("########## " + string(datetime('now')) + " - RUN " + runID + " OF PASS " + pass + " ##########")
            waitbar(runID/Niterations,h2,['Running Monte Carlo simulation: ' num2str(round(100*runID/Niterations)), '%']);
            run_output=fullfile(pth.folder, ['pass' num2str(pass) '_run' num2str(runID)]);
        
            % overwrite fixed user inputs with randomized inputs       
            for i=1:length(MonteCarlo_save.rand_fields)
                field=MonteCarlo_save.rand_fields(i);
                if ~ismember(field, lst.compositions) % parameters
                    parameters.(field) = MonteCarlo.(field).(strcat('array',num2str(pass)))(runID);
                else % compositions
                    rand_subfields = fieldnames(MonteCarlo.(field));
                    for j = length(rand_subfields)
                        subfield = string(rand_subfields(j));
                        compo.(field).(subfield) = MonteCarlo.(field).(subfield).(strcat('array',num2str(pass)))(runID);
                    end
                end
            end

        end

        % run MSM3
        [param_out, compo_out, add_outputs, exitFlag] = MSM3wrapper(parameters, compo, Kd, lst, pth, option, settings);

        calc_time=toc;
        fprintf('Calc. time = %0.2f sec. (tot. = %0.0f min.)\n',calc_time, minutes(datetime('now') - t0))

        % save Monte Carlo run outputs
        if option.MonteCarlo==1
            
            if exitFlag == 1
                
                MonteCarlo_save = MonteCarloOutput(MonteCarlo_save, targets, param_out, compo_out, add_outputs, passID, runID);
        
                % calculate Goodness & add to arrays
                GF=sum((abs(MonteCarlo_save.(str_output)(runID,:)-MonteCarlo_save.target))./MonteCarlo_save.error);
                Goodness=exp((-GF));
                for i=1:length(MonteCarlo_save.rand_fields)
                    field = MonteCarlo_save.rand_fields(i);
                    if (any(ismember(fieldnames(parameters), field)))
                        MonteCarlo.(field).(str_array)(runID,2)=Goodness;
                    else
                        rand_subfields=fieldnames(MonteCarlo.(field));
                        for j=1:length(rand_subfields)
                            subfield=string(rand_subfields(j));
                            MonteCarlo.(field).(subfield).(str_array)(runID,2)=Goodness;
                        end
                    end
                end
                
                % add Goodness to MonteCarlo_save
                MonteCarlo_save.(str_pass)(runID,MonteCarlo_save.Nrand_inputs+MonteCarlo_save.Ntarget_outputs+MonteCarlo_save.Npassive_outputs+1) = Goodness;
            
            end

        else

            % save run_output
            save_list_run=split('param_out compo_out'); % for run outputs
            save(run_output, save_list_run{:});

        end

    end

    %% Process pass
    if option.MonteCarlo==1
        
        plotMCoutputs

        % new distributions weighted by GOF
        for i=1:length(MonteCarlo_save.rand_fields)
            field = MonteCarlo_save.rand_fields(i); % for parameters
            if (any(ismember(fieldnames(parameters), field))) 
                if MonteCarlo.(field).adapt==1
                    MonteCarlo.(field).(strcat('array',num2str(pass+1)))=...
                        randsample(MonteCarlo.(field).(str_array)(:,1),settings.iterations,true,MonteCarlo.(field).(str_array)(:,2)./(sum(MonteCarlo.(field).(str_array)(:,2))));
                else
                    MonteCarlo.(field).(strcat('array',num2str(pass+1)))=MonteCarlo.(field).(str_array)(:,1);
                end
            else 
                rand_subfields=fieldnames(MonteCarlo.(field)); % for compo
                for j=1:length(rand_subfields)
                    subfield=string(rand_subfields(j));
                    if MonteCarlo.(field).(subfield).adapt==1
                        MonteCarlo.(field).(subfield).(strcat('array',num2str(pass+1)))=...
                            randsample(MonteCarlo.(field).(subfield).(str_array)(:,1),settings.iterations,true,MonteCarlo.(field).(subfield).(str_array)(:,2)./(sum(MonteCarlo.(field).(subfield).(str_array)(:,2))));
                    else
                        MonteCarlo.(field).(subfield).(strcat('array',num2str(pass+1)))=MonteCarlo.(field).(subfield).(str_array)(:,1);
                    end
                end
            end
        end
            
        % check termination criteria
        best_fit(pass+1,2)=max(MonteCarlo.(MonteCarlo_save.rand_fields(1)).(str_array)(:,2));
        [best_fit_i,~]=find(MonteCarlo.(MonteCarlo_save.rand_fields(1)).(str_array)(:,2)==best_fit(pass+1,2));
        best_fit(pass+1,1)=best_fit_i(1); % to handle cases where several iterations yielded same GOF
        
        if best_fit(pass+1,2)>=settings.target
            fprintf('Target goodness of fit reached. \nPass %d. Best goodness of fit =  %d.', pass, best_fit(pass+1,2));
            break
        end
        if pass>=1 && abs(best_fit(pass+1,2)-best_fit(pass,2)) < settings.stagnation_threshold
            stagnation_count = stagnation_count + 1;
        else
            stagnation_count = 0;
        end
        if stagnation_count > settings.stagnation_max
            fprintf('Maximum number of stagnating passes reached. \nPass %d. Best goodness of fit =  %d.', pass, best_fit(pass+1,2));
            break
        end
        if pass==settings.maxcycle
            fprintf('Maximum number of passes reached. \nPass %d. Best goodness of fit =  %d.', pass, best_fit(pass+1,2))
        end

        % add bestfit to MonteCarlo_save
        MonteCarlo_save.best_fit=best_fit;

    end
    
    pass=pass+1;

end

if option.MonteCarlo==1

    % save whole simulation output
    save_list=split('lst pth option settings parameters_init parameters_references compo_init compo_references Kd MonteCarlo');
    save(fullfile(pth.folder,'MonteCarlo_parameters.mat'),save_list{:});
    save(fullfile(pth.folder,'MonteCarlo_outputs.mat'), '-struct', 'MonteCarlo_save');
    close(h1);close(h2);

end
