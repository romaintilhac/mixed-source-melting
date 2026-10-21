function pth = userFolderPaths
% Set folder paths

    %% User changes
    %==============================================================

    % (no space in paths)

    pth.Alpha_install = "/home/romaintilhac/Documents/Alphamelts_1.9";     % path to the AlphaMELTS installation directory
    pth.Alpha_link = "/home/romaintilhac/bin";                                                      % path to the directory containing the AlphaMELTS executable

    pth.setup_name="default_MC";                                                                           % name of the subfolder where setup text files are located
    pth.folder_name="default_MC";                                                                          % name of the subfolder where model outputs will be stored
    pth.run_name="testrun";                                                                                       % for forward single runs only

    %==============================================================

    %% Allocations (NO user changes below)

    pth.Alpha_runID=fullfile(pth.Alpha_link, "run_alphamelts.command");
    pth.Alpha_defaultID=fullfile(pth.Alpha_install,"alphamelts_default_env.txt");

    pth.Alpha_batch=fullfile(pth.Alpha_install,"MSM3","batch"); 
    pth.Alpha_input=fullfile(pth.Alpha_install,"MSM3","input"); 
    pth.Alpha_output=fullfile(pth.Alpha_install,"MSM3","output"); 

    [~,pathLast] = fileparts(pwd);
    if ~contains(pathLast, "utils")
        pth.default=pwd;
    else
        cd ..
        pth.default=pwd;
    end
    
    pth.preload=fullfile(pth.default,"preload"); 
    pth.pMELTS_output1=fullfile(pth.preload,"stage-1 pMELTS");
    pth.pMELTS_output=fullfile(pth.preload,"stage-2 pMELTS");
    pth.MeltPX_output=fullfile(pth.preload,"stage-2 MeltPX");

    pth.utils = fullfile(pth.default,"utils");
    pth.setups = fullfile(pth.default,"setups");
    pth.resources = fullfile(pth.default,"resources"); 
    pth.outputs = fullfile(pth.default,"outputs");
    pth.folder = fullfile(pth.outputs, pth.folder_name);
    pth.run = fullfile(pth.folder, pth.run_name);

end