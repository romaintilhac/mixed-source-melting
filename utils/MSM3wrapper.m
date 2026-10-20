function [param_out, compo_out, add_outputs, exitFlag] = MSM3wrapper(param_in, compo_in, Kd, lst, pth, option, settings)
% Runs Stage 1 and Stage 2 models and merges outputs

%% Stage 1

param = param_in;
compo = compo_in;

param_out1=struct();
compo_out1=struct();

if option.stage1
   
    % Stage 1 mixing
    disp('=============== Stage 1 ===============')
    [param_out1, compo_out1] = Stage1model(param, compo, Kd, lst, pth, option, settings);

    fprintf('F: Peridotite = %.1f%% - Eclogite = %.1f%%\n', 100 * param_out1.F_Pdt1, 100 * param_out1.F_Ecl1)
    fprintf('SCSS: Peridotite = %.0f ppm - Eclogite = %.0f ppm\n', param_out1.SCSS_Pdt1, param_out1.SCSS_Ecl1)
    fprintf('S: Peridotite = %.0f ppm - Eclogite = %.0f ppm - Pyroxenite = %.0f ppm\n', param_out1.S_residuePdt1, param_out1.S_residueEcl1, param_out1.S_mixPyr1)
    fprintf('Pyroxenite: Se = %.0f ppb - d82Se = %.3f\n', compo_out1.CE_mixPyr1.Se, compo_out1.isoCE_mixPyr1.Se82);
    fprintf('Pyroxenite produced = %.1f%%\n', 100 * param_out1.PropPyr)

    [param, compo] = Stage1pass(param_in, param_out1, compo_in, compo_out1, lst, option, settings);

end

%% Stage 2

param_out2=struct();
compo_out2=struct();

if option.stage2

    % Stage 2 MSM3
    disp('=============== Stage 2 ===============')
    [param_out2, compo_out2, exitFlag] = Stage2model(param, compo, Kd, lst, pth, option, settings);

    if exitFlag == 1
        if param.PropEcl2>0
            fprintf('Eclogite contribution: %.1f%%\n', param_out2.Ecl_contrib*100);
        end
        if param.PropPyr>0
            fprintf('Pyroxenite contribution: %.1f%%\n', param_out2.Pyr_contrib*100);
        end
        fprintf('MORB Se: %.0f ppb - d82Se: %.3f\n', compo_out2.CE_finalmelt.Se, compo_out2.isoCE_finalmelt.Se82);
    end

end

%% Merge outputs

param_out = catstruct(param_out1, param_out2);
compo_out = catstruct(compo_out1, compo_out2);

% passive outputs
add_outputs=[];
if option.stage1==1; add_outputs=[add_outputs; fieldnames(param_out1)]; end
if option.stage2==1; add_outputs=[add_outputs; split("S_finalmelt Pyr_contrib Ecl_contrib")]; end
add_outputs = string(add_outputs)';

end      
