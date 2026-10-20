function [realP, realT, realF, ME_melt, ME_residue, mode, problem, exitFlag]=pMELTSexploration(P, T, F, C, litho, lst, pth, option)
%   Writes batch & input files, run pMELTS, collect and return outputs
  
% Modifications in alphamelts_default_env.txt:
%   ALPHAMELTS_MINP +0
%   ALPHAMELTS_SAVE_ALL true

struct=getTableMin; TableMin=struct.TableMin;

phasediagram = 1;   % 1 to match F, 0 to leave if free
Fmin = 0.01;        % minimum F
Fmax = 0.99;           % maximum F
limitCpx = 1;       % 1 to prevent multiple cpx, 0 to leave it free
limitOpx = 1;       % 1 to prevent multiple opx, 0 to leave it free

% exploration list
Tincr_list = [0 -100 100];           % 0 original T, T increments
guess_list = [0 1];             % 0 superliquidus, 1 subsolidus
monotonic_list = [0 1];         % 1 quick search, 0 monotonic special search
Fdecimal_list = [0 3 2 1];          % 0 original F, rounding decimal

%% Load paths and clear previous outputs

    Alpha_runID=pth.Alpha_runID;
    Alpha_defaultID=pth.Alpha_defaultID;
    Alpha_input=pth.Alpha_input;
    Alpha_batch=pth.Alpha_batch;
    Alpha_output=pth.Alpha_output;

    delete(fullfile(Alpha_output,'*')) %clear previous Alphamelts output files

%% pMELTS iterations

    exitFlag=0;
    while true  

        %%
        for F_i=1:length(Fdecimal_list)
            
            meltingD = F;
            if Fdecimal_list(F_i)>0
                meltingD=round(meltingD,Fdecimal_list(F_i));
                F_str=strcat("F rounded to ", string(meltingD), "\n");
                fprintf(F_str)
            end
            meltingD = min(max(meltingD, Fmin), Fmax);

            for monotonic_i = 1:length(monotonic_list)
            
                if monotonic_list(monotonic_i)==1
                    fprintf('+++ Monotonic special search +++\n');
                end
                monotonic = string(monotonic_list(monotonic_i));      

                %%
                for T_i=1:length(Tincr_list)
    
                    temperature=num2str(T+Tincr_list(T_i));
                    pressure=num2str(P*1e4); % convert to [bar]

                    %%
                    for guess_i=1:1:length(guess_list)
    
                        if guess_list(guess_i)==0
                            guess=["1"];
                            guess_str="suprasolidus";
                        elseif guess_list(guess_i)==1
                            guess=["0";"x"];
                            guess_str="subsolidus";
                        end

                        if ~(monotonic_i == 1 && F_i == 1 && T_i == 1 && guess_i == 1)
                            disp(litho + ": Trying with " + temperature + ", " + guess_str + " guess")
                        end
    
                        %% Write batch file
                    
                            Alpha_batchID= fullfile(Alpha_batch, strcat("batch_", litho, ".txt"));
                            Alpha_inputID=fullfile(Alpha_input, strcat('input_', litho, '.melts'));
                        
                            fileID = fopen(Alpha_batchID, 'w');
                            
                            if phasediagram==0
                                text=['-1'; % turn off outputs
                                '1';
                                Alpha_inputID;
                                '3';            % single (batch) calculation
                                guess;      
                                '0'];           % close
                            else %
                                text=['-1'; % turn off outputs
                                '1';
                                Alpha_inputID;   
                                '10';          % turn phase diagram mode on / off
                                '1';            % bracket reaction (includes follow phi, F, & aH2O contours)
                                '3';            % single (batch) calculation
                                guess;      
                                'liquid';    % track boundary of liquid
                                '1';            % 0 phi (melt fraction by volume), 1 F (melt fraction by mass)
                                meltingD;               % F increment
                                monotonic; % quick/monotonic search
                                '0'];           % close
                            end
                            
                            for n=1:length(text)
                                fprintf(fileID,'%s\n',text(n));
                            end
                            
                            fclose(fileID);
                        
                        %% Write input file
                    
                            fileID = fopen(Alpha_inputID, 'w');
                            
                            %title
                            text='Title: Mixed-source melting';
                            fprintf(fileID,'%s\n',text);
                            
                            %initial compositions
                            for n=1:lst.nME
                                ME_name=lst.ME_oxlist(n);
                                ME_value=num2str(C.(ME_name));
                                text=strcat('Initial Composition:', {' '}, ME_name, {' '}, ME_value);
                                fprintf(fileID,'%s\n',text);
                            end
                            
                            %limit phases
                            text=[""];
                            if limitCpx==1; text=[text "Limit coexisting: clinopyroxene"];end
                            if limitOpx==1; text=[text "Limit coexisting: orthopyroxene"];end
                            fprintf(fileID,'%s\n',text);
                            
                            %P-T-fO2 conditions
                            
                            text=['Initial Temperature: ' temperature];
                            fprintf(fileID,'%s\n',text);
                            
                            text=['Initial Pressure: ' pressure];
                            fprintf(fileID,'%s\n',text);
                            
                            text='Log fO2 Path: NONE';
                            fprintf(fileID,'%s\n',text);
                            
                            text='Log fO2 Delta: 0.0';
                            fprintf(fileID,'%s\n',text);

                            fclose(fileID);
                        
                        %% Run pMELTS & check outputs
                                        
                            problem=[Fdecimal_list(F_i) Tincr_list(T_i) guess_list(guess_i)];

                            % prepare command and run
                            command=Alpha_runID + ' -f ' + Alpha_defaultID + ' -p ' + Alpha_output + ' -b ' + Alpha_batchID;
                            if option.pMELTSnoprint==1
                                evalc('system(command)');
                            else
                                system(command)
                            end
                        
                            % check outputs
                            outputID= fullfile(Alpha_output, 'Solid_comp_tbl.txt');
                            fileID = fopen(outputID);
                    
                            if fileID~=-1
                                exitFlag=1;
                                break
                            end
            
                    end %end of guess iteration
                    if exitFlag; break; end
    
                end %end of T iteration
                if exitFlag; break; end
    
            end %end of monotonic iteration
            if exitFlag; break; end

        end %end of F iteration
        if exitFlag; break; end

        break

    end

%% Read outputs

if exitFlag==1

    % main parameters (P-T-F)

        outputID= fullfile(Alpha_output, 'System_main_tbl.txt');
        fileID = fopen(outputID);
        
        txt = textscan(fileID,'%s','delimiter','\n');
        output_header=split(string(txt{1}(4)));
        output_value=str2num(string(txt{1}(5)))';
        
        index = find(ismember(output_header,'F'));
        realF=output_value(index);
        
        index = find(ismember(output_header,'Temperature'));
        realT=output_value(index);
        
        index = find(ismember(output_header,'Pressure'));
        realP=output_value(index);

        fclose(fileID);

    % phase proportion

        outputID= fullfile(Alpha_output, 'Phase_mass_tbl.txt');
        fileID = fopen(outputID);
        
        txt = textscan(fileID,'%s','delimiter','\n');
        output_header=split(string(txt{1}(4)));
        output_value=str2num(string(txt{1}(5)))';
        
        %check liquid immiscibility
        index=find(contains(output_header,'liquid'));
        if length(index)>1
            if option.pMELTSnoprint==0
                warning(strcat('Liquid immiscibility has been encountered.'));
            end
            %keyboard
        end
        
        %oli-opx-cpx-plg-spl-grt
        mode=zeros(lst.nMineral,1);
        for n=1:lst.nMineral
            index = find(contains(output_header,TableMin.(lst.minerals(n))));
            if length(index)>1
                if option.pMELTSnoprint==0
                    warning(strcat('Several ', {' '}, TableMin.(lst.minerals(n)), {' '}, ' have been encountered.'))
                end
                %keyboard
                index=index(1);
            end
            if ~isempty(index)
                mode(n)=output_value(index);
            else
                mode(n)=0;
            end
        end
        
        % %manually add
        % index = find(contains(output_header,'whitlockite'));
        % if ~isempty(index)
        %     mode(lst.nMineral+1)=output_value(index);
        % else
        %     mode(lst.nMineral+1)=0;
        % end

        fclose(fileID);

    % solid compositions

        outputID= fullfile(Alpha_output, 'Solid_comp_tbl.txt');
        fileID = fopen(outputID);
        
        txt = textscan(fileID,'%s','delimiter','\n');
        output_header=split(string(txt{1}(4)));
        output_value=str2num(string(txt{1}(5)))';
        
        ME_residue=zeros(numel(lst.ME_oxlist),1);
        for n=1:numel(lst.ME_oxlist)
            index = find(contains(output_header,lst.ME_oxlist(n)));
            ME_residue(n)=output_value(index);
        end

        fclose(fileID);

    % liquid compositions

        outputID= fullfile(Alpha_output, 'Liquid_comp_tbl.txt');
        fileID = fopen(outputID);
        
        txt = textscan(fileID,'%s','delimiter','\n');
        output_header=split(string(txt{1}(4)));
        output_value=str2num(string(txt{1}(5)))';
        
        ME_melt=zeros(numel(lst.ME_oxlist),1);
        for n=1:lst.nME
            index = find(contains(output_header,lst.ME_oxlist(n)));
            ME_melt(n)=output_value(index);
        end

        fclose(fileID);
    
    fclose('all');

else

    exitFlag = 2;

    % fallback
    realP = NaN;
    realT = NaN;
    realF = NaN;
    ME_melt = NaN(lst.nME,1);
    ME_residue = NaN(numel(lst.ME_oxlist),1);
    mode = NaN(lst.nMineral,1);   

end
    
end