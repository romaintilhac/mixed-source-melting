function [lst, option, settings, parameters, parameters_references, compo, compo_references, Kd, MonteCarlo, targets] = initializeSetup(setup_name, lst, pth)
% Reads the setup text files and initializes the model
    
struct=getTableOxi; TableOxi=struct.TableOxi;
    
folder_name=fullfile('setups', setup_name);

[option,~,~] = importStructure(fullfile(folder_name,'options.txt'));
[settings,~,~] = importStructure(fullfile(folder_name,'settings.txt'));
[parameters, parameters_references,~] = importStructure(fullfile(folder_name,'parameters.txt'));
[input_name,input_references,~] = importStructure(fullfile(folder_name,'sources.txt'));

if option.MonteCarlo==1
    [MonteCarlo,~,~] = importStructure(fullfile(folder_name,'montecarlo_ranges.txt'));
    [targets,~,~] = importStructure(fullfile(folder_name,'montecarlo_targets.txt'));
else
    MonteCarlo=[];
    targets=[];
end
    
%% Load compositions

lst.compositions = [string(fieldnames(input_name))];
lst.compositions = [lst.compositions; split("CE_mixPyr1 isoCE_mixPyr1 CE_finalmelt isoCE_finalmelt")];

% read database
fields=string(fieldnames(input_name));
for i=1:length(fields)
    field=fields(i);
    if ~contains(field,'Kd')
        subfield=split(field,'_'); subfield=subfield(1);
        warning off
        
        import_compositions = readtable(fullfile(pth.resources, 'inputs.xlsx'), 'Sheet', subfield);
        warning on
        row_ref = find(import_compositions{:,2}==input_references.(field));
        row_name = find(import_compositions{:,1}==input_name.(field));
        row=intersect(row_ref,row_name);

        import_compositions=import_compositions(row,:);
        if subfield=='ME'
            item_list=lst.ME_oxlist;
        elseif subfield=='meltingmode'  || subfield=='melting'
            item_list=lst.minerals;
        else
            item_list=lst.(subfield);
        end
        compo.(field)=import_compositions(:, contains(import_compositions.Properties.VariableNames,item_list));
        compo_references.(strcat(field,'Ref'))=strcat(char(import_compositions.name), '-',char(import_compositions.reference));
    else
    end
end

% calculate normative compositions
[parameters.normX_Pdt, parameters.normY_Pdt]=normProjection(compo.ME_Pdt);
[parameters.normX_Ecl, parameters.normY_Ecl]=normProjection(compo.ME_Ecl);
[parameters.normX_Pyr, parameters.normY_Pyr]=normProjection(compo.ME_Pyr);

% read Kd (standalone implementation)
%Kd.TE = readtable(fullfile(pth.resources, 'inputs.xlsx'), 'Sheet', 'TE_Kd');
Kd.CE = readtable(fullfile(pth.resources, 'inputs.xlsx'), 'Sheet', 'CE_Kd');
    
% calculate Fe2O3 where missing if a non-zero rFe3 was provided
litho_list=split("Pdt Pyr Ecl");
for i=length(litho_list)
    litho=string(litho_list(i));
    if isfield(parameters, strcat('rFe3_',litho)) && parameters.(strcat('rFe3_',litho))~=0
        if ~any(strcmp('Fe2O3',compo.(strcat('ME_',litho)).Properties.VariableNames)) || compo.(strcat('ME_',litho)).Fe2O3==0 || isnan(compo.(strcat('ME_',litho)).Fe2O3)
            rFe3=parameters.(strcat('rFe3_',litho));
            FeOtot=compo.(strcat('ME_', litho)).FeO;
            Fetot=FeOtot/TableOxi.Fe.mol2;
            Fe2=Fetot*(1-rFe3);
            Fe3=Fetot*rFe3;
            FeO=Fe2*TableOxi.Fe.mol2;
            Fe2O3=Fe3*TableOxi.Fe.mol3;
            compo.(strcat('ME_', litho)).FeO=FeO;
            compo.(strcat('ME_', litho)).Fe2O3=Fe2O3;
        end
    end
end

% replace missing ME data by settings.missingME
for n=1:numel(lst.ME_oxlist)
    ME=lst.ME_oxlist(n); if ~any(strcmp(ME,compo.ME_Pdt.Properties.VariableNames)) || compo.ME_Pdt.(ME)==0 || isnan(compo.ME_Pdt.(ME)); compo.ME_Pdt.(ME)=settings.missingME;end
    ME=lst.ME_oxlist(n); if ~any(strcmp(ME,compo.ME_Pyr.Properties.VariableNames)) || compo.ME_Pyr.(ME)==0 || isnan(compo.ME_Pyr.(ME)); compo.ME_Pyr.(ME)=settings.missingME;end
    ME=lst.ME_oxlist(n); if ~any(strcmp(ME,compo.ME_Ecl.Properties.VariableNames)) || compo.ME_Ecl.(ME)==0 || isnan(compo.ME_Ecl.(ME)); compo.ME_Ecl.(ME)=settings.missingME;end
end
    
%% Allocate Monte Carlo

if option.MonteCarlo==1

    iterations=settings.iterations;
    iterations_preadapt=settings.iterations_preadapt;
    inputs=string(fieldnames(MonteCarlo));
    min=[]; max=[]; mean=[]; stdev=[]; pd=[];
       
    for i=1:length(inputs)
        input = inputs(i);
        if ~(any(ismember(lst.compositions, input))) % parameters
        
            distrib=MonteCarlo.(input).distrib;

            if isfield(MonteCarlo.(input),'min');min=MonteCarlo.(input).min;end
            if isfield(MonteCarlo.(input),'max');max=MonteCarlo.(input).max;end
            if isfield(MonteCarlo.(input),'mean');mean=MonteCarlo.(input).mean;end
            if isfield(MonteCarlo.(input),'stdev');stdev=MonteCarlo.(input).stdev;end
            if isfield(MonteCarlo.(input),'pd');pd=MonteCarlo.(input).pd;end
            
            % generate random array
            % The following distributions can be defined for each randomized parameter:
            % 0 = uniform (min-max),
            % 1 = normal (based on min-max),
            % 2 = normal (based on the mean and min-max or mean and stdev)
            % 3 = custom (pd must be provided)
            MonteCarlo.(input).array0=randomArray(iterations_preadapt,distrib,min,max,mean,stdev,pd);
            
            % round inputs if specified
            if ~isnan(MonteCarlo.(input).round)
                MonteCarlo.(input).array0=MonteCarlo.(input).round.*round((MonteCarlo.(input).array0)./MonteCarlo.(input).round);
            end
    
        else % compositions
    
            fields=fieldnames(MonteCarlo.(input));
    
            for j=1:length(fields)
                field=string(fields(j));
    
                distrib=MonteCarlo.(input).(field).distrib;
                % to be updated with proper handling of optional arguments
                if isfield(MonteCarlo.(input).(field),'min');min=MonteCarlo.(input).(field).min;end
                if isfield(MonteCarlo.(input).(field),'max');max=MonteCarlo.(input).(field).max;end
                if isfield(MonteCarlo.(input).(field),'mean');mean=MonteCarlo.(input).(field).mean;end
                if isfield(MonteCarlo.(input).(field),'stdev');stdev=MonteCarlo.(input).(field).stdev;end
                if isfield(MonteCarlo.(input).(field),'pd');pd=MonteCarlo.(input).(field).pd;end
                
                % generate random array
                if min==max
                    MonteCarlo.(input).(field).array0=repmat(min,iterations_preadapt,1);
                else
                MonteCarlo.(input).(field).array0=randomArray(iterations_preadapt,distrib,min,max,mean,stdev,pd);
                end
                if ~isnan(MonteCarlo.(input).(field).round)
                    MonteCarlo.(input).(field).array0=MonteCarlo.(input).(field).round.*round((MonteCarlo.(input).(field).array0)./MonteCarlo.(input).(field).round);
                end
    
            end
    
        end
        
    end

end

end 