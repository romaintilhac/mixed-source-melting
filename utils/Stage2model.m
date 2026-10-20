function  [param, compo, exitFlag] = Stage2model(param, compo, Kd, lst, pth, option, settings)
% Models the compositions of melts and residues from a three-lithology, mixed-mantle source in a 1D column simulating a triangular melting regime

struct_oxi=getTableOxi; TableOxi=struct_oxi.TableOxi; clear struct_oxi
struct_iso=getTableIso; TableIso=struct_iso.TableIso; clear struct_iso

%% Load parameters and compositions

unpackStruct(lst, "nMineral nME nTE nCE nisoTE nisoCE nSE minerals ME TE CE isoTE isoTE_index isoCE isoCE_index SE ME_oxlist", 0, 1);
unpackStruct(settings, "dP Prop_round F_thresh PropPyr_thresh fine_nodes coarse_nodes", 0, 0);
unpackStruct(param, "Tp PropPyr PropEcl2 Ptop Ur Xref PropCpx2 Pmax Pmin Pchange Pshift S_Pdt S_Pyr S_Ecl PropMSS_Pdt PropMSS_Pyr PropMSS_Ecl phi2 fixSCSS resSsil capSCSS_Ecl_P capSCSS_Ecl_F capSCSS_Ecl",0 , 0);
if nME>0; unpackStruct(compo, "ME_Pdt ME_Pyr ME_Ecl", 0, 1); end
if nTE>0; unpackStruct(compo, "TE_Pdt TE_Pyr TE_Ecl", 1, 1); end
if nCE>0; unpackStruct(compo, "CE_Pdt CE_Pyr CE_Ecl", 1, 1); end
if nisoTE>0; unpackStruct(compo, "isoTE_Pdt isoTE_Pyr isoTE_Ecl", 1, 1); end
if nisoCE>0; unpackStruct(compo, "isoCE_Pdt isoCE_Pyr isoCE_Ecl", 1, 1); end
unpackStruct(compo,"ME_PdtRef ME_PyrRef ME_EclRef SE_Pdt SE_Pyr SE_Ecl meltingmode_Pdt meltingmode_Pyr melting_Pdt melting_Pyr meltingmode_Ecl melting_Ecl", 0, 0);

param = struct(); compo = struct();

PropPdt=1-(PropPyr+PropEcl2);

%% Run Melt-PX

% hash
regime_string=[num2str(Tp) '-' num2str(Pmin) '-' num2str(Pmax) '-' num2str(dP) 'GPa_' num2str(PropCpx2) 'cpx'];
if PropEcl2>0; Ecl_string=['+' num2str(round(PropEcl2 / Prop_round) * Prop_round) '-' ME_EclRef]; else; Ecl_string=''; end
if PropPyr>PropPyr_thresh; Pyr_string=['+' num2str(round(PropPyr / Prop_round) * Prop_round) '-' ME_PyrRef]; else; Pyr_string=''; end
Pdt_string=[ME_PdtRef];
ouputIDcommon=[regime_string '_' Pdt_string Pyr_string Ecl_string '.mat'];
outputID= fullfile(pth.MeltPX_output, ouputIDcommon);

% run or load
if option.MeltPX==0 || fopen(outputID)==-1
    [Pressure,Temperature, MeltingD_Pdt,MeltingD_Pyr, MeltingD_Ecl, Tadiabat, Tsol_Pdt, Tsol_Pyr, Tsol_Ecl]=MeltPXdecompression3(Tp,Pmax, Pmin, dP, PropCpx2, PropPyr, PropEcl2, ME_Pyr, ME_Ecl);
    save(outputID,'Pressure','Temperature','Pressure','Temperature','MeltingD_Pdt','Tadiabat','MeltingD_Pdt','MeltingD_Pyr','MeltingD_Ecl', 'Tsol_Pdt', 'Tsol_Pyr', 'Tsol_Ecl');
    disp('Stage-2 Melt-PX calculations completed.')
else
    load(outputID);
    disp('Stage-2 Melt-PX outputs successfully loaded.')
end

%% Find solidus depths

% find deepest solidus
Pressure=(round(Pressure.*100))./100;
i0Pdt=find(MeltingD_Pdt>0,1); % find the point at which peridotite starts to melt
i0Pyr=find(MeltingD_Pyr>0,1); % find the point at which pyroxenite starts to melt
if PropEcl2>0
    i0Ecl=find(MeltingD_Ecl>0,1); % find the point at which eclogite starts to melt
    if i0Ecl<i0Pyr; i0=i0Ecl; else;  i0=i0Pyr; end
else
     i0=i0Pyr;
end
if i0Pdt<i0; i0=i0Pdt; end

% constrain to deepest solidus-i1
i1=1; %pressure increment below solidus (must be >=1 to allow calculation of compositions for first increment)
Pressure=Pressure(i0-i1:end);
Temperature=Temperature(i0-i1:end);
MeltingD_Pdt=MeltingD_Pdt(i0-i1:end)./100;
MeltingD_Pyr=MeltingD_Pyr(i0-i1:end)./100;
Tadiabat=Tadiabat(i0-i1:end);
Tsol_Pdt=Tsol_Pdt(i0-i1:end);
Tsol_Pyr=Tsol_Pyr(i0-i1:end);
Tsol_Ecl=Tsol_Ecl(i0-i1:end);
if PropEcl2>0; MeltingD_Ecl=MeltingD_Ecl(i0-i1:end)./100;end
ifull=length(Pressure); % full length of Melt-PX outputs above deepest solidus

% redefine i0 and itop in the cropped vectors
i0Pdt=find(MeltingD_Pdt>0,1); % find the point at which peridotite starts to melt
i0Pyr=find(MeltingD_Pyr>0,1); % find the point at which pyroxenite starts to melt
if PropEcl2>0 && ~isempty(find(MeltingD_Ecl>0,1)) % find the point at which eclogite starts to melt
    i0Ecl=find(MeltingD_Ecl>0,1);
else
    i0Ecl=ifull; %to handle refractory eclogite residue
end

itopPdt=find(MeltingD_Pdt>0,1, 'last'); % find the point at which peridotite stops melting
itopPyr=find(MeltingD_Pyr>0,1, 'last'); % find the point at which pyroxenite stops melting
if PropEcl2>0
    if ~isempty(find(MeltingD_Ecl>0,1, 'last'))
    itopEcl=find(MeltingD_Ecl>0,1, 'last');
    else
        itopEcl=ifull; % find the point at which eclogite stops melting
    end
end
% get solidus depths
P0Pdt=Pressure(i0Pdt); %depth of peridotite solidus
if PropEcl2>0; P0Ecl=Pressure(i0Ecl); end %depth of eclogite solidus
P0Pyr=Pressure(i0Pyr); %depth of pyroxenite solidus

%% Define nodes for pMELTS and SCSS calculations
    
ithreshold_Pdt=find(MeltingD_Pdt>F_thresh,1);
ithreshold_Pyr=find(MeltingD_Pyr>F_thresh,1);

node_list_Pdt=i0Pdt:fine_nodes:ithreshold_Pdt-fine_nodes;
node_list_Pdt=[node_list_Pdt i0Pdt+ithreshold_Pdt:coarse_nodes:itopPdt];
if PropPyr>0
    node_list_Pyr=i0Pyr:fine_nodes:ithreshold_Pyr-fine_nodes;
    node_list_Pyr=[node_list_Pyr ithreshold_Pyr:coarse_nodes:itopPyr];
end
if PropEcl2>0
    ithreshold_Ecl=find(MeltingD_Ecl>F_thresh,1);
    node_list_Ecl=i0Ecl:fine_nodes:ithreshold_Ecl-fine_nodes;
    node_list_Ecl=[node_list_Ecl ithreshold_Ecl:coarse_nodes:itopEcl];
end

node_list_Pdt=[node_list_Pdt itopPdt];
if PropPyr>0;node_list_Pyr=[node_list_Pyr itopPyr]; end
if PropEcl2>0; node_list_Ecl=[node_list_Ecl itopEcl]; end

% common nodes
node_Pdtbreak=find(ischange(MeltingD_Pdt,'linear','MaxNumChanges',3),2,'last');
if PropPyr>0; node_Pyrbreak=find(ischange(MeltingD_Pyr,'linear','MaxNumChanges',3),2,'last'); end
if PropEcl2>0; node_Eclbreak=find(ischange(MeltingD_Ecl,'linear','MaxNumChanges',3),2,'last'); end

for i=1:length(node_Pdtbreak)
    node=node_Pdtbreak(i);
    if node>i0Pdt&& node<itopPdt; node_Pdtbreak=[node_list_Pdt node];end
    if PropPyr>0 && node>i0Pyr && node<itopPyr; node_list_Pyr=[node_list_Pyr node];end
    if PropEcl2>0 && node>i0Ecl && node<itopEcl;node_list_Ecl=[node_list_Ecl node];end
end
if PropPyr>0
    for i=1:length(node_Pyrbreak)
        node=node_Pyrbreak(i);
        if node>i0Pyr && node<itopPyr; node_list_Pyr=[node_list_Pyr node];end
        if node>i0Pdt && node<itopPdt; node_list_Pdt=[node_list_Pdt node];end
        if PropEcl2>0 && node>i0Ecl && node<itopEcl;node_list_Ecl=[node_list_Ecl node];end
    end
end
if PropEcl2>0
    for i=1:length(node_Eclbreak)
        node=node_Eclbreak(i);
        if node>i0Ecl && node<itopEcl; node_list_Ecl=[node_list_Ecl node];end
        if node>i0Pdt && node<itopPdt; node_list_Pdt=[node_list_Pdt node];end
        if PropPyr>0 && node>i0Pyr && node<itopPyr;node_list_Pyr=[node_list_Pyr node];end
    end
end

if PropPyr>0 && i0Pdt>min(node_list_Pyr);node_list_Pyr=[node_list_Pyr i0Pdt];end
if PropEcl2>0 && i0Pdt>min(node_list_Ecl);node_list_Ecl=[node_list_Ecl i0Pdt]; end

if PropPyr>0
    if i0Pyr>min(node_list_Pdt);node_list_Pdt=[node_list_Pdt i0Pyr];end
    if PropEcl2>0 && i0Pyr>min(node_list_Ecl);node_list_Ecl=[node_list_Ecl i0Pyr];end
end
if PropEcl2>0
    if i0Ecl>min(node_list_Pdt);node_list_Pdt=[node_list_Pdt i0Ecl];end
    if PropPyr>0 && i0Ecl>min(node_list_Pyr);node_list_Pyr=[node_list_Pyr i0Ecl];end
end

node_list_Pdt=sort(unique(node_list_Pdt));
if PropPyr>0;node_list_Pyr=sort(unique(node_list_Pyr)); end
if PropEcl2>0; node_list_Ecl=sort(unique(node_list_Ecl)); end

%% Run pMELTS 

exitFlag = 1;

delete(fullfile(pth.Alpha_input,'*'))  %clear previous Alphamelts input files
delete(fullfile(pth.Alpha_batch,'*'))  %clear previous Alphamelts batch files

% hash
Pdt_string=['_' ME_PdtRef];
outputID_Pdt= fullfile(pth.pMELTS_output, ['Pdt_' ouputIDcommon]);
if PropPyr>0; outputID_Pyr= fullfile(pth.pMELTS_output, ['Pyr_' ouputIDcommon]); end
if PropEcl2>0; outputID_Ecl= fullfile(pth.pMELTS_output, ['Ecl_' ouputIDcommon]);end

% run or load
if option.pMELTS==0 || fopen(outputID_Pdt)==-1 || (PropPyr>0 && fopen(outputID_Pyr)==-1) || (PropEcl2>0 && fopen(outputID_Ecl)==-1) % if one is missing, all are rerun (to be changed?)
   
    % waitbar
    h3 = waitbar(0,'Setting up...', 'Position',[588 300 360 75],'windowstyle', 'modal');
    frames = java.awt.Frame.getFrames();frames(end).setAlwaysOnTop(1); 
    
    % pre-allocate nodes
    if option.stage2interpolation==1 % nodes only

        nodes_Pdt=length(node_list_Pdt); 
        realP_Pdt_nodes=zeros(nodes_Pdt,1);
        realT_Pdt_nodes=zeros(nodes_Pdt,1);
        realF_Pdt_nodes=zeros(nodes_Pdt,1);
        ME_meltPdt_nodes=zeros(nodes_Pdt,nME);
        ME_residuePdt_nodes=zeros(nodes_Pdt,nME);
        mode_Pdt_nodes=zeros(nodes_Pdt,nMineral);
        problem_Pdt_nodes=zeros(nodes_Pdt,3);
        exitFlag_Pdt = zeros(nodes_Pdt,1);
        j_Pdt=0;

        if PropPyr>0
            nodes_Pyr=length(node_list_Pyr);
            realP_Pyr_nodes=zeros(nodes_Pyr,1);
            realT_Pyr_nodes=zeros(nodes_Pyr,1);
            realF_Pyr_nodes=zeros(nodes_Pyr,1);
            ME_meltPyr_nodes=zeros(nodes_Pyr,nME);
            ME_residuePyr_nodes=zeros(nodes_Pyr,nME);
            mode_Pyr_nodes=zeros(nodes_Pyr,nMineral);
            problem_Pyr_nodes=zeros(nodes_Pyr,3);
            exitFlag_Pyr = zeros(nodes_Pdt,1);
            j_Pyr=0;
        end
        
        if PropEcl2>0
            nodes_Ecl=length(node_list_Ecl); 
            realP_Ecl_nodes=zeros(nodes_Ecl,1);
            realT_Ecl_nodes=zeros(nodes_Ecl,1);
            realF_Ecl_nodes=zeros(nodes_Ecl,1);
            ME_meltEcl_nodes=zeros(nodes_Ecl,nME);
            ME_residueEcl_nodes=zeros(nodes_Ecl,nME);
            mode_Ecl_nodes=zeros(nodes_Ecl,nMineral); 
            problem_Ecl_nodes=zeros(nodes_Ecl,3);
            exitFlag_Pdt = zeros(nodes_Ecl,1);
            j_Ecl=0;
        end

    else % all

        realP_Pdt=zeros(ifull,1);realT_Pdt=zeros(ifull,1);realF_Pdt=zeros(ifull,1);
        ME_meltPdt=zeros(ifull,nME);ME_residuePdt=zeros(ifull,nME);
        mode_Pdt=zeros(ifull,nMineral);
        problem_Pdt=zeros(ifull,1);
        exitFlag_Pdt = zeros(ifull,1);
        
        if PropPyr>0
            realP_Pyr=zeros(ifull,1);realT_Pyr=zeros(ifull,1);realF_Pyr=zeros(ifull,1);
            ME_meltPyr=zeros(ifull,nME);ME_residuePyr=zeros(ifull,nME);
            mode_Pyr=zeros(ifull,nMineral);       
            problem_Pyr=zeros(ifull,1);
            exitFlag_Pyr = zeros(ifull,1);
        end

        if PropEcl2>0
            realP_Ecl=zeros(ifull,1);realT_Ecl=zeros(ifull,1);realF_Ecl=zeros(ifull,1);
            ME_meltEcl=zeros(ifull,nME);ME_residueEcl=zeros(ifull,nME);
            mode_Ecl=zeros(ifull,nMineral);
            problem_Ecl=zeros(ifull,1);
            exitFlag_Ecl = zeros(ifull,1);
        end

    end
            
    disp('-- Stage-2 pMELTS calculations started --')
    if option.stage2interpolation==1 % compute nodes only
         for i=1:ifull
            if mod(round((i/ifull)*100),1)==0; waitbar(i/ifull,h3,['Running stage-2 pMELTS calculations: ' num2str(round(100*i/ifull)), '%']);end
            if ismember(i,node_list_Pdt)
                j_Pdt=j_Pdt+1;
                [realP_Pdt_nodes(j_Pdt), realT_Pdt_nodes(j_Pdt), realF_Pdt_nodes(j_Pdt), ME_meltPdt_nodes(j_Pdt,:),ME_residuePdt_nodes(j_Pdt,:), mode_Pdt_nodes(j_Pdt,:), problem_Pdt_nodes(j_Pdt,:), exitFlag_Pdt(j_Pdt)]=...
                    pMELTSexploration(Pressure(i),Temperature(i),MeltingD_Pdt(i),ME_Pdt, ['Pdt' num2str(i)], lst, pth, option);
            end
            if PropPyr>0 && ismember(i,node_list_Pyr)
                j_Pyr=j_Pyr+1;
                [realP_Pyr_nodes(j_Pyr), realT_Pyr_nodes(j_Pyr), realF_Pyr_nodes(j_Pyr), ME_meltPyr_nodes(j_Pyr,:),ME_residuePyr_nodes(j_Pyr,:), mode_Pyr_nodes(j_Pyr,:), problem_Pyr_nodes(j_Pyr,:), exitFlag_Pyr(j_Pyr)]=...
                    pMELTSexploration(Pressure(i),Temperature(i),MeltingD_Pyr(i),ME_Pyr, ['Pyr' num2str(i)], lst, pth, option);
            end
            if PropEcl2>0 && ismember(i,node_list_Ecl)
                j_Ecl=j_Ecl+1;
                [realP_Ecl_nodes(j_Ecl), realT_Ecl_nodes(j_Ecl), realF_Ecl_nodes(j_Ecl), ME_meltEcl_nodes(j_Ecl,:),ME_residueEcl_nodes(j_Ecl,:), mode_Ecl_nodes(j_Ecl,:), problem_Ecl_nodes(j_Ecl,:), exitFlag_Ecl(j_Ecl)]=...
                    pMELTSexploration(Pressure(i),Temperature(i),MeltingD_Ecl(i),ME_Ecl, ['Ecl' num2str(i)], lst, pth, option);
            end
         end
    else % compute all
        for i=1:ifull
            if mod(round((i/ifull)*100),1)==0; waitbar(i/ifull,h3,['Running stage-2 pMELTS calculations: ' num2str(round(100*i/ifull)), '%']);end
            if MeltingD_Pdt(i)>0
                [realP_Pdt(i), realT_Pdt(i), realF_Pdt(i), ME_meltPdt(i,:),ME_residuePdt(i,:), mode_Pdt(i,:), problem_Pdt(i), exitFlag_Pdt(i)]=...
                    pMELTSexploration(Pressure(i),Temperature(i),MeltingD_Pdt(i),ME_Pdt, ['Pdt' num2str(i)], lst, pth, option);
            end
            if PropPyr>0 && MeltingD_Pyr(i)>0
                [realP_Pyr(i), realT_Pyr(i), realF_Pyr(i), ME_meltPyr(i,:),ME_residuePyr(i,:), mode_Pyr(i,:), problem_Pyr(i), exitFlag_Pyr(i)]=...
                    pMELTSexploration(Pressure(i),Temperature(i),MeltingD_Pyr(i),ME_Pyr, ['Pyr' num2str(i)], lst, pth, option);
            end
            if PropEcl2>0 && MeltingD_Ecl(i)>0
                [realP_Ecl(i), realT_Ecl(i), realF_Ecl(i), ME_meltEcl(i,:),ME_residueEcl(i,:), mode_Ecl(i,:), problem_Ecl(i), exitFlag_Ecl(i)]=...
                    pMELTSexploration(Pressure(i),Temperature(i),MeltingD_Ecl(i),ME_Ecl, ['Ecl' num2str(i)], lst, pth, option);
            end
        end
    end

    % check pMELTs calculations status
    if any(exitFlag_Pdt == 2) || (PropPyr>0 && any(exitFlag_Pyr == 2)) || (PropEcl2>0 && any(exitFlag_Ecl == 2))
        exitFlag = 2;
        disp('!! Stage-2 pMELTS calculations failed !!')
    end
    close(h3)

    if exitFlag == 1

        % interpolate between nodes
        if option.stage2interpolation==1
    
            i_list=1:ifull;
            
            realP_Pdt=interp1(node_list_Pdt,realP_Pdt_nodes,i_list);
            realT_Pdt=interp1(node_list_Pdt,realT_Pdt_nodes,i_list);
            realF_Pdt=interp1(node_list_Pdt,realF_Pdt_nodes,i_list); 
            problem_Pdt=interp1(node_list_Pdt,problem_Pdt_nodes,i_list);
            ME_meltPdt = interp1(1:size(ME_meltPdt_nodes, 1), ME_meltPdt_nodes, linspace(1, size(ME_meltPdt_nodes, 1), ifull));
            ME_residuePdt = interp1(1:size(ME_residuePdt_nodes, 1), ME_residuePdt_nodes, linspace(1, size(ME_residuePdt_nodes, 1), ifull));
            mode_Pdt = interp1(1:size(mode_Pdt_nodes, 1), mode_Pdt_nodes, linspace(1, size(mode_Pdt_nodes, 1), ifull));
            
            if PropPyr>0 && i0Pyr~=ifull
                realP_Pyr=interp1(node_list_Pyr,realP_Pyr_nodes,i_list);
                realT_Pyr=interp1(node_list_Pyr,realT_Pyr_nodes,i_list);
                realF_Pyr=interp1(node_list_Pyr,realF_Pyr_nodes,i_list); 
                problem_Pyr=interp1(node_list_Pyr,problem_Pyr_nodes,i_list);
                ME_meltPyr = interp1(1:size(ME_meltPyr_nodes, 1), ME_meltPyr_nodes, linspace(1, size(ME_meltPyr_nodes, 1), ifull));
                ME_residuePyr = interp1(1:size(ME_residuePyr_nodes, 1), ME_residuePyr_nodes, linspace(1, size(ME_residuePyr_nodes, 1), ifull));
                mode_Pyr = interp1(1:size(mode_Pyr_nodes, 1), mode_Pyr_nodes, linspace(1, size(mode_Pyr_nodes, 1), ifull));
            else
                problem_Pyr=zeros(ifull,3);
                realP_Pyr=zeros(ifull,1);
                realT_Pyr=zeros(ifull,1);
                realF_Pyr=zeros(ifull,1);
                mode_Pyr=zeros(ifull,nMineral);
                ME_meltPyr=zeros(ifull, nME);
                ME_residuePyr=repmat(table2array(ME_Pyr),ifull,1);     
            end
    
            if PropEcl2>0 && i0Ecl~=ifull
                realP_Ecl=interp1(node_list_Ecl,realP_Ecl_nodes,i_list);
                realT_Ecl=interp1(node_list_Ecl,realT_Ecl_nodes,i_list);
                realF_Ecl=interp1(node_list_Ecl,realF_Ecl_nodes,i_list); 
                problem_Ecl=interp1(node_list_Ecl,problem_Ecl_nodes,i_list);
                ME_meltEcl = interp1(1:size(ME_meltEcl_nodes, 1), ME_meltEcl_nodes, linspace(1, size(ME_meltEcl_nodes, 1), ifull));
                ME_residueEcl = interp1(1:size(ME_residueEcl_nodes, 1), ME_residueEcl_nodes, linspace(1, size(ME_residueEcl_nodes, 1), ifull));
                mode_Ecl = interp1(1:size(mode_Ecl_nodes, 1), mode_Ecl_nodes, linspace(1, size(mode_Ecl_nodes, 1), ifull));
            elseif i0Ecl==ifull
                problem_Ecl=zeros(ifull,3);
                realP_Ecl=zeros(ifull,1);
                realT_Ecl=zeros(ifull,1);
                realF_Ecl=zeros(ifull,1);
                mode_Ecl=zeros(ifull,nMineral);
                ME_meltEcl=zeros(ifull, nME);
                ME_residueEcl=repmat(table2array(ME_Ecl),ifull,1);
            end
    
        end
       
        % build tables for saving
        ME_meltPdt=array2table(ME_meltPdt,VariableNames=ME_oxlist);
        ME_residuePdt=array2table(ME_residuePdt,VariableNames=ME_oxlist);
        mode_Pdt=array2table(mode_Pdt,VariableNames=minerals);
        ME_meltPyr=array2table(ME_meltPyr,VariableNames=ME_oxlist);
        ME_residuePyr=array2table(ME_residuePyr,VariableNames=ME_oxlist);
        mode_Pyr=array2table(mode_Pyr,VariableNames=minerals);
        if PropEcl2 > 0
            ME_meltEcl=array2table(ME_meltEcl,VariableNames=ME_oxlist);
            ME_residueEcl=array2table(ME_residueEcl,VariableNames=ME_oxlist);
            mode_Ecl=array2table(mode_Ecl,VariableNames=minerals);
        end
    
        save(outputID_Pdt,'realP_Pdt','realT_Pdt','realF_Pdt','ME_meltPdt','ME_residuePdt','mode_Pdt','problem_Pdt','node_list_Pdt');
        if PropPyr>0; save(outputID_Pyr,'realP_Pyr','realT_Pyr','realF_Pyr','ME_meltPyr','ME_residuePyr','mode_Pyr','problem_Pyr','node_list_Pyr'); end
        if PropEcl2>0; save(outputID_Ecl,'realP_Ecl','realT_Ecl','realF_Ecl','ME_meltEcl','ME_residueEcl','mode_Ecl','problem_Ecl','node_list_Ecl'); end
        disp('-- Stage-2 pMELTS calculations completed --')

    end

else

    load(outputID_Pdt);
    if PropPyr>0; load(outputID_Pyr); end
    if PropEcl2>0; load(outputID_Ecl); end
    disp('Stage-2 pMELTS outputs successfully loaded.')
    
end

%% Proceed
if exitFlag == 1
%% Crop top of the melting regime

%constrain Melt-PX outputs to itop
Ptop=(round(Ptop*100))/100;
if min(Pressure)<=Ptop
    itop=find(Pressure==Ptop); % find the top of melting column
    node_list_Pdt=node_list_Pdt(1:find(node_list_Pdt>itop,1)-1);
    node_list_Ecl=node_list_Ecl(1:find(node_list_Ecl>itop,1)-1);
    node_list_Pyr=node_list_Pyr(1:find(node_list_Pyr>itop,1)-1);
else
    itop=ifull;
end
P=Pressure(1:itop); 
T=Temperature(1:itop);
F_Pdt=MeltingD_Pdt(1:itop);
F_Ecl=MeltingD_Ecl(1:itop);
F_Pyr=MeltingD_Pyr(1:itop);

%constrain pMELTS outputs to itop
realP_Pdt=realP_Pdt(1:itop);
realT_Pdt=realT_Pdt(1:itop);
realF_Pdt=realF_Pdt(1:itop);
ME_meltPdt=ME_meltPdt(1:itop,:);
ME_residuePdt=ME_residuePdt(1:itop,:);
mode_Pdt=mode_Pdt(1:itop,:);
problem_Pdt=problem_Pdt(1:itop,3);

if PropPyr>0
    realP_Pyr=realP_Pyr(1:itop);
    realT_Pyr=realT_Pyr(1:itop);
    realF_Pyr=realF_Pyr(1:itop);
    ME_meltPyr=ME_meltPyr(1:itop,:);
    ME_residuePyr=ME_residuePyr(1:itop,:);
    mode_Pyr=mode_Pyr(1:itop,:);
    problem_Pyr=problem_Pyr(1:itop,3);
end

if PropEcl2>0
    realP_Ecl=realP_Ecl(1:itop);
    realT_Ecl=realT_Ecl(1:itop);
    realF_Ecl=realF_Ecl(1:itop);
    ME_meltEcl=ME_meltEcl(1:itop,:);
    ME_residueEcl=ME_residueEcl(1:itop,:);
    mode_Ecl=mode_Ecl(1:itop,:);
    problem_Ecl=problem_Ecl(1:itop,3);
end

%% Calculate volume, upwelling rate per depth increment

P_upwelling=P0Pdt+Pshift; % pressure at which relative upwelling velocity changes (GPa)
dmax=((max(P)-Ptop)*10^9)/(constants.rho_m*constants.grav);  % thickness (m) of the melting regime
dc=((Pchange-Ptop)*10^9)/(constants.rho_m*constants.grav); %thickness (m) of the triangular melting regime
dh=((dP)*10^9)/(constants.rho_m*constants.grav); % change of z per pressure increment 

v=zeros(itop,1);
Urate=ones(itop,1);
Xrefrate=ones(itop,1).*Xref;

for i=1:itop
    if P(i)>Pchange
        v(i)=(2*dh)/(dc+2*(dmax-dc));
    else
        if dmax>dc
            ichange=i-find(P==Pchange)+1; 
            v(i)=(2*dh*(dc-(ichange-0.5)*dh))/(dc*(dc+2*(dmax-dc))); 
        else
            v(i)=(2*dh*(dmax-(i-0.5)*dh))/(dmax^2); % to handle triangular melting regime only (at low Tp)
        end 
    end
    if P(i)>=P_upwelling
        Prange_i=P(i)-P_upwelling;
        Prange=P(1)-P_upwelling;
        Pratio=Prange_i/Prange;
        Urate(i)=1+Ur*(2*Pratio-Pratio^2);
    end
end

%% Calculate melt fractions from each depth increment

fPdt=zeros(itop,1);
fPdt_mass=zeros(itop,1);
FracPdtR=zeros(itop,1);
FracPdtsum=zeros(itop,1);
Pdt_contribR=zeros(itop,1);

if PropPyr>0
    fPyr=zeros(itop,1);
    fPyr_mass=zeros(itop,1);
    FracPyrR=zeros(itop,1);
    FracPyrsum=zeros(itop,1);
    Pyr_contribR=zeros(itop,1);
end

if PropEcl2>0
    fEcl=zeros(itop,1);
    fEcl_mass=zeros(itop,1);
    FracEclR=zeros(itop,1);
    FracEclsum=zeros(itop,1);
    Ecl_contribR=zeros(itop,1);
end

for i=1:itop
    if i==1
        fPdt(i)=F_Pdt(i);
        fPdt_mass(i)=F_Pdt(i);
        if PropPyr>0
            fPyr(i)=F_Pyr(i);
            fPyr_mass(i)=F_Pyr(i);
        end
        if PropEcl2>0
            fEcl(i)=F_Ecl(i);
            fEcl_mass(i)=F_Ecl(i);
        end
    else
        fPdt(i)=(F_Pdt(i)-F_Pdt(i-1))/(1-F_Pdt(i-1));
        fPdt_mass(i)=F_Pdt(i)-F_Pdt(i-1);
        if PropPyr>0
            fPyr(i)=(F_Pyr(i)-F_Pyr(i-1))/(1-F_Pyr(i-1));
            fPyr_mass(i)=F_Pyr(i)-F_Pyr(i-1);
        end
        if PropEcl2>0
            fEcl(i)=(F_Ecl(i)-F_Ecl(i-1))/(1-F_Ecl(i-1));
            fEcl_mass(i)=F_Ecl(i)-F_Ecl(i-1);
        end
    end
end

%replace negative values with zeros (to handle the last melt increment)
fPdt=max(fPdt,0);
fPdt_mass=max(fPdt_mass,0);
if PropPyr>0
    fPyr=max(fPyr,0);
    fPyr_mass=max(fPyr_mass,0);
end
if PropEcl2>0
    fEcl=max(fEcl,0);
    fEcl_mass=max(fEcl_mass,0);
end

%for final melt calculations
FracPdt=fPdt_mass.*v.*Urate;
if PropPyr>0; FracPyr=fPyr_mass.*v.*Urate; end
if PropEcl2>0; FracEcl=fEcl_mass.*v.*Urate; end

%for refertilization calculations only
for i=1:itop
    FracPdtsum(i)=sum(fPdt_mass(i0Pdt:i).*v(i0Pdt:i).*Urate(i0Pdt:i));
    FracPdtR(i)=FracPdt(i)/FracPdtsum(i);
    if PropPyr>0
        FracPyrsum(i)=sum(fPyr_mass(i0Pyr:i).*v(i0Pyr:i).*Urate(i0Pyr:i));
        FracPyrR(i)=FracPyr(i)/FracPyrsum(i);
    end
    if PropEcl2>0
        FracEclsum(i)=sum(fEcl_mass(i0Ecl:i).*v(i0Ecl:i).*Urate(i0Ecl:i));
        FracEclR(i)=FracEcl(i)/FracEclsum(i);
    end

    if option.ref==1

        sumR = sum(v(i0Pdt:i).*fPdt_mass(i0Pdt:i).*(PropPdt));
        if PropPyr>0; sumR = sumR + sum(v(i0Pyr:i).*fPyr_mass(i0Pyr:i).*(PropPyr)); end
        if PropEcl2>0; sumR = sumR + sum(v(i0Ecl:i).*fEcl_mass(i0Ecl:i).*(PropEcl2)); end
        
        Pdt_contribR(i)=sum(v(i0Pdt:i).*fPdt_mass(i0Pdt:i).*PropPdt)./ sumR;
        if PropPyr>0; Pyr_contribR(i)= sum(v(i0Pyr:i).*fPyr_mass(i0Pyr:i).*(PropPyr))./ sumR; end
        if PropEcl2>0; Ecl_contribR(i)=sum(v(i0Ecl:i).*fEcl_mass(i0Ecl:i).*(PropEcl2))./ sumR; end

    end

end

Pdt_contribR(isnan(Pdt_contribR))=0;
FracPdt=FracPdt./sum(FracPdt);
if PropPyr>0
    Pyr_contribR(isnan(Pyr_contribR))=0;
    FracPyr=FracPyr./sum(FracPyr);
end
if PropEcl2>0
    Ecl_contribR(isnan(Ecl_contribR))=0;
    FracEcl=FracEcl./sum(FracEcl);
end

%% Calculate melting reactions for TE

if nTE>0

tot=zeros(itop,1); 

for i=1:nMineral
    mineral=minerals(i);
    nPdt.(mineral)=meltingmode_Pdt.(mineral);
    rPdt.(mineral)=melting_Pdt.(mineral);
    if PropPyr>0
        nPyr.(mineral)=meltingmode_Pyr.(mineral);
        rPyr.(mineral)=melting_Pyr.(mineral);
    end
    if PropEcl2>0
        nEcl.(mineral)=meltingmode_Ecl.(mineral);
        rEcl.(mineral)=melting_Ecl.(mineral);
    end

    D.(mineral)=Kd.TE.(mineral)';
end

%peridotite
    PPdt=zeros(itop,1); DPdt=zeros(itop,nTE);
    %recalculate mode (Eq(2) Williams & Bizimis 2014)
    for i=1:nMineral
        mineral=minerals(i);
        n.(mineral)=nPdt.(mineral)-F_Pdt.*rPdt.(mineral);
        tot=tot+n.(mineral);
    end
    %recalculate to 100%
    for i=1:nMineral
        mineral=minerals(i);
        if n.(mineral)>0; n.(mineral)=n.(mineral)./tot;else; n.(mineral)=0;end
    end
    for i=1:nMineral
        mineral=minerals(i);
        PPdt=PPdt+D.(mineral).*rPdt.(mineral); %P of reaction    
        DPdt=DPdt+n.(mineral).*D.(mineral); %D0 of reaction
    end

%eclogite
    if PropEcl2>0
        PEcl=zeros(itop,1); DEcl=zeros(itop,nTE);
        %recalculate mode (Eq(2) Williams & Bizimis 2014)
        for i=1:nMineral
            mineral=minerals(i);
            n.(mineral)=nEcl.(mineral)-F_Ecl.*rEcl.(mineral);
            tot=tot+n.(mineral);
        end
        %recalculate to 100%
        for i=1:nMineral
            mineral=minerals(i);
            if n.(mineral)>0; n.(mineral)=n.(mineral)./tot;else; n.(mineral)=0;end
        end
        for i=1:nMineral
            mineral=minerals(i);
            PEcl=PEcl+D.(mineral).*rEcl.(mineral); %P of reaction    
            DEcl=n.(mineral).*D.(mineral); %D0 of reaction
        end
    end

%pyroxenite
    if PropPyr>0
        PPyr=zeros(itop,1);zeros(itop,nTE);
        %recalculate mode (Eq(2) Williams & Bizimis 2014)
        for i=1:nMineral
            mineral=minerals(i);
            n.(mineral)=nPyr.(mineral)-F_Pyr.*rPyr.(mineral);
            tot=tot+n.(mineral);
        end
        %recalculate to 100%
        for i=1:nMineral
            mineral=minerals(i);
            if n.(mineral)>0; n.(mineral)=n.(mineral)./tot;else; n.(mineral)=0;end
        end
        for i=1:nMineral
            mineral=minerals(i);
            PPyr=PPyr+D.(mineral).*rPyr.(mineral); %P of reaction    
            DPyr=n.(mineral).*D.(mineral); %D0 of reaction
        end
    end

end

%% Calculate SCSS and S contents

if option.stage2interpolation==1
    nodes_Pdt=length(node_list_Pdt); 
    SCSS_Pdt_nodes=zeros(nodes_Pdt,1);
    if PropPyr>0
        nodes_Pyr=length(node_list_Pyr);
        SCSS_Pyr_nodes=zeros(nodes_Pyr,1);
    end
    if PropEcl2>0
        nodes_Ecl=length(node_list_Ecl);
        SCSS_Ecl_nodes=zeros(nodes_Ecl,1);
    end
else
    SCSS_Pdt=zeros(itop,1);
    if PropPyr>0; SCSS_Pyr=zeros(itop,1); end
    if PropEcl2>0; SCSS_Ecl=zeros(itop,1); end
end
S_residuePdt=zeros(itop,1);S_meltPdt=zeros(itop,1);S_refmeltPdt=zeros(itop,1);
j_Pdt=0; 
if PropPyr>0
    S_residuePyr=zeros(itop,1);S_meltPyr=zeros(itop,1);S_refmeltPyr=zeros(itop,1);
    j_Pyr=0;
end
if PropEcl2>0
    S_residueEcl=zeros(itop,1);S_meltEcl=zeros(itop,1);S_refmeltEcl=zeros(itop,1);
    j_Ecl=0; 
end
S_refmelt=zeros(itop,1);

% SCSS
if option.SCSS==0
    SCSS_Pdt=ones(itop,1).*fixSCSS;
    if PropPyr>0; SCSS_Pyr=ones(itop,1).*fixSCSS; end
    if PropEcl2>0; SCSS_Ecl=ones(itop,1).*fixSCSS; end
else
    if option.stage2interpolation==1 % compute nodes only
        for i=1:itop
            if ismember(i,node_list_Pdt) 
                j_Pdt=j_Pdt+1;
                if i>=i0Pdt; [SCSS_Pdt_nodes(j_Pdt)]=calcSCSS(ME_meltPdt(i,:),SE_Pdt, P(i),T(i),option.SCSSref,option.SCSSideal, lst);end
            end
            if PropPyr>0 && ismember(i,node_list_Pyr) 
                j_Pyr=j_Pyr+1;
                 if i>=i0Pyr; [SCSS_Pyr_nodes(j_Pyr)]=calcSCSS(ME_meltPyr(i,:),SE_Pyr, P(i),T(i),option.SCSSref,option.SCSSideal, lst);end
            end
            if PropEcl2>0 && ismember(i,node_list_Ecl) 
                j_Ecl=j_Ecl+1;
                 if i>=i0Ecl; [SCSS_Ecl_nodes(j_Ecl)]=calcSCSS(ME_meltEcl(i,:),SE_Ecl, P(i),T(i),option.SCSSref,option.SCSSideal, lst);end
            end
        end
    else %compute all
        for i=1:itop
            if i>=i0Pdt; [SCSS_Pdt(i)]=calcSCSS(ME_meltPdt(i,:),SE_Pdt, P(i),T(i),option.SCSSref,option.SCSSideal, lst);end
            if PropEcl2>0 && i>=i0Ecl; [SCSS_Ecl(i)]=calcSCSS(ME_meltEcl(i,:),SE_Ecl, P(i),T(i),option.SCSSref,option.SCSSideal, lst);end
            if PropPyr>0 && i>=i0Pyr; [SCSS_Pyr(i)]=calcSCSS(ME_meltPyr(i,:),SE_Pyr, P(i),T(i),option.SCSSref,option.SCSSideal, lst);end
        end
    end
end

% interpolated between nodes
if option.stage2interpolation && option.SCSS
    i_list=1:itop;
    SCSS_Pdt=interp1(node_list_Pdt,SCSS_Pdt_nodes',i_list)';
    SCSS_Pdt(isnan(SCSS_Pdt))=0;
    if PropPyr>0
        SCSS_Pyr=interp1(node_list_Pyr,SCSS_Pyr_nodes',i_list)';
        SCSS_Pyr(isnan(SCSS_Pyr))=0;
    else
        SCSS_Pyr=zeros(itop,1);
    end
    if PropEcl2>0 && i0Ecl~=ifull
        SCSS_Ecl=interp1(node_list_Ecl,SCSS_Ecl_nodes',i_list)';
        SCSS_Ecl(isnan(SCSS_Ecl))=0;
    elseif i0Ecl==ifull
        SCSS_Ecl=zeros(itop,1);
    end
end

% cap SCSS to handle deep Ti-rich eclogite melt
if PropEcl2>0
    if capSCSS_Ecl > 0
        SCSS_Ecl(P > capSCSS_Ecl_P & F_Ecl < capSCSS_Ecl_F & F_Ecl > 0) = capSCSS_Ecl;
    end
end

for i=1:itop

    %peridotite
    if i<i0Pdt %below solidus
        S_residuePdt(i)=S_Pdt;
        S_meltPdt(i)=0;
        S_refmeltPdt(i)=0;
    elseif i==i0Pdt %at solidus
        %starting composition is previous residue
        S_residuePdt(i)=S_residuePdt(i-1);
        %saturate melt if sufficient S is available
        if S_residuePdt(i)>=fPdt(i)*SCSS_Pdt(i)
            S_meltPdt(i)=SCSS_Pdt(i);
            S_residuePdt(i)=(S_residuePdt(i)-S_meltPdt(i).*fPdt(i))./(1-fPdt(i));
        else
            S_residuePdt(i)=resSsil;
            S_meltPdt(i) = (S_residuePdt(i-1) - (1-fPdt(i))*S_residuePdt(i)) / fPdt(i);
        end
        S_refmeltPdt(i)=S_meltPdt(i);
    elseif i>i0Pdt %above solidus
        %calculate refertilized starting composition from previous residue
        S_refmeltPdt(i)=sum(S_meltPdt(i0Pdt:i-1).*FracPdtR(i0Pdt:i-1))/sum(FracPdtR(i0Pdt:i-1));
        if option.ref && option.refS
            S_residuePdt(i)=(1-Xrefrate(i-1))*S_residuePdt(i-1)+Xrefrate(i-1)*S_refmelt(i-1);
        else
            S_residuePdt(i)=S_residuePdt(i-1);
        end
        %saturate melt if sufficient S is available
        if S_residuePdt(i)>=fPdt(i)*SCSS_Pdt(i)
            S_meltPdt(i)=SCSS_Pdt(i);
            S_residuePdt(i)=(S_residuePdt(i)-S_meltPdt(i).*fPdt(i))./(1-fPdt(i));
        else
            S_residuePdt(i)=resSsil;
            S_meltPdt(i) = (S_residuePdt(i-1) - (1-fPdt(i))*S_residuePdt(i)) / fPdt(i);
        end
    end

    %pyroxenite
    if PropPyr>0
        if i<i0Pyr %below solidus
            S_residuePyr(i)=S_Pyr;
            S_meltPyr(i)=0;
            S_refmeltPyr(i)=0;
        elseif i==i0Pyr %at solidus
            %starting composition is previous residue
            S_residuePyr(i)=S_residuePyr(i-1);
            %saturate melt if sufficient S is available
            if S_residuePyr(i)>=fPyr(i)*SCSS_Pyr(i)
                S_meltPyr(i)=SCSS_Pyr(i);
                S_residuePyr(i)=(S_residuePyr(i)-S_meltPyr(i).*fPyr(i))./(1-fPyr(i));
            else
                S_residuePyr(i)=resSsil;
                S_meltPyr(i) = (S_residuePyr(i-1) - (1-fPyr(i))*S_residuePyr(i)) / fPyr(i);
            end
            S_refmeltPyr(i)=S_meltPyr(i);
        elseif i>i0Pyr %above solidus
            %calculate refertilized starting composition from previous residue
            S_refmeltPyr(i)=sum(S_meltPyr(i0Pyr:i-1).*FracPyrR(i0Pyr:i-1))/sum(FracPyrR(i0Pyr:i-1)); 
            if option.ref && option.refS
                S_residuePyr(i)=(1-Xrefrate(i-1))*S_residuePyr(i-1)+Xrefrate(i-1)*S_refmelt(i-1);
            else
                S_residuePyr(i)=S_residuePyr(i-1);
            end
            %saturate melt if sufficient S is available
            if S_residuePyr(i)>=fPyr(i)*SCSS_Pyr(i)
                S_meltPyr(i)=SCSS_Pyr(i);
                S_residuePyr(i)=(S_residuePyr(i)-S_meltPyr(i).*fPyr(i))./(1-fPyr(i));
            else
                S_residuePyr(i)=resSsil;
                S_meltPyr(i) = (S_residuePyr(i-1) - (1-fPyr(i))*S_residuePyr(i)) / fPyr(i);
            end
        end
    end
 
    %eclogite
    if PropEcl2>0
        if i<i0Ecl %below solidus
            S_residueEcl(i)=S_Ecl;
            S_meltEcl(i)=0;
            S_refmeltEcl(i)=0;
        elseif i==i0Ecl %at solidus
            %starting composition is previous residue
            S_residueEcl(i)=S_residueEcl(i-1);
            %saturate melt if sufficient S is available
            if S_residueEcl(i)>=fEcl(i)*SCSS_Ecl(i)
                S_meltEcl(i)=SCSS_Ecl(i);
                S_residueEcl(i)=(S_residueEcl(i)-S_meltEcl(i).*fEcl(i))./(1-fEcl(i));
            else
                S_residueEcl(i)=resSsil;
                S_meltEcl(i) = (S_residueEcl(i-1) - (1-fEcl(i))*S_residueEcl(i)) / fEcl(i);
            end
            S_refmeltEcl(i)=S_meltEcl(i);
        elseif i>i0Ecl %above solidus
            %calculate refertilized starting composition from previous residue
            S_refmeltEcl(i)=sum(S_meltEcl(i0Ecl:i-1).*FracEclR(i0Ecl:i-1))/sum(FracEclR(i0Ecl:i-1)); 
            if option.ref && option.refS
                S_residueEcl(i)=(1-Xrefrate(i-1))*S_residueEcl(i-1)+Xrefrate(i-1)*S_refmelt(i-1);
            else
                S_residueEcl(i)=S_residueEcl(i-1);
            end
            %saturate melt if sufficient S is available
            if S_residueEcl(i)>=fEcl(i)*SCSS_Ecl(i)
                S_meltEcl(i)=SCSS_Ecl(i);
                S_residueEcl(i)=(S_residueEcl(i)-S_meltEcl(i).*fEcl(i))./(1-fEcl(i));
            else
                S_residueEcl(i)=resSsil;
                S_meltEcl(i) = (S_residueEcl(i-1) - (1-fEcl(i))*S_residueEcl(i)) / fEcl(i);
            end
        end
    end

    if option.ref && option.refS
        S_refmelt(i) = S_refmeltPdt(i)*Pdt_contribR(i);
        if PropEcl2 > 0; S_refmelt(i) = S_refmelt(i) + S_refmeltEcl(i)*Ecl_contribR(i); end
        if PropPyr > 0; S_refmelt(i) = S_refmelt(i) + S_refmeltPyr(i)*Pyr_contribR(i); end
    end
end

%% Calculate sulfide/MSS mode and Se-Te partition coefficients

DPdtSul=zeros(itop,nCE);DPdtCESil=zeros(itop,nCE);DPdtCEMSS=zeros(itop,nCE);
if PropPyr>0; DPyrSul=zeros(itop,nCE);DPyrCESil=zeros(itop,nCE);DPyrCEMSS=zeros(itop,nCE); end
if PropEcl2>0; DEclSul=zeros(itop,nCE);DEclCESil=zeros(itop,nCE);DEclCEMSS=zeros(itop,nCE); end

m_Fe=TableOxi.Fe.mass;
m_Ni=TableOxi.Ni.mass;
m_Cu=TableOxi.Cu.mass;
m_O=TableOxi.O.mass;
m_S=TableOxi.S.mass;

if option.sul==1
    XsulFe=SE_Pdt.Fe/m_Fe; XsulNi=SE_Pdt.Ni/m_Ni;XsulCu=SE_Pdt.Cu/m_Cu;XsulO=SE_Pdt.O/m_O;XsulS=SE_Pdt.S/m_S;
    tot=XsulFe+XsulNi+XsulCu+XsulO+XsulS;
    SbulkSul_Pdt=XsulS/tot; % molar fraction
    SbulkSul_Pdt = SE_Pdt.S / (SE_Pdt.Fe + SE_Pdt.Ni + SE_Pdt.Cu + SE_Pdt.O + SE_Pdt.S); % mass fraction
    if PropPyr>0
        XsulFe=SE_Pyr.Fe/m_Fe; XsulNi=SE_Pyr.Ni/m_Ni;XsulCu=SE_Pyr.Cu/m_Cu;XsulO=SE_Pyr.O/m_O;XsulS=SE_Pyr.S/m_S;
        tot=XsulFe+XsulNi+XsulCu+XsulO+XsulS;
        SbulkSul_Pyr=XsulS/tot; % molar fraction
        SbulkSul_Pyr = SE_Pyr.S / (SE_Pyr.Fe + SE_Pyr.Ni + SE_Pyr.Cu + SE_Pyr.O + SE_Pyr.S); % mass fraction
    end
    if PropEcl2>0
        XsulFe=SE_Ecl.Fe/m_Fe; XsulNi=SE_Ecl.Ni/m_Ni;XsulCu=SE_Ecl.Cu/m_Cu;XsulO=SE_Ecl.O/m_O;XsulS=SE_Ecl.S/m_S;
        tot=XsulFe+XsulNi+XsulCu+XsulO+XsulS;
        SbulkSul_Ecl=XsulS/tot; % molar fraction
        SbulkSul_Ecl = SE_Ecl.S / (SE_Ecl.Fe + SE_Ecl.Ni + SE_Ecl.Cu + SE_Ecl.O + SE_Ecl.S); % mass fraction

    end
else
    SbulkSul_Pdt=m_S/(m_S+m_Fe);
    if PropPyr>0; SbulkSul_Pyr=m_S/(m_S+m_Fe); end
    if PropEcl2>0; SbulkSul_Ecl=m_S/(m_S+m_Fe); end
end

nPdtbulkSul = S_residuePdt ./ SbulkSul_Pdt ./ 1e6; % ppm to wt fraction
nPdtSil=1-nPdtbulkSul;
nPdtSul=nPdtbulkSul.*(1-PropMSS_Pdt);
nPdtMSS=nPdtbulkSul.*PropMSS_Pdt;
if PropPyr>0
    nPyrbulkSul = S_residuePyr ./ SbulkSul_Pyr ./ 1e6; % ppm to wt fraction
    nPyrSil=1-nPyrbulkSul;
    nPyrSul=nPyrbulkSul.*(1-PropMSS_Pyr);
    nPyrMSS=nPyrbulkSul.*PropMSS_Pyr;
end
if PropEcl2>0
    nEclbulkSul = S_residueEcl ./ SbulkSul_Ecl ./ 1e6; % ppm to wt fraction
    nEclSil=1-nEclbulkSul;
    nEclSul=nEclbulkSul.*(1-PropMSS_Ecl);
    nEclMSS=nEclbulkSul.*PropMSS_Ecl;
end

if option.SeTe==1 % Kd parametrization for Se and Te
    DPdtSul(:,1)=10.^(3.47-3.07*1e-2.*ME_meltPdt.FeO-9.13*1e-4.*ME_meltPdt.FeO.^2);
    DPdtSul(:,2)=10.^(4.3-2.58*1e-3.*ME_meltPdt.FeO.^2);
    if PropPyr>0
        DPyrSul(:,1)=10.^(3.47-3.07*1e-2.*ME_meltPyr.FeO-9.13*1e-4.*ME_meltPyr.FeO.^2);
        DPyrSul(:,2)=10.^(4.3-2.58*1e-3.*ME_meltPyr.FeO.^2);
    end
    if PropEcl2>0
        DEclSul(:,1)=10.^(3.47-3.07*1e-2.*ME_meltEcl.FeO-9.13*1e-4.*ME_meltEcl.FeO.^2);
        DEclSul(:,2)=10.^(4.3-2.58*1e-3.*ME_meltEcl.FeO.^2);
    end
else
    DPdtSul=ones(itop, 1) .* Kd.CE.Sul';
    if PropPyr>0; DPyrSul=ones(itop, 1) .* Kd.CE.Sul'; end
    if PropEcl2>0; DEclSul=ones(itop, 1) .* Kd.CE.Sul'; end
end

D_Sil=ones(itop, 1) .* Kd.CE.Sil';
D_MSS=ones(itop, 1) .* Kd.CE.MSS';
DPdtCE=nPdtSul.*DPdtSul+nPdtMSS.*D_MSS+nPdtSil.*D_Sil;
if PropPyr>0; DPyrCE=nPyrSul.*DPyrSul+nPyrMSS.*D_MSS+nPyrSil.*D_Sil; end
if PropEcl2>0; DEclCE=nEclSul.*DEclSul+nEclMSS.*D_MSS+nEclSil.*D_Sil; end

%% Calculate CE, TE and Se isotope compositions for each depth increment

if nTE>0
    TE_residuePdt=zeros(itop,nTE);TE_meltPdt=zeros(itop,nTE);TE_refmeltPdt=zeros(itop,nTE);
    if PropEcl2>0;TE_residueEcl=zeros(itop,nTE);TE_meltEcl=zeros(itop,nTE);TE_refmeltEcl=zeros(itop,nTE); end
    if PropPyr>0; TE_residuePyr=zeros(itop,nTE);TE_meltPyr=zeros(itop,nTE);TE_refmeltPyr=zeros(itop,nTE); end
    TE_refmelt=zeros(itop,nTE);
end
if nCE>0
    CE_residuePdt=zeros(itop,nCE);CE_meltPdt=zeros(itop,nCE);CE_refmeltPdt=zeros(itop,nCE);
    if PropEcl2>0; CE_residueEcl=zeros(itop,nCE);CE_meltEcl=zeros(itop,nCE);CE_refmeltEcl=zeros(itop,nCE); end
    if PropPyr>0; CE_residuePyr=zeros(itop,nCE);CE_meltPyr=zeros(itop,nCE);CE_refmeltPyr=zeros(itop,nCE); end
    CE_refmelt=zeros(itop,nCE);
end
if nisoTE>0
    isoTE_residuePdt=zeros(itop,nisoTE);isoTE_meltPdt=zeros(itop,nisoTE);isoTE_refmeltPdt=zeros(itop,nisoTE);
    if PropEcl2>0; isoTE_residueEcl=zeros(itop,nisoTE);isoTE_meltEcl=zeros(itop,nisoTE);isoTE_refmeltEcl=zeros(itop,nisoTE); end
    if PropPyr>0; isoTE_residuePyr=zeros(itop,nisoTE);isoTE_meltPyr=zeros(itop,nisoTE);isoTE_refmeltPyr=zeros(itop,nisoTE); end
    isoTE_refmelt=zeros(itop,nisoTE);
end
if nisoCE>0
    isoCE_residuePdt=zeros(itop,nisoCE);isoCE_meltPdt=zeros(itop,nisoCE);isoCE_refmeltPdt=zeros(itop,nisoCE);
    if PropEcl2>0; isoCE_residueEcl=zeros(itop,nisoCE);isoCE_meltEcl=zeros(itop,nisoCE);isoCE_refmeltEcl=zeros(itop,nisoCE); end
    if PropPyr>0; isoCE_residuePyr=zeros(itop,nisoCE);isoCE_meltPyr=zeros(itop,nisoCE);isoCE_refmeltPyr=zeros(itop,nisoCE); end
    isoCE_refmelt=zeros(itop,nisoCE);
end

for i=1:itop

    %peridotite

    if i<i0Pdt %below solidus

        if nTE>0; TE_residuePdt(i,:)=TE_Pdt(1,:); end
        if nisoTE>0; isoTE_residuePdt(i,:)=isoTE_Pdt(1,:); end
        if nCE>0; CE_residuePdt(i,:)=CE_Pdt(1,:); end
        if nisoCE>0, isoCE_residuePdt(i,:)=isoCE_Pdt(1,:); end

    elseif i==i0Pdt %at solidus

        %modal batch melting from previous residue for CE
         if nCE>0
            CE_residuePdt(i,:)=CE_residuePdt(i-1,:);
            [CE_meltPdt(i,:), ~]=meltingEquations(CE_residuePdt(i,:),fPdt(i),DPdtCE(i,:),phi2); 
            CE_refmeltPdt(i,:)=CE_meltPdt(i,:);
         end
        if nisoCE>0
            isoCE_residuePdt(i,:)=isoCE_residuePdt(i-1,:);
            isoCE_meltPdt(i,:)=isoCE_residuePdt(i-1,:);
            isoCE_refmeltPdt(i,:)=isoCE_meltPdt(i-1,:);
        end

        %non-modal batch melting from previous residue for TE
         if nTE >0
            TE_residuePdt(i,:)=TE_residuePdt(i-1,:);
            TE_meltPdt(i,:)=TE_residuePdt(i,:)./(DPdt(i,:)+fPdt(i).*(1-PPdt(1,:)));
            TE_residuePdt(i,:)=TE_residuePdt(i,:).*((DPdt(i,:)-PPdt(1,:).*fPdt(i))./(1-fPdt(i))).*(1./(DPdt(i,:)+fPdt(i).*(1-PPdt(1,:))));
            TE_refmeltPdt(i,:)=TE_meltPdt(i,:);
         end
         if nisoTE>0
            isoTE_residuePdt(i,:)=isoTE_residuePdt(i-1,:);
            isoTE_meltPdt(i,:)=isoTE_residuePdt(i-1,:);
            isoTE_refmeltPdt(i,:)=isoTE_meltPdt(i-1,:);
         end

    elseif i>i0Pdt %above solidus

        %modal batch melting from refertilized residue for CE
        if nCE>0; CE_refmeltPdt(i,:)=sum(CE_meltPdt(i0Pdt:i-1,:).*FracPdtR(i0Pdt:i-1),1)./sum(FracPdtR(i0Pdt:i-1)); end % the 1 argument forces the sum to work vertically(!) when the array is still 1xn
        if nisoCE>0
            isoCE_refmeltPdt(i,:)=sum(isoCE_meltPdt(i0Pdt:i-1,:).*CE_meltPdt(i0Pdt:i-1,isoCE_index(:)).*FracPdtR(i0Pdt:i-1), 1)...
                ./sum(FracPdtR(i0Pdt:i-1).*CE_meltPdt(i0Pdt:i-1,isoCE_index(:)),1);
        end

        if option.ref && option.refCE
            if nCE>0; CE_residuePdt(i,:)=(1-Xrefrate(i-1)).*CE_residuePdt(i-1,:)+Xrefrate(i-1).*CE_refmelt(i-1,:); end
            if nisoCE>0
                isoCE_residuePdt(i,:)=((1-Xrefrate(i-1)).*CE_residuePdt(i-1,isoCE_index(:)).*isoCE_residuePdt(i-1,:)+Xrefrate(i-1).*CE_refmelt(i-1,isoCE_index(:)).*isoCE_refmelt(i-1,:))...
                    ./((1-Xrefrate(i-1)).*CE_residuePdt(i-1,isoCE_index(:))+Xrefrate(i-1).*CE_refmelt(i-1,isoCE_index(:)));
            end
        else
            if nCE>0; CE_residuePdt(i,:)=CE_residuePdt(i-1,:); end
            if nisoCE>0; isoCE_residuePdt(i,:)=isoCE_residuePdt(i-1,:); end
        end

        if nCE>0; [CE_meltPdt(i,:), CE_residuePdt(i,:)]=meltingEquations(CE_residuePdt(i,:),fPdt(i),DPdtCE(i,:),phi2); end
        if nisoCE>0; isoCE_meltPdt(i,:)=isoCE_residuePdt(i-1,:); end

        %non-modal batch melting from refertilized residue for TE
        if nTE>0; TE_refmeltPdt(i,:)=sum(TE_meltPdt(i0Pdt:i-1,:).*FracPdtR(i0Pdt:i-1),1)./sum(FracPdtR(i0Pdt:i-1)); end % the 1 argument forces the sum to work vertically(!) when the array is still 1xn
        if nisoTE>0
            isoTE_refmeltPdt(i,:)=sum(isoTE_meltPdt(i0Pdt:i-1,:).*TE_meltPdt(i0Pdt:i-1,isoTE_index(:)).*FracPdtR(i0Pdt:i-1),1)...
                ./sum(FracPdtR(i0Pdt:i-1).*TE_meltPdt(i0Pdt:i-1,isoTE_index(:)),1);
        end
        if option.ref && option.refTE
            if nTE>0; TE_residuePdt(i,:)=(1-Xrefrate(i-1)).*TE_residuePdt(i-1,:)+Xrefrate(i-1).*TE_refmelt(i-1,:); end
            if nisoTE>0
                isoTE_residuePdt(i,:)=((1-Xrefrate(i-1)).*TE_residuePdt(i-1,isoTE_index(:)).*isoTE_residuePdt(i-1,:)+Xrefrate(i-1).*TE_refmelt(i-1,isoTE_index(:)).*isoTE_refmelt(i-1,:))...
                ./((1-Xrefrate(i-1)).*TE_residuePdt(i-1,isoTE_index(:))+Xrefrate(i-1).*TE_refmelt(i-1,isoTE_index(:)));
            end
        else
            if nTE>0; TE_residuePdt(i,:)=TE_residuePdt(i-1,:); end
            if nisoTE>0;  isoTE_residuePdt(i,:)=isoTE_residuePdt(i-1,:); end
        end
        if nTE>0
            TE_meltPdt(i,:)=TE_residuePdt(i,:)./(DPdt(i,:)+fPdt(i).*(1-PPdt(1,:)));
            TE_residuePdt(i,:)=TE_residuePdt(i,:).*((DPdt(i,:)-PPdt(1,:).*fPdt(i))./(1-fPdt(i))).*(1./(DPdt(i,:)+fPdt(i).*(1-PPdt(1,:))));
        end
        if nisoTE>0; isoTE_meltPdt(i,:)=isoTE_residuePdt(i-1,:); end
    end

    %pyroxenite
    if PropPyr>0
        if i<i0Pyr %below solidus
    
            if nTE>0; TE_residuePyr(i,:)=TE_Pyr(1,:); end
            if nisoTE>0; isoTE_residuePyr(i,:)=isoTE_Pyr(1,:); end
            if nCE>0; CE_residuePyr(i,:)=CE_Pyr(1,:); end
            if nisoCE>0, isoCE_residuePyr(i,:)=isoCE_Pyr(1,:); end
    
        elseif i==i0Pyr %at solidus
    
            %modal batch melting from previous residue for CE
             if nCE>0
                CE_residuePyr(i,:)=CE_residuePyr(i-1,:);
                [CE_meltPyr(i,:), ~]=meltingEquations(CE_residuePyr(i,:),fPyr(i),DPyrCE(i,:),phi2); 
                CE_refmeltPyr(i,:)=CE_meltPyr(i,:);
             end
            if nisoCE>0
                isoCE_residuePyr(i,:)=isoCE_residuePyr(i-1,:);
                isoCE_meltPyr(i,:)=isoCE_residuePyr(i-1,:);
                isoCE_refmeltPyr(i,:)=isoCE_meltPyr(i-1,:);
            end
    
            %non-modal batch melting from previous residue for TE
             if nTE >0
                TE_residuePyr(i,:)=TE_residuePyr(i-1,:);
                TE_meltPyr(i,:)=TE_residuePyr(i,:)./(DPyr(i,:)+fPyr(i).*(1-PPyr(1,:)));
                TE_residuePyr(i,:)=TE_residuePyr(i,:).*((DPyr(i,:)-PPyr(1,:).*fPyr(i))./(1-fPyr(i))).*(1./(DPyr(i,:)+fPyr(i).*(1-PPyr(1,:))));
                TE_refmeltPyr(i,:)=TE_meltPyr(i,:);
             end
             if nisoTE>0
                isoTE_residuePyr(i,:)=isoTE_residuePyr(i-1,:);
                isoTE_meltPyr(i,:)=isoTE_residuePyr(i-1,:);
                isoTE_refmeltPyr(i,:)=isoTE_meltPyr(i-1,:);
             end
    
        elseif i>i0Pyr %above solidus
    
            %modal batch melting from refertilized residue for CE
            if nCE>0; CE_refmeltPyr(i,:)=sum(CE_meltPyr(i0Pyr:i-1,:).*FracPyrR(i0Pyr:i-1),1)./sum(FracPyrR(i0Pyr:i-1)); end % the 1 argument forces the sum to work vertically(!) when the array is still 1xn
            if nisoCE>0
                isoCE_refmeltPyr(i,:)=sum(isoCE_meltPyr(i0Pyr:i-1,:).*CE_meltPyr(i0Pyr:i-1,isoCE_index(:)).*FracPyrR(i0Pyr:i-1),1)...
                    ./sum(FracPyrR(i0Pyr:i-1).*CE_meltPyr(i0Pyr:i-1,isoCE_index(:)),1);
            end
            if option.ref && option.refCE
                if nCE>0; CE_residuePyr(i,:)=(1-Xrefrate(i-1)).*CE_residuePyr(i-1,:)+Xrefrate(i-1).*CE_refmelt(i-1,:); end
                if nisoCE>0
                    isoCE_residuePyr(i,:)=((1-Xrefrate(i-1)).*CE_residuePyr(i-1,isoCE_index(:)).*isoCE_residuePyr(i-1,:)+Xrefrate(i-1).*CE_refmelt(i-1,isoCE_index(:)).*isoCE_refmelt(i-1,:))...
                        ./((1-Xrefrate(i-1)).*CE_residuePyr(i-1,isoCE_index(:))+Xrefrate(i-1).*CE_refmelt(i-1,isoCE_index(:)));
                end
            else
                if nCE>0; CE_residuePyr(i,:)=CE_residuePyr(i-1,:); end
                if nisoCE>0; isoCE_residuePyr(i,:)=isoCE_residuePyr(i-1,:); end
            end
            if nCE>0; [CE_meltPyr(i,:), CE_residuePyr(i,:)]=meltingEquations(CE_residuePyr(i,:),fPyr(i),DPyrCE(i,:),phi2); end
            if nisoCE>0; isoCE_meltPyr(i,:)=isoCE_residuePyr(i-1,:); end
    
            %non-modal batch melting from refertilized residue for TE
            if nTE>0; TE_refmeltPyr(i,:)=sum(TE_meltPyr(i0Pyr:i-1,:).*FracPyrR(i0Pyr:i-1),1)./sum(FracPyrR(i0Pyr:i-1)); end % the 1 argument forces the sum to work vertically(!) when the array is still 1xn
            if nisoTE>0
                isoTE_refmeltPyr(i,:)=sum(isoTE_meltPyr(i0Pyr:i-1,:).*TE_meltPyr(i0Pyr:i-1,isoTE_index(:)).*FracPyrR(i0Pyr:i-1),1)...
                    ./sum(FracPyrR(i0Pyr:i-1).*TE_meltPyr(i0Pyr:i-1,isoTE_index(:)),1);
            end
            if option.ref && option.refTE
                if nTE>0; TE_residuePyr(i,:)=(1-Xrefrate(i-1)).*TE_residuePyr(i-1,:)+Xrefrate(i-1).*TE_refmelt(i-1,:); end
                if nisoTE>0
                    isoTE_residuePyr(i,:)=((1-Xrefrate(i-1)).*TE_residuePyr(i-1,isoTE_index(:)).*isoTE_residuePyr(i-1,:)+Xrefrate(i-1).*TE_refmelt(i-1,isoTE_index(:)).*isoTE_refmelt(i-1,:))...
                    ./((1-Xrefrate(i-1)).*TE_residuePyr(i-1,isoTE_index(:))+Xrefrate(i-1).*TE_refmelt(i-1,isoTE_index(:)));
                end
            else
                if nTE>0; TE_residuePyr(i,:)=TE_residuePyr(i-1,:); end
                if nisoTE>0;  isoTE_residuePyr(i,:)=isoTE_residuePyr(i-1,:); end
            end
            if nTE>0
                TE_meltPyr(i,:)=TE_residuePyr(i,:)./(DPyr(i,:)+fPyr(i).*(1-PPyr(1,:)));
                TE_residuePyr(i,:)=TE_residuePyr(i,:).*((DPyr(i,:)-PPyr(1,:).*fPyr(i))./(1-fPyr(i))).*(1./(DPyr(i,:)+fPyr(i).*(1-PPyr(1,:))));
            end
            if nisoTE>0; isoTE_meltPyr(i,:)=isoTE_residuePyr(i-1,:); end
        end
    end

    %eclogite
    if PropEcl2>0
        if i<i0Ecl %below solidus

            if nTE>0; TE_residueEcl(i,:)=TE_Ecl(1,:); end
            if nisoTE>0; isoTE_residueEcl(i,:)=isoTE_Ecl(1,:); end
            if nCE>0; CE_residueEcl(i,:)=CE_Ecl(1,:); end
            if nisoCE>0, isoCE_residueEcl(i,:)=isoCE_Ecl(1,:); end

        elseif i==i0Ecl %at solidus

            %modal batch melting from previous residue for CE
             if nCE>0
                CE_residueEcl(i,:)=CE_residueEcl(i-1,:);
                [CE_meltEcl(i,:), ~]=meltingEquations(CE_residueEcl(i,:),fEcl(i),DEclCE(i,:),phi2); 
                CE_refmeltEcl(i,:)=CE_meltEcl(i,:);
             end
            if nisoCE>0
                isoCE_residueEcl(i,:)=isoCE_residueEcl(i-1,:);
                isoCE_meltEcl(i,:)=isoCE_residueEcl(i-1,:);
                isoCE_refmeltEcl(i,:)=isoCE_meltEcl(i-1,:);
            end

            %non-modal batch melting from previous residue for TE
             if nTE >0
                TE_residueEcl(i,:)=TE_residueEcl(i-1,:);
                TE_meltEcl(i,:)=TE_residueEcl(i,:)./(DEcl(i,:)+fEcl(i).*(1-PEcl(1,:)));
                TE_residueEcl(i,:)=TE_residueEcl(i,:).*((DEcl(i,:)-PEcl(1,:).*fEcl(i))./(1-fEcl(i))).*(1./(DEcl(i,:)+fEcl(i).*(1-PEcl(1,:))));
                TE_refmeltEcl(i,:)=TE_meltEcl(i,:);
             end
             if nisoTE>0
                isoTE_residueEcl(i,:)=isoTE_residueEcl(i-1,:);
                isoTE_meltEcl(i,:)=isoTE_residueEcl(i-1,:);
                isoTE_refmeltEcl(i,:)=isoTE_meltEcl(i-1,:);
             end

        elseif i>i0Ecl %above solidus

            %modal batch melting from refertilized residue for CE
            if nCE>0; CE_refmeltEcl(i,:)=sum(CE_meltEcl(i0Ecl:i-1,:).*FracEclR(i0Ecl:i-1),1)./sum(FracEclR(i0Ecl:i-1)); end % the 1 argument forces the sum to work vertically(!) when the array is still 1xn
            if nisoCE>0
                isoCE_refmeltEcl(i,:)=sum(isoCE_meltEcl(i0Ecl:i-1,:).*CE_meltEcl(i0Ecl:i-1,isoCE_index(:)).*FracEclR(i0Ecl:i-1),1)...
                    ./sum(FracEclR(i0Ecl:i-1).*CE_meltEcl(i0Ecl:i-1,isoCE_index(:)),1);
            end
            if option.ref && option.refCE
                if nCE>0; CE_residueEcl(i,:)=(1-Xrefrate(i-1)).*CE_residueEcl(i-1,:)+Xrefrate(i-1).*CE_refmelt(i-1,:); end
                if nisoCE>0
                    isoCE_residueEcl(i,:)=((1-Xrefrate(i-1)).*CE_residueEcl(i-1,isoCE_index(:)).*isoCE_residueEcl(i-1,:)+Xrefrate(i-1).*CE_refmelt(i-1,isoCE_index(:)).*isoCE_refmelt(i-1,:))...  
                        ./((1-Xrefrate(i-1)).*CE_residueEcl(i-1,isoCE_index(:))+Xrefrate(i-1).*CE_refmelt(i-1,isoCE_index(:)));
                end
            else
                if nCE>0; CE_residueEcl(i,:)=CE_residueEcl(i-1,:); end
                if nisoCE>0; isoCE_residueEcl(i,:)=isoCE_residueEcl(i-1,:); end
            end
            if nCE>0; [CE_meltEcl(i,:), CE_residueEcl(i,:)]=meltingEquations(CE_residueEcl(i,:),fEcl(i),DEclCE(i,:),phi2); end
            if nisoCE>0; isoCE_meltEcl(i,:)=isoCE_residueEcl(i-1,:); end

            %non-modal batch melting from refertilized residue for TE
            if nTE>0; TE_refmeltEcl(i,:)=sum(TE_meltEcl(i0Ecl:i-1,:).*FracEclR(i0Ecl:i-1),1)./sum(FracEclR(i0Ecl:i-1)); end % the 1 argument forces the sum to work vertically(!) when the array is still 1xn
            if nisoTE>0
                isoTE_refmeltEcl(i,:)=sum(isoTE_meltEcl(i0Ecl:i-1,:).*TE_meltEcl(i0Ecl:i-1,isoTE_index(:)).*FracEclR(i0Ecl:i-1),1)...
                    ./sum(FracEclR(i0Ecl:i-1).*TE_meltEcl(i0Ecl:i-1,isoTE_index(:)),1);
            end
            if option.ref && option.refTE
                if nTE>0; TE_residueEcl(i,:)=(1-Xrefrate(i-1)).*TE_residueEcl(i-1,:)+Xrefrate(i-1).*TE_refmelt(i-1,:); end
                if nisoTE>0
                    isoTE_residueEcl(i,:)=((1-Xrefrate(i-1)).*TE_residueEcl(i-1,isoTE_index(:)).*isoTE_residueEcl(i-1,:)+Xrefrate(i-1).*TE_refmelt(i-1,isoTE_index(:)).*isoTE_refmelt(i-1,:))...
                    ./((1-Xrefrate(i-1)).*TE_residueEcl(i-1,isoTE_index(:))+Xrefrate(i-1).*TE_refmelt(i-1,isoTE_index(:)));
                end
            else
                if nTE>0; TE_residueEcl(i,:)=TE_residueEcl(i-1,:); end
                if nisoTE>0;  isoTE_residueEcl(i,:)=isoTE_residueEcl(i-1,:); end
            end
            if nTE>0
                TE_meltEcl(i,:)=TE_residueEcl(i,:)./(DEcl(i,:)+fEcl(i).*(1-PEcl(1,:)));
                TE_residueEcl(i,:)=TE_residueEcl(i,:).*((DEcl(i,:)-PEcl(1,:).*fEcl(i))./(1-fEcl(i))).*(1./(DEcl(i,:)+fEcl(i).*(1-PEcl(1,:))));
            end
            if nisoTE>0; isoTE_meltEcl(i,:)=isoTE_residueEcl(i-1,:); end
        end
    end

    % aggregate all melts for refertilization
        
    % TE
    if nTE>0
        TE_refmelt(i,:) = TE_refmeltPdt(i,:).*Pdt_contribR(i);
        if PropEcl2>0; TE_refmelt(i,:) = TE_refmelt(i,:) + TE_refmeltEcl(i,:).*Ecl_contribR(i); end
        if PropPyr>0; TE_refmelt(i,:) = TE_refmelt(i,:) + TE_refmeltPyr(i,:).*Pyr_contribR(i); end
    end
    % CE
    if nCE>0
        CE_refmelt(i,:) = CE_refmeltPdt(i,:).*Pdt_contribR(i);
        if PropEcl2>0; CE_refmelt(i,:) = CE_refmelt(i,:) + CE_refmeltEcl(i,:).*Ecl_contribR(i); end
        if PropPyr>0; CE_refmelt(i,:) = CE_refmelt(i,:) + CE_refmeltPyr(i,:).*Pyr_contribR(i); end
    end
    
    % isoTE
    if nisoTE>0

        num_TE = isoTE_refmeltPdt(i,:).*TE_refmeltPdt(i,isoTE_index(:)).*Pdt_contribR(i);
        if PropPyr>0; num_TE = num_TE + isoTE_refmeltPyr(i,:).*TE_refmeltPyr(i,isoTE_index(:)).*Pyr_contribR(i); end
        if PropEcl2>0; num_TE = num_TE + isoTE_refmeltEcl(i,:).*TE_refmeltEcl(i,isoTE_index(:)).*Ecl_contribR(i); end

        den_TE = TE_refmeltPdt(i,isoTE_index(:)).*Pdt_contribR(i);
        if PropPyr>0; den_TE = den_TE + TE_refmeltPyr(i,isoTE_index(:)).*Pyr_contribR(i); end
        if PropEcl2>0; den_TE = den_TE + TE_refmeltEcl(i,isoTE_index(:)).*Ecl_contribR(i); end

        isoTE_refmelt(i,:)= num_TE ./ den_TE;

    end
    
    % isoCE
    if nisoCE>0

        num_CE = isoCE_refmeltPdt(i,:).*CE_refmeltPdt(i,isoCE_index(:)).*Pdt_contribR(i);
        if PropPyr>0; num_CE = num_CE + isoCE_refmeltPyr(i,:).*CE_refmeltPyr(i,isoCE_index(:)).*Pyr_contribR(i); end
        if PropEcl2>0; num_CE = num_CE + isoCE_refmeltEcl(i,:).*CE_refmeltEcl(i,isoCE_index(:)).*Ecl_contribR(i); end

        den_CE = CE_refmeltPdt(i,isoCE_index(:)).*Pdt_contribR(i);
        if PropPyr>0; den_CE = den_CE + CE_refmeltPyr(i,isoCE_index(:)).*Pyr_contribR(i); end
        if PropEcl2>0; den_CE = den_CE + CE_refmeltEcl(i,isoCE_index(:)).*Ecl_contribR(i); end

        isoCE_refmelt(i,:)= num_CE ./ den_CE;

    end

end

%% Calculate aggregated melt compositions

%pre-allocate

ME_aggmeltPdt=zeros(itop,nME);
if PropPyr>0; ME_aggmeltPyr=zeros(itop,nME); end
if PropEcl2>0; ME_aggmeltEcl=zeros(itop,nME); end

ME_meltPdtarray=table2array(ME_meltPdt);
if PropPyr>0;ME_meltPyrarray=table2array(ME_meltPyr); end
if PropEcl2>0;  ME_meltEclarray=table2array(ME_meltEcl); end

if nTE>0
    TE_meltPdt(1:i0Pdt-1,:)=missing; 
    TE_aggmeltPdt=zeros(itop,nTE);
    TE_finalmeltPdt=zeros(nTE);
    if PropPyr>0
        TE_meltPyr(1:i0Pyr-1,:)=missing;
        TE_aggmeltPyr=zeros(itop,nTE);  
        TE_finalmeltPyr=zeros(nTE);
    end
    if PropEcl2>0
        TE_meltEcl(1:i0Ecl-1,:)=missing;
        TE_aggmeltEcl=zeros(itop,nTE);
    end
end

if nCE>0
    CE_meltPdt(1:i0Pdt-1,:)=missing; 
    CE_aggmeltPdt=zeros(itop,nCE);
    CE_finalmeltPdt=zeros(nCE);
    if PropPyr>0
        CE_meltPyr(1:i0Pyr-1,:)=missing;
        CE_aggmeltPyr=zeros(itop,nCE);
        CE_finalmeltPyr=zeros(nCE);
    end
    if PropEcl2>0
        CE_meltEcl(1:i0Ecl-1,:)=missing;
        CE_aggmeltEcl=zeros(itop,nCE);
    end
end

if nisoTE>0
    isoTE_meltPdt(1:i0Pdt-1,:)=missing;
    isoTE_aggmeltPdt=zeros(itop,nisoTE);
    isoTE_finalmeltPdt=zeros(nisoTE);
    if PropPyr>0
        isoTE_meltPyr(1:i0Pyr-1,:)=missing;
        isoTE_aggmeltPyr=zeros(itop,nisoTE);
        isoTE_finalmeltPyr=zeros(nisoTE);
    end
    if PropEcl2>0
        isoTE_meltEcl(1:i0Ecl-1,:)=missing;
        isoTE_aggmeltEcl=zeros(itop,nisoTE);
        isoTE_finalmeltEcl=zeros(nisoTE);
    end
end

if nisoCE>0
    isoCE_meltPdt(1:i0Pdt-1,:)=missing;
    isoCE_aggmeltPdt=zeros(itop,nisoCE);
    isoCE_finalmeltPdt=zeros(nisoCE);
    if PropPyr>0
        isoCE_meltPyr(1:i0Pyr-1,:)=missing;
        isoCE_aggmeltPyr=zeros(itop,nisoCE);
        isoCE_finalmeltPyr=zeros(nisoCE);
    end
    if PropEcl2>0
        isoCE_meltEcl(1:i0Ecl-1,:)=missing;
        isoCE_aggmeltEcl=zeros(itop,nisoCE);
        isoCE_finalmeltEcl=zeros(nisoCE);
    end
end

S_aggmeltPdt=zeros(itop,1);
if PropPyr>0; S_aggmeltPyr=zeros(itop,1); end
if PropEcl2>0; S_aggmeltEcl=zeros(itop,1); end

%aggmelt increments
for i=1:itop

    ME_aggmeltPdt(i,:)=ME_meltPdtarray(i,:).*FracPdt(i);
    if nisoTE>0; TE_aggmeltPdt(i,:)=TE_meltPdt(i,:).*FracPdt(i); end
    if nisoCE>0; CE_aggmeltPdt(i,:)=CE_meltPdt(i,:).*FracPdt(i); end
    if nisoTE>0; isoTE_aggmeltPdt(i,:)=isoTE_meltPdt(i,:)*TE_meltPdt(i,isoTE_index(:))*FracPdt(i); end
    if nisoCE>0; isoCE_aggmeltPdt(i,:)=isoCE_meltPdt(i,:)*CE_meltPdt(i,isoCE_index)*FracPdt(i); end
    S_aggmeltPdt(i)=S_meltPdt(i)*FracPdt(i);
    
    if PropPyr>0
        ME_aggmeltPyr(i,:)=ME_meltPyrarray(i,:)*FracPyr(i);    
        if nisoTE>0; TE_aggmeltPyr(i,:)=TE_meltPyr(i,:).*FracPyr(i); end
        if nisoCE>0; CE_aggmeltPyr(i,:)=CE_meltPyr(i,:).*FracPyr(i); end
        if nisoTE>0; isoTE_aggmeltPyr(i,:)=isoTE_meltPyr(i,:)*TE_meltPyr(i,isoTE_index(:))*FracPyr(i); end
        if nisoCE>0; isoCE_aggmeltPyr(i,:)=isoCE_meltPyr(i,:)*CE_meltPyr(i,isoCE_index)*FracPyr(i); end
        S_aggmeltPyr(i)=S_meltPyr(i)*FracPyr(i);
    end

    if PropEcl2>0
        ME_aggmeltEcl(i,:)=ME_meltEclarray(i,:)*FracEcl(i);    
        if nisoTE>0; TE_aggmeltEcl(i,:)=TE_meltEcl(i,:).*FracEcl(i); end
        if nisoCE>0; CE_aggmeltEcl(i,:)=CE_meltEcl(i,:).*FracEcl(i); end
        if nisoTE>0; isoTE_aggmeltEcl(i,:)=isoTE_meltEcl(i,:)*TE_meltEcl(i,isoTE_index(:))*FracEcl(i); end
        if nisoCE>0; isoCE_aggmeltEcl(i,:)=isoCE_meltEcl(i,:)*CE_meltEcl(i,isoCE_index)*FracEcl(i); end
        S_aggmeltEcl(i)=S_meltEcl(i)*FracEcl(i);
    end

end

% calculate sum
ME_finalmeltPdt=sum(ME_aggmeltPdt(i0Pdt:end,:));
if PropPyr>0; ME_finalmeltPyr=sum(ME_aggmeltPyr(i0Pyr:end,:)); end
if PropEcl2>0; ME_finalmeltEcl=sum(ME_aggmeltEcl(i0Ecl:end,:)); end
S_finalmeltPdt=sum(S_aggmeltPdt(i0Pdt:end,:));
if PropPyr>0; S_finalmeltPyr=sum(S_aggmeltPyr(i0Pyr:end,:)); end
if PropEcl2>0; S_finalmeltEcl=sum(S_aggmeltEcl(i0Ecl:end,:)); end

%isnan=0
if nTE>0
    TE_meltPdt(isnan(TE_meltPdt))=0;
    TE_aggmeltPdt(isnan(TE_aggmeltPdt))=0;
    TE_finalmeltPdt=sum(TE_aggmeltPdt(i0Pdt:end,:));
    if PropPyr>0
        TE_meltPyr(isnan(TE_meltPyr))=0;
        TE_aggmeltPyr(isnan(TE_aggmeltPyr))=0;
        TE_finalmeltPyr=sum(TE_aggmeltPyr(i0Pyr:end,:));
    end
    if PropEcl2>0
        TE_meltEcl(isnan(TE_meltEcl))=0;
        TE_aggmeltEcl(isnan(TE_aggmeltEcl))=0;
        TE_finalmeltEcl=sum(TE_aggmeltEcl(i0Ecl:end,:));
    end
end

if nCE>0
    CE_meltPdt(isnan(CE_meltPdt))=0;
    CE_aggmeltPdt(isnan(CE_aggmeltPdt))=0;
    CE_finalmeltPdt=sum(CE_aggmeltPdt(i0Pdt:end,:));
    if PropPyr>0
        CE_meltPyr(isnan(CE_meltPyr))=0;
        CE_aggmeltPyr(isnan(CE_aggmeltPyr))=0;
        CE_finalmeltPyr=sum(CE_aggmeltPyr(i0Pyr:end,:));
    end
    if PropEcl2>0
        CE_meltEcl(isnan(CE_meltEcl))=0;
        CE_aggmeltEcl(isnan(CE_aggmeltEcl))=0;
        CE_finalmeltEcl=sum(CE_aggmeltEcl(i0Ecl:end,:));
    end
end

if nisoTE>0
    isoTE_aggmeltPdt(isnan(isoTE_aggmeltPdt))=0;
    isoTE_finalmeltPdt(:)=sum(isoTE_aggmeltPdt(i0Pdt:end,:))./sum(TE_meltPdt(i0Pdt:end,isoTE_index(:)).*FracPdt(i0Pdt:end));
    if PropPyr>0
        isoTE_aggmeltPyr(isnan(isoTE_aggmeltPyr))=0;
        isoTE_finalmeltPyr(:)=sum(isoTE_aggmeltPyr(i0Pyr:end,:))./sum(TE_meltPyr(i0Pyr:end,isoTE_index(:)).*FracPyr(i0Pyr:end));
    end
    if PropEcl2>0
        isoTE_aggmeltEcl(isnan(isoTE_aggmeltEcl))=0;
        isoTE_finalmeltEcl(:)=sum(isoTE_aggmeltEcl(i0Ecl:end,:))./sum(TE_meltEcl(i0Ecl:end,isoTE_index(:)).*FracEcl(i0Ecl:end));
    end
end

if nisoCE>0
    isoCE_aggmeltPdt(isnan(isoCE_aggmeltPdt))=0;
    isoCE_finalmeltPdt(:)=sum(isoCE_aggmeltPdt(i0Pdt:end,:))./sum(CE_meltPdt(i0Pdt:end,isoCE_index(:)).*FracPdt(i0Pdt:end));
    if PropPyr>0
        isoCE_aggmeltPyr(isnan(isoCE_aggmeltPyr))=0;
        isoCE_finalmeltPyr(:)=sum(isoCE_aggmeltPyr(i0Pyr:end,:))./sum(CE_meltPyr(i0Pyr:end,isoCE_index(:)).*FracPyr(i0Pyr:end));
    end
    if PropEcl2>0
        isoCE_aggmeltEcl(isnan(isoCE_aggmeltEcl))=0;
        isoCE_finalmeltEcl(:)=sum(isoCE_aggmeltEcl(i0Ecl:end,:))./sum(CE_meltEcl(i0Ecl:end,isoCE_index(:)).*FracEcl(i0Ecl:end));
    end
end

%% Calculate final melt compositions
    
    hasEcl = 0; hasPyr = 0;
    if PropEcl2 > 0 && ~any(isnan(ME_finalmeltEcl)); hasEcl = 1; end
    if PropPyr > 0 && ~any(isnan(ME_finalmeltPyr)); hasPyr = 1; end
    
    finalsum = sum(Urate .* v .* fPdt_mass .* PropPdt);
    if hasPyr; finalsum = finalsum + sum(Urate .* v .* fPyr_mass .* PropPyr); end
    if hasEcl; finalsum = finalsum + sum(Urate .* v .* fEcl_mass .* PropEcl2); end
    
    Pdt_contrib = sum(Urate .* v .* fPdt_mass .* PropPdt) / finalsum;
    if hasPyr
        Pyr_contrib = sum(Urate .* v .* fPyr_mass .* PropPyr) / finalsum;
    else
        Pyr_contrib = 0;
    end
    if hasEcl
        Ecl_contrib = sum(Urate .* v .* fEcl_mass * PropEcl2) / finalsum;
    else
        Ecl_contrib = 0;
    end
    
    
    ME_finalmelt = Pdt_contrib .* ME_finalmeltPdt;
    if hasPyr;  ME_finalmelt = ME_finalmelt + Pyr_contrib .* ME_finalmeltPyr; end
    if hasEcl; ME_finalmelt = ME_finalmelt + Ecl_contrib .* ME_finalmeltEcl; end
    
    if nTE > 0
        TE_finalmelt = Pdt_contrib .* TE_finalmeltPdt;
        if hasPyr;  TE_finalmelt = TE_finalmelt + Pyr_contrib .* TE_finalmeltPyr; end
        if hasEcl; TE_finalmelt = TE_finalmelt + Ecl_contrib .* TE_finalmeltEcl; end
    end
    
    if nCE > 0
        CE_finalmelt = Pdt_contrib .* CE_finalmeltPdt;
        if hasPyr;  CE_finalmelt = CE_finalmelt + Pyr_contrib .* CE_finalmeltPyr; end
        if hasEcl; CE_finalmelt = CE_finalmelt + Ecl_contrib .* CE_finalmeltEcl; end
    end
    
    if nisoTE > 0
        num = Pdt_contrib .* TE_finalmeltPdt(isoTE_index(:)) .* isoTE_finalmeltPdt;
        den = Pdt_contrib .* TE_finalmeltPdt(isoTE_index(:));
    
        if hasPyr
            num = num + Pyr_contrib .* TE_finalmeltPyr(isoTE_index(:)) .* isoTE_finalmeltPyr;
            den = den + Pyr_contrib .* TE_finalmeltPyr(isoTE_index(:));
        end
    
        if hasEcl
            num = num + Ecl_contrib .* TE_finalmeltEcl(isoTE_index(:)) .* isoTE_finalmeltEcl;
            den = den + Ecl_contrib .* TE_finalmeltEcl(isoTE_index(:));
        end
    
        isoTE_finalmelt = num ./ den;
    end
    
    
    if nisoCE > 0
        num = Pdt_contrib .* CE_finalmeltPdt(isoCE_index(:)) .* isoCE_finalmeltPdt;
        den = Pdt_contrib .* CE_finalmeltPdt(isoCE_index(:));
    
        if hasPyr
            num = num + Pyr_contrib .* CE_finalmeltPyr(isoCE_index(:)) .* isoCE_finalmeltPyr;
            den = den + Pyr_contrib .* CE_finalmeltPyr(isoCE_index(:));
        end
    
        if hasEcl
            num = num + Ecl_contrib .* CE_finalmeltEcl(isoCE_index(:)) .* isoCE_finalmeltEcl;
            den = den + Ecl_contrib .* CE_finalmeltEcl(isoCE_index(:));
        end
    
        isoCE_finalmelt = num ./ den;
    end
    
    S_finalmelt = Pdt_contrib .* S_finalmeltPdt;
    if hasPyr; S_finalmelt = S_finalmelt + Pyr_contrib .* S_finalmeltPyr; end
    if hasEcl; S_finalmelt = S_finalmelt + Ecl_contrib .* S_finalmeltEcl; end
    
%% Store output structures

    % parameters

    parameters_fields = split("v Xrefrate Urate P Pressure Temperature Tadiabat Pyr_contrib Ecl_contrib S_finalmelt");

    parameters_fields = [parameters_fields; split("nPdtSul DPdtCE S_residuePdt S_meltPdt SCSS_Pdt F_Pdt fPdt_mass fPdt FracPdt Pdt_contribR MeltingD_Pdt Tsol_Pdt realP_Pdt realT_Pdt realF_Pdt node_list_Pdt")];
    if PropPyr > 0; parameters_fields = [parameters_fields; split("nPyrSul DPyrCE S_residuePyr S_meltPyr SCSS_Pyr F_Pyr fPyr_mass fPyr FracPyr Pyr_contrib Pyr_contribR MeltingD_Pyr Tsol_Pyr realP_Pyr realT_Pyr realF_Pyr node_list_Pyr")]; end    
    if PropEcl2 > 0; parameters_fields = [parameters_fields; split("nEclSul DEclCE S_residueEcl S_meltEcl SCSS_Ecl F_Ecl fEcl_mass fEcl FracEcl Ecl_contrib Ecl_contribR MeltingD_Ecl Tsol_Ecl realP_Ecl realT_Ecl realF_Ecl node_list_Ecl")]; end
    
    param = addFields(param, parameters_fields, 1);

    % compositions

    if nME > 0
        compo = addFields(compo, "ME_finalmelt ME_residuePdt ME_meltPdt ME_finalmeltPdt", 1); % array
        if PropPyr > 0; compo = addFields(compo, "ME_residuePyr ME_meltPyr ME_finalmeltPyr", 1); end
        if PropEcl2 > 0; compo = addFields(compo, "ME_residueEcl ME_meltEcl ME_finalmeltEcl", 1); end
    end
    
    if nTE > 0
        compo = addFields(compo, "TE_finalmelt TE_residuePdt TE_meltPdt TE_finalmeltPdt", 0, TE); % table
        if PropPyr > 0; compo = addFields(compo, "TE_residuePyr TE_meltPyr TE_finalmeltPyr", 0, TE); end
        if PropEcl2 > 0; compo = addFields(compo, "TE_residueEcl TE_meltEcl TE_finalmeltEcl", 0, TE); end
    end
    
    if nCE > 0
        compo = addFields(compo, "CE_finalmelt CE_residuePdt CE_meltPdt CE_finalmeltPdt", 0, CE); % table
        if PropPyr > 0; compo = addFields(compo, "CE_residuePyr CE_meltPyr CE_finalmeltPyr", 0, CE); end
        if PropEcl2 > 0; compo = addFields(compo, "CE_residueEcl CE_meltEcl CE_finalmeltEcl", 0, CE); end
    end
    
    if nisoTE > 0
        compo = addFields(compo, "isoTE_finalmelt isoTE_residuePdt isoTE_meltPdt isoTE_finalmeltPdt", 0, isoTE); % table
        if PropPyr > 0; compo = addFields(compo, "isoTE_residuePyr isoTE_meltPyr isoTE_finalmeltPyr", 0, isoTE); end
        if PropEcl2 > 0; compo = addFields(compo, "isoTE_residueEcl isoTE_meltEcl isoTE_finalmeltEcl", 0, isoTE); end
    end
    
    if nisoCE > 0
        compo = addFields(compo, "isoCE_finalmelt isoCE_residuePdt isoCE_meltPdt isoCE_finalmeltPdt", 0, isoCE); % table
        if PropPyr > 0; compo = addFields(compo, "isoCE_residuePyr isoCE_meltPyr isoCE_finalmeltPyr", 0, isoCE); end
        if PropEcl2 > 0; compo = addFields(compo, "isoCE_residueEcl isoCE_meltEcl isoCE_finalmeltEcl", 0, isoCE); end
    end
    
    compo.mode_Pdt = mode_Pdt;
    if PropPyr > 0; compo.mode_Pyr = mode_Pyr; end
    if PropEcl2 > 0; compo.mode_Ecl = mode_Ecl; end
    
    compo = orderfields(compo);
    param = orderfields(param);

%% Plots

    if option.MonteCarlo==0
        plotStage2Inputs;
        plotStage2Outputs
    end

end

fclose('all');

end