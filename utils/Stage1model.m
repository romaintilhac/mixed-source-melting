function [param_out, compo_out] = Stage1model(param, compo, Kd, lst, pth, option, settings)
% Models solid-state mixing of peridotite and eclogite or mixing of eclogite melt with peridotite melt and residue

%% Load parameters and compositions

struct_oxi=getTableOxi; TableOxi=struct_oxi.TableOxi; clear struct_oxi

unpackStruct(lst, "nMineral nME nTE nCE nisoTE nisoCE nSE minerals ME TE CE isoTE isoTE_index isoCE isoCE_index SE ME_oxlist", 0, 1)
unpackStruct(settings, "F_round", 0, 0);
if nME>0; unpackStruct(compo, "ME_Pdt ME_Ecl", 0); end
if nTE>0; unpackStruct(compo, "TE_Pdt TE_Ecl", 1); end % convert to array
if nCE>0; unpackStruct(compo, "CE_Pdt CE_Ecl", 1); end % convert to array  
if nisoTE>0; unpackStruct(compo, "isoTE_Pdt isoTE_Ecl", 1); end % convert to array
if nisoCE>0; unpackStruct(compo, "isoCE_Pdt isoCE_Ecl", 1); end % convert to array
unpackStruct(compo, "ME_PdtRef ME_EclRef SE_Pdt SE_Ecl meltingmode_Pdt melting_Pdt meltingmode_Ecl melting_Ecl", 0, 0) 

% handle stage1MeltPX
if option.stage1MeltPX==1
    param_fields=split("P1 T1_excess PropCpx1 PropEcl1 R1 phi1 S_Pdt S_Ecl PropMSS_Pdt1 PropMSS_Ecl1 fixSCSS resSsil");
else
    param_fields=split("P1 T1 F_Ecl1 F_Pdt1 PropCpx1 PropEcl1 R1 phi1 S_Pdt S_Ecl PropMSS_Pdt1 PropMSS_Ecl1 fixSCSS resSsil");
end
unpackStruct(param, param_fields, 0, 0)

%% Switch scenario

switch option.stage1reaction  
    
    case 1 %% melting-mixing

        %% Run Melt-PX (optional)
        
        if option.stage1MeltPX==1
    
            % calculate T from eclogite solidus
            Tsol_Ecl = eclogiteSolidus(P1, PropCpx1, ME_Ecl);
            T1 = Tsol_Ecl + T1_excess;
    
            % run Melt-PX
            [MeltingD_Pdt,MeltingD_Ecl]=MeltPXisobaric2(P1, T1, PropCpx1, ME_Ecl);
            disp('Stage-1 Melt-PX run successfully completed.')
            
            F_Pdt1=round(MeltingD_Pdt/100/F_round) * F_round;
            F_Ecl1=round(MeltingD_Ecl/100/F_round) * F_round;
    
        end
        
        %% Run pMELTS 
        
        if F_Ecl1==0; error('Eclogite must be above solidus for stage 1'); end
        
        delete(fullfile(pth.Alpha_input,'*')) 
        delete(fullfile(pth.Alpha_batch,'*')) 
        
        % hash
        common_string=[num2str(round(P1,1)) 'GPa'];
        Pdt_string=['F' num2str(F_Pdt1) '-' ME_PdtRef];
        Ecl_string=['F' num2str(F_Ecl1) '-' ME_EclRef];
        outputID_Pdt= fullfile(pth.pMELTS_output1, ['Pdt_' common_string '_' Pdt_string '.mat']);
        outputID_Ecl= fullfile(pth.pMELTS_output1, ['Ecl_' common_string '_' Ecl_string '.mat']);
        residueEcl_string=[common_string '_' Ecl_string];
        Pyr_string=[common_string '_' Pdt_string '_' num2str(round(PropEcl1,2))  '_' Ecl_string];
        
        if option.pMELTS==0 || fopen(outputID_Pdt)==-1 || fopen(outputID_Ecl)==-1
        
            disp('-------------------------- Stage-1 pMELTS calculations started  --------------------------')
            
            if F_Pdt1>0
                [realP_Pdt, realT1_Pdt, realF_Pdt, ME_meltPdt,ME_residuePdt, mode_Pdt, problem_Pdt]=pMELTSexploration(P1,T1,F_Pdt1,ME_Pdt,"Pdt", lst, pth, option);
                % build tables for saving
                ME_meltPdt=array2table(ME_meltPdt',VariableNames=ME_oxlist);
                ME_residuePdt=array2table(ME_residuePdt',VariableNames=ME_oxlist);
                mode_Pdt=array2table(mode_Pdt',VariableNames=minerals);
                save(outputID_Pdt,'realP_Pdt','realT1_Pdt','realF_Pdt','ME_meltPdt','ME_residuePdt','mode_Pdt','problem_Pdt');
            else
                ME_residuePdt=ME_Pdt;
                ME_meltPdt=array2table(zeros(1,nME),VariableNames=ME_oxlist);
                % handling previous version where realT1 wasn't saved
                if option.pMELTS==0
                    realT1_Pdt = nan;
                else
                    realT1_Pdt = T1;
                end
            end
        
            if F_Ecl1<1
               [realP_Ecl, realT1_Ecl, realF_Ecl, ME_meltEcl,ME_residueEcl, mode_Ecl, problem_Ecl]=pMELTSexploration(P1,T1,F_Ecl1,ME_Ecl,"Ecl", lst, pth, option);
               % build tables for saving 
                ME_meltEcl=array2table(ME_meltEcl',VariableNames=ME_oxlist);
                ME_residueEcl=array2table(ME_residueEcl',VariableNames=ME_oxlist);
                mode_Ecl=array2table(mode_Ecl',VariableNames=minerals);
                save(outputID_Ecl,'realP_Ecl','realT1_Ecl','realF_Ecl','ME_meltEcl','ME_residueEcl','mode_Ecl','problem_Ecl');
            else
                ME_meltEcl=ME_Ecl;
                % handling previous version where realT1 wasn't saved
                if option.pMELTS==0
                    realT1_Ecl = nan;
                else
                    realT1_Ecl = T1;
                end 
            end
            disp('-------------------------- Stage-1 pMELTS calculations completed  --------------------------')
        
        else
        
            load(outputID_Pdt); load(outputID_Ecl);
            disp('Stage-1 pMELTS outputs successfully loaded.')
            % reorder in case ME list has changed
            if F_Pdt1>0; ME_meltPdt=ME_meltPdt(:,ME_oxlist); end
            ME_residuePdt=ME_residuePdt(:,ME_oxlist);
            ME_meltEcl=ME_meltEcl(:,ME_oxlist);
            if F_Ecl1<1; ME_residueEcl=ME_residueEcl(:,ME_oxlist); end
    
            % handling previous version where realT1 wasn't saved
            if exist('realT_Pdt','var')
                realT1_Pdt = realT_Pdt;
                clear realT_Pdt
            end
            if exist('realT_Ecl','var')
                realT1_Ecl = realT_Ecl;
                clear realT_Ecl
            end
       
        end
        
        %% Calculate melting reactions for TE
        
        for i=1:nMineral
            mineral=minerals(i);
            nPdt.(mineral)=meltingmode_Pdt.(mineral);
            nEcl.(mineral)=meltingmode_Ecl.(mineral);
            rPdt.(mineral)=melting_Pdt.(mineral);
            rEcl.(mineral)=melting_Ecl.(mineral);
            D.(mineral)=Kd.TE.(mineral)';
        end
    
        %peridotite 
            if F_Pdt1>0
                %recalculate mode (Eq(2) Williams & Bizimis 2014)
                tot=0;
                for i=1:nMineral
                    mineral=minerals(i);
                    n.(mineral)=nPdt.(mineral)-F_Pdt1.*rPdt.(mineral);
                    tot=tot+n.(mineral);
                end
                %recalculate to 100%
                for i=1:nMineral
                    mineral=minerals(i);
                    if n.(mineral)>0; n.(mineral)=n.(mineral)./tot;else; n.(mineral)=0;end
                end
               PPdt=0;DPdt=0;
                for i=1:nMineral
                    mineral=minerals(i);
                    PPdt=PPdt+D.(mineral).*rPdt.(mineral); %P1 of reaction    
                    DPdt=n.(mineral).*D.(mineral); %D0 of reaction
                end
            end
        
        %eclogite
            %recalculate mode (Eq(2) Williams & Bizimis 2014)
            tot=0;
            for i=1:nMineral
                mineral=minerals(i);
                n.(mineral)=nEcl.(mineral)-F_Ecl1.*rEcl.(mineral);
                tot=tot+n.(mineral);
            end
            %recalculate to 100%
            for i=1:nMineral
                mineral=minerals(i);
                if n.(mineral)>0; n.(mineral)=n.(mineral)./tot;else; n.(mineral)=0;end
            end
           PEcl=0;DEcl=0;
            for i=1:nMineral
                mineral=minerals(i);
                PEcl=PEcl+D.(mineral).*rEcl.(mineral); %P1 of reaction    
                DEcl=n.(mineral).*D.(mineral); %D0 of reaction
            end
        
        %% Calculate SCSS and S contents
        
        SCSS_Pdt=0;
        SCSS_Ecl=0;
        
        % SCSS
        if option.SCSS==0
            SCSS_Pdt=fixSCSS;
            SCSS_Ecl=fixSCSS;
        else
            if F_Pdt1>0; [SCSS_Pdt]=calcSCSS(ME_meltPdt,SE_Pdt,P1,realT1_Pdt,option.SCSSref,option.SCSSideal, lst);end
            if F_Ecl1>0 && F_Ecl1<1;[SCSS_Ecl]=calcSCSS(ME_meltEcl,SE_Ecl, P1,realT1_Ecl,option.SCSSref,option.SCSSideal, lst);end
        end
        
        %peridotite
        if F_Pdt1==0 %below solidus
            S_residuePdt=S_Pdt;
            S_meltPdt=0;
        elseif F_Pdt1>0 %above solidus
            %saturate melt if sufficient S is available
            if S_Pdt>=F_Pdt1*SCSS_Pdt
                S_meltPdt=SCSS_Pdt;
                S_residuePdt=(S_Pdt-S_meltPdt*F_Pdt1)/(1-F_Pdt1);
            else
                S_residuePdt = resSsil;
                S_meltPdt = (S_Pdt - (1-F_Pdt1)*S_residuePdt) / F_Pdt1; % fixed mass balance (edit 2026)
            end
        end
        
        %eclogite
        if F_Ecl1<1
            if F_Ecl1==0 %below solidus
                S_residueEcl=S_Ecl;
                S_meltEcl=0;
            elseif F_Ecl1>0 %above solidus
                %saturate melt if sufficient S is available
                if S_Ecl>=F_Ecl1*SCSS_Ecl % saturated case
                    S_meltEcl=SCSS_Ecl;
                    S_residueEcl=(S_Ecl-S_meltEcl*F_Ecl1)/(1-F_Ecl1);
                else % undersaturated case
                    S_residueEcl = resSsil;
                    S_meltEcl = (S_Ecl - (1-F_Ecl1)*S_residueEcl) / F_Ecl1; % fixed mass balance (edit 2026)
                end
            end
        else
            S_meltEcl=S_Ecl;
            S_residueEcl=0;
        end
        
        %% Calculate sulfide/MSS mode and CE partition coefficients
        
        D_Sil=Kd.CE.Sil';
        D_MSS=Kd.CE.MSS';
        
        m_Fe=TableOxi.Fe.mass;
        m_Ni=TableOxi.Ni.mass;
        m_Cu=TableOxi.Cu.mass;
        m_O=TableOxi.O.mass;
        m_S=TableOxi.S.mass;
        
        %peridotite
        if F_Pdt1>0
            if option.sul==1
                XsulFe=SE_Pdt.Fe/m_Fe; XsulNi=SE_Pdt.Ni/m_Ni;XsulCu=SE_Pdt.Cu/m_Cu;XsulO=SE_Pdt.O/m_O;XsulS=SE_Pdt.S/m_S;
                tot=XsulFe+XsulNi+XsulCu+XsulO+XsulS;
                SbulkSul_Pdt=XsulS/tot; % molar fraction
                SbulkSul_Pdt = SE_Pdt.S / (SE_Pdt.Fe + SE_Pdt.Ni + SE_Pdt.Cu + SE_Pdt.O + SE_Pdt.S); % mass fraction (edit 2026)
            else
                SbulkSul_Pdt=m_S/(m_S+m_Fe);
            end
            
            nPdt.bulkSul = S_residuePdt ./ SbulkSul_Pdt ./ 1e6; % ppm to wt fraction fixed (edit 2026)
    
            nPdt.Sil=1-nPdt.bulkSul;
            nPdt.Sul=nPdt.bulkSul.*(1-PropMSS_Pdt1);
            nPdt.MSS=nPdt.bulkSul.*PropMSS_Pdt1;
            if option.SeTe==1 % custom for Se and Te (edit 2026)
                DPdtSul(1)=10.^(3.47-3.07*1e-2.*ME_meltPdt.FeO-9.13*1e-4.*ME_meltPdt.FeO.^2);
                DPdtSul(2)=10.^(4.3-2.58*1e-3.*ME_meltPdt.FeO.^2);
            else
                DPdtSul=Kd.CE.Sul';
            end
            DPdtCE=nPdt.Sul.*DPdtSul+nPdt.MSS.*D_MSS+nPdt.Sil.*D_Sil;
        end
        
        %eclogite
        if F_Ecl1<1
            if option.sul==1
                XsulFe=SE_Ecl.Fe/m_Fe; XsulNi=SE_Ecl.Ni/m_Ni;XsulCu=SE_Ecl.Cu/m_Cu;XsulO=SE_Ecl.O/m_O;XsulS=SE_Ecl.S/m_S;
                tot=XsulFe+XsulNi+XsulCu+XsulO+XsulS;
                SbulkSul_Ecl=XsulS/tot;
                SbulkSul_Ecl = SE_Ecl.S / (SE_Ecl.Fe + SE_Ecl.Ni + SE_Ecl.Cu + SE_Ecl.O + SE_Ecl.S); % mass fraction (edit 2026)
            else
                SbulkSul_Ecl=m_S/(m_S+m_Fe); 
            end
    
            nEcl.bulkSul = S_residueEcl ./ SbulkSul_Ecl ./ 1e6; % ppm to wt fraction fixed (edit 2026)
            
            nEcl.Sil=1-nEcl.bulkSul;
            nEcl.Sul=nEcl.bulkSul.*(1-PropMSS_Ecl1);
            nEcl.MSS=nEcl.bulkSul.*PropMSS_Ecl1;
            if option.SeTe==1 % custom for Se and Te (edit 2026)
                DEclSul(1)=10.^(3.47-3.07*1e-2.*ME_meltEcl.FeO-9.13*1e-4.*ME_meltEcl.FeO.^2);
                DEclSul(2)=10.^(4.3-2.58*1e-3.*ME_meltEcl.FeO.^2);
            else
                DEclSul=Kd.CE.Sul';
            end
            DEclCE=nEcl.Sul.*DEclSul+nEcl.MSS.*D_MSS+nEcl.Sil.*D_Sil;
        end
        
        %% Calculate TE, CE, isoTE and isoCE
        
        %peridotite
        if F_Pdt1 ==0 %below solidus
            if nTE>0
                TE_residuePdt(:)=TE_Pdt;
                TE_meltPdt(:)=zeros(1,nTE);
            end
            if nCE>0
                CE_meltPdt=zeros(1,nCE);
                CE_residuePdt(:)=CE_Pdt;
            end
            if nisoCE>0
                isoCE_meltPdt=zeros(1,nisoCE);
                isoCE_residuePdt(:)=isoCE_Pdt;
            end
            if nisoTE>0
                isoTE_meltPdt=zeros(1,nisoTE);
                isoTE_residuePdt(:)=tisoTE_Pdt;
            end
        elseif F_Pdt1>0  %above solidus
            %modal batch melting for CE
            if nCE>0; [CE_meltPdt, CE_residuePdt]=meltingEquations(CE_Pdt,F_Pdt1,DPdtCE,phi1); end
            if nisoCE
                isoCE_meltPdt=isoCE_Pdt;
                isoCE_residuePdt=isoCE_Pdt;
            end
            %non-modal batch melting for TE
            if nTE>0
                TE_meltPdt=TE_Pdt./(DPdt+F_Pdt1.*(1-PPdt));
                TE_residuePdt=TE_Pdt.*((DPdt-PPdt.*F_Pdt1)./(1-F_Pdt1)).*(1./(DPdt+F_Pdt1.*(1-PPdt)));
            end
            if nisoTE>0
                isoTE_meltPdt=isoTE_Pdt;
                isoTE_residuePdt=isoTE_Pdt;
            end
        end
        
        %eclogite
        if F_Ecl1==1 %above liquidus
            if nTE>0
                TE_residueEcl(:)=zeros(1,nTE);
                TE_meltEcl=TE_Ecl;
            end
            if nisoTE>0
                isoTE_residueEcl=isoTE_Ecl;
                isoTE_meltEcl=isoTE_Ecl;
            end
            if nCE>0
                CE_residueEcl(:)=zeros(1,nCE);
                CE_meltEcl=CE_Ecl;
            end
            if nisoCE>0
                isoCE_residueEcl=isoCE_Ecl;
                isoCE_meltEcl=isoCE_Ecl;
            end
        elseif F_Ecl1 ==0 %below solidus
            if nTE>0
                TE_residueEcl(:)=TE_Ecl;
                TE_meltEcl(:)=zeros(1,nTE);
            end
            if nCE>0
                CE_meltEcl=zeros(1,nCE);
                CE_residueEcl(:)=CE_Ecl;
            end
            if nisoCE>0
                isoCE_meltEcl=zeros(1,nisoCE);
                isoCE_residueEcl(:)=isoCE_Ecl;
            end
            if nisoTE>0
                isoTE_meltEcl=zeros(1,nisoTE);
                isoTE_residueEcl(:)=isoTE_Ecl;
            end
        elseif F_Ecl1>0  %above solidus
            %modal batch melting for CE
            if nCE>0; [CE_meltEcl, CE_residueEcl]=meltingEquations(CE_Ecl,F_Ecl1,DEclCE,phi1); end
            if nisoCE
                isoCE_meltEcl=isoCE_Ecl;
                isoCE_residueEcl=isoCE_Ecl;
            end
            %non-modal batch melting for TE
            if nTE>0
                TE_meltEcl=TE_Ecl./(DEcl+F_Ecl1.*(1-PEcl));
                TE_residueEcl=TE_Ecl.*((DEcl-PEcl.*F_Ecl1)./(1-F_Ecl1)).*(1./(DEcl+F_Ecl1.*(1-PEcl)));
            end
            if nisoTE>0
                isoTE_meltEcl=isoTE_Ecl;
                isoTE_residueEcl=isoTE_Ecl;
            end
        end
        
        %% Mixing proportions
       
        PropPdt = 1 - PropEcl1;
        
        PropEcl_melt = PropEcl1 * F_Ecl1;
        PropPdt_reacted = R1 * PropEcl_melt;

        PropPdt_residue = PropPdt_reacted * (1 - F_Pdt1);
        PropPdt_melt = PropPdt_reacted * F_Pdt1;      
   
        PropPyr = PropEcl_melt + PropPdt_reacted;
        PropEcl2 = PropEcl1 * (1- F_Ecl1);

        sum = PropEcl_melt + PropPdt_reacted;
        
        PropEcl_melt = PropEcl_melt / sum;
        PropPdt_residue = PropPdt_residue / sum;
        PropPdt_melt = PropPdt_melt / sum;

        %% Mixing

        ME_Pdtarray=table2array(ME_Pdt);
        ME_Eclarray=table2array(ME_Ecl);
        ME_meltPdtarray=table2array(ME_meltPdt);
        ME_residuePdtarray=table2array(ME_residuePdt);
        ME_meltEclarray=table2array(ME_meltEcl);
        ME_residueEclarray=table2array(ME_residueEcl);

        ME_mixPyr = PropEcl_melt.*ME_meltEclarray + ...
                    PropPdt_melt.*ME_meltPdtarray + ...
                    PropPdt_residue.*ME_residuePdtarray;
        
        if nCE>0
            CE_mixPyr = PropEcl_melt.*CE_meltEcl + ...
                        PropPdt_melt.*CE_meltPdt + ...
                        PropPdt_residue.*CE_residuePdt;
        end
        
        S_mixPyr = PropEcl_melt*S_meltEcl + ...
                   PropPdt_melt*S_meltPdt + ...
                   PropPdt_residue*S_residuePdt;
        
        if nTE>0
            TE_mixPyr = PropEcl_melt.*TE_meltEcl + ...
                        PropPdt_melt.*TE_meltPdt + ...
                        PropPdt_residue.*TE_residuePdt;
        end
        
        if nisoTE>0
            isoTE_mixPyr = (PropEcl_melt.*isoTE_meltEcl.*TE_meltEcl(isoTE_index) + ...
                         PropPdt_melt.*isoTE_meltPdt.*TE_meltPdt(isoTE_index) + ...
                         PropPdt_residue.*isoTE_residuePdt.*TE_residuePdt(isoTE_index)) ./ ...
                        (PropEcl_melt.*TE_meltEcl(isoTE_index) + ...
                         PropPdt_melt.*TE_meltPdt(isoTE_index) + ...
                         PropPdt_residue.*TE_residuePdt(isoTE_index));
        end
        
        if nisoCE>0
            isoCE_mixPyr = (PropEcl_melt.*isoCE_meltEcl.*CE_meltEcl(isoCE_index) + ...
                         PropPdt_melt.*isoCE_meltPdt.*CE_meltPdt(isoCE_index) + ...
                         PropPdt_residue.*isoCE_residuePdt.*CE_residuePdt(isoCE_index)) ./ ...
                        (PropEcl_melt.*CE_meltEcl(isoCE_index) + ...
                         PropPdt_melt.*CE_meltPdt(isoCE_index) + ...
                         PropPdt_residue.*CE_residuePdt(isoCE_index));
        end

    case 0 % solid-state mixing
    
        ME_mixPyr=PropEcl1.*ME_Eclarray+(1-PropEcl1).*ME_Pdtarray;
        if nTE>0; TE_mixPyr=PropEcl1.*TE_Ecl+(1-PropEcl1).*TE_Pdt;end
        if nCE>0; CE_mixPyr=PropEcl1.*CE_Ecl+(1-PropEcl1).*CE_Pdt;end
        S_mixPyr=PropEcl1*S_Ecl+(1-PropEcl1)*S_Pdt;
        if nisoTE>0; isoTE_mixPyr=((1-PropEcl1).*isoTE_Pdt.*TE_Pdt(isoTE_index)+PropEcl1.*isoTE_Ecl.*TE_Ecl(isoTE_index))./((1-PropEcl1).*TE_Pdt(isoTE_index)+PropEcl1.*TE_Ecl(isoTE_index));end
        if nisoCE>0;isoCE_mixPyr=((1-PropEcl1).*isoCE_Pdt.*CE_Pdt(isoCE_index)+PropEcl1.*isoCE_Ecl.*CE_Ecl(isoCE_index))./((1-PropEcl1).*CE_Pdt(isoCE_index)+PropEcl1.*CE_Ecl(isoCE_index));end

end

%% Process and store output structures

param_out.PropPyr = PropPyr;
param_out.PropEcl2 = PropEcl2;

param_out.realT1_Pdt = realT1_Pdt;
param_out.realT1_Ecl = realT1_Ecl;
if option.stage1MeltPX
    param_out.T1 = T1;
    param_out.F_Pdt1 = F_Pdt1;
    param_out.F_Ecl1 = F_Ecl1;
end

param_out.SCSS_Ecl1 = SCSS_Ecl;
param_out.SCSS_Pdt1 = SCSS_Pdt;

param_out.PropEcl_melt = PropEcl_melt;
param_out.PropPdt_melt = PropPdt_melt;
param_out.PropPdt_residue = PropPdt_residue;

param_out.S_mixPyr1 = S_mixPyr;
param_out.S_residueEcl1 = S_residueEcl;
param_out.S_residuePdt1 = S_residuePdt;

if nME>0
    ME_mixPyr = array2table(ME_mixPyr, VariableNames=ME_oxlist); 
    [normX_mixPyr1,normY_mixPyr1]=normProjection(ME_mixPyr);
    param_out.normX_mixPyr1 = normX_mixPyr1;
    param_out.normY_mixPyr1 = normY_mixPyr1;
    compo_out.ME_mixPyr1 = ME_mixPyr;
    compo_out.ME_mixPyr1Ref = Pyr_string;
    compo_out.ME_residueEcl1 = ME_residueEcl;
    compo_out.ME_residueEcl1Ref = residueEcl_string;
end

if nTE>0
    compo_out.TE_mixPyr1 = array2table(TE_mixPyr, VariableNames=TE);
    compo_out.TE_residueEcl1 = array2table(TE_residueEcl, VariableNames=TE);
end

if nCE>0
    compo_out.CE_mixPyr1 = array2table(CE_mixPyr, VariableNames=CE);
    compo_out.CE_residueEcl1 = array2table(CE_residueEcl, VariableNames=CE);
end

if nisoTE>0
    compo_out.isoTE_mixPyr1 = array2table(isoTE_mixPyr, VariableNames=isoTE);
    compo_out.isoTE_residueEcl1 = array2table(isoTE_residueEcl, VariableNames=isoTE);
end

if nisoCE>0
    compo_out.isoCE_mixPyr1 = array2table(isoCE_mixPyr, VariableNames=isoCE);
    compo_out.isoCE_residueEcl1 = array2table(isoCE_residueEcl, VariableNames=isoCE);
end

compo_out = orderfields(compo_out);
param_out = orderfields(param_out);

end

