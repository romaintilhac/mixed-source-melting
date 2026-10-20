function [param, compo] = Stage1pass(param, param_out1, compo, compo_out1, lst, option, settings)
% Passes the outputs of Stage1model to the main param and compo structures

unpackStruct(settings, "Prop_round missingME pristineEcl_thresh", 0, 0);

if option.stage1Prop==1

    param.PropPyr = round(param_out1.PropPyr / Prop_round) * Prop_round;
    param.PropEcl2 = round(param_out1.PropEcl2 / Prop_round) * Prop_round;

elseif F_Ecl1 == 1

    param.PropEcl2 = 0;

end

%% Pyroxenite
    
    if lst.nTE>0
        compo.TE_Pyr=compo_out1.TE_mixPyr1;
    end
    if lst.nCE>0
        compo.CE_Pyr=compo_out1.CE_mixPyr1;
    end
    if lst.nisoTE>0
        compo.isoTE_Pyr=compo_out1.isoTE_mixPyr1;
    end
    if lst.nisoCE>0
        compo.isoCE_Pyr=compo_out1.isoCE_mixPyr1;
    end

    % ME
    if option.stage1PyrME==1
       
        for n=1:numel(lst.ME_oxlist)
    
            ME=lst.ME_oxlist(n);

            if compo_out1.ME_mixPyr1.(ME)<missingME
                compo_out1.ME_mixPyr1.(ME)=missingME;
            end
    
        end
    
        compo.ME_Pyr=compo_out1.ME_mixPyr1;
        compo.ME_PyrRef=compo_out1.ME_mixPyr1Ref;
    
    end  

    % S
    if option.stage1PyrS==1
        param.S_Pyr=param_out1.S_mixPyr1;
    end

%% Eclogite

    % ME
    if option.stage1EclresidueME==1
    
        for n=1:numel(lst.ME_oxlist)
    
            ME=lst.ME_oxlist(n);
    
            if compo_out1.ME_residueEcl1.(ME)<missingME
                compo_out1.ME_residueEcl1.(ME)=missingME;
            end
    
        end

        if param_out1.F_Ecl1 > pristineEcl_thresh
            compo.ME_Ecl=compo_out1.ME_residueEcl1;
            compo.ME_EclRef=compo_out1.ME_residueEcl1Ref;
        end
    
    end

    % TE
    if option.stage1EclresidueTE==1
        if lst.nTE>0; compo.TE_Ecl=compo_out1.TE_residueEcl1; end
        if lst.nisoTE>0; compo.isoTE_Ecl=compo_out1.isoTE_residueEcl1; end
    end

    % CE
    if option.stage1EclresidueCE==1
        if lst.nCE>0; compo.CE_Ecl=compo_out1.CE_residueEcl1; end
        if lst.nisoCE>0; compo.isoCE_Ecl=compo_out1.isoCE_residueEcl1; end
    end

    % S
    if option.stage1EclresidueS==1
        param.S_Ecl=param_out1.S_residueEcl1;
    end

end