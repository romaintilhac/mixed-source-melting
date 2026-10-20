function lst = userModelComponents
% Set model components

    %% User changes
    %==============================================================
    
    lst.minerals=split("Oli Opx Cpx Plg Spl Grt");                       % list minerals to be included
    lst.ME=split("Si Ti Al Cr Fe3 Fe2 Mg Ni Mn Ca Na K P");     % list lithophile major elements to be included
    lst.TE=split("");                                                                             % list lithophile trace elements to be included
    lst.isoTE=split("");                                                                        % list lithophile isotopes to be included
    lst.SE=split("Fe Ni Cu O S");                                                      % list chalcophile major elements to be included
    lst.CE=split("Se Te");                                                                   % list chalcophile trace elements to be included
    lst.isoCE=split("Se82");                                                              % list chalcophile isotopes to be included

    %==============================================================

    %% Allocations (NO user change below)

    struct=getTableOxi; TableOxi=struct.TableOxi;
    struct=getTableIso; TableIso=struct.TableIso;

    if ~(lst.minerals=="");lst.nMineral=length(lst.minerals);else;lst.nMineral=0;end
    if ~(lst.ME==""); lst.nME=length(lst.ME);else;lst.nME=0;end
    if ~(lst.TE==""); lst.nTE=length(lst.TE);else;lst.nTE=0;end
    if ~(lst.SE==""); lst.nSE=length(lst.SE);else;lst.nSE=0;end
    if ~(lst.CE==""); lst.nCE=length(lst.CE);else;lst.nCE=0;end 
    if ~(lst.isoTE==""); lst.nisoTE=length(lst.isoTE);else;lst.nisoTE=0;end 
    if ~(lst.isoCE==""); lst.nisoCE=length(lst.isoCE);else;lst.nisoCE=0;end 
    
    if lst.nisoTE >=1
        for i=1:lst.nisoTE
            try 
                elem=TableIso.(lst.isoTE(i)).elem;
                lst.isoTE_index(i)=find(lst.TE==elem);
            catch
                warning(append(elem, ' is required for ', lst.isoTE(i), ' but is missing from TE_list'))
            end
        end
    else
        lst.isoTE_index="";
    end

    if lst.nisoCE >=1
        for i=1:lst.nisoCE
            try 
                elem=TableIso.(lst.isoCE(i)).elem;
                lst.isoCE_index(i)=find(lst.CE==elem);
            catch
                warning(append(elem, ' is required for ', lst.isoCE(i), ' but is missing from CE_list'))
            end
        end
    else
        lst.isoCE_index="";
    end
    
    if lst.nME>=1
        for i=1:lst.nME
            ato = regexprep(lst.ME(i), '[0-9]', '');
            oxn=string(regexp(lst.ME(i), '\d+', 'match'));
            lst.atomass(i)=TableOxi.(ato).mass;
            if isempty(oxn) && ~isfield(TableOxi.(ato), 'default') % single oxidation state
                lst.moxmass(i)=TableOxi.(ato).mox;
                lst.molmass(i)=TableOxi.(ato).mol;
                lst.charge(i)=TableOxi.(ato).charge;  
                lst.catnum(i)=TableOxi.(ato).catn;
                lst.ME_oxlist(i,1)=string(TableOxi.(ato).ox);
            elseif ~isempty(oxn) || isfield(TableOxi.(ato), 'default') % multiple oxidation states
                if isfield(TableOxi.(ato), 'default'); oxn=string(regexp(TableOxi.(ato).default, '\d+', 'match'));end % default oxidation state
                lst.moxmass(i)=TableOxi.(ato).(strcat('mox',oxn));
                lst.molmass(i)=TableOxi.(ato).(strcat('mol',oxn));
                lst.charge(i)=TableOxi.(ato).(strcat('charge',oxn));  
                lst.catnum(i)=TableOxi.(ato).(strcat('catn',oxn));
                lst.ME_oxlist(i,1)=string(TableOxi.(ato).(strcat('ox',oxn)));
            end   
        end
    end

    if lst.nSE >=1
        for i=1:lst.nSE
            ato = regexprep(lst.SE(i), '[0-9]', '');
            lst.atomass_Sul(i)=TableOxi.(ato).mass;
        end
    end

end
