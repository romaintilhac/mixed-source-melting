function S = addFields(S, vars, asArray, varNames)
% Adds fields to a struct S from the caller workspace

    % ---- parse vars ----
    if isstring(vars)
        vars = split(strtrim(vars));
    elseif ischar(vars)
        vars = strsplit(strtrim(vars));
    else
        vars = cellstr(vars);
    end

    % ============================================================
    % CASE 1: RAW STORAGE (varNames NOT USED)
    % ============================================================
    if asArray == 1
        for i = 1:numel(vars)
            name = vars{i};

            if ~evalin('caller', sprintf('exist(''%s'',''var'')', name))
                warning('Variable "%s" not found.', name);
                continue;
            end

            S.(name) = evalin('caller', name);
        end
        return;
    end

    % ============================================================
    % CASE 2: ARRAY -> TABLE (varNames REQUIRED)
    % ============================================================

    % ---- pack ----
    for i = 1:numel(vars)
        name = vars{i};

        if ~evalin('caller', sprintf('exist(''%s'',''var'')', name))
            warning('Variable "%s" not found.', name);
            continue;
        end

        val = evalin('caller', name);

        if istable(val)
            S.(name) = val;
        else
            S.(name) = array2table(val, 'VariableNames', varNames);
        end
    end
    
end