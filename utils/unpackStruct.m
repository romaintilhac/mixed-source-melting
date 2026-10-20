function unpackStruct(S, vars, asArray, keep)
% Unpack struct fields into caller workspace

% vars     : field names to unpack
% asArray  : true -> convert table values using table2array
% keep     : true -> keep original struct (default = true)

% default: keep original struct
if nargin < 4 || isempty(keep)
    keep = true;
end

if nargin < 2 || isempty(vars)
    f = fieldnames(S);
else
    if ischar(vars) || (isstring(vars) && isscalar(vars))
        f = strsplit(char(vars));
    else
        f = cellstr(vars);
    end

    f = f(isfield(S, f));
end

if nargin < 3 || isempty(asArray)
    asArray = false;
end

for k = 1:numel(f)
    val = S.(f{k});

    % convert table -> array if requested
    if asArray && istable(val)
        val = table2array(val);
    end

    assignin('caller', f{k}, val);
end

% optionally clear original struct variable
if ~keep
    callerName = inputname(1);
    if ~isempty(callerName)
        evalin('caller', ['clear ', callerName]);
    end
end

end