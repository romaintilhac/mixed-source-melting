function struct = getTableIso()
% Create a persistent structure from TableIsotope.mat

    persistent persistentStructure;

    if isempty(persistentStructure)
        % Initialize the structure if it's empty
        persistentStructure = load('resources/TableIsotope.mat');
    end

    struct = persistentStructure;

end
