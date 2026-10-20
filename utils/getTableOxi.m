function struct = getTableOxi()
% Create a persistent structure from TableOxide.mat

    persistent persistentStructure;

    if isempty(persistentStructure)
        % Initialize the structure if it's empty
        persistentStructure = load('resources/TableOxide.mat');
    end

    struct = persistentStructure;

end
