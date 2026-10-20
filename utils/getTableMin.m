function struct = getTableMin()
% Create a persistent structure from TableMineral.mat

    persistent persistentStructure;

    if isempty(persistentStructure)
        % Initialize the structure if it's empty
        persistentStructure = load('resources/TableMineral.mat');
    end

    struct = persistentStructure;

end
