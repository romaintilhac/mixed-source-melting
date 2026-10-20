function buildTableIsotope(pth)
% Builds TableIsotope.mat

fileID=fullfile(pth.resources, 'TableIsotope.mat');

TableIso.Nd143.label='143Nd/144Nd';
TableIso.Nd143.elem='Nd';

TableIso.Se82.label='d82Se';
TableIso.Se82.elem='Se';

TableIso.Te169.label='d169Te';
TableIso.Te169.elem='Te';

save(fileID, 'TableIso')