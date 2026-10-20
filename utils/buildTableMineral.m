function buildTableMineral(pth)
% Builds TableMineral.mat

fileID=fullfile(pth.resources, 'TableMineral.mat');

TableMin.Oli='olivine';
TableMin.Opx='orthopyroxene';
TableMin.Cpx='clinopyroxene';
TableMin.Plg='plagioclase';
TableMin.Spl='spinel';
TableMin.Grt='garnet';

save(fileID, 'TableMin')