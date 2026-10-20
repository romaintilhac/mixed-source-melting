function buildTableOxide(pth)
% Builds TableOxide.mat from oxides.csv

fileID=fullfile(pth.resources, 'TableOxide.mat');

import = readtable(fullfile(pth.resources,'oxides.csv'));

import_size =  size(import);
import_nrows =  import_size(1);
import_ncols= import_size(2);

for row = 1:import_nrows
    elem=string(import.Element(row));
    mass=import.Mass(row);
    TableOxi.(elem).mass=mass;
end   

form_list=split("R2O RO R2O3 RO2 R2O5 RO3");
catn_list=[2 1 2 1 2 1];
oxn_list=[1 1 3 2 5 3];
charge_list=[1 2 3 4 5 6];

for row = 1:import_nrows
    elem=string(import.Element(row));

    for i=1:length(form_list)
        form=form_list(i);
        if ~isnan(import.(form)(row))
            catn=catn_list(i); oxn=oxn_list(i);
            charge=charge_list(i);
            ox=erase(strcat(elem,string(catn),'O',string(oxn)),'1');
            if import.(form)(row)==0 %single oxidation state
                TableOxi.(elem).catn=catn;
                TableOxi.(elem).oxn=oxn;
                TableOxi.(elem).charge=charge;
                TableOxi.(elem).ox=ox;
                TableOxi.(elem).mol=TableOxi.(elem).catn*TableOxi.(elem).mass+TableOxi.(elem).oxn*TableOxi.O.mass;
                TableOxi.(elem).mox=TableOxi.(elem).mass+(TableOxi.(elem).charge/2)*TableOxi.O.mass;
            elseif import.(form)(row)>=1 %multiple oxidation state
                TableOxi.(elem).(strcat('catn',string(charge)))=catn;
                TableOxi.(elem).(strcat('oxn',string(charge)))=oxn;
                TableOxi.(elem).(strcat('charge',string(charge)))=charge;
                TableOxi.(elem).(strcat('ox',string(charge)))=ox;
                TableOxi.(elem).(strcat('mol',string(charge)))=TableOxi.(elem).(strcat('catn',string(i)))*TableOxi.(elem).mass+TableOxi.(elem).(strcat('oxn',string(i)))*TableOxi.O.mass;
                TableOxi.(elem).(strcat('mox',string(charge)))=TableOxi.(elem).mass+(TableOxi.(elem).(strcat('charge',string(i)))/2)*TableOxi.O.mass;
                if import.(form)(row)==2 % default oxidation state
                    TableOxi.(elem).default=strcat(elem,string(charge));
                end
            end
        end
    end
end   

save(fileID,'TableOxi')