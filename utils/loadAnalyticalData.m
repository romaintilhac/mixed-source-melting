% Load analytical data

pth=userFolderPaths

import=readtable(fullfile(pth.resources, 'data.xlsx'));

import_Pdt=import(contains(import.Rock_type,'Peridotite','IgnoreCase',true),:);
import_Ecl=import(contains(import.Rock_type,'Eclogite','IgnoreCase',true),:);
import_Pyr=import(contains(import.Rock_type,'Pyroxenite','IgnoreCase',true),:);
import_MORB=import(contains(import.Rock_type,'MORB','IgnoreCase',true),:);

data.Se82_Pdt=import_Pdt.Se82;
data.Se82_Ecl=import_Ecl.Se82;
data.Se82_Pyr=import_Pyr.Se82;
data.Se82_MORB=import_MORB.Se82;
data.Se_Pdt=import_Pdt.Se;
data.Se_Ecl=import_Ecl.Se;
data.Se_Pyr=import_Pyr.Se;
data.Se_MORB=import_MORB.Se;

% mean
Se_Ecl=mean(data.Se_Ecl);
Se_Pdt=mean(data.Se_Pdt);
Se_Pyr=mean(data.Se_Pyr);
Se_MORB=mean(data.Se_MORB);
Se82_Ecl=mean(data.Se82_Ecl);
Se82_Pdt=mean(data.Se82_Pdt);
Se82_Pyr=mean(data.Se82_Pyr);
Se82_MORB=mean(data.Se82_MORB);

% fit distributions
Se82_Pdt_pd=fitdist(data.Se82_Pdt,"Normal"); 
Se82_Ecl_pd=fitdist(data.Se82_Ecl,"Normal");
Se82_Pyr_pd=fitdist(data.Se82_Pyr,"Normal");
Se_Pdt_pd=fitdist(data.Se_Pdt,"Normal"); 
Se_Ecl_pd=fitdist(data.Se_Ecl,"Normal");
Se_Pyr_pd=fitdist(data.Se_Pyr,"Normal");