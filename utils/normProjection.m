function [normX, normY] = normProjection(ME)
% Converts major-element compositions into x-y coordinates for plotting in the CaTs-Fo-Qz ternary projection

struct=getTableOxi; TableOxi=struct.TableOxi;

for i= 1: width(ME)
    ME.(i)(isnan(ME.(i))) = 0;
end   

SiO2=ME.SiO2;
TiO2=ME.TiO2;
Al2O3=ME.Al2O3;
Cr2O3=ME.Cr2O3;
FeO=ME.FeO;
Fe2O3=ME.Fe2O3;
MnO=ME.MnO;
MgO=ME.MgO;
CaO=ME.CaO;
Na2O=ME.Na2O;
K2O=ME.K2O;
NiO=ME.NiO;
P2O5=ME.P2O5;

FeOT=FeO+TableOxi.Fe.mol2*Fe2O3/TableOxi.Fe.mol3;
tot=SiO2+TiO2+Al2O3+Cr2O3+FeOT+MnO+MgO+CaO+Na2O+K2O+P2O5;

SiO2=100*SiO2/tot;
TiO2=100*TiO2/tot;
Al2O3=100*Al2O3/tot;
Cr2O3=100*Cr2O3/tot;
FeO=100*FeO/tot;
Fe2O3=100*Fe2O3/tot;
MnO=100*MnO/tot;
MgO=100*MgO/tot;
CaO=100*CaO/tot;
Na2O=100*Na2O/tot;
K2O=100*K2O/tot;
NiO=100*NiO/tot;
P2O5=100*P2O5/tot;

normCaTs=TiO2+Al2O3+Cr2O3+Na2O+K2O;
normFo=0.5.*(Al2O3+Cr2O3+FeOT+MnO+MgO)-0.5.*(CaO+Na2O+K2O)+1.75.*P2O5;
normQz=SiO2+TiO2+0.5.*(Al2O3+Cr2O3)-0.5.*(FeOT+MnO+MgO)-1.5.*CaO-4.5.*(Na2O+K2O)+5.25.*P2O5;

normCaTs=normCaTs.*(TableOxi.Ca.mass+2*TableOxi.Al.mass+TableOxi.Si.mass+6*TableOxi.O.mass);
normFo=normFo.*(2*TableOxi.Mg.mass+TableOxi.Si.mass+4*TableOxi.O.mass);
normQz=normQz*(TableOxi.Si.mass+2*TableOxi.O.mass);

tot=normCaTs+normFo+normQz;
normCaTs=100.*normCaTs./tot; %apex A
normFo=100.*normFo./tot; %apex B
normQz=100.*normQz./tot; %apex C

normX=normQz+normCaTs./2; %X=C+A/2
normY=normCaTs; %Y=A

end