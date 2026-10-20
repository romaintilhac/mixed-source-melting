function Tsol_Ecl = eclogiteSolidus(P, PropCpx, ME_Ecl)
% Estimates eclogite solidus temperature based on Lambart et al. 2016

    struct = getTableOxi;
    TableOxi = struct.TableOxi;

    %% Convert eclogite composition to mol %
    sum_Ecl = ME_Ecl.SiO2 + ME_Ecl.TiO2 + ME_Ecl.Al2O3 + ME_Ecl.Cr2O3 + ...
              ME_Ecl.FeO + ME_Ecl.MnO + ME_Ecl.MgO + ME_Ecl.CaO + ...
              ME_Ecl.Na2O + ME_Ecl.K2O;

    ox = {'Si','Ti','Al','Cr','Fe','Mn','Mg','Ca','Na','K'};
    wt = [ME_Ecl.SiO2 ME_Ecl.TiO2 ME_Ecl.Al2O3 ME_Ecl.Cr2O3 ME_Ecl.FeO ...
          ME_Ecl.MnO ME_Ecl.MgO ME_Ecl.CaO ME_Ecl.Na2O ME_Ecl.K2O] / sum_Ecl * 100;

    mol = wt ./ [TableOxi.Si.mol TableOxi.Ti.mol TableOxi.Al.mol ...
                 TableOxi.Cr.mol3 TableOxi.Fe.mol2 TableOxi.Mn.mol2 ...
                 TableOxi.Mg.mol TableOxi.Ca.mol TableOxi.Na.mol TableOxi.K.mol];

    mol = mol / sum(mol) * 100;

    Al2O3 = mol(3);
    CaO   = mol(8);
    alk   = mol(9) + mol(10);
    MgNo  = mol(7) / (mol(5) + mol(7));

    %% Lambart et al. (2016) parameters
    a=1.4; b=-7.08; c=-12.5; d=150; e=25.5; f=1105;
    aa=-14; bb=177; cc=219; dd=-17.2; ee=10.4; ff=1000;
    Aprim=-2.63; Bprim=-1.58; Cprim=26.9; Dprim=-7.5; Eprim=4.23;
    AA=2.86; BB=18;

    %% Reference temperatures
    T5_Ecl = (a*alk + b)*P^2 + (c*alk + d*MgNo + e)*P + f;
    Tcpxout_Ecl = aa*P^2 + bb*P + cc*MgNo + dd*Al2O3 + ee*CaO + ff;

    %% Quadratic melt law coefficients
    A = (Aprim + Bprim*MgNo)*P^2 + (Cprim + Dprim*MgNo)*P + Eprim;
    B = AA*CaO + BB;

    %% Solidus (F = 0)
    TprimF0 = (-B + sqrt(B^2 - 4*A*5)) / (2*A);
    Tsol_Ecl = T5_Ecl + TprimF0 * (Tcpxout_Ecl - T5_Ecl);

end
