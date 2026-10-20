function [MeltingD_Pdt,MeltingD_Ecl]=MeltPXisobaric2(P, T, PropCpx, ME_Ecl)
% Computes the degrees of melting of peridotite and eclogite in an isobaric situation following Lambart et al. 2016

    struct=getTableOxi; TableOxi=struct.TableOxi;
    
    %% Convert eclogite composition from wt % to mol%
    
    sum_Ecl=ME_Ecl.SiO2+ME_Ecl.TiO2+ME_Ecl.Al2O3+ME_Ecl.Cr2O3+ME_Ecl.FeO+ME_Ecl.MnO+ME_Ecl.MgO+ME_Ecl.CaO+ME_Ecl.Na2O+ME_Ecl.K2O;
    
    SiO2_Ecl=100*ME_Ecl.SiO2/sum_Ecl;
    TiO2_Ecl=100*ME_Ecl.TiO2/sum_Ecl;
    Al2O3_Ecl=100*ME_Ecl.Al2O3/sum_Ecl;
    Cr2O3_Ecl=100*ME_Ecl.Cr2O3/sum_Ecl;
    FeO_Ecl=100*ME_Ecl.FeO/sum_Ecl;
    MnO_Ecl=100*ME_Ecl.MnO/sum_Ecl;
    MgO_Ecl=100*ME_Ecl.MgO/sum_Ecl;
    CaO_Ecl=100*ME_Ecl.CaO/sum_Ecl;
    Na2O_Ecl=100*ME_Ecl.Na2O/sum_Ecl;
    K2O_Ecl=100*ME_Ecl.K2O/sum_Ecl;
    
    SiO2_Ecl=SiO2_Ecl/TableOxi.Si.mol;
    TiO2_Ecl=TiO2_Ecl/TableOxi.Ti.mol;
    Al2O3_Ecl=Al2O3_Ecl/TableOxi.Al.mol;
    Cr2O3_Ecl=Cr2O3_Ecl/TableOxi.Cr.mol3;
    FeO_Ecl=FeO_Ecl/TableOxi.Fe.mol2;
    MnO_Ecl=MnO_Ecl/TableOxi.Mn.mol2;
    MgO_Ecl=MgO_Ecl/TableOxi.Mg.mol;
    CaO_Ecl=CaO_Ecl/TableOxi.Ca.mol;
    Na2O_Ecl=Na2O_Ecl/TableOxi.Na.mol;
    K2O_Ecl=K2O_Ecl/TableOxi.K.mol;
    
    sum_Ecl=SiO2_Ecl+TiO2_Ecl+Al2O3_Ecl+Cr2O3_Ecl+FeO_Ecl+MnO_Ecl+MgO_Ecl+CaO_Ecl+Na2O_Ecl+K2O_Ecl;
    
    SiO2_Ecl=SiO2_Ecl*100/sum_Ecl;
    TiO2_Ecl=TiO2_Ecl*100/sum_Ecl;
    Al2O3_Ecl=Al2O3_Ecl*100/sum_Ecl;
    Cr2O3_Ecl=Cr2O3_Ecl*100/sum_Ecl;
    FeO_Ecl=FeO_Ecl*100/sum_Ecl;
    MnO_Ecl=MnO_Ecl*100/sum_Ecl;
    MgO_Ecl=MgO_Ecl*100/sum_Ecl;
    CaO_Ecl=CaO_Ecl*100/sum_Ecl;
    Na2O_Ecl=Na2O_Ecl*100/sum_Ecl;
    K2O_Ecl=K2O_Ecl*100/sum_Ecl;
    alkali_Ecl=Na2O_Ecl+K2O_Ecl;
    MgNo_Ecl=MgO_Ecl/(FeO_Ecl+MgO_Ecl);
    
    %% Melt-PX parameters (Lambart et al. 2016)
    
    a=1.4;b=-7.08;c=-12.5;d=150;e=25.5;f=1105;
    aa=-14;bb=177;cc=219;dd=-17.2;ee=10.4;ff=1000;
    Aprim=-2.63;Bprim=-1.58;Cprim=26.9;Dprim=-7.5;Eprim=4.23;AA=2.86;BB=18;
    
    %% Eclogite parametrization
    
    % T 5% pyroxenite
    T5_Ecl=(a*(alkali_Ecl)+b)*(P^2)+(c*(alkali_Ecl)+d*MgNo_Ecl+e)*P+f;
    
    % T cpx-out pyroxenite
    Tcpxout_Ecl=aa*(P^2)+bb*P+cc*MgNo_Ecl+dd*Al2O3_Ecl+ee*CaO_Ecl+ff;
    
    %T' pyroxenite
    A=(Aprim+Bprim*MgNo_Ecl)*(P^2)+(Cprim+Dprim*MgNo_Ecl)*P+Eprim;
    B=AA*CaO_Ecl+BB;
    Tprim_Ecl=(T-T5_Ecl)/(Tcpxout_Ecl-T5_Ecl);
    F_Ecl=A*Tprim_Ecl^2+B*Tprim_Ecl+5;
    
    Tpp_Ecl=(T+1-T5_Ecl)/(Tcpxout_Ecl-T5_Ecl);
    FTpp_Ecl=A*Tpp_Ecl+B*Tpp_Ecl+5;
    
    if Tprim_Ecl<0 || B^2-4*A*5>0; intersect0=true; else; intersect0=false; end
    if intersect0==true; TprimF0_Ecl=(-B+(B^2-4*A*5)^0.5)/(2*A); else; error('F cannot be calculated'); end
    if FTpp_Ecl>0 && intersect0==true && Tprim_Ecl<TprimF0_Ecl; Tprimvsracine=true; else; Tprimvsracine=false; end
    if Tpp_Ecl<F_Ecl && intersect0==false; final=true; else; final=false; end
    
    if T5_Ecl<Tcpxout_Ecl; cpxfree=false; Fcpxfree_Ecl=F_Ecl; else; cpxfree=true; error('Cpx-free assemblage'); end
    if F_Ecl<100; F_Ecl=Fcpxfree_Ecl; else; F_Ecl=100;end
    if F_Ecl>0; F_Ecl=F_Ecl; else; F_Ecl=0;end
    if Tprimvsracine==true;F5_Ecl=0;else; F5_Ecl=F_Ecl; end
    if final==false;Ffinal_Ecl=F5_Ecl;else; error('F cannot be calculated'); end
    if cpxfree==true; F_Ecl=Fcpxfree_Ecl; else; F_Ecl=Ffinal_Ecl; end
    
    %% Peridotite parametrization (Katz 2003)
    
    %peridotite solidus temperature
    Tsol_Pdt=1085.7+132.9*P-5.1*(P^2);
    
    % lherzolite liquidus temperature
    Tliq_Lhz=1475+80*P-3.2*P^2;
    
    % Fcpx-out
    Fcpxout_Pdt=PropCpx/(0.5+0.08*P);
    
    % Tcpx-out peridotite
    Tcpxout_Pdt=(Fcpxout_Pdt.^(1/1.5))*(Tliq_Lhz-Tsol_Pdt)+Tsol_Pdt;
    
    % peridotite liquidus temperature
    Tliq_Pdt=1780+45*P-2*P^2;
    
    %F calculations
    TprimPdt=(T-Tsol_Pdt)/(Tliq_Lhz-Tsol_Pdt);
    
    if TprimPdt<0
        Fcpx_Pdt=0;
    else
        Fcpx_Pdt=(TprimPdt^1.5)*100;
    end
    
    if T<Tcpxout_Pdt
        Fopx_Pdt=0;
    else
        Fopx_Pdt=((Fcpxout_Pdt+(1-Fcpxout_Pdt)*((T-Tcpxout_Pdt)/(Tliq_Pdt-Tcpxout_Pdt)))^1.5)*100;
    end
    
    if Fopx_Pdt>0
        F_Pdt=Fopx_Pdt;
    else
        F_Pdt=Fcpx_Pdt;
    end
    
    MeltingD_Pdt=F_Pdt; 
    MeltingD_Ecl=F_Ecl; 

end