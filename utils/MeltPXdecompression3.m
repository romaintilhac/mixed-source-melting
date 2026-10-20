function [Pressure,Temperature,MeltingD_Pdt,MeltingD_Pyr,MeltingD_Ecl,Tadiabat,Tsol_Pdt,Tsol_Pyr, Tsol_Ecl]=MeltPXdecompression3(Tp, Pmax, Pmin, Pincr, PropCpx, PropPyr, PropEcl, ME_Pyr, ME_Ecl)
% Computes the degrees of melting of peridotite, pyroxenite and eclogite in a decompression situation following Lambart et al. 2016

    struct=getTableOxi; TableOxi=struct.TableOxi;
    
    PropPdt=1-(PropPyr+PropEcl);
    
    %% Pre-defined parameters
    
    % constants
    rho_Pdt=constants.rho_m; % peridotite density
    rho_Pyr=constants.rho_p; % pyroxenite density
    rho_Ecl=constants.rho_p; % eclogite density
    TK=constants.TK;
    PGPa=constants.PGPa;
    
    % peridotite (Lambart et al. 2016, Table S2)
    Cp_Pdt=1190; % heat capacity
    alpha_Pdt=0.00003; % thermal expansivity
    S_Pdt=395; % entropy of fusion
    
    % pyroxenite (taken as MIX1G, Lambart et al. 2016, Table S2)
    Cp_Pyr=1090; % heat capacity
    alpha_Pyr=0.00003; % thermal expansivity
    S_Pyr=360; % entropy of fusion
    defaultdFdTcpxout_Pyr=0.3; % default isobaric productivtiy after cpx exhaustion [%/ºC]
    
    % eclogite (taken as G2, Lambart et al. 2016, Table S2)
    Cp_Ecl=1045; % heat capacity
    alpha_Ecl=0.00003; % thermal expansivity
    S_Ecl=350; % entropy of fusion 
    defaultdFdTcpxout_Ecl=0.3; % default isobaric productivtiy after cpx exhaustion [%/ºC]
    
    %% Set solid adiabat from bulk thermal expansivity
    
    Nincr=(Pmax-Pmin)/Pincr;
    P=linspace(Pmax,Pmin,Nincr);
    P=P';
    
    T_SA=Tp+(PropPdt.*((alpha_Pdt.*(Tp+TK))./(Cp_Pdt.*rho_Pdt).*PGPa).*P...
        +PropPyr.*((alpha_Pyr.*(Tp+TK))./(Cp_Pyr.*rho_Pyr).*PGPa).*P...
        +PropEcl.*((alpha_Ecl.*(Tp+TK))./(Cp_Ecl*rho_Ecl).*PGPa).*P);
    
    %% Convert pyroxenite composition from wt % to mol%
    
    sum_Pyr=ME_Pyr.SiO2+ME_Pyr.TiO2+ME_Pyr.Al2O3+ME_Pyr.Cr2O3+ME_Pyr.FeO+ME_Pyr.MnO+ME_Pyr.MgO+ME_Pyr.CaO+ME_Pyr.Na2O+ME_Pyr.K2O;
    
    SiO2_Pyr=100*ME_Pyr.SiO2/sum_Pyr;
    TiO2_Pyr=100*ME_Pyr.TiO2/sum_Pyr;
    Al2O3_Pyr=100*ME_Pyr.Al2O3/sum_Pyr;
    Cr2O3_Pyr=100*ME_Pyr.Cr2O3/sum_Pyr;
    FeO_Pyr=100*ME_Pyr.FeO/sum_Pyr;
    MnO_Pyr=100*ME_Pyr.MnO/sum_Pyr;
    MgO_Pyr=100*ME_Pyr.MgO/sum_Pyr;
    CaO_Pyr=100*ME_Pyr.CaO/sum_Pyr;
    Na2O_Pyr=100*ME_Pyr.Na2O/sum_Pyr;
    K2O_Pyr=100*ME_Pyr.K2O/sum_Pyr;
    
    SiO2_Pyr=SiO2_Pyr/TableOxi.Si.mol;
    TiO2_Pyr=TiO2_Pyr/TableOxi.Ti.mol;
    Al2O3_Pyr=Al2O3_Pyr/TableOxi.Al.mol;
    Cr2O3_Pyr=Cr2O3_Pyr/TableOxi.Cr.mol3;
    FeO_Pyr=FeO_Pyr/TableOxi.Fe.mol2;
    MnO_Pyr=MnO_Pyr/TableOxi.Mn.mol2;
    MgO_Pyr=MgO_Pyr/TableOxi.Mg.mol;
    CaO_Pyr=CaO_Pyr/TableOxi.Ca.mol;
    Na2O_Pyr=Na2O_Pyr/TableOxi.Na.mol;
    K2O_Pyr=K2O_Pyr/TableOxi.K.mol;
    
    sum_Pyr=SiO2_Pyr+TiO2_Pyr+Al2O3_Pyr+Cr2O3_Pyr+FeO_Pyr+MnO_Pyr+MgO_Pyr+CaO_Pyr+Na2O_Pyr+K2O_Pyr;
    
    SiO2_Pyr=SiO2_Pyr*100/sum_Pyr;
    TiO2_Pyr=TiO2_Pyr*100/sum_Pyr;
    Al2O3_Pyr=Al2O3_Pyr*100/sum_Pyr;
    Cr2O3_Pyr=Cr2O3_Pyr*100/sum_Pyr;
    FeO_Pyr=FeO_Pyr*100/sum_Pyr;
    MnO_Pyr=MnO_Pyr*100/sum_Pyr;
    MgO_Pyr=MgO_Pyr*100/sum_Pyr;
    CaO_Pyr=CaO_Pyr*100/sum_Pyr;
    Na2O_Pyr=Na2O_Pyr*100/sum_Pyr;
    K2O_Pyr=K2O_Pyr*100/sum_Pyr;
    alkali_Pyr=Na2O_Pyr+K2O_Pyr;
    MgNo_Pyr=MgO_Pyr/(FeO_Pyr+MgO_Pyr);
    
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
    A=-2.63;B=-1.58;C=26.9;D=-7.5;E=4.23;AA=2.86;BB=18;
    
    %% Pyroxenite parametrization
    
    % T 5% pyroxenite
    T5_Pyr=(a.*(alkali_Pyr)+b).*(P(1:end-1).^2)+(c.*(alkali_Pyr)+d.*MgNo_Pyr+e).*P(1:end-1)+f;
    
    % T cpx-out pyroxenite
    Tcpxout_Pyr=aa.*(P(1:end-1).^2)+bb.*P(1:end-1)+cc.*MgNo_Pyr+dd.*Al2O3_Pyr+ee.*CaO_Pyr+ff;
    
    %T' pyroxenite
    Tprim_Pyr=(-(AA.*CaO_Pyr+BB)+((AA.*CaO_Pyr+BB).^2-4.*5.*((A+B.*MgNo_Pyr).*(P(1:end-1).^2)+(C+D.*MgNo_Pyr).*P(1:end-1)+E)).^0.5)./(2.*((A+B.*MgNo_Pyr).*(P(1:end-1).^2)+(C+D.*MgNo_Pyr).*P(1:end-1)+E));
    
    %T_sol pyroxenite
    Tsol_Pyr=Tprim_Pyr.*(Tcpxout_Pyr-T5_Pyr)+T5_Pyr;
    
    %% Eclogite parametrization
    
    % T 5% pyroxenite
    T5_Ecl=(a.*(alkali_Ecl)+b).*(P(1:end-1).^2)+(c.*(alkali_Ecl)+d.*MgNo_Ecl+e).*P(1:end-1)+f;
    
    % T cpx-out pyroxenite
    Tcpxout_Ecl=aa.*(P(1:end-1).^2)+bb.*P(1:end-1)+cc.*MgNo_Ecl+dd.*Al2O3_Ecl+ee.*CaO_Ecl+ff;
    
    %T' pyroxenite
    Tprim_Ecl=(-(AA.*CaO_Ecl+BB)+((AA.*CaO_Ecl+BB).^2-4.*5.*((A+B.*MgNo_Ecl).*(P(1:end-1).^2)+(C+D.*MgNo_Ecl).*P(1:end-1)+E)).^0.5)./(2.*((A+B.*MgNo_Ecl).*(P(1:end-1).^2)+(C+D.*MgNo_Ecl).*P(1:end-1)+E));
    
    %T_sol pyroxenite
    Tsol_Ecl=Tprim_Ecl.*(Tcpxout_Ecl-T5_Ecl)+T5_Ecl;
    
    %% Peridotite parametrization (Katz 2003)
    
    %peridotite solidus temperature
    Tsol_Pdt=1085.7+132.9.*P(1:end-1)-5.1.*(P(1:end-1).^2);
    
    % lherzolite liquidus temperature
    Tliq_Lhz=1475+80.*P(1:end-1)-3.2.*P(1:end-1).^2;
    
    % Fcpx-out
    Fcpxout_Pdt=PropCpx./(0.5+0.08.*P(1:end-1));
    
    % Tcpx-out peridotite
    Tcpxout_Pdt=(Fcpxout_Pdt.^(1./1.5)).*(Tliq_Lhz-Tsol_Pdt)+Tsol_Pdt;
    
    % peridotite liquidus temperature
    Tliq_Pdt=1780+45.*P(1:end-1)-2.*P(1:end-1).^2;
    
    %% Along adiabat
    
    % relative T between T5 and Tcpxout along the solid adiabat
        
        % pyroxenite
        Tprim_Pyr_SA=zeros(length(P)-1,1);
        for i=1:length(P)-1
            if Tcpxout_Pyr(i)<T5_Pyr(i)
                Tprim_Pyr_SA(i)=-1;
            else
                Tprim_Pyr_SA(i)=(T_SA(i)-T5_Pyr(i))/(Tcpxout_Pyr(i)-T5_Pyr(i));
            end
        end
        
        % eclogite
        Tprim_Ecl_SA=zeros(length(P)-1,1);
        for i=1:length(P)-1
            if Tcpxout_Ecl(i)<T5_Ecl(i)
                Tprim_Ecl_SA(i)=-1;
            else
                Tprim_Ecl_SA(i)=(T_SA(i)-T5_Ecl(i))/(Tcpxout_Ecl(i)-T5_Ecl(i));
            end
        end
        
        % peridotite
        TprimPdt_SA=(T_SA(1:end-1)-Tsol_Pdt)./(Tliq_Lhz-Tsol_Pdt);
    
    % melting degree along solid adiabat
    
        % pyroxenite
        F_Pyr_SA=zeros(length(P)-1,1);
        for i=1:length(P)-1
            if Tprim_Pyr_SA(i)<Tprim_Pyr(i)
                F_Pyr_SA(i)=0;
            else
                F_Pyr_SA(i)=((A+B*MgNo_Pyr)*(P(i)^2)+(C+D*MgNo_Pyr)*P(i)+E)*(Tprim_Pyr_SA(i)^2)+(AA*CaO_Pyr+BB)*Tprim_Pyr_SA(i)+5;
            end
        end
        F_Pyr_SA(end+1)=F_Pyr_SA(end);
    
        % eclogite
        F_Ecl_SA=zeros(length(P)-1,1);
        for i=1:length(P)-1
            if Tprim_Ecl_SA(i)<Tprim_Ecl(i)
                F_Ecl_SA(i)=0;
            else
                F_Ecl_SA(i)=((A+B*MgNo_Ecl)*(P(i)^2)+(C+D*MgNo_Ecl)*P(i)+E)*(Tprim_Ecl_SA(i)^2)+(AA*CaO_Ecl+BB)*Tprim_Ecl_SA(i)+5;
            end
        end
        F_Ecl_SA(end+1)=F_Ecl_SA(end);
    
        % peridotite (before and after cpx exhaustion)
    
        Fcpx_Pdt_SA=zeros(length(P)-1,1);
        for i=1:length(P)-1
            if TprimPdt_SA(i)<0
                Fcpx_Pdt_SA(i)=0;
            else
                Fcpx_Pdt_SA(i)=(TprimPdt_SA(i)^1.5)*100;
            end
        end
        
        Fopx_Pdt_SA=zeros(length(P)-1,1);
        for i=1:length(P)-1
            if T_SA(i)<Tcpxout_Pdt(i)
                Fopx_Pdt_SA(i)=0;
            else
                Fopx_Pdt_SA(i)=((Fcpxout_Pdt(i)+(1-Fcpxout_Pdt(i))*((T_SA(i)-Tcpxout_Pdt(i))/(Tliq_Pdt(i)-Tcpxout_Pdt(i))))^1.5)*100;
            end
        end
    
        F_Pdt_SA=zeros(length(P)-1,1);
        for i=1:length(P)-1
            if Fopx_Pdt_SA(i)>0
                F_Pdt_SA(i)=Fopx_Pdt_SA(i);
            else
                F_Pdt_SA(i)=Fcpx_Pdt_SA(i);
            end
        end
    
    %% Solidus slopes (dTs/dP)s
    
    % pyroxenite
    dTsdP_s_Pyr=zeros(length(P)-1,1);
    for i=2:length(P)-1
        dTsdP_s_Pyr(i)=(Tsol_Pyr(i)-Tsol_Pyr(i-1))/(P(i)-P(i-1));
    end
    dTsdP_s_Pyr(1)=dTsdP_s_Pyr(2);
    
    % eclogite
    dTsdP_s_Ecl=zeros(length(P)-1,1);
    for i=2:length(P)-1
        dTsdP_s_Ecl(i)=(Tsol_Ecl(i)-Tsol_Ecl(i-1))/(P(i)-P(i-1));
    end
    dTsdP_s_Ecl(1)=dTsdP_s_Ecl(2);
    
    % peridotite
    dTsdP_s_Pdt=zeros(length(P)-1,1);
    for i=2:length(P)-1
        dTsdP_s_Pdt(i)=(Tsol_Pdt(i)-Tsol_Pdt(i-1))/(P(i)-P(i-1));
    end
    dTsdP_s_Pdt(1)=dTsdP_s_Pdt(2);
    
    %% Initialize all variables for the melting column
    
    T=zeros(length(P)-1,1);
    
    Tprim_PyrT=zeros(length(P)-1,1);
    FTprim_Pyr=zeros(length(P)-1,1);
    Tpp1_Pyr=zeros(length(P)-1,1);
    FTpp1_Pyr=zeros(length(P)-1,1);
    dFdT_P_Pyr=zeros(length(P)-1,1);
    F_Pyr=zeros(length(P)-1,1);
    moddFdT_P_Pyr=zeros(length(P)-1,1);
    dTdP_F_Pyr=zeros(length(P)-1,1);
    dTdP_F_PyrT=zeros(length(P)-1,1);
    dFdP1_Pyr=zeros(length(P)-1,1);
    dFdP2_Pyr=zeros(length(P)-1,1);
    dFdP3_Pyr=zeros(length(P)-1,1);
    dFdP4_Pyr=zeros(length(P)-1,1);
    dFdP_S_Pyr=zeros(length(P)-1,1);
    
    Tprim_EclT=zeros(length(P)-1,1);
    FTprim_Ecl=zeros(length(P)-1,1);
    Tpp1_Ecl=zeros(length(P)-1,1);
    FTpp1_Ecl=zeros(length(P)-1,1);
    dFdT_P_Ecl=zeros(length(P)-1,1);
    F_Ecl=zeros(length(P)-1,1);
    moddFdT_P_Ecl=zeros(length(P)-1,1);
    dTdP_F_Ecl=zeros(length(P)-1,1);
    dTdP_F_EclT=zeros(length(P)-1,1);
    dFdP1_Ecl=zeros(length(P)-1,1);
    dFdP2_Ecl=zeros(length(P)-1,1);
    dFdP3_Ecl=zeros(length(P)-1,1);
    dFdP4_Ecl=zeros(length(P)-1,1);
    dFdP_S_Ecl=zeros(length(P)-1,1);
    
    Tprim_PdtT=zeros(length(P)-1,1);
    F0_cpxTprim_Pdt=zeros(length(P)-1,1);
    Tpp1_Pdt=zeros(length(P)-1,1);
    F_cpxTpp_Pdt=zeros(length(P)-1,1);
    F0Tprim_Pdt=zeros(length(P)-1,1);
    FTpp_Pdt=zeros(length(P)-1,1);
    dFdT_P_Pdt=zeros(length(P)-1,1);
    F_Pdt=zeros(length(P)-1,1);
    dTdP_F_Pdt=zeros(length(P)-1,1);
    dTdP_F_PdtT=zeros(length(P)-1,1);
    dFdP1_Pdt=zeros(length(P)-1,1);
    dFdP2_Pdt=zeros(length(P)-1,1);
    dFdP3_Pdt=zeros(length(P)-1,1);
    dFdP_S_Pdt=zeros(length(P)-1,1);
    
    TpM_Pyr=zeros(length(P)-2,1);
    TM_Pyr=zeros(length(P)-2,1);
    
    TpM_Ecl=zeros(length(P)-2,1);
    TM_Ecl=zeros(length(P)-2,1);
    
    TpM_Pdt=zeros(length(P)-2,1);
    TM_Pdt=zeros(length(P)-2,1);
    
    DTDP=zeros(length(P)-1,1);
    dTdP_S_final=zeros(length(P)-1,1);
    
    for i=1:length(P)-1
        if i==1
        %% i==1
        
            T(i)=T_SA(i);
    
            % T'
    
                % pyroxenite
                if T5_Pyr(i)>Tcpxout_Pyr(i) 
                    Tprim_PyrT(i)=-1;
                else
                    Tprim_PyrT(i)=(T(i)-T5_Pyr(i))/(Tcpxout_Pyr(i)-T5_Pyr(i));
                end
              
                % eclogite
                if T5_Ecl(i)>Tcpxout_Ecl(i) 
                    Tprim_EclT(i)=-1;
                else
                    Tprim_EclT(i)=(T(i)-T5_Ecl(i))/(Tcpxout_Ecl(i)-T5_Ecl(i));
                end
        
                % peridotite
                Tprim_PdtT(i)=(T(i)-Tsol_Pdt(i))/(Tliq_Lhz(i)-Tsol_Pdt(i));
    
            % F(T') Eq. 5 & 6 into Eq. 1 (Lambart et al. 2016)
    
                % pyroxenite
                if Tprim_PyrT(i)<Tprim_Pyr(i)
                    FTprim_Pyr(i)=0;
                else
                    if (((A+B*MgNo_Pyr)*(P(i)^2)+(C+D*MgNo_Pyr)*P(i)+E)*(Tprim_PyrT(i)^2)+(AA*CaO_Pyr+BB)*Tprim_PyrT(i)+5)<0
                        FTprim_Pyr(i)=0;
                    else
                        FTprim_Pyr(i)=(((A+B*MgNo_Pyr)*(P(i)^2)+(C+D*MgNo_Pyr)*P(i)+E)*(Tprim_PyrT(i)^2)+(AA*CaO_Pyr+BB)*Tprim_PyrT(i)+5);
                    end
                end
    
                % eclogite
                if Tprim_EclT(i)<Tprim_Ecl(i)
                    FTprim_Ecl(i)=0;
                else
                    if (((A+B*MgNo_Ecl)*(P(i)^2)+(C+D*MgNo_Ecl)*P(i)+E)*(Tprim_EclT(i)^2)+(AA*CaO_Ecl+BB)*Tprim_EclT(i)+5)<0
                        FTprim_Ecl(i)=0;
                    else
                        FTprim_Ecl(i)=(((A+B*MgNo_Ecl)*(P(i)^2)+(C+D*MgNo_Ecl)*P(i)+E)*(Tprim_EclT(i)^2)+(AA*CaO_Ecl+BB)*Tprim_EclT(i)+5);
                    end
                end
    
            % T''(T+1)
    
                % pyroxenite
                if T5_Pyr(i)>Tcpxout_Pyr(i) 
                    Tpp1_Pyr(i)=-1;
                else
                    Tpp1_Pyr(i)=Tprim_PyrT(i)+(1/(Tcpxout_Pyr(i)-T5_Pyr(i)));
                end
       
                % eclogite
                if T5_Ecl(i)>Tcpxout_Ecl(i) 
                    Tpp1_Ecl(i)=-1;
                else
                    Tpp1_Ecl(i)=Tprim_EclT(i)+(1/(Tcpxout_Ecl(i)-T5_Ecl(i)));
                end
    
                % peridotite
                Tpp1_Pdt(i)=(T(i)+1-Tsol_Pdt(i))/(Tliq_Lhz(i)-Tsol_Pdt(i));
         
            % F(T'')
    
                % pyroxenite
                if Tpp1_Pyr(i)<Tprim_Pyr(i)
                    FTpp1_Pyr(i)=0;
                else
                    if (((A+B*MgNo_Pyr)*(P(i)^2)+(C+D*MgNo_Pyr)*P(i)+E)*(Tpp1_Pyr(i)^2)+(AA*CaO_Pyr+BB)*Tpp1_Pyr(i)+5)<0
                        FTpp1_Pyr(i)=0;
                    else
                        FTpp1_Pyr(i)=(((A+B*MgNo_Pyr)*(P(i)^2)+(C+D*MgNo_Pyr)*P(i)+E)*(Tpp1_Pyr(i)^2)+(AA*CaO_Pyr+BB)*Tpp1_Pyr(i)+5);
                    end
                end
    
                % eclogite
                if Tpp1_Ecl(i)<Tprim_Ecl(i)
                    FTpp1_Ecl(i)=0;
                else
                    if (((A+B*MgNo_Ecl)*(P(i)^2)+(C+D*MgNo_Ecl)*P(i)+E)*(Tpp1_Ecl(i)^2)+(AA*CaO_Ecl+BB)*Tpp1_Ecl(i)+5)<0
                        FTpp1_Ecl(i)=0;
                    else
                        FTpp1_Ecl(i)=(((A+B*MgNo_Ecl)*(P(i)^2)+(C+D*MgNo_Ecl)*P(i)+E)*(Tpp1_Ecl(i)^2)+(AA*CaO_Ecl+BB)*Tpp1_Ecl(i)+5);
                    end
                end
                
                % peridotite
    
                    % F0-cpx(T') for peridotite
                    if Tprim_PdtT(i)<0
                        F0_cpxTprim_Pdt(i)=0;
                    else
                        F0_cpxTprim_Pdt(i)=(Tprim_PdtT(i)^(1.5))*100;
                    end
            
                    % F-cpx(T'') peridotite
                    if Tpp1_Pdt(i)<0
                        F_cpxTpp_Pdt(i)=0;
                    else
                        F_cpxTpp_Pdt(i)=(Tpp1_Pdt(i)^(1.5))*100;
                    end
                    
                    % F0(T') peridotite
                    if F0_cpxTprim_Pdt(i)<Fcpxout_Pdt(i)*100
                        F0Tprim_Pdt(i)=F0_cpxTprim_Pdt(i);
                    else
                        F0Tprim_Pdt(i)=(Fcpxout_Pdt(i)+(1-Fcpxout_Pdt(i))*((T(i)-Tcpxout_Pdt(i))/(Tliq_Pdt(i)-Tcpxout_Pdt(i)))^1.5)*100;
                    end
                   
                    % F(T'') peridotite
                    if F_cpxTpp_Pdt(i)<Fcpxout_Pdt(i)*100
                        FTpp_Pdt(i)=F_cpxTpp_Pdt(i);
                    else
                        FTpp_Pdt(i)=(Fcpxout_Pdt(i)+(1-Fcpxout_Pdt(i))*((T(i)-Tcpxout_Pdt(i))/(Tliq_Pdt(i)-Tcpxout_Pdt(i)))^1.5);
                    end
                
            % melt fraction
    
                % pyroxenite
                F_Pyr(i)=F_Pyr_SA(i);
     
                % eclogite
                F_Ecl(i)=F_Ecl_SA(i);
    
                % peridotite
                F_Pdt(i)=F0Tprim_Pdt(i);
    
            % dF/dT_P
    
                % pyroxenite
                if T(i)>Tcpxout_Pyr(i)
                    dFdT_P_Pyr(i)=defaultdFdTcpxout_Pyr;
                else
                    dFdT_P_Pyr(i)=FTpp1_Pyr(i)-F_Pyr(i);
                end
    
                if dFdT_P_Pyr(i)<0
                    moddFdT_P_Pyr(i)=0.0001;
                else
                    moddFdT_P_Pyr(i)=dFdT_P_Pyr(i);
                end
    
                % eclogite
                if T(i)>Tcpxout_Ecl(i)
                    dFdT_P_Ecl(i)=defaultdFdTcpxout_Ecl;
                else
                    dFdT_P_Ecl(i)=FTpp1_Ecl(i)-F_Ecl(i);
                end
    
                if dFdT_P_Ecl(i)<0
                    moddFdT_P_Ecl(i)=0.0001;
                else
                    moddFdT_P_Ecl(i)=dFdT_P_Ecl(i);
                end
    
                % peridotite
                dFdT_P_Pdt(i)=FTpp_Pdt(i)-F_Pdt(i);   
    
            % dTdP_F
    
                % pyroxenite
                dTdP_F_Pyr(i)=dTsdP_s_Pyr(i);
                if F_Pyr(i)<=0
                    dTdP_F_PyrT(i)=dTsdP_s_Pyr(i);
                else
                    dTdP_F_PyrT(i)=dTdP_F_Pyr(i);
                end
        
                % eclogite
                dTdP_F_Ecl(i)=dTsdP_s_Ecl(i);
                if F_Ecl(i)<=0
                    dTdP_F_EclT(i)=dTsdP_s_Ecl(i);
                else
                    dTdP_F_EclT(i)=dTdP_F_Ecl(i);
                end
    
                % peridotite
                dTdP_F_Pdt(i)=dTsdP_s_Pdt(i);
                if F_Pdt(i)<=0
                    dTdP_F_PdtT(i)=dTsdP_s_Pdt(i);
                else
                    dTdP_F_PdtT(i)=dTdP_F_Pdt(i);
                end
            
            % thermal equilibrium (Phipps Morgan 2001 Eq. 26-29)
        
                % pyroxenite
                dFdP1_Pyr(i)=0;
                dFdP2_Pyr(i)=0;
                dFdP3_Pyr(i)=0;
                dFdP4_Pyr(i)=0;
                dFdP_S_Pyr(i)=0;
              
                % eclogite
                dFdP1_Ecl(i)=0;
                dFdP2_Ecl(i)=0;
                dFdP3_Ecl(i)=0;
                dFdP4_Ecl(i)=0;
                dFdP_S_Ecl(i)=0;
        
                % peridotite
                dFdP1_Pdt(i)=0;
                dFdP2_Pdt(i)=0;
                dFdP3_Pdt(i)=0;
                dFdP_S_Pdt(i)=0;
            
            % change in temperature with pressure
            if (F_Pdt(i)==0) && (F_Pyr(i)==0) && (F_Ecl(i)==0)
                DTDP(i)=(T_SA(i+1)-T_SA(i))/(P(i+1)-P(i));
            else
                if (F_Pyr(i)==0) && (F_Ecl(i)==0)
                    DTDP(i)=dTdP_F_PdtT(i)+(1/dFdT_P_Pdt(i))*dFdP3_Pdt(i);
                else
                    if F_Ecl(i)==0
                        DTDP(i)=dTdP_F_PyrT(i)+(1/moddFdT_P_Pyr(i))*dFdP4_Pyr(i);
                    else
                        DTDP(i)=dTdP_F_EclT(i)+(1/moddFdT_P_Ecl(i))*dFdP4_Ecl(i);
                    end
                end
            end
            
            % dT/dP_S
            if (FTprim_Pyr(i)>=100) && (FTprim_Ecl(i)>=100)
                dTdP_S_final(i)=dTdP_F_PdtT(i)+(1/dFdT_P_Pdt(i))*dFdP3_Pdt(i);
            else
                dTdP_S_final(i)=DTDP(i);
            end
            
        else
        %% i>1
            
            % temperature of solid mantle
            if (Tsol_Pyr(i)>T_SA(i)) && (Tsol_Ecl(i)>T_SA(i)) && (Tsol_Pdt(i)>T_SA(i))
                T(i)=T_SA(i);
            else
                T(i)=T(i-1)+dTdP_S_final(i-1)*(P(i)-P(i-1));
            end
    
            % T' 
    
                % pyroxenite
                if T5_Pyr(i)>Tcpxout_Pyr(i) 
                    Tprim_PyrT(i)=-1;
                else
                    Tprim_PyrT(i)=(T(i)-T5_Pyr(i))/(Tcpxout_Pyr(i)-T5_Pyr(i));
                end
    
                % pyroxenite
                if T5_Ecl(i)>Tcpxout_Ecl(i) 
                    Tprim_EclT(i)=-1;
                else
                    Tprim_EclT(i)=(T(i)-T5_Ecl(i))/(Tcpxout_Ecl(i)-T5_Ecl(i));
                end          
            
                % peridotite
                Tprim_PdtT(i)=(T(i)-Tsol_Pdt(i))/(Tliq_Lhz(i)-Tsol_Pdt(i));
    
            % F(T') Eq. 5 & 6 into Eq. 1 (Lambart et al. 2016)
    
                % pyroxenite
                if Tprim_PyrT(i)<Tprim_Pyr(i) 
                    FTprim_Pyr(i)=0;
                else
                    if (((A+B*MgNo_Pyr)*(P(i)^2)+(C+D*MgNo_Pyr)*P(i)+E)*(Tprim_PyrT(i)^2)+(AA*CaO_Pyr+BB)*Tprim_PyrT(i)+5)<0
                        FTprim_Pyr(i)=0;
                    else
                        FTprim_Pyr(i)=(((A+B*MgNo_Pyr)*(P(i)^2)+(C+D*MgNo_Pyr)*P(i)+E)*(Tprim_PyrT(i)^2)+(AA*CaO_Pyr+BB)*Tprim_PyrT(i)+5);
                    end
                end
    
                % eclogite
                if Tprim_EclT(i)<Tprim_Ecl(i) 
                    FTprim_Ecl(i)=0;
                else
                    if (((A+B*MgNo_Ecl)*(P(i)^2)+(C+D*MgNo_Ecl)*P(i)+E)*(Tprim_EclT(i)^2)+(AA*CaO_Ecl+BB)*Tprim_EclT(i)+5)<0
                        FTprim_Ecl(i)=0;
                    else
                        FTprim_Ecl(i)=(((A+B*MgNo_Ecl)*(P(i)^2)+(C+D*MgNo_Ecl)*P(i)+E)*(Tprim_EclT(i)^2)+(AA*CaO_Ecl+BB)*Tprim_EclT(i)+5);
                    end
                end
           
            % T''(T+1)
    
                % pyroxenite
                if T5_Pyr(i)>Tcpxout_Pyr(i) 
                    Tpp1_Pyr(i)=-1;
                else
                    Tpp1_Pyr(i)=Tprim_PyrT(i)+(1/(Tcpxout_Pyr(i)-T5_Pyr(i)));
                end
    
                % eclogite
                if T5_Ecl(i)>Tcpxout_Ecl(i) 
                    Tpp1_Ecl(i)=-1;
                else
                    Tpp1_Ecl(i)=Tprim_EclT(i)+(1/(Tcpxout_Ecl(i)-T5_Ecl(i)));
                end
    
                % peridotite
                Tpp1_Pdt(i)=(T(i)+1-Tsol_Pdt(i))/(Tliq_Lhz(i)-Tsol_Pdt(i));
    
            % F(T'')
    
                % pyroxenite
                if Tpp1_Pyr(i)<Tprim_Pyr(i)
                    FTpp1_Pyr(i)=0;
                else
                    if (((A+B*MgNo_Pyr)*(P(i)^2)+(C+D*MgNo_Pyr)*P(i)+E)*(Tpp1_Pyr(i)^2)+(AA*CaO_Pyr+BB)*Tpp1_Pyr(i)+5)<0
                        FTpp1_Pyr(i)=0;
                    else
                        FTpp1_Pyr(i)=(((A+B*MgNo_Pyr)*(P(i)^2)+(C+D*MgNo_Pyr)*P(i)+E)*(Tpp1_Pyr(i)^2)+(AA*CaO_Pyr+BB)*Tpp1_Pyr(i)+5);
                    end
                end
    
                % eclogite
                if Tpp1_Ecl(i)<Tprim_Ecl(i)
                    FTpp1_Ecl(i)=0;
                else
                    if (((A+B*MgNo_Ecl)*(P(i)^2)+(C+D*MgNo_Ecl)*P(i)+E)*(Tpp1_Ecl(i)^2)+(AA*CaO_Ecl+BB)*Tpp1_Ecl(i)+5)<0
                        FTpp1_Ecl(i)=0;
                    else
                        FTpp1_Ecl(i)=(((A+B*MgNo_Ecl)*(P(i)^2)+(C+D*MgNo_Ecl)*P(i)+E)*(Tpp1_Ecl(i)^2)+(AA*CaO_Ecl+BB)*Tpp1_Ecl(i)+5);
                    end
                end
    
                % peridotite
    
                    % F0-CPX (T') for peridotite
                    if Tprim_PdtT(i)<0 
                        F0_cpxTprim_Pdt(i)=0;
                    else
                        F0_cpxTprim_Pdt(i)=(Tprim_PdtT(i)^(1.5))*100;
                    end
                                    
                    % F-cpx(T'') peridotite
                    if Tpp1_Pdt(i)<0
                        F_cpxTpp_Pdt(i)=0;
                    else
                        F_cpxTpp_Pdt(i)=(Tpp1_Pdt(i)^(1.5))*100;
                    end
                    
                    % F0(T') peridotite
                    if F0_cpxTprim_Pdt(i)<Fcpxout_Pdt(i)*100
                        F0Tprim_Pdt(i)=F0_cpxTprim_Pdt(i);
                    else
                        F0Tprim_Pdt(i)=(Fcpxout_Pdt(i)+((1-Fcpxout_Pdt(i))*((T(i)-Tcpxout_Pdt(i))/(Tliq_Pdt(i)-Tcpxout_Pdt(i))))^1.5)*100;
                    end
            
                    % F(T'') peridotite
                    if F_cpxTpp_Pdt(i)<Fcpxout_Pdt(i)*100
                        FTpp_Pdt(i)=F_cpxTpp_Pdt(i);
                    else
                        FTpp_Pdt(i)=(Fcpxout_Pdt(i)+((1-Fcpxout_Pdt(i))*((T(i)+1-Tcpxout_Pdt(i))/(Tliq_Pdt(i)-Tcpxout_Pdt(i))))^1.5)*100;
                    end
            
                % F_Pyr
    
                    % pyroxenite
                    if FTprim_Pyr(i-1)>0
                        F_Pyr(i)=F_Pyr(i-1)+dFdP4_Pyr(i-1)*(P(i)-P(i-1));
                    else
                        F_Pyr(i)=FTprim_Pyr(i);
                    end
       
                    % eclogite
                    if FTprim_Ecl(i-1)>0
                        F_Ecl(i)=F_Ecl(i-1)+dFdP4_Ecl(i-1)*(P(i)-P(i-1));
                    else
                        F_Ecl(i)=FTprim_Ecl(i);
                    end
    
                    % peridotite
                    if F0Tprim_Pdt(i-1)<=0
                        F_Pdt(i)=F0Tprim_Pdt(i);
                    else
                        F_Pdt(i)=F_Pdt(i-1)+dFdP3_Pdt(i-1)*(P(i)-P(i-1));
                    end
            
            % dF/dT_P
    
                % pyroxenite
                if T(i)>Tcpxout_Pyr(i)
                    dFdT_P_Pyr(i)=defaultdFdTcpxout_Pyr;
                else
                    dFdT_P_Pyr(i)=FTpp1_Pyr(i)-F_Pyr(i);
                end
                
                if dFdT_P_Pyr(i)<0
                    moddFdT_P_Pyr(i)=0.0001;
                else
                    moddFdT_P_Pyr(i)=dFdT_P_Pyr(i);
                end
           
                % eclogite
                if T(i)>Tcpxout_Ecl(i)
                    dFdT_P_Ecl(i)=defaultdFdTcpxout_Ecl;
                else
                    dFdT_P_Ecl(i)=FTpp1_Ecl(i)-F_Ecl(i);
                end
                
                if dFdT_P_Ecl(i)<0
                    moddFdT_P_Ecl(i)=0.0001;
                else
                    moddFdT_P_Ecl(i)=dFdT_P_Ecl(i);
                end
    
                % peridotite
                dFdT_P_Pdt(i)=FTpp_Pdt(i)-F_Pdt(i);
            
            % dTdP_F
    
                % pyroxenite
    
                    % calculate Tm and T'm first
                    TpM_Pyr(i-1)=(-(AA*CaO_Pyr+BB)+(((AA*CaO_Pyr+BB)^2)-4*(5-FTprim_Pyr(i))*((A+B*MgNo_Pyr)*(P(i-1)^2)+(C+D*MgNo_Pyr)*P(i-1)+E))^0.5)/(2*((A+B*MgNo_Pyr)*(P(i-1)^2)+(C+D*MgNo_Pyr)*P(i-1)+E));
                    TM_Pyr(i-1)=TpM_Pyr(i-1)*(Tcpxout_Pyr(i-1)-T5_Pyr(i-1))+T5_Pyr(i-1);
                    
                    dTdP_F_Pyr(i)=dTsdP_s_Pyr(i-1);
                    if F_Pyr(i)<=0
                        dTdP_F_PyrT(i)=dTsdP_s_Pyr(i);
                    else
                        dTdP_F_PyrT(i)=dTdP_F_Pyr(i);
                    end
    
                % eclogite
    
                    % calculate Tm and T'm first
                    TpM_Ecl(i-1)=(-(AA*CaO_Ecl+BB)+(((AA*CaO_Ecl+BB)^2)-4*(5-FTprim_Ecl(i))*((A+B*MgNo_Ecl)*(P(i-1)^2)+(C+D*MgNo_Ecl)*P(i-1)+E))^0.5)/(2*((A+B*MgNo_Ecl)*(P(i-1)^2)+(C+D*MgNo_Ecl)*P(i-1)+E));
                    TM_Ecl(i-1)=TpM_Ecl(i-1)*(Tcpxout_Ecl(i-1)-T5_Ecl(i-1))+T5_Ecl(i-1);
                    
                    dTdP_F_Ecl(i)=dTsdP_s_Ecl(i-1);
                    if F_Ecl(i)<=0
                        dTdP_F_EclT(i)=dTsdP_s_Ecl(i);
                    else
                        dTdP_F_EclT(i)=dTdP_F_Ecl(i);
                    end
    
                % peridotite
        
                    % calculate Tm and T'm first
                    TpM_Pdt(i-1)=(F0Tprim_Pdt(i)/100)^(1/1.5);
                    if F_Pdt(i-1)<(Fcpxout_Pdt(i-1)*100)
                        TM_Pdt(i-1)=TpM_Pdt(i-1)*(Tliq_Lhz(i-1)-Tsol_Pdt(i-1))+Tsol_Pdt(i-1);
                    else
                        TM_Pdt(i-1)=TpM_Pdt(i-1)*(Tliq_Pdt(i-1)-Tcpxout_Pdt(i-1))+Tcpxout_Pdt(i-1);
                    end
            
                    if TpM_Pdt(i-1)<Tprim_PdtT(i-1)
                        dTdP_F_Pdt(i)=dTdP_F_Pdt(i-1);
                    else
                        dTdP_F_Pdt(i)=(T(i)-TM_Pdt(i-1))/(P(i)-P(i-1));
                    end
                    if F_Pdt(i)<=0
                        dTdP_F_PdtT(i)=dTsdP_s_Pdt(i);
                    else
                        dTdP_F_PdtT(i)=dTdP_F_Pdt(i);
                    end
                        
            % thermal equilibrium (Phipps Morgan 2001 Eq. 26-29)
       
                % peridotite
        
                    if (FTprim_Pyr(i)<=0) || (F_Pyr(i)>=100) && (FTprim_Pyr(i)<=0) || (F_Pyr(i)>=100) % 1 lithology (original Pdt)
                        dFdP1_Pdt(i)=((dTdP_F_PdtT(i)-(alpha_Pdt*(Tp+TK)/(PropPyr*rho_Pyr*Cp_Pyr+PropPdt*rho_Pdt*Cp_Pdt)*PGPa))/(PropPdt*(Tp+TK)*S_Pdt/(PropPyr*Cp_Pyr+PropPdt*Cp_Pdt)+(100/dFdT_P_Pdt(i))))*(-100);
                    elseif (FTprim_Ecl(i)<=0) || (F_Ecl(i)>=100) % 2 lithologies (here is original Pdt-Pyr)
                        dFdP1_Pdt(i)=(dTdP_F_PdtT(i)-(PropPdt*alpha_Pdt+PropPyr*alpha_Pyr)*((Tp+TK)/(PropPyr*rho_Pyr*Cp_Pyr+PropPdt*rho_Pdt*Cp_Pdt)*PGPa)+PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)*(dTdP_F_PdtT(i)-dTdP_F_PyrT(i))/(100/moddFdT_P_Pyr(i)))/(PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)+PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)*(100/dFdT_P_Pdt(i))/(100/moddFdT_P_Pyr(i))+(100/dFdT_P_Pdt(i)))*(-100);
                    elseif (FTprim_Pyr(i)<=0) || (F_Pyr(i)>=100) % 2 lithologies (original replaced by Pdt-Ecl)
                        dFdP1_Pdt(i)=(dTdP_F_PdtT(i)-(PropPdt*alpha_Pdt+PropEcl*alpha_Ecl)*((Tp+TK)/(PropEcl*rho_Ecl*Cp_Ecl+PropPdt*rho_Pdt*Cp_Pdt)*PGPa)+PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)*(dTdP_F_PdtT(i)-dTdP_F_EclT(i))/(100/moddFdT_P_Ecl(i)))/(PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)+PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)*(100/dFdT_P_Pdt(i))/(100/moddFdT_P_Ecl(i))+(100/dFdT_P_Pdt(i)))*(-100);
                    else % 3 lithologies
                        dFdP1_Pdt(i)=...
                            (...
                            dTdP_F_PdtT(i)-... % n1 (current=1)
                            (PropPdt*alpha_Pdt+PropPyr*alpha_Pyr+PropEcl*alpha_Ecl)*... % nTa
                            ((Tp+TK)/(PropPdt*rho_Pdt*Cp_Pdt+PropPyr*rho_Pyr*Cp_Pyr+PropEcl*rho_Ecl*Cp_Ecl)*PGPa)... %nTb
                            +PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)*... %n2a (2)
                            (dTdP_F_PdtT(i)-dTdP_F_PyrT(i))/(100/moddFdT_P_Pyr(i))... %n2b (current=1 - 2)
                            +PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)*... %n3a (3)
                            (dTdP_F_PdtT(i)-dTdP_F_EclT(i))/(100/moddFdT_P_Ecl(i))... %n3b (current=1 - 3)
                            )...
                            /...
                            (...
                            PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)... %d1a (current=1)
                            +PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)*... &d2a (2)
                            (100/dFdT_P_Pdt(i))/(100/moddFdT_P_Pyr(i))... &d2b (current=1 / 2)
                            +PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)*... &d3a (3)
                            (100/dFdT_P_Pdt(i))/(100/moddFdT_P_Ecl(i))... &d3b (current=1 / 3)
                            +(100/dFdT_P_Pdt(i))... %d1b
                            )*(-100); 
                    end
                
                if F_Pdt(i)<=0
                    dFdP2_Pdt(i)=0;
                else
                    dFdP2_Pdt(i)=dFdP1_Pdt(i);
                end
                
                if dFdT_P_Pdt(i)<0 || dFdP2_Pdt(i)>=0
                    dFdP3_Pdt(i)=0;
                else
                    dFdP3_Pdt(i)=dFdP2_Pdt(i);
                end
                
                dFdP_S_Pdt(i)=-dFdP3_Pdt(i);
    
                % pyroxenite
        
                    if dFdP_S_Pdt(i)<=0 && (FTprim_Ecl(i)<=0) || (F_Ecl(i)>=100) % 1 lithology (original Pyr)
                        dFdP1_Pyr(i)=((dTdP_F_PyrT(i)-((alpha_Pyr*(Tp+TK)/(PropPyr*rho_Pyr*Cp_Pyr+PropPdt*rho_Pdt*Cp_Pdt))*PGPa))/(PropPyr*(Tp+TK)*S_Pyr/(PropPyr*Cp_Pyr+PropPdt*Cp_Pdt)+(100/moddFdT_P_Pyr(i))))*(-100);
                    elseif (FTprim_Ecl(i)<=0) || (F_Ecl(i)>=100) % 2 lithologies (original Pdt-Pyr)
                        dFdP1_Pyr(i)=((dTdP_F_PyrT(i)-(PropPyr*alpha_Pyr+PropPdt*alpha_Pdt)*((Tp+TK)/(PropPyr*rho_Pyr*Cp_Pyr+PropPdt*rho_Pdt*Cp_Pdt))*PGPa+PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)*(dTdP_F_PyrT(i)-dTdP_F_PdtT(i))/(100/dFdT_P_Pdt(i)))/(PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)+PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)*(100/moddFdT_P_Pyr(i))/(100/dFdT_P_Pdt(i))+(100/moddFdT_P_Pyr(i))))*(-100);
                    elseif dFdP_S_Pdt(i)<=0 % 2 lithologies (original replaced by Pyr-Ecl)
                        dFdP1_Pyr(i)=((dTdP_F_PyrT(i)-(PropPyr*alpha_Pyr+PropEcl*alpha_Ecl)*((Tp+TK)/(PropPyr*rho_Pyr*Cp_Pyr+PropEcl*rho_Ecl*Cp_Ecl))*PGPa+PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)*(dTdP_F_PyrT(i)-dTdP_F_EclT(i))/(100/dFdT_P_Ecl(i)))/(PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)+PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)*(100/moddFdT_P_Pyr(i))/(100/dFdT_P_Ecl(i))+(100/moddFdT_P_Pyr(i))))*(-100);
                    else % 3 lithologies
                        dFdP1_Pyr(i)=...
                            (...
                            (dTdP_F_PyrT(i)-... % n1 (current=2)
                            (PropPyr*alpha_Pyr+PropPdt*alpha_Pdt+PropPdt*alpha_Ecl)*... % nTa
                            ((Tp+TK)/(PropPyr*rho_Pyr*Cp_Pyr+PropPdt*rho_Pdt*Cp_Pdt+PropEcl*rho_Ecl*Cp_Ecl)*PGPa)... %nTb
                            +PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)*... %n2a (1)
                            (dTdP_F_PyrT(i)-dTdP_F_PdtT(i))/(100/dFdT_P_Pdt(i))... %n2b (current=2 - 1)
                            +PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)*... %n3a (3)
                            (dTdP_F_PyrT(i)-dTdP_F_EclT(i))/(100/dFdT_P_Ecl(i))... %n3b (current=2 - 3)
                            )...
                            /...
                            (...
                            PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)... %d1a (current=2)
                            +PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)*... &d2a (1)
                            (100/moddFdT_P_Pyr(i))/(100/dFdT_P_Pdt(i))... &d2b (current=2 / 1) 
                            +PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)*... &d3a (3)
                            (100/moddFdT_P_Pyr(i))/(100/dFdT_P_Ecl(i))... &d2b (current=2 / 3) 
                            +(100/moddFdT_P_Pyr(i)))... %d1b (current=2)
                            )*(-100);
                    end
                
                if moddFdT_P_Pyr(i)==0
                    dFdP2_Pyr(i)=0;
                else 
                    if F_Pyr_SA(i+1)==0
                        dFdP2_Pyr(i)=0;
                    else
                        dFdP2_Pyr(i)=dFdP1_Pyr(i);
                    end
                end
                
                if dFdP2_Pyr(i)>=0
                    dFdP3_Pyr(i)=0;
                else
                    dFdP3_Pyr(i)=dFdP2_Pyr(i);
                end
                
                if (F_Pyr(i)>100)
                    dFdP4_Pyr(i)=0;
                else
                    dFdP4_Pyr(i)=dFdP3_Pyr(i);
                end
                
                dFdP_S_Pyr(i)=-dFdP4_Pyr(i);
           
                % eclogite
        
                    if dFdP_S_Pdt(i)<=0 && (FTprim_Pyr(i)<=0) || (F_Pyr(i)>=100) % 1 lithology (original replaced by Ecl)
                        dFdP1_Ecl(i)=((dTdP_F_EclT(i)-((alpha_Ecl*(Tp+TK)/(PropPyr*rho_Ecl*Cp_Ecl+PropPdt*rho_Pdt*Cp_Pdt))*PGPa))/(PropPyr*(Tp+TK)*S_Ecl/(PropPyr*Cp_Ecl+PropPdt*Cp_Pdt)+(100/moddFdT_P_Ecl(i))))*(-100);
                    elseif (FTprim_Pyr(i)<=0) || (F_Pyr(i)>=100) % 2 lithologies (original replaced by Pdt-Ecl)
                        dFdP1_Ecl(i)=((dTdP_F_EclT(i)-(PropEcl*alpha_Ecl+PropPdt*alpha_Pdt)*((Tp+TK)/(PropEcl*rho_Ecl*Cp_Ecl+PropPdt*rho_Pdt*Cp_Pdt))*PGPa+PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)*(dTdP_F_EclT(i)-dTdP_F_PdtT(i))/(100/dFdT_P_Pdt(i)))/(PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)+PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)*(100/moddFdT_P_Ecl(i))/(100/dFdT_P_Pdt(i))+(100/moddFdT_P_Ecl(i))))*(-100);
                    elseif dFdP_S_Pdt(i)<=0 % 2 lithologies (original replaced by Pyr-Ecl)
                        dFdP1_Ecl(i)=((dTdP_F_EclT(i)-(PropEcl*alpha_Ecl+PropPyr*alpha_Pyr)*((Tp+TK)/(PropEcl*rho_Ecl*Cp_Ecl+PropPyr*rho_Pyr*Cp_Pyr))*PGPa+PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)*(dTdP_F_EclT(i)-dTdP_F_PyrT(i))/(100/dFdT_P_Pyr(i)))/(PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)+PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)*(100/moddFdT_P_Ecl(i))/(100/dFdT_P_Pyr(i))+(100/moddFdT_P_Ecl(i))))*(-100);
                    else % 3 lithologies
                        dFdP1_Ecl(i)=...
                            (...
                            (dTdP_F_EclT(i)-... % n1 (current=3)
                            (PropPyr*alpha_Pyr+PropPdt*alpha_Pdt+PropPdt*alpha_Ecl)*... % nTa
                            ((Tp+TK)/(PropPyr*rho_Pyr*Cp_Pyr+PropPdt*rho_Pdt*Cp_Pdt+PropEcl*rho_Ecl*Cp_Ecl)*PGPa)... %nTb
                            +PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)*... %n2a (1)
                            (dTdP_F_EclT(i)-dTdP_F_PdtT(i))/(100/dFdT_P_Pdt(i))... %n2b (current=3 - 1)
                            +PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)*... %n3a (2)
                            (dTdP_F_EclT(i)-dTdP_F_PyrT(i))/(100/moddFdT_P_Pyr(i))... %n3b (current=3 - 2)
                            )...
                            /...
                            (...
                            PropEcl*(Tp+TK)*(S_Ecl/Cp_Ecl)+... %d1a (current=3)
                            +PropPdt*(Tp+TK)*(S_Pdt/Cp_Pdt)*... &d2a (1)
                            (100/moddFdT_P_Ecl(i))/(100/dFdT_P_Pdt(i))... &d2b (current=3 / 1) 
                            +PropPyr*(Tp+TK)*(S_Pyr/Cp_Pyr)*... &d3a (2)
                            (100/moddFdT_P_Ecl(i))/(100/moddFdT_P_Pyr(i))... &d2b (current=3 / 2) 
                            +(100/moddFdT_P_Ecl(i)))... %d1b (current=3)
                            )*(-100);
                    end
                
                if moddFdT_P_Ecl(i)==0
                    dFdP2_Ecl(i)=0;
                else 
                    if F_Ecl_SA(i+1)==0
                        dFdP2_Ecl(i)=0;
                    else
                        dFdP2_Ecl(i)=dFdP1_Ecl(i);
                    end
                end
                
                if dFdP2_Ecl(i)>=0
                    dFdP3_Ecl(i)=0;
                else
                    dFdP3_Ecl(i)=dFdP2_Ecl(i);
                end
                
                if (F_Ecl(i)>100)
                    dFdP4_Ecl(i)=0;
                else
                    dFdP4_Ecl(i)=dFdP3_Ecl(i);
                end
                
                dFdP_S_Ecl(i)=-dFdP4_Ecl(i);
           
            % change in temperature with pressure
            if (F_Pdt(i)==0) && (F_Pyr(i)==0) && (F_Ecl(i)==0)
                DTDP(i)=(T_SA(i+1)-T_SA(i))/(P(i+1)-P(i));
            else
                if (F_Pyr(i)==0) && (F_Ecl(i)==0)
                    DTDP(i)=dTdP_F_PdtT(i)+(1/dFdT_P_Pdt(i))*dFdP3_Pdt(i);
                else
                    if F_Ecl(i)==0
                        DTDP(i)=dTdP_F_PyrT(i)+(1/moddFdT_P_Pyr(i))*dFdP4_Pyr(i);
                    else
                        DTDP(i)=dTdP_F_EclT(i)+(1/moddFdT_P_Ecl(i))*dFdP4_Ecl(i);
                    end
                end
            end
    
            % dT/dP_S
            if (FTprim_Pyr(i)>=100) && (FTprim_Ecl(i)>=100)
                dTdP_S_final(i)=dTdP_F_PdtT(i)+(1/dFdT_P_Pdt(i))*dFdP3_Pdt(i);
            else
                dTdP_S_final(i)=DTDP(i);
            end
    
        end
    
    end
    
    Temperature=T;
    Pressure=P(1:end-1);
    MeltingD_Pdt=F_Pdt;
    MeltingD_Pyr=F_Pyr;
    MeltingD_Ecl=F_Ecl;
    Tadiabat=T_SA(1:end-1);
    
end