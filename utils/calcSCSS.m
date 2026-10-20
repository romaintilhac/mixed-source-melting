function [SCSS]=calcSCSS(ME_melt,ME_sul,P,T,ref,ideal, lst)
% Calculates Sulfur Concentration at Sulfide Saturation (SCSS)

    % If ref=0 (Ding et al. 2017), only ME_melt and T are required to calculated sulfide-free SCSS.
    % If ref=1 (Smythe et al. 2017), ME_sul, P and ideal (0/1) are also required to calculated SCSS in the presence of sulfide.

    struct=getTableOxi; TableOxi=struct.TableOxi;

    T=T+constants.TK;
    
    ME_melt=table2array(ME_melt);
    
    if height(ME_melt)>1; error('SCSS requires a single set of compositions');end
    
    if ref==0 % Ding et al. 2017
    
        if nargin>3
            warning('SCSS calculation based on Ding et al. 2017 is only based on melt compositions and T. Other arguments provided will be ignored.')
        end
    
        ME_melt=ME_melt./lst.molmass; tot=sum(ME_melt);
        ME_melt=100.*ME_melt./tot; tot=sum(ME_melt);
        ME_melt=ME_melt.*lst.catnum; tot=sum(ME_melt);
        ME_melt=ME_melt./tot;
        ME_melt=array2table(ME_melt,VariableNames=lst.ME);
        
        SCSS=exp(12.10023817-4951.220517/T+(4.02527185*ME_melt.Ti+4.173632609.*ME_melt.Ca-3.643865073.*ME_melt.Si-3.936000202.*ME_melt.Al+5.574892678.*(ME_melt.Fe2+ME_melt.Fe3))-40.67763841*(ME_melt.Fe2+ME_melt.Fe3)*ME_melt.Ti-273.4844764*(1/T));
    
    elseif ref==1 % Smythe et al. 2017
    
        if nargin<6
            error('SCSS calculation based on Smythe et al. 2017 requires melt and sulfide compositions, P, T and a type of solution (ideal/non-ideal).')
        end
    
        Oxi=(TableOxi.O.mass/2)*sum(ME_melt./lst.moxmass.*lst.charge);
        ME_melt=ME_melt./lst.moxmass.*lst.atomass; %Silicate Weight Percent Oxide
        tot=sum(ME_melt)+Oxi;
        Oxi=Oxi/TableOxi.O.mass;
        ME_melt=ME_melt./lst.atomass; %Silicate Weight Percent Elemental
        tot=sum(ME_melt)+Oxi;
        Oxi=Oxi/tot*100;
        ME_melt=ME_melt/tot*100; %Silicate Atomic Percent
        tot=sum(ME_melt);
        ME_melt=ME_melt./tot; %Silicate Mole Fraction Oxide
        
        ME_melt=array2table(ME_melt,VariableNames=lst.ME);
        
        %Sulfide Mole Fraction
        XsulFe=ME_sul.Fe/TableOxi.Fe.mass; 
        XsulNi=ME_sul.Ni/TableOxi.Ni.mass;
        XsulCu=ME_sul.Cu/TableOxi.Cu.mass;
        XsulO=ME_sul.O/TableOxi.O.mass;
        XsulS=ME_sul.S/TableOxi.S.mass;
        tot=XsulFe+XsulNi+XsulCu+XsulO+XsulS;
        
        %Sulfide Atomic Percent
        XsulFe=XsulFe/tot*100; 
        XsulNi=XsulNi/tot*100;
        XsulCu=XsulCu/tot*100;
        XsulO=XsulO/tot*100;
        XsulS=XsulS/tot*100;
        
        XsulFeS=XsulFe/(XsulFe+XsulNi+XsulCu);
        XsulNiS=XsulNi/(XsulFe+XsulNi+XsulCu);
        XsulCuS=XsulCu/(XsulFe+XsulNi+XsulCu);
        
        if ideal==1 %ideal
            B=9.087; C=-269.4;
            ASi=-27561; ATi=-11220; AAl=-18450; AMg=-13970; ACa=-7831; AFe=-34274; ANa=-13247; AK=-29015; AH=-17495; ASiFe=116568;
            sumXA=(ASi*ME_melt.Si+ATi*ME_melt.Ti+AAl*ME_melt.Al+AMg*ME_melt.Mg+ACa*ME_melt.Ca+AFe*ME_melt.Fe2+ANa*ME_melt.Na+AK*ME_melt.K+ASiFe*(ME_melt.Si*ME_melt.Fe2));
            SCSS=exp((122175-80.28*T+8.474*T*log(T))/(8.314*T)+B+sumXA/T+log(XsulFeS)-log(ME_melt.Fe2)+C*(P/T));
        elseif ideal==0 %non-ideal
            B=9.352; C=-264.85; E=546;
            ASi=-27996; ATi=-10715; AAl=-19000; AMg=-14512; ACa=-8832; AFe=-34895; ANa=-13713; AK=-28584; AH=-17766; ASiFe=117816;
            sumXA=(ASi*ME_melt.Si+ATi*ME_melt.Ti+AAl*ME_melt.Al+AMg*ME_melt.Mg+ACa*ME_melt.Ca+AFe*ME_melt.Fe2+ANa*ME_melt.Na+AK*ME_melt.K+ASiFe*(ME_melt.Si*ME_melt.Fe2));
            SCSS=exp((122175-80.28*T+8.474*T*log(T))/(8.314*T)+B+sumXA/T+log(XsulFeS)-log(ME_melt.Fe2)+C*(P/T)+((XsulCuS^2+XsulCuS*XsulNiS)/T)*E);
        end
    
    end

end