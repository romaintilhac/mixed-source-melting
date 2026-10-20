function [CL, CS] = meltingEquations(C0,F,D,varargin)
% Calculates liquid and solid compositions based melting model equations

    % batch=0
    % fractional=1;
    % dynamic=2;
    
    % check if phi is provided
    if numel(varargin) > 1
        error('Too many optional arguments. Provide only one optional argument.');
    elseif numel(varargin) == 1
        phi = varargin{1};
    elseif isnan(varargin{1})
        phi = []; % Default value if the optional argument isnan
    else
        phi = []; % Default value if the optional argument is not provided
    end
    
    if phi==1
        CL=(C0)./(F+D.*(1-F)); %liquid
        CS=CL.*D; % residual solid
    elseif phi==0
        Cl=(C0./D).*(1-F).^((1./D)-1); %instantaneous liquid
        CL=(C0./F).*(1-(1-F).^(1./D)); %accumulated liquid
        CS=C0.*(1-F).^((1./D)-1); %residual solid
    else
        phi=F*phi;
        X=(F-phi)./(1-phi);
        Cl=(C0./(phi+(1-phi).*D)).*(1-X).^((1./(phi+(1-phi).*D))-1); %instantaneous liquid
        CL=(C0./X).*(1-(1-X).^(1./(phi+(1-phi)*D))); %accumulated liquid
        Cs=D.*Cl; %residual solid
        CS=phi.*Cl+(1-phi).*Cs; %total residual solid
    end

end