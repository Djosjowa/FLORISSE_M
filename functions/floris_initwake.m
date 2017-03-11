function [ wake ] = floris_initwake( model,turbine,wake,turb_type )
% This function computes the coefficients that determine wake behaviour.
% The initial deflection and diameter of the wake are also computed

    % Calculate ke, the basic expansion coefficient
    wake.Ke = model.Ke + model.KeCorrCT*(turbine.Ct-model.baselineCT);

    % Calculate mU, the zone multiplier for different wake zones
    if model.useaUbU
        wake.mU = model.MU/cosd(model.aU+model.bU*turbine.YawWF);
    else
        wake.mU = model.MU;
    end

    % Calculate initial wake deflection due to blade rotation etc.
%     wake.zetaInit = 0.5*(cosd(turbine.YawWF).^2)*sind(turbine.YawWF)*turbine.Ct; % Eq. 8
    % The originial matlab-FLORIS omits the cos(Y)^2 ??
    wake.zetaInit = 0.5*sind(turbine.YawWF)*turbine.Ct; % Eq. 8
    if model.useWakeAngle
        wake.zetaInit = wake.zetaInit + deg2rad(model.kd);
    end;

    % Calculate initial wake diameter
    if model.adjustInitialWakeDiamToYaw
        wake.wakeDiameterInit = turb_type.rotorDiameter*cosd(turbine.YawWF);
    else
        wake.wakeDiameterInit = turb_type.rotorDiameter;
    end
end