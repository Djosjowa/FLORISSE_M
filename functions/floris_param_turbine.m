function [ turb ] = floris_turbine( name )
switch lower(name)
    case 'nrel5mw'
        turb.rotorDiameter           = 126.4;
        turb.rotorArea               = pi*turb.rotorDiameter(1)*turb.rotorDiameter(1)/4.0;
        turb.generator_efficiency    = 0.944;
        turb.hub_height              = 90.0;
        
    otherwise
        error(['Turbine parameters with name "' turb.name '" not defined']);
end;
end