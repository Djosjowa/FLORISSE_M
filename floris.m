clear all; clc;
timer.script = tic;
addpath('functions');

%% Script settings
plotLayout    = true;     % plot farm layout w.r.t. inertial and wind frame
plot2DFlowfield = true ; % 2DflowFieldvisualisation in wind-aligned frame
plot3DFlowfield = false ;  % 3DflowFieldvisualisation in wind-aligned frame

if (plot2DFlowfield || plot3DFlowfield)
   % resz is not used when only 2Dflowfield is plotted
   flowField.resx   = 5;     % resolution in x-axis in meters (windframe)
   flowField.resy   = 5;     % resolution in y-axis in meters (windframe)
   flowField.resz   = 10;     % resolution in z-axis in meters (windframe)
   flowField.fixYaw  = true;  % Account for yaw in near-turbine region in 2Dplot
   % TODO: implement fixyaw for 3d plot
end

%% Simulation setup
model = floris_param_model('default');          % Import model settings
turbType  = floris_param_turbine('nrel5mw');    % Import turbine settings

% Wind turbine locations in inertial frame 
LocIF =   [300,    100.0,  turbType.hub_height
           300,    300.0,  turbType.hub_height
           300,    500.0,  turbType.hub_height
           1000,   100.0,  turbType.hub_height
           1000,   300.0,  turbType.hub_height
           1000,   500.0,  turbType.hub_height
           1600,   100.0,  turbType.hub_height
           1600,   300.0,  turbType.hub_height
           1600,   500.0,  turbType.hub_height];

% Turbine operation settings in wind frame
% Yaw misalignment with flow (counterclockwise, wind frame)
% Axial induction control setting (used only if model.axialIndProvided == true)
turbines = struct(  'Tilt',num2cell([0 0 0 15 15 15 0 0 0].'), ...
                    'YawWF',num2cell([-27 10 -30 -30 -20 -15 0 10 0].'), ...
                    'axialInd',1/3,'windSpeed',[],'Cp',[],'Ct',[], ...
                    'power',[],'downstream',[]);
% TODO: implement effects of turbine tilt
                
wakes = struct( 'Ke',num2cell(zeros(1,length(turbines))),'mU',{[]}, ...
                'zetaInit',[],'wakeDiameterInit',[],'centerLine',[], ...
                'diameters',[],'OverlapAreaRel',[]);

% Atmospheric settings
site.uInfIf   = 12;       % x-direction flow speed inertial frame (m/s)
site.vInfIf   = 4;        % y-direction flow speed inertial frame (m/s)
site.rho      = 1.1716;   % Atmospheric air density (kg/m3)

%% Internal code of FLORIS
% Determine wind farm layout in wind-aligned frame. Note that the
% turbines are renumbered in the order of appearance w.r.t wind direction
[site,turbines,wtRows] = floris_frame(site,turbines,LocIF);
% The first row of turbines has the freestream as inflow windspeed
[turbines(wtRows{1}).windSpeed] = deal(site.uInfWf); 

% Setup flowField visualisation grid if neccesary
if (plot2DFlowfield || plot3DFlowfield)
    % flowField.dims holds the X, Y and Z windframe dimensions in which
    % the turbines exist
    flowField.dims = max([turbines.LocWF],[],2);
    
    % The X, Y and Z variables form a 3D or 2D mesh
    if plot3DFlowfield
        [flowField.X,flowField.Y,flowField.Z] = meshgrid(...
        -200 : flowField.resx : flowField.dims(1)+500,...
        -200 : flowField.resy : flowField.dims(2)+200,...
        0 : flowField.resz : 200);
    else
        [flowField.X,flowField.Y,flowField.Z] = meshgrid(...
        -200 : flowField.resx : flowField.dims(1)+1000,...
        -200 : flowField.resy : flowField.dims(2)+200,...
        turbType.hub_height);
    end
    
    % initialize the flowfield as freestream in the U direction
    flowField.U  = site.uInfWf*ones(size(flowField.X));
    flowField.V  = zeros(size(flowField.X));
    flowField.W  = zeros(size(flowField.X));
end;

% Start the core model. Without any visualization this is all that runs, It
% computes the power produced at all turbines given the flow and
% turbine settings
timer.core = tic;
for turbirow = 1:length(wtRows) % for first to last row of turbines
    for turbNum = wtRows{turbirow} % for each turbine in this row
        
        % Determine Cp, Ct, axialInduction and power for a turbine
        turbines(turbNum) = floris_cpctpower(model,site.rho,turbType,turbines(turbNum));
        
        % calculate ke, mU, and initial wake deflection & diameter
        wakes(turbNum) = floris_initwake( model,turbines(turbNum),wakes(turbNum),turbType );
        
        % Compute  the X locations of  the downstream turbines rows
        wakes(turbNum).centerLine(1,:) = arrayfun(@(x) x.LocWF(1), turbines(cellfun(@(x) x(1),wtRows(turbirow+1:end)))).';
        wakes(turbNum).centerLine(1,end+1) = turbines(end).LocWF(1)+300;

        % Compute the wake centerLines and diameters at those X locations
        wakes(turbNum) = floris_wakeCenterLine_and_diameter(...
             turbType.rotorDiameter, model, turbines(turbNum), wakes(turbNum));

        % Calculate overlap of this turbine on downstream turbines
        wakes(turbNum) = floris_overlap( (turbirow+1):length(wtRows),wtRows,wakes(turbNum),turbines,turbType );
    end
    
    % If this is not the last turbine row compute the windspeeds at the next row
    if turbirow < length(wtRows)
        % Pass all the upstream turbines and wakes including the next
        % downstream row to the function: wt_rows{1:turbirow+1}
        % Return only the downstream turbine row: wt_rows{turbirow+1}.
        turbines(wtRows{turbirow+1}) = floris_compute_windspeed(...
            turbines([wtRows{1:turbirow+1}]),wakes([wtRows{1:turbirow}]),site,turbType,wtRows,turbirow);
    end;
end;
disp(['TIMER: core operations: ' num2str(toc(timer.core)) ' s.']);

%% Plot the layout and flowfield visualization
% Plot a map with the turbine layout and wake centerLines
if plotLayout
    figure;
    plot_layout( wtRows,site,turbType,turbines,wakes );
end

if (plot2DFlowfield || plot3DFlowfield)
    % Compute the flowfield velocity at every voxel(3D) or pixel(2D)
    [wakes,flowField] = floris_compute_flowfield(site,model,turbType,flowField,turbines,wakes);
end

% Plot the flowfield as a cutthourgh at hubHeigth
if plot2DFlowfield
    figure;
    plot_2d_field( flowField,site,turbines,turbType )
end;

% Plot the 3D flowfield as
if plot3DFlowfield
    figure;

    q = quiver3(flowField.X, flowField.Y, flowField.Z, flowField.U, flowField.V, flowField.W, ...
        1.8,'linewidth',2.5,'ShowArrowHead','off');
    quiverMagCol(q,gca);
    axis equal;
    set(gca,'view',[-55 35]);
    xlabel('x-direction (m)');
    ylabel('y-direction (m)');
    colorbar;
    caxis([floor(min(flowField.U(:))) ceil(site.uInfWf)])

%     [X, Y, Z] = meshgrid(...
%     -200 : flowField.resx : flowField.dims(1)+1000,...
%     -200 : flowField.resy : flowField.dims(2)+200,...
%     0 : flowField.resz : 200);
%     q=coneplot(flowField.X, flowField.Y, flowField.Z, flowField.U, flowField.V, flowField.W,X,Y,Z,.7,flowField.U);
%     set(q,'EdgeColor','none');
%     alpha(q, .2)
%     axis equal;

end

disp(['TIMER: script: ' num2str(toc(timer.script)) ' s.']);