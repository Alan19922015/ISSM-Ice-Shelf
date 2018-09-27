 close all; clear; clc;

%% general script settings (1)
plot_data_sets = false;
md.damage.isdamage=1;
plot_grounding_line = true;
plot_friction_coefficient = true;
stress_balance = 'SSA';  % SSA or HO

%% add to path (2)
addpath('/home/lidaan/Documents/Project/examples/exercise/test/filch/bin/');    % my scripts

% if any(steps==1) 
    
% create output folders
if ~exist('models', 'file')
    mkdir models
end

if ~exist('figures', 'file')
    mkdir figures
end

%% Read input data (3)

% read Bedmap2 data if not loaded
if ~exist('bm2', 'var')
    bm2 = read_bedmap2();
end

% read ShenQiang velocity if not loaded
if ~exist('rignot', 'var')
    rignot = read_daan_velocity();
end

%plot_velocity(rignot, 1);

% visualize input data
if plot_data_sets
    close all
    plot_bed(bm2)
    plot_surface(bm2)
    plot_grounded(bm2)
    plot_velocity(rignot, 1)
end


%% Define model limits (4)
% 1. Press 'add a contour (closed)'
% 2. Click around area of interest (no need to close polygon)
% 3. Press <Enter>
% 4. Close exptool dialog box by pressing 'Quit' button
% 5. Close figure
domain = 'DomainOutline.exp';
if ~exist(domain, 'file')
    if ~plot_data_sets
        plot_velocity(rignot, 1)
    end
    exptool(domain)
end


%% Mesh generation (5)

md=triangle(model,'DomainOutline.exp',2000);

% interpolate velocities onto coarse mesh
vx_obs = InterpFromGridToMesh( ...
    rignot.x, rignot.y, ...
    rignot.vx, ...
    md.mesh.x, md.mesh.y, ...
    0);

vy_obs = InterpFromGridToMesh( ...
    rignot.x, rignot.y, ...
    rignot.vy, ...
    md.mesh.x, md.mesh.y, ...
    0);

vel_obs = sqrt(vx_obs.^2 + vy_obs.^2);

% md.initialization.vx = vx_obs;
% md.initialization.vy = vy_obs;
% md.initialization.vel = vel_obs;


%refine mesh using surface velocities as metric

md=bamg(md,'hmin',1000,'hmax',2000,'gradation',1.7,'field',vel_obs,'err',8);
%[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);
plotmodel(md,'data','mesh');
saveas(gcf, 'figures/model_mesh')

save models/model_mesh_generation md;

% end


%% Apply masks for grounded/floating ice (6)

% if any(steps==2)
    
md = loadmodel('models/model_mesh_generation');

% interpolate onto our mesh vertices
groundedice = double(InterpFromGridToMesh(...
    bm2.x', bm2.y', bm2.grounded, ...
    md.mesh.x, md.mesh.y, 0));

% fill in the md.mask structure
% ice is grounded for mask equal one
md.mask.groundedice_levelset = groundedice;
clear groundedice

% interpolate onto our mesh vertices
ice = double(InterpFromGridToMesh(...
    bm2.x', bm2.y', bm2.grounded, ...
    md.mesh.x, md.mesh.y, 0));
ice(ice > 1.0e-5) = -1;   % make ice shelves count as ice

% ice is present when negative
%md.mask.ice_levelset = -1 * ones(md.mesh.numberofvertices, 1);% all is ice
md.mask.ice_levelset = ice;
clear ice

if plot_grounding_line
    plotmodel(md, ...
        'data', md.mask.groundedice_levelset, ...
        'title', 'grounded/floating', ...
        'data', md.mask.ice_levelset, ...
        'title', 'ice/no-ice');
    saveas(gcf, 'figures/model_grounding_line')
    saveas(gcf, 'figures/model_grounding_line.pdf')
end

% Save model
save models/model_set_mask md;

% end


%% Parameterization (7)
md = loadmodel('models/model_set_mask');
md = parameterize(md, 'model_params.m');

% define stress balance

if strcmp(stress_balance, 'SSA')
    md = setflowequation(md, 'SSA', 'all');
    
    % use SIA for slow ice
    %md = setflowequation(md, 'SSA', md.inversion.vel_obs > 100, ...
    %  'fill', 'SIA');
    
elseif strcmp(stress_balance, 'HO')
    n_layers = 3;
    md = extrude(md, n_layers, 0.9);
    md = setflowequation(md, 'HO', 'all');
    clear n_layers;
end


    
% Save model
save models/model_parameterization md;


%% Find stress balance (control method) (8)

md = loadmodel('models/model_parameterization');

% Control general
md.inversion.iscontrol = 1;
md.inversion.maxsteps = 20;
md.inversion.maxiter = 40;
md.inversion.dxmin = 0.1;
md.inversion.gttol = 1.0e-4;
md.verbose=verbose('solution', true, 'control', true);
%md.inversion.step_threshold=0.99*ones(md.inversion.nsteps,1);
%md.inversion.maxiter_per_step=5*ones(md.inversion.nsteps,1);

% Cost functions
md.inversion.cost_functions = [101 502];  %[101 502]
md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,2);
md.inversion.cost_functions_coefficients(:,1)     = 1; %ones(md.mesh.numberofvertices,1); % 40
md.inversion.cost_functions_coefficients(:,2) = 1e-23; %10^-10*ones(md.mesh.numberofvertices,1);
	
% Controls
md.inversion.control_parameters = {'MaterialsRheologyBbar'};
%%md.inversion.gradient_scaling(1:md.inversion.nsteps)=30;
%md.inversion.control_parameters = {'FrictionCoefficient'};
%md.inversion.min_parameters = 1 * ones(md.mesh.numberofvertices, 1);



% Min/max allowed values of FrictionCoefficient or MaterialsRheologyBbaro 
md.inversion.min_parameters    = cuffey(273)*ones(md.mesh.numberofvertices,1);
md.inversion.max_parameters    = cuffey(200)*ones(md.mesh.numberofvertices,1);
    
% md.inversion.min_parameters = 0 * ones(md.mesh.numberofvertices, 1);
% md.inversion.max_parameters = 200 * ones(md.mesh.numberofvertices, 1);

% Additional parameters
% md.stressbalance.restol = 0.01;
% md.stressbalance.reltol = 0.1;
% md.stressbalance.abstol = NaN;

% Solve
md.toolkits = toolkits;
md.cluster = generic('name', oshostname, 'np', 4);
%md = solve(md, StressbalanceSolutionEnum);
md = solve(md,'Stressbalance');

% Update model friction fields accordingly
% md.friction.coefficient = ...
%     md.results.StressbalanceSolution.FrictionCoefficient;

plotmodel(md,...
		'data',md.results.StressbalanceSolution.MaterialsRheologyBbar,'title','B inversion',...
		'data',md.materials.rheology_B,'title','B observation');
saveas(gcf, 'figures/Rheology_B');

if plot_friction_coefficient
    plotmodel(md, 'data', md.friction.coefficient, ...
        'FontSize#all', 12, ...
        'colormap#all', 'parula')
    saveas(gcf, 'figures/model_friction')
    saveas(gcf, 'figures/model_friction.pdf')
end

% Save model
save models/model_control_drag md;

%% Temperature from SeaRise data (9)

plotmodel(md,'data',md.initialization.temperature,'title','Observed temperature field');
saveas(gcf, 'figures/Observed temperature field');


%% Calculate stress balance and basal drag
% Find sliding exponents
s = averaging(md, 1 ./ md.friction.p, 0);
r = averaging(md, md.friction.q ./ md.friction.p, 0);

% Compute horizontal basal velocity [m/a]
type='basal_drag';
if strcmpi(type, 'basal_drag')
	ub = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2) / ...
        md.constants.yts;
elseif strcmpi(type, 'basal_dragx')
	ub = md.initialization.vx / md.constants.yts;
elseif strcmpi(type, 'basal_dragy')
	ub = md.initialization.vy / md.constants.yts;
end

% Compute basal drag in Pa
basal_drag = (max(md.constants.g * ...
    (md.materials.rho_ice * md.geometry.thickness + ...
    md.materials.rho_water * md.geometry.base), 0)...
    ).^r .* (md.friction.coefficient).^2 .* ub.^s;

% Compute basal shear heat production rate per square meter [J/(a*m^2)]
basal_shear_heating_rate = basal_drag .* ub;

% Compute vertical heat diffusion
vertical_heat_diffusion = 0;%...
    %md.materials.thermalconductivity * dT_dz;

% Find basal melt rate
basal_melt_rate = ...
    (md.basalforcings.geothermalflux + ...
    basal_shear_heating_rate + ...
    vertical_heat_diffusion) / ...
    (md.materials.latentheat * ... % latent heat in J/kg
    md.materials.rho_ice);

clear s r;

%%  observed and modeled velocity field


% plotmodel(md,...
% 		'data',md.initialization.vel,'title','Observed velocity',...
% 		'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity',...
% 		'data',md.geometry.base,'title','Bed elevation',...
% 		'data',md.results.StressbalanceSolution.FrictionCoefficient,'title','Friction Coefficient',...
% 		'colorbar#all','on','colorbartitle#1-2','(m/yr)',...
% 		'caxis#1-2',([1.5,2000]),...
% 		'colorbartitle#3','(m)', 'log#1-2',10);
    
plotmodel(md,...
		'data',md.initialization.vel,'title','Observed velocity',...
		'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity',...
		'data',md.geometry.base,'title','Bed elevation',...
        'colorbar#all','on','colorbartitle#1-2','(m/yr)',...
		'caxis#1-2',([1.5,2000]),...
        'colorbartitle#3','(m)')
% 		'colorbartitle#3','(m)', 'log#1-2',10);
saveas(gcf, 'figures/Observation-Modeled velocity field');

plotmodel(md,'data',md.initialization.vel,'title','Observed velocity',...
             'data',md.geometry.thickness,'title','Ice Thickness');
saveas(gcf,'figures/Ice Thickness');
%plotmodel(md,'data',md.results.StressbalanceSolution.Vel-md.initialization.vel,'title','compared velocity')

%% Calculate stress feild

%md=mechanicalproperties(md,md.initialization.vx,md.initialization.vy);
    
md=mechanicalproperties(md,md.inversion.vx_obs,md.inversion.vy_obs);

plotmodel(md,'data',md.results.deviatoricstress.effectivevalue,'title','Stress field');
saveas(gcf, 'figures/Stress field')
plotmodel(md,'data',md.results.deviatoricstress.xx,'title','Stress-x field');
saveas(gcf, 'figures/Stress-x field')
plotmodel(md,'data',md.results.deviatoricstress.yy,'title','Stress-y field');
saveas(gcf, 'figures/Stress-y field')
plotmodel(md,'data',md.results.deviatoricstress.xy,'title','Stress-xy field');
saveas(gcf, 'figures/Stress-xy field')



%% Calculate damage

 md.damage.D=damagefrominversion(md);
 plotmodel(md,'data',md.damage.D,'title','damage field');
 
 saveas(gcf, 'figures/Damage field');
 
 backstress=backstressfrominversion(md,'smoothing',2,'coordsys','longitudinal','tempmask',true);
    
 %backstress=calcbackstress(md,'smoothing',2,'coordsys','longitudinal','tempmask',true);
    
 plotmodel(md,'data',backstress,'title','Backstress field');
 saveas(gcf, 'figures/Backstress field');


%% Calculate hydrology
effective_pressure = zeros(length(md.geometry.thickness), 1);

water_pressure = ...
    md.materials.rho_ice * md.constants.g * md.geometry.thickness - ...
    effective_pressure;

hydro_potential = ...
    md.materials.rho_water * md.constants.g * md.geometry.base + ...
    water_pressure;

%[hydro_potential_gradient_x, hydro_potential_gradient_y] = ...
%    gradient(hydro_potential);
%hydro_potential_gradient_norm = ...
%    sqrt(hydro_potential_gradient_x.^2 + hydro_potential_gradient_y.^2);

%% Post visualization
% plotmodel(md, 'nlines', 3, 'ncols', 3, ...
%     'unit#all', 'km', 'axis#all', 'equal', ...
%     'xlim#all', [min(md.mesh.x) max(md.mesh.x)] / 10^3, ...
%     'ylim#all', [min(md.mesh.y) max(md.mesh.y)] / 10^3, ...
%     'FontSize#all', 12, ...
%     'colormap#all', 'parula', ...
%     'data', md.initialization.vel, ...
%     'title', 'Observed velocity', ...
%     'data', md.results.StressbalanceSolution.Vel, ...
%     'title', 'Modeled Velocity', ...
%     'colorbar#all', 'on', 'colorbartitle#1-2', '[m/yr]', ...
%     'caxis#1-2', ([1.5, 4000]), ...
%     'log#1-2', 10, ...
%     'data', md.geometry.base, ...
%     'title', 'Bed elevation', ...
%     'colorbartitle#3', '[m]', ...
% %     'data', md.results.StressbalanceSolution.FrictionCoefficient, ...
% %     'title', 'Friction Coefficient', ...
%     'title', 'SMB', 'data', md.smb.mass_balance, ...
%     'title', 'Geothermal heatflux', 'data', md.basalforcings.geothermalflux, ...
%     'title', 'Hydropotential', 'data', hydro_potential, ...
%     'title', 'Driving stress [kPa]', 'data', 'driving_stress', ...
%     'title', 'Basal drag [kPa]', 'data', 'basal_drag', ...
%     'title', 'Basal melt rate [m/a]', 'data', 'basal_drag' ...
%     );
% 
% saveas(gcf, 'figures/model_combined')
% saveas(gcf, 'figures/model_combined.pdf')

% plotmodel(md, 'nlines', 1, 'ncols', 2, ...
%     'unit#all', 'km', 'axis#all', 'equal', ...
%     'xlim#all', [min(md.mesh.x) max(md.mesh.x)] / 10^3, ...
%     'ylim#all', [min(md.mesh.y) max(md.mesh.y)] / 10^3, ...
%     'FontSize#all', 12, ...
%     'colormap#all', 'parula', ...
%     'title', 'Basal drag [kPa]', 'data', 'basal_drag' ...
%     );


plotmodel(md, 'nlines', 1, 'ncols', 1, ...
    'unit#all', 'km', 'axis#all', 'equal', ...
    'xlim#all', [min(md.mesh.x) max(md.mesh.x)] / 10^3, ...
    'ylim#all', [min(md.mesh.y) max(md.mesh.y)] / 10^3, ...
    'FontSize#all', 12, ...
    'log', 10, ...
    'colormap#all', 'parula', ...
    'title', 'Basal melt rate [m/a]', 'data', basal_melt_rate ...
    );

%% Cleanup time
clear plot_data_sets ...
    plot_meshes ...
    plot_geometry ...
    plot_grounding_line ...
    plot_friction_coefficient ...
    stress_balance


save models/model_filch md;



