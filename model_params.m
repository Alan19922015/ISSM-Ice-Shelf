% Parameters to change/Try
friction_coefficient = 10; % default [10]
Temp_change          =  0;  % default [0 K]

%Name and Coordinate system
md.miscellaneous.name = 'siple';
md.mesh.epsg = 3031;

%NetCdf Loading
disp('   Loading Bedmap2 data from NetCDF');
bm2 = read_bedmap2();

disp('   Loading SeaRISE data from NetCDF');
ncdata = '../data/Antarctica_5km_withshelves_v0.75.nc';
sr.x     = ncread(ncdata, 'x1');  % These coordinates differ from bedmap/rignot
sr.y     = ncread(ncdata, 'y1');  % These coordinates differ from bedmap/rignot
%sr.usrf  = ncread(ncdata, 'usrf')';
%sr.topg  = ncread(ncdata, 'topg')';
sr.temp  = ncread(ncdata, 'presartm')';
sr.smb   = ncread(ncdata, 'presprcp')';
sr.gflux = ncread(ncdata, 'bheatflx_fox')'; % fox or shapiro
%sr.gflux = ncread(ncdata, 'bheatflx_shapiro')'; % fox or shapiro


disp('   Loading surface velocities');
rignot = read_daan_velocity();

% Geometry
disp('   Interpolating surface and ice base');

md.geometry.base    = InterpFromGridToMesh( ...
    bm2.x', bm2.y', bm2.bed, ...
    md.mesh.x, md.mesh.y, 0);
md.geometry.surface = InterpFromGridToMesh( ...
    bm2.x', bm2.y', bm2.surface, ...
    md.mesh.x, md.mesh.y, 0);
%clear usrf topg;

disp('   Constructing thickness');
md.geometry.thickness = md.geometry.surface - md.geometry.base;

% ensure hydrostatic equilibrium on ice shelf: 
di = md.materials.rho_ice / md.materials.rho_water;

%Get the node numbers of floating nodes
pos = find(md.mask.groundedice_levelset < 0); 

%apply a flotation criterion on the precedingly defined nodes and
%redefine base and thickness accordingly
md.geometry.thickness(pos) = 1 / (1 - di) * md.geometry.surface(pos);% rho_water*(thickness-surface)
md.geometry.base(pos) = md.geometry.surface(pos) ...                 % =rho_ice*thickness
    - md.geometry.thickness(pos);
md.geometry.hydrostatic_ratio = ones(md.mesh.numberofvertices, 1);

%Set min thickness to 1 meter
pos0=find(md.geometry.thickness <= 0);
md.geometry.thickness(pos0) = 1;
md.geometry.surface = md.geometry.thickness + md.geometry.base;

%Initialization parameters
disp('   Interpolating temperatures');
md.initialization.temperature = InterpFromGridToMesh( ...
    linspace(-1*max(sr.y), abs(min(sr.y)), length(sr.y))', ...
    sr.x, ...
    rot90(sr.temp, 3), ...
    md.mesh.x, md.mesh.y, 0) + 273.15 + Temp_change;
% md.initialization.temperature = InterpFromGridToMesh( ...
%     sr.x, sr.y, ...
%     rot90(sr.temp, 3), ...
%     md.mesh.x, md.mesh.y, 0) + 273.15 + Temp_change;
clear temp;

disp('   Set observed velocities')
vx_obs=InterpFromGridToMesh( ...
    rignot.x, rignot.y, rignot.vx, ...
    md.mesh.x, md.mesh.y, 0);
vy_obs=InterpFromGridToMesh( ...
    rignot.x, rignot.y, rignot.vy, ...
    md.mesh.x, md.mesh.y, 0);
clear velx vely;

vel_obs = sqrt(vx_obs.^2 + vy_obs.^2);
md.initialization.vx = vx_obs;
md.initialization.vy = vy_obs;
md.initialization.vz = zeros(md.mesh.numberofvertices, 1);
md.initialization.vel = vel_obs;

disp('   Set Pressure');
md.initialization.pressure = ...
    md.materials.rho_ice * md.constants.g * md.geometry.thickness;

disp('   Construct ice rheological properties');
md.materials.rheology_n = 3 * ones(md.mesh.numberofelements, 1);
md.materials.rheology_B = paterson(md.initialization.temperature);
md.damage.D=zeros(md.mesh.numberofvertices,1);

%Forcings
disp('   Interpolating surface mass balance'); % is the SeaRISE SMB projection different?
md.smb.mass_balance = InterpFromGridToMesh( ...
    linspace(-1*max(sr.y), abs(min(sr.y)), length(sr.y))', ...
    sr.x, ...
    rot90(sr.smb, 3), ...
    md.mesh.x, md.mesh.y, 0);
md.smb.mass_balance = ...
    md.smb.mass_balance * md.materials.rho_water / md.materials.rho_ice;
clear smb;

disp('   Set geothermal heat flux');
md.basalforcings.geothermalflux = InterpFromGridToMesh( ...
    linspace(-1*max(sr.y), abs(min(sr.y)), length(sr.y))', ...
    sr.x, ...
    rot90(sr.gflux, 3), ...
    md.mesh.x, md.mesh.y, 0);
clear gflux;

%Friction and inversion set up
disp('   Construct basal friction parameters');
md.friction.coefficient = ...
    friction_coefficient * ones(md.mesh.numberofvertices, 1);
md.friction.p = ones(md.mesh.numberofelements, 1);
md.friction.q = ones(md.mesh.numberofelements, 1);

%no friction applied on floating ice
pos = find(md.mask.groundedice_levelset < 0);
md.friction.coefficient(pos) = 0;

md.inversion = m1qn3inversion();
md.inversion.vx_obs = vx_obs;
md.inversion.vy_obs = vy_obs;
md.inversion.vel_obs = vel_obs;

% exclude areas without data from inversion
md.inversion.cost_functions_coefficients(...
    md.inversion.vel_obs < 1, 1:2) = 0;  % choose to accord with the judgement as lines

disp('   Set boundary conditions');
md = SetMarineIceSheetBC(md);
md.basalforcings.floatingice_melting_rate = ...
    zeros(md.mesh.numberofvertices, 1);
md.basalforcings.groundedice_melting_rate = ...
    zeros(md.mesh.numberofvertices, 1);

%impose observed temperature on surface
md.thermal.spctemperature = md.initialization.temperature; 
md.masstransport.spcthickness = NaN * ones(md.mesh.numberofvertices, 1);
