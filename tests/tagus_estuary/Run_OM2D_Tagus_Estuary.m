% This script runs OceanMesh2D for the Tagus Estuary test case of OMesh2KRAKEN
%
% (c) 2022 Nuno Monteiro and Tiago Oliveira, University of Aveiro
%
%

clear all, close all, clc

% Add paths to call OceanMesh2D routines, load data, etc.
addpath(genpath('..\..\OceanMesh2D'));
addpath(genpath('..\..\m_map'));
addpath(genpath('..\..\data'));
addpath(genpath('..\..\utils'));
addpath(genpath('..\..\atWin10_2020_11_4\Matlab\ReadWrite'));  % to use the write_env.m function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 1 - Input options %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 1.1 Inputs for OceanMesh2D %%%%%%%%%

d_limits = [-9.41 -9.24 38.6 38.72; -9.36 -9.27 38.64 38.7]; % boundary box (or mesh domains) limits
N = size(d_limits,1); % number of domains
min_el_all = [250 100]; % minimum element size (max resolution)
max_el_all = [800 400]; % maximum element size (min resolution)
max_el_ns_all = [400 250]; % maximum element size nearshore            
R_all = [2 2]; % nr of elements to solve the feature size of coastlines
grade_all = [0.32 0.27]; % mesh grade
dems = ["EMODNET_DTM_SouthLisbon.nc"; "EMODNET_DTM_SouthLisbon.nc"]; % bathymetry data filenames
coastlines = ["CNTR_RG_01M_2020_4326"; "CNTR_RG_01M_2020_4326"]; % coastline filenames
alpha_all = [0 25]; % alpha parameter to solve the slope mesh size function
coord_source = [-9.317, 38.665]; % sound source coords in lat/lon

%%%%% 1.2 Inputs for KRAKEN %%%%%%%%%%

% Struct array with input parameters to be written in all ENV files
params_env.freq = 1000; % frequency (Hz)
params_env.NMedia = 2; % number of media
params_env.top_opt = "CVW ."; % top option
params_env.bot_opt = "V"; % bottom option
params_env.sigma_int = [0 0]; % RMS roughness at surface and bottom interfaces
params_env.N_layers = [10e3 10e3]; % number of mesh points in the water and bottom layers
params_env.Rdepth = 150; % max depth of receiver (m)
params_env.sigma_bot = 0; % RMS roughness at the bottom
params_env.Sdepth = 10; % source depth (m)
params_env.C_bot = [1700 1700]; % bottom layer sound speeds
params_env.C_lim = [1400 15000]; % min and max phase sound speed (m/s)
params_env.C_high = 15000; % max phase sound speed (m/s)
params_env.deltassp = 1; % increment of sound speed profile's depth (m)
params_env.sndspd_bot =  1700; % soundspeed at bottom layer

%%%%% 1.3 Inputs for FIELD3D %%%%%%%%%%

% Struct array with input parameters to be written in the FLP file
params_flp.title = "Tagus Estuary" % title of run
params_flp.calc = "STD"; % type of calculation
params_flp.Nm = 999; % number of modes 
coord_or_km = [0 0]; % coords of the origin (km)
params_flp.coord_or = coord_or_km;                      
params_flp.Sxy_lim = [coord_or_km(1) ...
    coord_or_km(1); coord_or_km(2) coord_or_km(2)]; % 2D range of source(s) (m)
params_flp.Ns = [1 1]; % number of source coords in x and y
params_flp.Sz_lim = [10 10]; % source depth range (m)
params_flp.NSz = 1; % number of source depths
params_flp.Rz_lim = [0 150]; % receivers' depth range (m)
params_flp.NRz = 151; % number of receiver depths
params_flp.Rr_lim = [0 10]; % receivers' radial range (km)
dr = 30; % receiver range increment (km)
rR=0:dr:10e3; % receiver depths (m)
NRr = length(rR); % number of receiver depths
params_flp.NRr = NRr; % number of receiver distances
params_flp.Rtheta_lim = [-180 180]; % min receiver radial (deg)
params_flp.NRtheta = 361; % number of receiver radials
params_flp.GBtheta_lim = [-180 180]; % angle range of Gaussian beams (deg)
params_flp.NGBtheta = 361; % number of Gaussian beam angles
params_flp.GBstep = 30; % radial step size of Gaussian beams (m)
params_flp.GBsteps = 334; % number of Gaussian beam radial steps
params_flp.epsilon_mult = 0.3; % epsilon multiplier for Gaussian beam initial conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 2 - Input processing %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Temperature-Salinity data extraction and arrangement
fname_TSdata = 'cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m_1647891936780.nc';
lon_TS=ncread(fname_TSdata,'longitude');
lat_TS=ncread(fname_TSdata,'latitude');
t_TS=ncread(fname_TSdata,'time')./24+datenum(1950,1,1)
idx_t = 1; % Extracting the first instant (daily mean data for the first day only)
t_TS=t_TS(idx_t); 
T=ncread(fname_TSdata,'thetao');
T=squeeze(T(:,:,:,idx_t));
S=ncread(fname_TSdata,'so');
S=squeeze(S(:,:,:,idx_t));
z_TS=ncread(fname_TSdata,'depth');

% Concatenating the original first depth level values of T and S to the zero level (if first depth is not zero meters)
if z_TS(1) ~= 0      
    z_TS=[0; z_TS];
    S=cat(3,S(:,:,1),S); 
    T=cat(3,T(:,:,1),T);
end

% Create a struct with the temperature and salinity data, to be read in OM2D's "write" function
TS_data.S = S;
TS_data.T = T;
TS_data.lon = lon_TS;
TS_data.lat = lat_TS;
TS_data.z = z_TS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 3 - Run OceanMesh2D %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 3.1 geodata and edgefx  %%%%%%%%%

for i=1:N
    min_el = min_el_all(i);
    max_el = max_el_all(i);
    max_el_ns = max_el_ns_all(i);
    R = R_all(i); grade=grade_all(i);
    alpha = alpha_all(i);
    coastline = convertStringsToChars(coastlines(i));
    dem = convertStringsToChars(dems(i));
    bbox = [d_limits(i,1:2);d_limits(i,3:4)];
    gdat{i} = geodata('shp',coastline,'dem',dem,'bbox',bbox,'h0',min_el);
    fh{i} = edgefx('geodata',gdat{i},'fs',R,'max_el_ns',max_el_ns,'max_el',max_el,'g',grade,'slp',alpha);
end

%%%%% 3.2 meshgen  %%%%%%%%%%%%%%%%%%%%

mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'proj','mercator','pfix',coord_source,'plot_on',0);
mshopts = mshopts.build;
m = mshopts.grd;

%%%%% 3.3 msh  %%%%%%%%%%%%%%%%%%%%%%%%

m = interp(m,gdat);
m = make_bc(m,'auto',gdat{1}); 
tri = m.t;
lon = m.p(:,1);
lat = m.p(:,2);
z = m.b;

% Define Lat/Lon coords of the new referential's origin in km for KRAKEN
[lo,la,idx] = find_nearest_point(lon, lat, coord_source(1), coord_source(2)); % Source coordinates and index in the mesh
coord_or = [lo la]; % Coords of source (in degrees)

% Other corrections to the mesh bathymetry
z(isnan(z)) = 0;  % Convert nan nodes into zeros
z(z<0) = 0;       % Convert negative numbers (land) into zeros
m.b = z;          % Save correction to unstructured mesh struct variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 4 - Write mesh data %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_env = 1;   % Flag to write/store ENV files (0 - no; 1 - yes)
flag_flp = 1;   % Flag to write/store FLP files (0 - no; 1 - yes)

fname = 'tagus_estuary';
inp_env.fname = fname;        % prefix of ENV filenames
inp_flp.fname = fname;        % FLP filename

% Struct array with ENV inputs (as struct fields)
inp_env.flag = flag_env;        % flag for ENV file creation
inp_env.params = params_env;    
inp_env.TS_data = TS_data; % not required field: "TS_data", with a struct with the temperature and salinity data for 
                           % sound speed calculation and storage in vertical profiles at ENV files.
                           % (alternatively, the user can create a 'SS_constant' field and attribute a number as a constant 
                           % value of sound speed. Or even if the user doesn't use 'TS_data' or 'SS_constant' fields, the 
                           % writekraken3d.m program attributes 1500 m/s value in all points.

% Struct array with all required FLP inputs (as struct fields)
inp_flp.coord_or = coord_or; % origin coordinates in Lat/Lon
inp_flp.flag = flag_flp; % flag for FLP file creation (to run KRAKEN+FIELD3D): 0 - disable FLP file creation, 1 - enable
inp_flp.params = params_flp; % parameters of the FLP file

% Write data to ENV and FLP files for KRAKEN3D
write(m,fname,'kraken3d','nob',1,'inp_env',inp_env,'inp_flp',inp_flp);

% Create a directory for input/output data of KRAKEN, if it doesn't exist
if ~exist('data_kraken','dir')
    mkdir('data_kraken')
else
    delete('data_kraken\*') % Remove everything from the KRAKEN data directory
end

movefile('*.env','data_kraken\') % Move all envs to the input/ouput kraken directory

if ~exist('plots\','dir')
    mkdir('plots\') % Create a directory for plotting if it doesnt exist
end

% Create a Struct with mesh data (bathymetry, nodes' position, etc.)
mesh_data.z=z;
mesh_data.lon=lon;
mesh_data.lat=lat;
mesh_data.tri=tri;
mesh_data.pfix=coord_or;
save('mesh_data.mat','-struct','mesh_data');    % Save the MATLAB struct with mesh info (to be used in Plotting_Tagus_Estuary.m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 5 - Plotting mesh  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

view_tri = 1; % Flag to view mesh triangles (under 'plot_mesh_bathy')
view_cb = 1; % Flag to view colorbar (under 'plot_mesh_bathy')
caxis_lims = [0 50]; % depth limits, in meters

% Plot the entire mesh (Tagus Estuary)
figure
lonlat_lims = d_limits(1,:);
plot_mesh_bathy(tri, z, lon, lat, lonlat_lims, caxis_lims, view_tri, view_cb)
plot3(lo,la,1e6,'ro','MarkerFaceColor','r')
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'g-','LineWidth',2.5)
lonlat_lims = d_limits(2,:);
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
saveas(gcf,['plots/' fname '_mesh_box1.png'])

% Plot mesh in boundary box 2 (Tagus Estuary's inlet)
figure
lonlat_lims = d_limits(2,:);
view_tri = 1;
plot_mesh_bathy(tri, z, lon, lat, lonlat_lims, caxis_lims, view_tri, view_cb)
plot3(lo,la,1e6,'ro','MarkerFaceColor','r')
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
saveas(gcf,['plots/' fname_env '_mesh_box2.png'])
