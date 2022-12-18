% This script runs OceanMesh2D for the Azores test case of OMesh2KRAKEN
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
addpath(genpath('..\..\atWin10_2020_11_4\Matlab\ReadWrite'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 1 - Input options %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 1.1 Inputs for OceanMesh2D %%%%%%%%%

d_limits = [-32 -24.5 36 40.5; -31.5 -30.7 39 40;
            -29 -26.7 38 39.4; -26.2 -24.8 36.7 38.1]; % boundary box (or mesh domains) limits
N = size(d_limits,1); % number of domains
min_el_all = [10e3 3e3 3e3 3e3]; % minimum element size (max resolution)
max_el_all = [22e3 15e3 15e3 15e3]; % maximum element size (min resolution)
max_el_ns_all = [12e3 7e3 7e3 7e3]; % maximum element size nearshore            
R_all = [2 3 3 3]; % nr of elements to solve the feature size of coastlines
grade_all = [0.3 0.2 0.2 0.2]; % mesh grade
dems = ["gebco_azores_2020.nc"; "EMODNET_DTM_Azores_1.nc";
        "EMODNET_DTM_Azores_2.nc"; "EMODNET_DTM_Azores_3.nc"]; % digital elevation model filenames, defining the input bathymetry
coastlines = ["GSHHS_h_L1"; "GSHHS_h_L1"; "GSHHS_h_L1"; "GSHHS_h_L1"]; % coastline filenames
alpha_all = [0 25 25 25]; % alpha parameter to solve the slope mesh size function                                                        % filename of output mesh file (FLP) and prefix name of ENV files
coord_source = [-29 38]; % sound source coords in lat/lon


%%%%% 1.2 Inputs for KRAKEN %%%%%%%%%%

% Struct array with input parameters to be written in all ENV files
params_env.freq = 25; % frequency (Hz)
params_env.NMedia = 2; % number of media
params_env.top_opt = "CVW ."; % top option
params_env.bot_opt = "V"; % bottom option
params_env.sigma_int = [0 0]; % RMS roughness at surface and bottom interfaces
params_env.N_layers = [10e3 10e3]; % number of mesh points in the water and bottom layers
params_env.Rdepth = 6e3; % max depth of receiver (m)
params_env.sigma_bot = 0; % RMS roughness at the bottom
params_env.Sdepth = 10; % source depth (m)
params_env.C_bot = [1700 1700]; % bottom layer sound speeds
params_env.C_lim = [1400 15000]; % min and max phase sound speed (m/s)
params_env.C_high = 15000; % max phase sound speed (m/s)
params_env.deltassp = 10; % increment of sound speed profile's depth (m)
params_env.sndspd_bot =  1700; % soundspeed at bottom layer

%%%%% 1.3 Inputs for FIELD3D %%%%%%%%%%

% Struct array with input parameters to be written in the FLP file
title_flp = "Azores archipelago"            % title of run 
calc_flp = "STD";                           % type of calculation
Nm = 999;                                   % number of modes
Nsx = 1;                                    % number of source coords in x
Nsy = 1;                                    % number of source coords in y
coord_or_km = [0 0];                       % coords of the origin (km)
NSz = 1;                                    % number of source depths   
Sz = 10;                                    % source depths (m)
Rzmin = 0;                                  % min receiver depth (m)
Rzmax = Rdepth;                             % max receiver depth (m)
dz = deltassp;                              % receiver depth discretization
zR=Rzmin:dz:Rzmax;                          % receiver depths
NRz = length(zR);                           % number of receiver depths
Rrmin = 0;                                  % min receiver range (km)
Rrmax = 350;                                % max receiver range (km)
dr = 1;                                     % receiver range discretization (km)
rR=Rrmin:dr:Rrmax;                          % receiver depths (km)
NRr = length(rR);                           % number of receiver depths
Rthetamin = -180;                           % min receiver radial (deg)
Rthetamax = 180;                            % max receiver radial (deg)
dtheta = 1;                                 % receiver radial discretization (deg)
thetaR=Rthetamin:dtheta:Rthetamax;          % receiver radials (deg)
NRtheta = length(thetaR);                   % number of receiver radials
GBthetamin = -180;                          % min Gaussian beam radial (deg)
GBthetamax = 180;                           % max Gaussian beam radial (deg)
dthetaGB = 1;                               % Gaussian beam radial discretization (deg)
thetaGB=GBthetamin:dthetaGB:GBthetamax;     % Gaussian beam radials (deg)
NGBtheta = length(thetaGB);                 % number of Gaussian beam radials
GBstep = 1;                                 % step size of Gaussian beams (m)
GBsteps = 1000;                             % number of Gaussian beam steps
epsilon_mult = 0.3;                         % epsilon multiplier for Gaussian beam initial conditions
fname_flp='tagus_estuary'; 

params_flp.title = "Azores archipelago" % title of run
params_flp.calc = "STD"; % type of calculation
params_flp.Nm = 999; % number of modes 
coord_or_km = [0 0]; % coords of the origin (km)
params_flp.coord_or = coord_or_km;                      
params_flp.Sxy_lim = [coord_or_km(1) ...
    coord_or_km(1); coord_or_km(2) coord_or_km(2)]; % 2D range of source(s) (m)
params_flp.Ns = [1 1]; % number of source coords in x and y
params_flp.Sz_lim = [10 10]; % source depth range (m)
params_flp.NSz = 1; % number of source depths
dz = 10;                           
zR = Rzmin:dz:Rzmax;
NRz = length(zR);
params_flp.Rz_lim = [zR(1) zR(end)]; % receivers' depth range (m)
params_flp.NRz = 151; % number of receiver depths
dr = 1; % receiver range increment (km)
rR=0:dr:350;
NRr = length(rR);
params_flp.Rr_lim = [rR(1) rR(end)]; % receivers' radial range (km)
params_flp.NRr = NRr; % number of receiver radial distances
params_flp.Rtheta_lim = [-180 180]; % min receiver angle (deg)
params_flp.NRtheta = 361; % number of receiver angles
params_flp.GBtheta_lim = [-180 180]; % angle range of Gaussian beams (deg)
params_flp.NGBtheta = 361; % number of Gaussian beam angles
params_flp.GBstep = 1e3; % radial step size of Gaussian beams (m)
params_flp.GBsteps = 350; % number of Gaussian beam radial steps
params_flp.epsilon_mult = 0.3; % epsilon multiplier for Gaussian beam initial conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 2 - Input processing %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Temperature-Salinity data extraction and arrangement
fname_TSdata = 'global-analysis-forecast-phy-001-024-monthly_1656015476392.nc';
lon_TS=ncread(fname_TSdata,'longitude');
lat_TS=ncread(fname_TSdata,'latitude');
t_TS=ncread(fname_TSdata,'time')./24+datenum(1950,1,1)
idx_t = 1; % month index (data will be only extracted for January)
t_TS=t_TS(idx_t); % extract the first instant
T=ncread(fname_TSdata,'thetao');
T=squeeze(T(:,:,:,idx_t));
S=ncread(fname_TSdata,'so');
S=squeeze(S(:,:,:,idx_t));
z_TS=ncread(fname_TSdata,'depth');

% Concatenating the previous first depth level values of T and S to the zero level:
if z_TS(1)~=0
    z_TS=[0; z_TS];
    S=cat(3,S(:,:,month),S); 
    T=cat(3,T(:,:,month),T);
end

% Create a struct out of the T-S input data to be read in OM2D's "write" function
TS_data.S=S;
TS_data.T=T;
TS_data.lon=lon_TS;
TS_data.lat=lat_TS;
TS_data.z=z_TS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 3 - Run OceanMesh2D %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 3.1 geodata and edgefx %%%%%%%%%%

for i=1:N
    min_el=min_el_all(i);
    max_el=max_el_all(i);
    max_el_ns=max_el_ns_all(i);
    R=R_all(i);
    grade=grade_all(i);
    alpha=alpha_all(i);
    coastline=convertStringsToChars(coastlines(i));
    dem=convertStringsToChars(dems(i));
    bbox=[d_limits(i,1:2);d_limits(i,3:4)];
    gdat{i}=geodata('shp',coastline,'dem',dem,'bbox',bbox,'h0',min_el);
    fh{i}=edgefx('geodata',gdat{i},'fs',R,'max_el_ns',max_el_ns,'max_el',max_el,'g',grade,'slp',alpha);
end

%%%%% 3.2 meshgen %%%%%%%%%%%%%%%%%%%%%

mshopts=meshgen('ef',fh,'bou',gdat,'plot_on',0,'proj','mercator','plot_on',0);
mshopts=mshopts.build;
m=mshopts.grd;

%%%%% 3.3 msh %%%%%%%%%%%%%%%%%%%%%%%%%

m=interp(m,gdat);
tri=m.t;
lon=m.p(:,1);lat=m.p(:,2);
z=m.b;

% Define Lat/Lon coords of the new referential's origin in km for KRAKEN
[lo,la,idx] = find_nearest_point(lon, lat, coord_source(1), coord_source(2)); % source coordinates and index in the mesh
coord_or = [lo la]; % Coords of source (in degrees)

% Other corrections to the bathymetry
z(isnan(z))=0;  % Convert nan nodes into zeros
z(z<0)=0;       % Convert negative numbers (land) into zeros
m.b=z;          %  Save correction to unstructured mesh struct variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 4 - Write mesh data %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_env=1;    % Flag to write/store ENV files (0 - no; 1 - yes)
flag_flp=0;    % Flag to write/store FLP files (0 - no; 1 - yes)

fname = 'azores';
inp_env.fname = fname; % prefix of ENV filenames
inp_flp.fname = fname; % FLP filename

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
mesh_data.z = z;
mesh_data.lon = lon;
mesh_data.lat = lat;
mesh_data.tri = tri;
mesh_data.pfix = coord_or;
save('mesh_data.mat','-struct','mesh_data');    % Save the MATLAB struct with mesh info (to be used in Plotting_Tagus_Estuary.m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 5 - Plotting mesh %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

view_tri=1;   % Flag to view mesh triangles
view_cb=1;    % Flag to view colorbar

% Plot the entire mesh (Azores archipelago)
figure
caxis_lims=[0 3e3]; % depth limits, in meters
lonlat_lims=d_limits(1,:);
plot_mesh_bathy(tri, z, lon, lat, lonlat_lims, caxis_lims,view_tri,view_cb)
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'g-','LineWidth',2.5)
lonlat_lims=d_limits(2,:);
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
lonlat_lims=d_limits(3,:);
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
lonlat_lims=d_limits(4,:);
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
saveas(gcf,['plots/' fname_env '_mesh_box1.png'])

% Plot mesh in boundary box 2 (Western Group of Azores)
figure
caxis_lims=[0 1e3]
view_cb=1;
lonlat_lims=d_limits(2,:);
view_tri=1;
plot_mesh_bathy(tri, z, lon, lat, lonlat_lims, caxis_lims, view_tri, view_cb)
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
saveas(gcf,['plots/' fname_env '_mesh_box2.png'])

% Plot mesh in boundary box 3 (Central Group of Azores)
figure
caxis_lims=[0 1e3]
view_cb=1;
lonlat_lims=d_limits(3,:);
view_tri=1;
plot_mesh_bathy(tri, z, lon, lat, lonlat_lims, caxis_lims, view_tri, view_cb)
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
saveas(gcf,['plots/' fname_env '_mesh_box3.png'])

% Plot mesh in boundary box 4 (Eastern Group of Azores)
figure
caxis_lims=[0 1500]
view_cb=1;
lonlat_lims=d_limits(4,:);
view_tri=1;
plot_mesh_bathy(tri, z, lon, lat, lonlat_lims, caxis_lims, view_tri, view_cb)
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
saveas(gcf,['plots/' fname_env '_mesh_box4.png'])
