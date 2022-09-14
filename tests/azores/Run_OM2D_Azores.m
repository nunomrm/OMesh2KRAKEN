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

d_limits=[-32 -24.5 36 40.5; -31.5 -30.7 39 40;
            -29 -26.7 38 39.4; -26.2 -24.8 36.7 38.1];                  % boundary box (or mesh domains) limits
N=size(d_limits,1);                                                     % number of domains
min_el_all=[10e3 3e3 3e3 3e3];                                          % minimum element size (max resolution)
max_el_all=[22e3 15e3 15e3 15e3];                                       % maximum element size (min resolution)
max_el_ns_all=[12e3 7e3 7e3 7e3];                                       % maximum element size nearshore            
R_all=[2 3 3 3];                                                        % nr of elements to solve the feature size of coastlines
grade_all=[0.3 0.2 0.2 0.2];                                            % mesh grade
dems=["gebco_azores_2020.nc"; "EMODNET_DTM_Azores_1.nc";
        "EMODNET_DTM_Azores_2.nc"; "EMODNET_DTM_Azores_3.nc"];          % digital elevation model filenames, defining the input bathymetry
coastlines=["GSHHS_h_L1"; "GSHHS_h_L1"; "GSHHS_h_L1"; "GSHHS_h_L1"];    % coastline filenames
alpha_all = [0 25 25 25];                                               % alpha parameter to solve the slope mesh size function                                                        % filename of output mesh file (FLP) and prefix name of ENV files
coord_source=[-29 38];                                                  % sound source coords in lat/lon


%%%%% 1.2 Inputs for KRAKEN %%%%%%%%%%

freq = 25;                          % frequency (Hz)
NMedia = 2;                         % number of media
top_opt = "CVW .";                  % top option
bot_opt = "V";                      % bottom option
sigma_int1 = 0;                     % RMS roughness at the surface interface
sigma_int2 = 0;                     % RMS roughness at the bottom interface
N_layer1 = 10e3;                    % number of mesh points used in the internal discretization (water column)
N_layer2 = 10e3;                    % number of mesh points used in the internal discretization (bottom layer)
Rdepth = 6e3;                       % max depth of receiver (m)
sigma_bot = 0;                      % RMS roughness at the bottom
Sdepth = 10;                        % source depth (m)
C_low = 1400;                       % min phase sound speed (m/s)
C_high = 15000;                     % max phase sound speed (m/s)
deltassp = 10;                      % increment of the SSP depths
sndspd_bot = 1700;                  % soundspeed at bottom layer
fname_env='azores'; 

%%%%% 1.3 Inputs for FIELD3D %%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% end of Input options %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 2 - Input processing %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cell array containing parameters for the ENV file provided by the user in 1.2
params_env = {freq NMedia top_opt sigma_int1 sigma_int2 N_layer1 N_layer2 ...
    Rdepth Sdepth sigma_bot C_low C_high deltassp sndspd_bot};

% Cell array containing parameters for the FLP file provided by user in 1.3
params_flp={title_flp calc_flp Nm Nsx coord_or_km(1) Nsy coord_or_km(2) NSz Sz ...
    NRz Rzmin Rzmax NRr Rrmin Rrmax NRtheta Rthetamin Rthetamax NGBtheta ...
    GBthetamin GBthetamax GBstep GBsteps epsilon_mult};

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
% concatenating the previous first depth level values of T and S to the zero level:
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
%%%%%%%% end of Input processing %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 3 - Run OceanMesh2D %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 3.1 geodata and edgefx %%%%%%%%%%

for i=1:N
    min_el=min_el_all(i); max_el=max_el_all(i); max_el_ns=max_el_ns_all(i);
    R=R_all(i); grade=grade_all(i); alpha=alpha_all(i);
    coastline=convertStringsToChars(coastlines(i)); dem=convertStringsToChars(dems(i));
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
coord_or = [lo la]; % Coords of source (in degrees), which will be the origin of the referential in km coords

% other corrections to the bathymetry
z(isnan(z))=0;  % Convert nan nodes into zeros
z(z<0)=0;       % Convert negative numbers (land) into zeros
m.b=z;          % Apply correction to mesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% end of Run OceanMesh2D %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 4 - Write mesh data %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_env=1;    % Flag to write/store ENV files (0 - no; 1 - yes)
flag_flp=0;    % Flag to write/store FLP files (0 - no; 1 - yes)

% Struct array with all ENV inputs
inp_env.flag = flag_env;
inp_env.params = params_env;
inp_env.TS_data = TS_data;

% Struct array with all FLP inputs (not needed, except for the flag)
inp_flp.coord_or = coord_or;
inp_flp.flag = flag_flp;
inp_flp.params = params_flp;
fnames = {fname_env fname_flp};

% Write data to ENV files for KRAKEN3D
write(m,fnames,'kraken3d','nob',1,'inp_env',inp_env,'inp_flp',inp_flp);

% Create a directory for input/output data of KRAKEN, if it doesnt exist
if ~exist('data_kraken','dir')
    mkdir('data_kraken')
else
    delete('data_kraken\*')    % Remove everything from the KRAKEN data directory
end

% Move all envs to the input/ouput kraken directory
movefile('*.env','data_kraken\')

% Create a directory for plotting if it doesnt exist
if ~exist('plots\','dir')
    mkdir('plots\')
end

% Create a Struct with mesh data (bathymetry, nodes' position, etc.)
mesh_data.z=z;
mesh_data.lon=lon;
mesh_data.lat=lat;
mesh_data.tri=tri;
save('mesh_data.mat','-struct','mesh_data');    % Save the MATLAB struct with mesh info (to be used in Plotting_Azores.m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% end of Write mesh data %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% end of Plotting mesh %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
