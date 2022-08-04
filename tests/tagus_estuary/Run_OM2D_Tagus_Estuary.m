% This script runs the Tagus Estuary example for OMesh2kraken
clear all, close all, clc

% addpaths
addpath(genpath('..\..\OceanMesh2D'));
addpath(genpath('..\..\m_map'));
addpath(genpath('..\..\data'));
addpath(genpath('..\..\utils'));
addpath(genpath('..\..\atWin10_2020_11_4\Matlab\ReadWrite'));  % so that write_env.m function is called during writekraken3d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 1 - Input options %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 1.1 Inputs for OceanMesh2D %%%%%%%%%

d_limits=[-9.41 -9.24 38.6 38.72; -9.36 -9.27 38.64 38.7];              % boundary box (or mesh domains) limits
N=size(d_limits,1);                                                     % number of domains
min_el_all=[250 100];                                                   % minimum element size (max resolution)
max_el_all=[800 400];                                                   % maximum element size (min resolution)
max_el_ns_all=[400 250];                                                % maximum element size nearshore            
R_all=[2 2];                                                            % nr of elements to solve the feature size of coastlines
grade_all=[0.32 0.27];                                                  % mesh grade
dems=["EMODNET_DTM_SouthLisbon.nc"; "EMODNET_DTM_SouthLisbon.nc"];      % digital elevation model filenames, defining the input bathymetry
coastlines=["CNTR_RG_01M_2020_4326"; "CNTR_RG_01M_2020_4326"];          % coastline filenames
alpha_all = [0 25];                                                     % alpha parameter to solve the slope mesh size function
fname='tagus_estuary';                                                  % filename of output mesh file (FLP) and prefix name of ENV files
coord_source=[-9.317, 38.665];                                          % sound source coords in lat/lon

%%%%% 1.2 Inputs for KRAKEN %%%%%%%%%%

% definition of parameters
freq = 1000;                        % frequency (Hz)
NMedia = 2;                         % number of media
top_opt = "CVW .";                  % top option
bot_opt = "V";                      % bottom option
sigma_int1 = 0;                     % RMS roughness at the surface interface
sigma_int2 = 0;                     % RMS roughness at the bottom interface
N_layer1 = 10e3;                    % number of mesh points used in the internal discretization (water column)
N_layer2 = 10e3;                    % number of mesh points used in the internal discretization (bottom layer)
Rdepth = 150;                       % max depth of receiver (m)
sigma_bot = 0;                      % RMS roughness at the bottom
Sdepth = 10;                        % source depth (m)
C_low = 1400;                       % min phase sound speed (m/s)
C_high = 15000;                     % max phase sound speed (m/s)
deltassp = 1;                       % increment of the SSP depths
sndspd_bot = 1700;                  % soundspeed at bottom layer

%%%%% 1.3 Inputs for FIELD3D %%%%%%%%%%

title_flp = "Tagus Estuary"                 % title of run 
calc_flp = "STD";                           % type of calculation
Nm = 999;                                   % number of modes
Nsx = 1;                                    % number of source coords in x
Nsy = 1;                                    % number of source coords in y
coords_or_km = [0 0];                       % coords of the origin (km)
NSz = 1;                                    % number of source depths
Sz = 10;                                    % source depths (m)
Rzmin = 0;                                  % min receiver depth (m)
Rzmax = Rdepth;                             % max receiver depth (m)
dz = 1;                                     % receiver depth discretization
zR=Rzmin:dz:Rzmax;                          % receiver depths
NRz = length(zR);                           % number of receiver depths
Rrmin = 0;                                  % min receiver range (km)
Rrmax = 10;                                 % max receiver range (km)
dr = 30/1e3;                                % receiver range discretization (km)
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
params_flp={title_flp calc_flp Nm Nsx coords_or_km(1) Nsy coords_or_km(2) NSz Sz ...
    NRz Rzmin Rzmax NRr Rrmin Rrmax NRtheta Rthetamin Rthetamax NGBtheta ...
    GBthetamin GBthetamax GBstep GBsteps epsilon_mult};

% Temperature-Salinity data extraction and arrangement
fname_TSdata = 'cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m_1647891936780.nc';
lon_TS=ncread(fname_TSdata,'longitude');
lat_TS=ncread(fname_TSdata,'latitude');
t_TS=ncread(fname_TSdata,'time')./24+datenum(1950,1,1)
t_TS=t_TS(1); % extract the first instant
T=ncread(fname_TSdata,'thetao');
T=squeeze(T(:,:,:,1));
S=ncread(fname_TSdata,'so');
S=squeeze(S(:,:,:,1));
z_TS=ncread(fname_TSdata,'depth');
% concatenating the previous first depth level values of T and S to the zero level:
if z_TS(1)~=0      
    z_TS=[0; z_TS];
    S=cat(3,S(:,:,1),S); 
    T=cat(3,T(:,:,1),T);
end
% Create a struct out of the T-S input data to be read in the
% writekraken3d.m function:
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

%%%%% 3.1 geodata and edgefx  %%%%%%%%%

for i=1:N
    min_el=min_el_all(i); max_el=max_el_all(i); max_el_ns=max_el_ns_all(i);
    R=R_all(i); grade=grade_all(i); alpha=alpha_all(i);
    coastline=convertStringsToChars(coastlines(i)); dem=convertStringsToChars(dems(i));
    bbox=[d_limits(i,1:2);d_limits(i,3:4)];
    gdat{i}=geodata('shp',coastline,'dem',dem,'bbox',bbox,'h0',min_el);
    fh{i}=edgefx('geodata',gdat{i},'fs',R,'max_el_ns',max_el_ns,'max_el',max_el,'g',grade,'slp',alpha);
end

%%%%% 3.2 meshgen  %%%%%%%%%%%%%%%%%%%%

mshopts=meshgen('ef',fh,'bou',gdat,'plot_on',1,'proj','mercator','pfix',coord_source,'plot_on',0);
mshopts=mshopts.build;
m=mshopts.grd;

%%%%% 3.3 msh  %%%%%%%%%%%%%%%%%%%%%%%%

m=interp(m,gdat);
m=make_bc(m,'auto',gdat{1}); 
tri=m.t;
lon=m.p(:,1);lat=m.p(:,2);
z=m.b;

% define Lat/Lon coords of origin the new referential in km for KRAKEN
[lo,la,idx] = find_nearest_point(lon, lat, coord_source(1), coord_source(2)); % source coordinates and index in the mesh
coords_or = [lo la]; % coords of source (in degrees), which will be the origin of the referential in km coords

% other corrections to the bathymetry
z(isnan(z))=0;  % convert nan nodes into zeros
z(z<0)=0;       % convert negative numbers (land) into zeros
m.b=z;          % apply correction to mesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% end of Run OceanMesh2D %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 4 - Write mesh data %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write(m,fname,'kraken3d','nob',1,'fname_env',fname,'TS_data',TS_data, ...
    'params_env',params_env, 'coords_or',coords_or,'params_flp',params_flp);

if ~exist('io_files_kraken','dir')
    mkdir('io_files_kraken')                        % create a directory for input/output files of kraken if it doesnt exist
else
    delete('io_files_kraken\*')                     % remove everything from the input/ouput kraken directory
end

movefile('*.env','io_files_kraken\')                % move all envs to the input/ouput kraken directory

if ~exist('plots\','dir')
    mkdir('plots\')                                 % create a directory for plotting if it doesnt exist
end

grid_info.z=z;
grid_info.lon=lon;
grid_info.lat=lat;
grid_info.tri=tri;
grid_info.pfix=coords_or;
save('grid_info.mat','-struct','grid_info');    % save a MATLAB struct with mesh info (to be used in Plotting_Tagus_Estuary.m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% end of Write mesh data %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 5 - Plotting mesh (bathym.) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
caxis_lims=[0 50]; % meters
view_tri=1;view_cb=1;
lonlat_lims=d_limits(1,:);
plot_mesh_bathy(tri, z, lon, lat, lonlat_lims, caxis_lims,view_tri,view_cb)
plot3(lo,la,1e6,'ro','MarkerFaceColor','r')
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'g-','LineWidth',2.5)
lonlat_lims=d_limits(2,:);
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
saveas(gcf,['plots/' fname '_mesh_box1.png'])

figure
view_cb=1;
lonlat_lims=d_limits(2,:);
view_tri=1;
plot_mesh_bathy(tri, z, lon, lat, lonlat_lims, caxis_lims, view_tri, view_cb)
plot3(lo,la,1e6,'ro','MarkerFaceColor','r')
plot3([lonlat_lims(1) lonlat_lims(1) lonlat_lims(2) lonlat_lims(2) lonlat_lims(1)], ...
    [lonlat_lims(3) lonlat_lims(4) lonlat_lims(4) lonlat_lims(3) lonlat_lims(3)], ...
    ones(1,5).*1e6,'r-','LineWidth',2.5)
saveas(gcf,['plots/' fname '_mesh_box2.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% end of Plotting mesh %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%