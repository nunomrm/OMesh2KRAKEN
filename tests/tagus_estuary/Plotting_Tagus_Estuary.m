% This script has the purpose to display the plots of TL results of the Tagus Estuary test case of OMesh2KRAKEN
%
% (c) Nuno Monteiro and Tiago Oliveira, University of Aveiro
% September 2022
%

clear all, close all, clc

% addpaths
addpath(genpath('..\..\atWin10_2020_11_4\Matlab\ReadWrite'));   % to call read_shd.m in Acoustics Toolbox
addpath(genpath('..\..\utils'));
addpath(genpath('..\..\data'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Pre-processing %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'tagus_estuary'; % Filename of shade file (.shd)
tl_caxis=[20 120]; % Sound transmission loss (TL) limits, in dB
[ PlotTitle, ~, freqVec, ~, ~, Pos, pressure ] = read_shd([filename '.shd']); % Read shade file (function from Acoustics Toolbox)
tl = -20.0 * log10(double(abs(pressure)));      % Calculate TL
theta = Pos.theta;                              % Extract angles
dist = Pos.r.r;                                 % Extract radial distances
z = Pos.r.z;                                    % Extract depths
depths=[1 15];                                  % Select depths, in meters, for TL surface plots
S=shaperead('CNTR_RG_01M_2020_4326.shp');       % Read coastline shapefile

% Read mesh data (nodes, elements, depth)
mesh=load('mesh_data.mat');                       
z_m=mesh.z;
lon_m=mesh.lon;
lat_m=mesh.lat;
tri=mesh.tri;
pfix=mesh.pfix;

% Definition of plot/coastline limits and coordinates of source
x1=min(lon_m); x2=max(lon_m); y1=min(lat_m); y2=max(lat_m);
xy_lims=[x1 x2 y1 y2];
coord_src=[pfix(1) pfix(2)];
lo_s=[S.X]; la_s=[S.Y];
dom=[-9.45 -9.2 38.6 38.72];
ff=find(lo_s>xy_lims(1)-0.1 & lo_s<xy_lims(2)+0.1 & la_s>xy_lims(3)-0.1 & la_s<xy_lims(4)+0.1);
lo_s=[lo_s(ff) xy_lims(2)+0.2 xy_lims(2)+0.2 lo_s(1)];
la_s=[la_s(ff) xy_lims(4)+0.1 xy_lims(3)-0.1 la_s(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of pre-processing %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



N=length(depths)*length(dist);
for i=1:length(depths)
    iz=find(Pos.r.z==depths(i));
    tl_z=transpose(squeeze(tl(:,:,iz,:)));
    k=1;
    tl_n=zeros(N,1);
    for ii=1:length(dist)
        for jj=1:length(theta)
            tl_n(k) = tl_z(ii,jj);
            k=k+1;
        end
    end
    figure
    plot2D_TL_slice(tl_n,coord_src,dist,theta,tl_caxis,xy_lims)
    c=colorbar;
    c.Label.String = 'Transmission Loss (dB)';
    c.FontSize=14;
    hold on
    fill3(lo_s,la_s,ones(1,length(lo_s)).*1e6,'k')
    xlabel('Longitude','FontSize',15)
    ylabel('Latitude','FontSize',15)
    a=get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',14);
    saveas(gcf,['plots/TLfield_z_' num2str(depths(i)) 'm.png'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of plotting TL %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

