% This script has the purpose to display plots of Group Velocity field (for January) of of the Azores test case of OMesh2KRAKEN
%
% (c) 2022 Nuno Monteiro and Tiago Oliveira, University of Aveiro
% 
%

clear all, close all, clc

dirname='data_kraken/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Pre-processing %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_env = 'azores';               % prefix name of all ENV files
files = dir ([dirname '*.prt']);
N=length({files.name});             % number of ENV files
S=shaperead('..\..\data\shorelines\CNTR_RG_01M_2020_4326.shp');
lo_s=[S.X]; la_s=[S.Y];

% Read mesh data (nodes, elements, depth)
mesh=load('mesh_data.mat');
z_m=mesh.z;
lon_m=mesh.lon;
lat_m=mesh.lat;
tri=mesh.tri;
xy_lims=[min(lon_m) max(lon_m) min(lat_m) max(lat_m)];

% Initialize Cg matrix for first 10 modes
N_mod = zeros(N,1); % array of the number of modes of all prt files
N_prt_mod = zeros(N,1); % array of the printed nr of modes of all prt files
m=10; % number of modes to extract
Cg = zeros(N,m); % matrix of group velocities obtained by KRAKEN (first m modes at each point are extracted)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plotting Group Velocity (Cg) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N
    fname = [dirname fname_env '_' sprintf( '%06d', i ) '.prt'];
    fid = fopen(fname);
    lines = textread(fname,'%s','delimiter','\n');
    i_m = find(~cellfun(@isempty,strfind(lines,'Number of modes')));
    N_mod(i) = str2num(char(regexp(char(lines(i_m)),'\d*','Match'))); % number of modes
    i_0 = find(~cellfun(@isempty,strfind(lines,'I    k (1/m)')))+1;
    i_f = find(~cellfun(@isempty,strfind(lines,'_')))-1;
    N_prt_mod(i) = length(i_0:i_f(end)); % number of printed modes
    lin = lines(i_0:i_f(end)); % lines for cg extraction
    for j = 1:m
        vals = split(char(lin(j))); % array of values in each column
        Cg(i,j) = str2num(char(vals(5))); % extract fifth column and store in the Cg matrix
    end
    fclose(fid);
end

trisurf(tri,lon_m,lat_m,Cg(:,1))
shading interp, view(2)
colormap(turbo(20))
caxis([1500 1530])
hold on
fill3(lo_s,la_s,ones(1,length(lo_s)).*1e6,'k')
xlim([xy_lims(1) xy_lims(2)])
ylim([xy_lims(3) xy_lims(4)])
saveas(gcf,['plots/Cg_mode_1.png'])
