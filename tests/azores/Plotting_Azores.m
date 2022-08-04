clear all, close all, clc

dirname='io_files_kraken/';
files = dir ([dirname '*.prt']);

S=shaperead('..\..\data\shorelines\CNTR_RG_01M_2020_4326.shp');
lo_s=[S.X]; la_s=[S.Y];

N=length({files.name});

N_mod=zeros(N,1);                   % array of the number of modes of all prt files
N_prt_mod=zeros(N,1);               % array of the printed nr of modes of all prt files
m=10;                               % number of modes to extract
Cg = zeros(N,m);                    % matrix of group velocities obtained by KRAKEN (first m modes at each point are extracted)

mesh=load('grid_info.mat');
z_m=mesh.z;
lon_m=mesh.lon;
lat_m=mesh.lat;
tri=mesh.tri;
xy_lims=[min(lon_m) max(lon_m) min(lat_m) max(lat_m)];

lo_s=[S.X]; la_s=[S.Y];
ff=find(lo_s>xy_lims(1)-0.1 & lo_s<xy_lims(2)+0.1 & la_s>xy_lims(3)-0.1 & la_s<xy_lims(4)+0.1);
lo_s=lo_s(ff);
la_s=la_s(ff);

for i=1:N
    fname = [dirname 'azores_' sprintf( '%06d', i ) '.prt'];
    fid = fopen(fname);
    lines=textread(fname,'%s','delimiter','\n');
    i_m=find(~cellfun(@isempty,strfind(lines,'Number of modes')));
    N_mod(i)=str2num(char(regexp(char(lines(i_m)),'\d*','Match')));     % number of modes
    i_0=find(~cellfun(@isempty,strfind(lines,'I    k (1/m)')))+1;
    i_f=find(~cellfun(@isempty,strfind(lines,'_')))-1;
    N_prt_mod(i)=length(i_0:i_f(end));                                  % number of printed modes
    lin=lines(i_0:i_f(end));                                            % lines for cg extraction
    for j=1:m
        vals=split(char(lin(j)));       % array of values in each column
        Cg(i,j)=str2num(char(vals(5))); % extract fifth column and store in the Cg matrix
    end
    fclose(fid);
end

figure
ff=find(Cg(:,1)>1690);
Cg(ff,1)=NaN;
trisurf(tri,lon_m,lat_m,Cg(:,1))
hold on
plot3(lo_s,la_s,ones(1,length(lo_s)).*1e6,'k.')
shading interp, view(2)
colorbar
colormap(jet)
caxis([1470 1550])
xlim([xy_lims(1) xy_lims(2)])
ylim([xy_lims(3) xy_lims(4)])
saveas(gcf,['plots/Cg_mode_1.png'])