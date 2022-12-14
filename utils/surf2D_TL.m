function surf2D_TL(tl, coord_src, dist, theta, tl_caxis, xy_lims)
%
% (c) 2022 Nuno Monteiro and Tiago Oliveira, University of Aveiro
% (Function for OMesh2KRAKEN - https://github.com/nunomrm/OMesh2KRAKEN/)
%
% Plots top-view transmission loss (dB) surface field for a given depth
% in lon/lat coordinates (deg)
% 
% INPUTS:
% tl - transmission loss data
% coord_src - coordinates of sound source in (longitude, latitude)

tej = flipud( jet( 256 ) );  % 'jet' colormap reversed

% convert polar coords to Lat/Lon
N=length(dist)*length(theta);
lat=zeros(N,1);
lon=zeros(N,1);
k=1;
for i=1:length(dist)
    for j=1:length(theta)
        d_x=km2deg(dist(i)*cosd(theta(j))/1e3);
        d_y=km2deg(dist(i)*sind(theta(j))/1e3);
        lat(k) = coord_src(2)+d_y;
        lon(k) = coord_src(1)+d_x;
        k=k+1;
    end
end
dt = delaunayTriangulation(lon,lat);
tri = dt.ConnectivityList;
xi = dt.Points(:,1);
yi = dt.Points(:,2);
F = scatteredInterpolant(lon,lat,tl);
tli = F(xi,yi) ;
trisurf(tri,xi,yi,tli)
shading interp, view(2), colormap( tej ) 
caxis([tl_caxis(1) tl_caxis(2)])
xlim([xy_lims(1) xy_lims(2)])
ylim([xy_lims(3) xy_lims(4)])
end
