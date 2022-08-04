function plot_mesh_bathy(tri, z, lon, lat,lonlat_lims, caxis_lims, view_tri, view_cb)
%
%
% (c) 2022 Nuno Monteiro
%
% Plot unstructured mesh bathymetry field
%

trisurf(tri,lon,lat,z)

hold on, shading interp, grid off
if view_tri==1
    trimesh(tri,lon,lat,z, 'facecolor', 'none', 'edgecolor', 'k')
end
caxis([caxis_lims(1) caxis_lims(2)]);
if view_cb==1
    c=colorbar;
    c.Label.String = 'Depth (m)';
end
map = brewermap(15,'Blues'); 
colormap(map)

xlim([lonlat_lims(1) lonlat_lims(2)])
ylim([lonlat_lims(3) lonlat_lims(4)])
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
view(0,90)
ax = gca;
ax.FontSize = 14; 
end

