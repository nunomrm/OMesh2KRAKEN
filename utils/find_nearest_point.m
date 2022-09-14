function [lo,la,idx] = find_nearest_point(lon, lat, lon_p, lat_p)
%
% (c) 2022 Nuno Monteiro, University of Aveiro
%
% Finds the nearest geographical point between a given (Lon,Lat) location and a set of 
% N points in (lon,lat) coordinates. This function is originally intended for 2D finite 
% element meshes where point coordinates are stored in Nx1 arrays
%
%%% INPUTS
% lon - longitudes of the input set of points in a Nx1 array
% lat - latitudes of the input set of points in a Nx1 array
% lon_p - longitude of input point
% lat_p - latitude of input point
%
%%% OUTPUTS:
% lo - longitude of the nearest found point of the input array
% la - latitude of the nearest found point of the input array
% idx - index of the first dimension of the input Nx2 array (or mesh index, if it is a 
% finite element mesh)
%
%
anom=1e10;
idx=1;
for i =1:length(lon)
    if anom>abs(lon(i)-lon_p)+abs(lat(i)-lat_p)
        anom=abs(lon(i)-lon_p)+abs(lat(i)-lat_p);
        lo=lon(i);
        la=lat(i);
        idx=i;
    end
end
end
