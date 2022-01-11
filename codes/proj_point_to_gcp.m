% Feb. 1, 2017, Yunfeng Chen, project a point onto the great circle path.
% this code could be used in case you want to project the station location
% onto a seismic profile
% convert geographic coordinate to cartesian coordinates
function [latp,lonp,d] = proj_point_to_gcp(lat0,lon0,lat1,lon1,lat,lon)
% input: lat0,lon0 starting point of GCP
%        lat1,lon1 end point of GCP
%        lat, lon  the target point needs to be projected
% output: latp,lonp the coordinates of point on the GCP after projection
%         d, distance between the point and GCP
R = 6371;
x = R*cosd(lon)*cosd(lat);
y = R*sind(lon)*cosd(lat);
z = R*sind(lat);
X = [x y z];
% define the starting and end points that define the GCP
x0 = cosd(lon0)*cosd(lat0);
y0 = sind(lon0)*cosd(lat0);
z0 = sind(lat0);
x1 = cosd(lon1)*cosd(lat1);
y1 = sind(lon1)*cosd(lat1);
z1 = sind(lat1);
% calculate the cross product of X0 and X1, two vectors correspond to
% starting and end points
X0 = [x0 y0 z0];
X1 = [x1 y1 z1];
V = cross(X0,X1);
% calculate normalized V vector
U = V/sqrt(sum(V.^2));
% obtain the distance between the given point and the plain defined by the
% GCP, which is the dot product of vector X and U
d = X * U';
% project the point to the GCP plane and extend it radially outward to the
% earth's surface
Xp = X - d * U;
Xp = Xp/sqrt(sum(Xp.^2));
xp = Xp(1);
yp = Xp(2);
zp = Xp(3);
% convert the coordinates back to lat lon
latp =  asind(zp/1);
if (xp > 0)
    lonp = atand(yp/xp);
elseif (yp > 0)
    lonp = atand(yp/xp) + 180;
else
    lonp = atand(yp/xp) - 180;
end
% figure
% plot(lon,lat,'*')
% hold on;
% plot(lon1,lat1,'*')
% plot(lon0,lat0,'*')
% plot([lon0 lon1],[lat0 lat1]);
% plot(lonp,latp,'g*')
