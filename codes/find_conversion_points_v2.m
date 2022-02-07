% calculate the PmS conversion point of RF for a given velocity model
% Dec. 18, 2016, Yunfeng Chen, Global Seismology Group, University of
% Alberta. Modified after processRFmaster
% also calcuate the location of conversion points, which requires the Vs
% profile
% Jan. 22, 2017, Yunfeng Chen, consider the spherical earth for seismic ray
% traveling at deeper depth
% May 13, 2018, Y.C., CSRIO, add comments to the code, note the ray path is
% computed using 1D reference model this is because the ray-parameter was
% estimated using this model.  Rewrite the code to compute RayMatrix and
% Midpoints matrix more efficiently.
% Input: 
% p: ray parameter in s/km
% backaz: back-azimuth
% dz: depth increment
% zmax: maximum depth
% vp: P-wave velocity
% vs: S-wave velocity
% lat: station latitude
% lon: station longitude
% model_type: flat or spherical earth, default is flat 
% Output:
% cp: struct contains location of conversion point
% RayMatrix0: a nz*nx*7 ray piercing point matrix where the columns are
% RRF Amplitudes, depth (these two varaibles will be saved to matrix later),
% Ps - lat,lon, and Pp - lat,lon and index of receiver function (from 1 to nx)
% MidPoints0: a (nz-1)*nx*7 matrix where the columns are Ps - lon,lat,raylength
% and Pp - lon,lat,raylength, index of receiver function (from 1 to nx).
function [cp, RayMatrix0, MidPoints0] = find_conversion_points_v2( p, backaz, dz, zmax, z, vp, vs, lat, lon, model_type)
global dha dhb
%- Code:
% model_type: 'flat' or 'spherical'
if nargin == 9
    model_type = 'flat';
end

EPS = 1e-6;
R = 6371; % radius of earth
nx = numel( p );
% force rayp to be row vector
if( nx > 1 )
  if( size(p,1) > size(p,2) ) p = p'; end
end

% get the depths
zpos = (0.0:dz:zmax)';
nz = numel(zpos);
% deal with discontinuities in the vel model
idisc = find( z(1:end-1) == z(2:end) );
z(idisc) = z(idisc) - EPS;

r = repmat(R - zpos(1:end-1),1,nx);

% interpolate the vel model between each layer
vp = interp1( z, vp, zpos(1:end-1)+0.5*dz, 'linear','extrap');
vs = interp1( z, vs, zpos(1:end-1)+0.5*dz, 'linear','extrap');

% repeat matrices
p = repmat( p, nz-1, 1 );
% convert rayp to s/rad
if strcmp(model_type,'spherical')
    p = p*R;
end

backaz = repmat( backaz, nz-1, 1 );
vp = repmat( vp, 1, nx );
vs = repmat( vs, 1, nx );
% get associated vertical slowness
switch model_type
    case 'flat'
        % get horizontal position
        qa = sqrt(1./vp.^2 - p.^2);
        qb = sqrt(1./vs.^2 - p.^2);
        dha = p.*dz./qa;
        dhb = p.*dz./qb;
        dla=dz./(qa.*vp);
        dlb=dz./(qb.*vs);
    case 'spherical'
        qa = sqrt((r./vp).^2 - p.^2);
        qb = sqrt((r./vs).^2 - p.^2);
        dha = p.*dz./(r.*qa);
        dhb = p.*dz./(r.*qb);
        dla=r.*dz./(qa.*vp);
        dlb=r.*dz./(qb.*vs);
        % convert rad to degree
        dha = dha .* (360/(2*pi));
        dhb = dhb .* (360/(2*pi));
        % convert to km
        dha = dha .* (2*pi*R)/360;
        dhb = dhb .* (2*pi*R)/360;
end

% search for imaginary values in the distance calculations
for n = 1:size(dha,2)
    StopIndex = find(imag(dha(:,n)),1);
    if ~isempty(StopIndex)
        dha(StopIndex:end,n) = NaN * ones(length(dha(:,n))-StopIndex+1,1);
        dhb(StopIndex:end,n) = NaN * ones(length(dhb(:,n))-StopIndex+1,1);
    end
end

hposa = cumsum( dha, 1 );   % horizontal position of P phase piercing point
eposa = [zeros(1,nx) ; hposa.*sind(backaz)]; % easting of P phase piercing point
nposa = [zeros(1,nx) ; hposa.*cosd(backaz)]; % northing of P phase piercing point
la = dla;  % length of raypath segment of P phase

hposb = cumsum( dhb, 1 );   % horizontal position of Ps phase piercing point
eposb = [zeros(1,nx) ; hposb.*sind(backaz)]; % easting of Ps phase piercing point
nposb = [zeros(1,nx) ; hposb.*cosd(backaz)]; % northing of Ps phase piercing point
lb = dlb;  % length of raypath segment of P phase

% calculate the lat lon of the conversion point
arclena = km2deg(hposa);
arclenb = km2deg(hposb);
az = backaz;

lat = lat*ones(size(dha));   % latitude of station
lon = lon*ones(size(dha));   % longitude of station
[lata,lona] = reckon(lat,lon,arclena,az);
[latb,lonb] = reckon(lat,lon,arclenb,az);

% output conversion point information as a structure array
cp.lata = lata;
cp.lona = lona;
cp.latb = latb;
cp.lonb = lonb;
cp.zpos = zpos(2:end);
cp.hposa = hposa;
cp.nposa = nposa;
cp.eposa = eposa;
cp.la=la;
cp.hposb = hposb;
cp.nposb = nposb;
cp.eposb = eposb;
cp.lb=lb;
%% Get the RayMatrix and MidPoints matrices
nz=length(zpos);
RayMatrix0=NaN*ones(length(zpos),nx,7);
RayMatrix0(2:end,:,3)=[cp.latb];
RayMatrix0(2:end,:,4)=[cp.lonb];
RayMatrix0(2:end,:,5)=[cp.lata];
RayMatrix0(2:end,:,6)=[cp.lona];
% now deal with the first point, which is the station location
RayMatrix0(1,:,3)=lat(1,:);
RayMatrix0(1,:,4)=lon(1,:);
RayMatrix0(1,:,5)=lat(1,:);
RayMatrix0(1,:,6)=lon(1,:);
RayMatrix0(1:end,:,7)=ones(nz,1)*[1:nx];
% compute the mid points, which is roughly the average of each ray segment
MidPoints0=0.*ones(length(zpos)-1,nx,7);
MidPoints0(:,:,3)=[cp.lb];  % length of raypath segment of Ps
MidPoints0(:,:,6)=[cp.la];  % length of raypath segment of P
for n=1:nz-1
   MidPoints0(n,:,1)=(RayMatrix0(n,:,3)+RayMatrix0(n+1,:,3))/2;  % Ps midpoint lat
   MidPoints0(n,:,2)=(RayMatrix0(n,:,4)+RayMatrix0(n+1,:,4))/2;  % Ps midpoint lon
   MidPoints0(n,:,4)=(RayMatrix0(n,:,5)+RayMatrix0(n+1,:,5))/2;  % P midpoint lon
   MidPoints0(n,:,5)=(RayMatrix0(n,:,6)+RayMatrix0(n+1,:,6))/2;  % P midpoint lon
   MidPoints0(n,:,7)=RayMatrix0(n,:,7);  % Index of RF
end
%% compare with funclab result
%{
DepthMax = zmax;
DepthIncrement = dz;
YAxisRange = 0:(DepthIncrement/2):DepthMax;
Vp = interp1(zpos(1:end-1),vp(:,1),YAxisRange,'linear','extrap')';
Vs = interp1(zpos(1:end-1),vs(:,1),YAxisRange,'linear','extrap')';
Depths = YAxisRange';
Record = 0;
dz = [0; diff(Depths)];
for n = 1:nx
    Record = Record + 1;
    RayTracingMatrix(:,Record,:) = NaN*ones(length(YAxisRange),8);
    % 1D ray tracing
    %----------------------------------------------------------------------
    slat = lat(1,1);
    slon = lon(1,1);
    baz = backaz(1,n);
    rayp = p(1,n);
    R = 6371 - Depths;
    x_s = cumsum((dz./R) ./ sqrt((1./(rayp^2.* (R./Vs).^-2)) - 1));
    raylength_s = (dz.*R) ./  (sqrt(((R./Vs).^2) - (rayp^2)).* Vs);
    x_p = cumsum((dz./R) ./ sqrt((1./(rayp^2.* (R./Vp).^-2)) - 1));
    raylength_p = (dz.*R) ./  (sqrt(((R./Vp).^2) - (rayp^2)).* Vp);

    % Search for imaginary values in the distance calculations of the
    % P-wave
    %----------------------------------------------------------------------
    StopIndex = find(imag(x_p),1);
    if ~isempty(StopIndex)
        x_p(StopIndex:end) = NaN * ones(length(x_p)-StopIndex+1,1);
        x_s(StopIndex:end) = NaN * ones(length(x_s)-StopIndex+1,1);
    end
    % convert distance in km to distance in degrees along a Great Circle
    % path
    %----------------------------------------------------------------------
    gcarc_dist_s = rad2deg(x_s);
    gcarc_dist_p = rad2deg(x_p);

    % Find latitude of piercing points using law of cosines
    %----------------------------------------------------------------------
    pplat_s = asind((sind(slat) .* cosd(gcarc_dist_s)) + (cosd(slat) .* sind(gcarc_dist_s) .* cosd(baz)));
    pplat_p = asind((sind(slat) .* cosd(gcarc_dist_p)) + (cosd(slat) .* sind(gcarc_dist_p) .* cosd(baz)));
    % Do trig calculations
    %----------------------------------------------------------------------
    c1 = cosd(gcarc_dist_s);
    c2 = cosd(90 - slat);
    c3 = cosd(90 - pplat_s);
    c4 = sind(gcarc_dist_s);
    c5 = sind(baz);
    c6 = cosd(pplat_s);
    c7 = cosd(gcarc_dist_p);
    c8 = cosd(90 - pplat_p);
    c9 = sind(gcarc_dist_p);
    c10 = cosd(pplat_p);
    for k = 1:length(gcarc_dist_s)
        % Find longitude of piercing points using law of sines
        %------------------------------------------------------------------
        if isnan(gcarc_dist_s(k))
            pplon_p = NaN;
            pplon_s = NaN;
            raylength_p(k) = NaN;
            raylength_s(k) = NaN;
        else
            if ( c1(k) >= (c2 .* c3(k)))
                pplon_s = slon + asind(c4(k) .* c5 ./ c6(k));
            else
                pplon_s = slon + asind(c4(k) .* c5 ./ c6(k)) + 180;
            end
            if ( c7(k) >= (c2 .* c8(k)))
                pplon_p = slon + asind(c9(k) .* c5 ./ c10(k));
            else
                pplon_p = slon + asind(c9(k) .* c5 ./ c10(k)) + 180;
            end
        end
        RayTracingMatrix(k,Record,2) = pplat_s(k);
        RayTracingMatrix(k,Record,3) = pplon_s;
        RayTracingMatrix(k,Record,4) = raylength_s(k);
        RayTracingMatrix(k,Record,5) = pplat_p(k);
        RayTracingMatrix(k,Record,6) = pplon_p;
        RayTracingMatrix(k,Record,7) = raylength_p(k);
        RayTracingMatrix(k,Record,8) = n;
    end
end
% lontmp=RayTracingMatrix(:,:,3);
% lattmp=RayTracingMatrix(:,:,2);
% scatter(lontmp(:),lattmp(:))
% lontmp1=cp.lonb;
% lattmp1=cp.latb;
% hold on;
% scatter(lontmp1(:),lattmp1(:),'r')

% Make the MidPoint Matrix where the columns are S - lat,lon,raylength
% and P - lat,lon,raylength.  Each of these rays are referenced to the
% bottom depth of the segment.
count = 0;
% create MidPoints matrix
MidPoints = zeros(nz-1,nx,7);
for n = 2:2:size(RayTracingMatrix,1)-1
    count = count + 1;
    MidPoints(count,:,:) = RayTracingMatrix(n,:,(2:end));
    seglength_s = RayTracingMatrix(n,:,4) + RayTracingMatrix(n+1,:,4);
    seglength_p = RayTracingMatrix(n,:,7) + RayTracingMatrix(n+1,:,7);
    MidPoints(count,:,3) = seglength_s;
    MidPoints(count,:,6) = seglength_p;   
end
% Make the ray piercing point matrix where the columns are RRF,TRF
% Amplitudes, S - lat,lon, and P - lat,lon.
count = 0;
RayMatrix = zeros(nz,nx,7);
for n = 1:2:size(RayTracingMatrix,1)
    count = count + 1;
    RayMatrix(count,:,:) = RayTracingMatrix(n,:,[1 1 2 3 5 6 8]);
end
% CCPMatFile = 'CCPData.mat';
% save(CCPMatFile,'RayMatrix','MidPoints');
figure;
set(gcf,'Position',[50 50 1200 500])
subplot(131)
plot(RayMatrix(2,:,4),RayMatrix(2,:,3),'ko'); hold on;
plot(RayMatrix0(2,:,4),RayMatrix0(2,:,3),'rx'); hold on;
xlabel('Longitude (deg)')
ylabel('Latitude (deg)');
legend({'My function','Funclab'})
subplot(132)
plot(MidPoints(1,:,2),MidPoints(1,:,1),'ko'); hold on;
plot(MidPoints0(1,:,2),MidPoints0(1,:,1),'rx'); hold on;
xlabel('Longitude (deg)')
ylabel('Latitude (deg)');
subplot(133)
plot(zpos(2:end),MidPoints(:,1,3),'ko'); hold on;
plot(zpos(2:end),MidPoints0(:,1,3),'rx'); hold on;
xlabel('Depth (km)')
ylabel('Ray segment length (km)');
%}
return;