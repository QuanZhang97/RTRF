%{
  Real data test for Radon transform based receiver function preprocessing
  
  written by Quan Zhang et al., 2022
  
  Copyright: Zhejiang University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  
  Reference: 
  
%}

clear; clc; close all;
addpath ./data
addpath(genpath('./codes'))
[z, rho, vp, vs, qk, qm] = ak135( 'cont' );
zmax = 800;
dz = 0.5;
tic;
%% Step 1: Load data
datadir = './';
gauss=2.5;
load(fullfile(datadir,['event_1.mat']));
load(fullfile(datadir,['event_6~16.mat']));
log=lognew;
sta=unique({log.sta});
% CCPData
stalist = sta;
for n = 1:length(stalist)
    index = find(strcmp(sta{n},{log.sta}));
    
    slat(n) = log(index(1)).slat;
    slon(n) = log(index(1)).slon;
end
%% Step 2: define the CCP grid
disp('Create CCP bin starts')
LatMin = min(slat)-0.2;
LatMax = max(slat)+0.2;
LonMin = min(slon)-0.2;
LonMax = max(slon)+0.2;
BinSpacing = 20;
BinSize = 50;
% call function to set up grid nodes
CCPGrid = ccp_setup_grid(LatMin,LatMax,LonMin,LonMax,BinSpacing,BinSize);
centers = [cell2mat(CCPGrid(:,2)) cell2mat(CCPGrid(:,1))];
radii = km2deg(cell2mat(CCPGrid(:,3)));
viscircles(centers,radii); hold on;
% plot station location
scatter(slon,slat,100,'^','filled')
disp('Create CCP bin completes')
%% read in 3D model, use crust 1.0 model
% The speed of extracting crust 1.0 model could be improved by using 
% getcrust function from seizmo
disp('Create velocity model starts')
X = [];
Y = [];
Z = [];
VP = [];
VS = [];
knode = 0;
latall = linspace(min([CCPGrid{:,1}]),max([CCPGrid{:,1}]),5);
lonall = linspace(min([CCPGrid{:,2}]),max([CCPGrid{:,2}]),5);
% plot it on the CCP map
[lontemp,lattemp] = meshgrid(lonall,latall);
% scatter(lontemp(:),lattemp(:),'ko','filled')
for i = 1:length(latall)
    for j = 1:length(lonall)
        knode = knode + 1;
        disp(['Now extracting velocity model at node #',num2str(knode)]);
        lat = latall(i);
        lon = lonall(j);
        % note that m0 starts at 0 km depth. The elevation infomration is
        % missing, to include topography set if_topo to 1
        if_topo = 0;
        [m0,nsedi] = obtain_crust1_v2(lat,lon,[],if_topo);
        % define the interface in the crust 1.0 model
        m1 = [];
        for l = 1:size(m0,1)
            if l == 1
                m1(l,:) = m0(l,:);
            else
                m1(2*(l-1),:) = [m0(l,1) m0(l-1,2:4)];
                m1(2*(l-1)+1,:) = m0(l,:);
            end
        end
        dmoho(knode) = m0(end,1);
        dmax = 1000;
        % find the depth beneath the moho from prem
%         prem_model = prem;
%         z_prem = prem_model.depth;
%         vp_prem = prem_model.vp;
%         vs_prem = prem_model.vs;
%         rho_prem = prem_model.rho;
%         keepz = z_prem > dmoho(knode) & z_prem < 1000;
        keepz=z>dmoho(knode) & z<dmax;
        % deal with the upper mantle, use the PREM model
        m1 = [m1; z(keepz) vp(keepz) vs(keepz) rho(keepz)];
        EPS = 1e-6;
        ztemp = m1(:,1);
        idisc = find( ztemp(1:end-1) == ztemp(2:end) );
        ztemp(idisc) = ztemp(idisc) - EPS;
        zpos = 0:dz:zmax;
        vptemp = interp1( ztemp, m1(:,2), zpos(1:end-1)+0.5*dz, 'linear','extrap');
        vstemp = interp1( ztemp, m1(:,3), zpos(1:end-1)+0.5*dz, 'linear','extrap');
        Z = [Z; zpos(1:end-1)'];
        X = [X; ones(size(vptemp'))*lon];
        Y = [Y; ones(size(vptemp'))*lat];
        VP = [VP; vptemp'];
        VS = [VS; vstemp'];
        % plot(vptemp,zpos(1:end-1)); hold on;
        % ylim([-5 100])
        % axis ij;
    end
end
% plot the model
% imagesc(1:knode,zpos(1:end-1),reshape(VS,length(zpos)-1,knode));
% interpolate 
Fvp = scatteredInterpolant(X,Y,Z,VP);
Fvs = scatteredInterpolant(X,Y,Z,VS);
disp('Create velocity model completes')
%% Step 3: 1D raytracing and time to depth conversion
disp('Migration starts')
start_index = 1;
cp_lat = [];
cp_lon = [];
for n = 1:length(stalist)
    disp(['Now processing station: ', stalist{n}]);
    seis = {};
    time = {};
%     load([datadir,file(nsta(n)).name]);
    keep=strcmp(sta{n},{log.sta}) & [log.snr]>=2;
    if sum(keep) == 0
        continue
    end
    good_RFs=log(keep);
    seis = {good_RFs.itr};
    nseis=length(seis);
    time = {good_RFs.ittax};
    lat = good_RFs(1).slat;
    lon = good_RFs(1).slon;
    backaz = [good_RFs.baz];
    p = [good_RFs.rayp]; 
    if mean(p)>1
        % deg to km
        p=km2deg(p);
    end
    % 1D ray tracing
    [cp, RayMatrix, MidPoints] = find_conversion_points_v2(p, backaz, dz, zmax, z, vp, vs, lat, lon, 'spherical');
    % save the conversion point location for GMT plot
    cp_lat = [cp_lat; cp.latb(80,:)'];
    cp_lon = [cp_lon; cp.lonb(80,:)'];
    RayDepths = 1*dz:dz:zmax;
    RayDepths = RayDepths(:);
    % corret for heterogenity
    [TimeCorrections, Tpds3D, Tpds1D] = correct_RFs(MidPoints, RayDepths, Fvp, Fvs, z, vp, vs);
    % reindex the RFs
    end_index = start_index + length(good_RFs) - 1;
    index = start_index:1:end_index;
    index_matrix = repmat(index,size(RayMatrix,1),1);
    RayMatrix(:,:,7) = index_matrix;
    % update the index
    start_index = end_index + 1;
    % map RFs to depth
    [timeout, seisout, depth0] = migrate_RFs( time, seis, p, dz, zmax, z, vp, vs, TimeCorrections);
    % normalization
    seisout_norm = normalize_RFs(seisout);
    % read in migrated RFs and depth0
    RayMatrix(:,:,1) = cell2mat(seisout_norm);
    RayMatrix(:,:,2) = repmat(depth0',1,nseis);
    % save the Raymatrix
    CCPMatFile = ['CCPData',stalist{n},'.mat'];
    save(CCPMatFile,'RayMatrix','MidPoints');
end
disp('Migration completes')
%% Step 4: find the RFs fall within each bin
disp('Assign RF starts')
RayMatrix_all = [];
for n = 1:length(stalist)
    filename=['CCPData',stalist{n},'.mat'];
    if exist(filename,'file')
        load(filename,'RayMatrix');
    else
        continue
    end
    % save matrix for individual station into a big matrix
    RayMatrix_all = [RayMatrix_all, RayMatrix];
%     nptsTotal = size(CCPGrid,1) * size(RayMatrix,1);
%     nptsUpdate = floor(nptsTotal / 100);
    kk = 0;
    for m = 1:size(CCPGrid,1) % Each Bin
        for k = 1:size(RayMatrix,1) % Each Depth
            lons = RayMatrix(k,:,4) * pi/180;
            lats = RayMatrix(k,:,3) * pi/180;
            tlon = CCPGrid{m,2} * pi/180;
            tlat = CCPGrid{m,1} * pi/180;
            dlon = lons - tlon;
            dlat = lats - tlat;
            a = (sin(dlat/2)).^2 + cos(lats) .* cos(tlat) .* (sin(dlon/2)).^2;
            angles = 2 .* atan2(sqrt(a),sqrt(1-a));
            dist = 6371 * angles;
            % It would be really easy to adjust the size of your bin as you
            % go to deeper depths.
            %--------------------------------------------------------------
            Indx = find(dist <= CCPGrid{m,3}); % This ONLY takes RFs falling into the bin
            Temp{k,1} = RayMatrix(k,Indx,7); % Record IDs
            TempDist{k,1} = dist(Indx);
            if ~isempty(Indx)
                disp('')
            end
        end
        CCPGrid{m,4} = [CCPGrid{m,4} Temp];
        CCPGrid{m,5} = [CCPGrid{m,5} TempDist];
        disp(['Processed bin: ',num2str(m),'/',num2str(size(CCPGrid,1))])
    end
    
end
disp('Assign RF completes')
%%
% plot the results (you may want to turn plotting off when processing a
% large number of stations)
figure()
% set(gcf,'Position',[100 100 1200 800])
set(gcf,'Position',[0 0 1600 800],'Color','w')
for n=1:1:size(RayMatrix_all,2)
    lattmp=RayMatrix_all(:,n,3);
    lontmp=RayMatrix_all(:,n,4);
    ztmp=RayMatrix_all(:,n,2);
    atmp=RayMatrix_all(:,n,1);
%     scatter(lattmp,ztmp,10,atmp,'filled'); hold on;
    scatter3(lontmp,lattmp,ztmp,30,atmp,'filled'); hold on;
end
viscircles(centers,radii); hold on;
% plot station location
scatter(slon,slat,200,'b^','filled')
set(gca,'Zdir','reverse')
xlabel(['Longitude', char(176)]);
ylabel(['Latitude', char(176)]);
zlabel('Depth (km)');
set(gca,'FontSize',16)
colormap(seismic(3))
caxis([-0.5 0.5])
colorbar
xlim([82 88])
ylim([28 35])
zlim([0 200])
view([140,20])
export_fig(gcf,['ray-path_multistations_gauss',num2str(gauss),'_src.png'])
% saveas(gcf, 'ray-path_multistations.png','png')
% interpolate 
ilat=29:0.05:34;
idep=0:0.5:200;
[LAT,DEP]=meshgrid(ilat,idep);
x=RayMatrix_all(:,:,3);
y=RayMatrix_all(:,:,2);
v=RayMatrix_all(:,:,1);
keep=~isnan(v(:)) & ~isinf(v(:));
AMP=griddata(x(keep),y(keep),v(keep),LAT,DEP);
K = (1/50)*ones(5,10);
AMP_smooth = conv2(AMP,K,'same');
% AMP_smooth = medfilt2(AMP);
% AMP_smooth = wiener2(AMP,[5 5]);
figure;
set(gcf,'Position',[0 0 1600 1600],'Color','w')
subplot(211)
% imagesc(ilat,idep,AMP);
scatter(x(:),y(:),5,v(:),'filled');
xlabel(['Latitude', char(176)]);
ylabel('Depth (km)');
ylim([0 200])
xlim([29 34])
axis ij;
set(gca,'FontSize',16)
colormap(seismic(3))
caxis([-0.5 0.5])
colorbar
subplot(212)
imagesc(ilat,idep,AMP_smooth);
xlabel(['Latitude', char(176)]);
ylabel('Depth (km)');
set(gca,'FontSize',16)
colormap(seismic(3))
caxis([-0.5 0.5])
colorbar
export_fig(gcf,['project_to_latitude_gauss',num2str(gauss),'_src.png'])
%% Step 5: ccp stacking
disp('CCP stacking starts')
for i = 1:size(CCPGrid,1) % For each bin
    disp(['Now stacking CCP bin #',num2str(i)]);
    RRFAmps = [];
    Weights = [];
    MaxHits(i) = 0;
    for k = 1:length(CCPGrid{i,4})
        temp(k) = numel([CCPGrid{i,4}{k,:}]);
    end
    MaxHits(i) = max(temp);
    clear temp
    if MaxHits(i) > 0
        % Build 2 matrices with MaxHits # columns and length(DepthAxis) rows
        Rtemp = NaN * ones(length(CCPGrid{i,5}),MaxHits(i));
        Wtemp = NaN * ones(length(CCPGrid{i,5}),MaxHits(i));
        % Add the value in the RayMatrix using the RF indices to the
        % matrices we just built.
        for k = 1:length(CCPGrid{i,4})
            Ids = cell2mat(CCPGrid{i,4}(k,:));
            %                 Dist = CCPGrid{n,4}{k};
            if ~isempty(Ids)
                for l = 1:length(Ids)
                    temp = find(RayMatrix_all(k,:,7) == Ids(l));
                    if ~isempty(temp)   % 06/12/2011: ADDED IF STATEMENT TO CHECK IF THE RECORD IDS EXIST IN THIS MATRIX
                        Rtemp(k,l) = RayMatrix_all(k,temp,1);
                        
                        % Ttemp(k,l) = RayMatrix(k,temp,2);
                        %                         Wtemp(k,l) = exp(-(Dist(l)).^2./(2.*CCPGrid{n,3}.^2));
                        Wtemp(k,l) = 1; % weight of the matrix
                    end
                end
            end
        end
        RRFAmps = [RRFAmps Rtemp];
        Weights = [Weights Wtemp];
        DepthAxis = depth0;
        ResampleNumber = 50;
        RandomSelectionCount = 0.8; %resample 80% of the original traces
        % apply bootstrap to each bin
        [RRFBootMean,RRFBootSigma,Peaks410,Peaks660,Amps410,Amps660] = ...
            ccp_bootstrap(RRFAmps,Weights,DepthAxis,ResampleNumber,round(size(RRFAmps,2)*RandomSelectionCount));
    else
        RRFBootMean = NaN * ones(length(CCPGrid{i,5}),1);
        RRFBootSigma = NaN * ones(length(CCPGrid{i,5}),1);
    end
    % save the bin
    BootAxis = depth0;
    LatCCP = CCPGrid{i,1};
    LonCCP = CCPGrid{i,2};
    CCPBinLead = ['CCPData' filesep 'Bin_'];

    CCPPath = [CCPBinLead, num2str(sprintf('%0.0f',i)), '.mat'];
    save([CCPPath],'RRFBootMean','RRFBootSigma','BootAxis','LatCCP','LonCCP');
end
disp('CCP stacking completes')
% latb = cp.latb;
% lonb = cp.lonb;
% scatter(lonb(:),latb(:),30,'b'); hold on;
% plat_s = Raymatrix(:,:,2);
% plon_s = Raymatrix(:,:,3);
% scatter(plon_s(:),plat_s(:),10,'r')
%% rough plot
rf = {};
nbin  = size(CCPGrid,1);
lat = [];
lon = [];
dep = [];
k = 0;
dep0 = BootAxis;
figure(1)
for n = 1:nbin
    % check if bin exists
    binfile = ['CCPData/Bin_',num2str(n),'.mat'];
    if exist(binfile,'file')
        k = k + 1;
        load(binfile);
        rf{k} = RRFBootMean;
        lat(k) = LatCCP;
        lon(k) = LonCCP;
    end
    scatter3(ones(size(rf{k}))*lon(k), ones(size(rf{k}))*lat(k), -dep0,30, rf{k}, 'filled' ); hold on;
end
rf1 = cell2mat(rf);
zlim([-400 0]);
% figure;
% imagesc(rf1); hold on;
% wigb(rf1,1,1:nbin,dep0)
% ylim([0 200]);
% colormap('jet')
% caxis([-0.1 0.1])
%% create the interpolation volume
k = 0;
V = [];
for i = 1:length(lat)
    for j = 1:length(depth0)
        k = k + 1;
        V(k,:) = [lon(i) lat(i) depth0(j) rf1(j,i)];
    end
end
F = scatteredInterpolant(V(:,1),V(:,2),V(:,3),V(:,4));
%% create 2D cross-section
lon1 = 85.80837236+0.; 
lat1 = 29.2673;
lon2 = 83.87268109+0.;
lat2 = 33.9673;
nlatlon = 100;
[latp,lonp] = gcwaypts(lat1,lon1,lat2,lon2,nlatlon);
[deg0,az0]= distance(lat1,lon1,latp,lonp);
% degree to distance
dist0 = deg0*2*pi*6371/360;
% plot the bin and profile location
% viscircles(centers,radii); hold on;
% plot station location
% select the ccp bin
keepm = [];
for m = 1:length(CCPGrid)
    lons = lonp * pi/180;
    lats = latp * pi/180;
    tlon = CCPGrid{m,2} * pi/180;
    tlat = CCPGrid{m,1} * pi/180;
    dlon = lons - tlon;
    dlat = lats - tlat;
    a = (sin(dlat/2)).^2 + cos(lats) .* cos(tlat) .* (sin(dlon/2)).^2;
    angles = 2 .* atan2(sqrt(a),sqrt(1-a));
    dist = 6371 * angles;
    Indx = find(dist <= CCPGrid{m,3});
    if ~isempty(Indx)
        keepm = [keepm m];
    end
end

% It would be really easy to adjust the size of your bin as you
% go to deeper depths.
%--------------------------------------------------------------
% obtain the staion index contribute to this cross-section
keep_sta_indx = [];
for m = 1:length(keepm)
    tempIndx = CCPGrid{keepm(m),4};
    % station loop
    for j = 1:size(tempIndx,2)
        if_keep = 0;
        % depth loop
        for i = 1:20
            if_keep = ~isempty(tempIndx{i,j});
            if if_keep == 1
                break;
            end
        end
        if if_keep == 1
           keep_sta_indx = [keep_sta_indx, j]; 
        end
    end
end
keep_sta_indx = unique(keep_sta_indx);
% plot the geometry of RFs
LAT_profile = repmat(latp',length(depth0),1);
LON_profile = repmat(lonp',length(depth0),1);
DEP_profile = repmat(depth0',1,length(latp));
DIST_profile = repmat(dist0',length(depth0),1);
Vprofile = F(LON_profile,LAT_profile,DEP_profile);
Vprofile(isnan(Vprofile)) = NaN; 
%%
figure
set(gcf,'Position',[0 0 1600 800],'Color','w')
scatter(slon,slat,200,'^','filled'); hold on;
scatter(lonp,latp,50,'ko','filled');
xlabel(['Longitude' ,char(176)])
ylabel(['Latitude' ,char(176)])
zlabel('Depth (km)')
set(gca,'FontSize',16)

viscircles(centers(keepm,:),radii(keepm)); hold on;
% viscircles(centers,radii); hold on;
surface(LON_profile,LAT_profile,-DEP_profile,Vprofile,'EdgeColor','none')
zlim([-200 0])
colormap(seismic(3))
caxis([-0.2 0.2])
colorbar
grid on; box on;
view(140,20)
export_fig(gcf,['ccp_3D_gauss',num2str(gauss),'_src.png'])
% saveas(gcf, 'example_ccp.png','png')
toc;

%%

figure
set(gcf,'Position',[0 0 800 400],'Color','w')
h=imagesc(LAT_profile(1,:),DEP_profile(:,1),Vprofile);
ylim([0 200])
xlabel(['Latitude' ,char(176)])
ylabel('Depth (km)')
set(gca,'FontSize',16)
set(gca,'YDir','reverse')
colormap(seismic(3))
% set(h,'alphadata',~isnan(Vprofile))
caxis([-0.2 0.2])
colorbar
