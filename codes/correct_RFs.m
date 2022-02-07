% Jan. 24, 2017, Yunfeng Chen, remigrate the RFs based on some regional
% earth model (e.g., CRUST 1.0), modified after ccp_correction code,
% the code also calculate the amount of corrections (travel time difference)
% between 1D and 3D velocity models
function [TimeCorrections, Tpds3D, Tpds1D] = correct_RFs( MidPoints, RayDepths, Fvp, Fvs, z, vp, vs)
% Input:
% MidPoints: a nz*nseis*7 matrix where the columns are S - lon,lat,raylength
% and P - lon,lat,raylength, index of ray (from 1 to nseis)
% RayDepths: Depths axis from dz to zmax
% Fvp: struct for Scattered data interpolation (P wave)
% Fvs: struct for Scattered data interpolation (S wave)
% z: depth defined in 1D model
% vp: P-wave velocity defined in 1D model
% vs: S-wave velocity defined in 1D model
% Output:
% TimeCorrections: the amount of time correction
% Tpds3D: Ps travel time in 3D model
% Tpds1D: Ps travel time in 1D model
for m = 1:size(MidPoints,2)
    % interpolate for P velocities
    xi = MidPoints(:,m,5);
    yi = MidPoints(:,m,4);
    zi = RayDepths;
    StopIndex = find(isnan(xi),1);
    if ~isempty(StopIndex)
        xi = xi(1:StopIndex-1); yi = yi(1:StopIndex-1); zi = zi(1:StopIndex-1);
        dvp = Fvp(xi,yi,zi);
        dvp = [dvp;(NaN * ones(size(MidPoints,1)-StopIndex+1,1))];
    else
        dvp = Fvp(xi,yi,zi);
    end
    % interpolate for S velocities
    xi = MidPoints(:,m,2); yi = MidPoints(:,m,1); zi = RayDepths;
    StopIndex = find(isnan(xi),1);
    if ~isempty(StopIndex)
        xi = xi(1:StopIndex-1); yi = yi(1:StopIndex-1); zi = zi(1:StopIndex-1);
        dvs = Fvs(xi,yi,zi);
        dvs = [dvs;(NaN * ones(size(MidPoints,1)-StopIndex+1,1))];
    else
        dvs = Fvs(xi,yi,zi);
    end
    EPS = 1e-6;
    
    idisc = find( z(1:end-1) == z(2:end) );
    z(idisc) = z(idisc) - EPS;
    CVp = interp1(z,vp,RayDepths);
    CVs = interp1(z,vs,RayDepths);
    % a quick comparison of 1D and 3D models
%     figure;
%     set(gcf,'Position',[100 100 600 800],'Color','w')
%     plot(CVp,RayDepths,'k'); hold on; plot(dvp,RayDepths,'k--')
%     plot(CVs,RayDepths,'b'); plot(dvs,RayDepths,'b--')
%     axis ij
%     xlabel('Velocity (km/sec)'); ylabel('Depth (km)')
%     set(gca,'fontsize',14);
%     legend({'AK135 Vp','GYPSUMS Vp','AK135 Vp','GYPSUMS Vp'})
% 	export_fig('model_1D_3D_comparison.png')
    
    dls = MidPoints(:,m,3); dlp = MidPoints(:,m,6);
    Tpds3D = dls./dvs - dlp./dvp;
    
    Tpds1D = dls./CVs - dlp./CVp;
    Temp = (dls./dvs - dls./CVs) - (dlp./dvp - dlp./CVp);
    Temp(isnan(Temp)) = 0;
    TimeCorrections(:,m) = cumsum([0; Temp]);    
end
end

