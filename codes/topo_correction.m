% Jan. 24, 2017, Yunfeng Chen, calculate the amount of topography
% correction assuming vertical incident ray path
% Feb. 2nd, 2017, consider non-vertical incidence for a reference RF 
function [TopoCorrections] = topo_correction(m)
% input m: model from crust 1.0
z = m(:,1);
vp = m(:,2);
vs = m(:,3);
EPS = 1e-6;
idisc = find( z(1:end-1) == z(2:end) );
z(idisc) = z(idisc) - EPS;
% find the velocity at 0 km depth
% vp0km = interp1(z,vp,0);
% vs0km = interp1(z,vs,0);
ztopo = [z(z<0); 0];
thi_topo = diff(ztopo);
vp_topo = vp(z<0);
vs_topo = vs(z<0);
p0 = 0.06*ones(size(vp_topo));
Temp = thi_topo.*sqrt(1./vs_topo.^2 - p0.^2) - thi_topo.*sqrt(1./vp_topo.^2 - p0.^2);
% Temp = thi_topo./vs_topo - thi_topo./vp_topo;
TopoCorrections = -sum(Temp);
end