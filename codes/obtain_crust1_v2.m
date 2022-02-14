% Calculate the average crustal velocity model using crust 1.0 model
% FEB 13 2022, Originally written by Yunfeng Chen, modified by Quan Zhang
function varargout = obtain_crust1_v2(varargin)
lat = varargin{1};
lon = varargin{2};
H = varargin{3};
if_topo = varargin{4};
if nargin < 3
    H = 0;
end
if nargin < 4
    if_topo = 0;
end
% lat = 54.5;
% lon = -115.5;
nla=180;
nlo=360;
nd=9;
mdldir='./crust1.0/';
vp=load([mdldir,'crust1.vp']);
vs=load([mdldir,'crust1.vs']);
rho=load([mdldir,'crust1.rho']);
bnd=load([mdldir,'crust1.bnds']);
vp1=zeros(nd,nla,nlo);
vs1=zeros(nd,nla,nlo);
rho1=zeros(nd,nla,nlo);
bnd1=zeros(nd,nla,nlo);
nline=1;
for j = 1:nla
   for i = 1:nlo
      for k = 1:nd
         vp1(k,j,i)=vp(nline,k);
         vs1(k,j,i)=vs(nline,k);
         rho1(k,j,i)=rho(nline,k);
         bnd1(k,j,i)=bnd(nline,k);
      end
      nline=nline+1;
   end
end
ilat = int16(floor(90.-lat))+1;
ilon = int16(floor(180.+lon))+1;
vptmp=vp1(:,ilat,ilon);
vstmp=vs1(:,ilat,ilon);
rhotmp=rho1(:,ilat,ilon);
ztmp=bnd1(:,ilat,ilon);
% remove the first two layers since they are water and ice, start with
% the third layer
% sediments is from layer 3 to layer 5, and crust is from layer 6 to 8
% remember converts the crustal depth to thickness (total thickness), such
% that the depth starts at 0 km
if if_topo == 0
    model = [-(ztmp(3:end)-ztmp(1)) vptmp(3:end) vstmp(3:end) rhotmp(3:end)];
else
    model = [-ztmp(3:end) vptmp(3:end) vstmp(3:end) rhotmp(3:end)];
end
% remove the layer with 0 velcoity or density (non-existing layer) as well
% as 0 thickness
thick = diff(model(:,1));
% deal with the upper mantle, assuming its a half space
thick(end+1) = 999;
knock_out = model(:,2) == 0 |  model(:,3) == 0 | model(:,4) == 0 | thick == 0;
% find the number sedimentary layers that exist
nsedi = sum(vptmp(3:5)>0 & thick(1:3) >0);
model(knock_out,:) = [];
% set the Moho depth to H
if H ~= 0
    model(end,1) = H;
end
varargout{1} = model;
varargout{2} = nsedi;
end
% visualize the velocity model
% plot_velocity_model(model(:,1),model(:,2),model(:,3))
% model = [2.2 5.7504  3.3200  2.6101; 9.8 6.6338  3.8300  2.8928; 
%     24.0 6.4259  3.7100  2.8263; 17.0 6.7723  3.9100  2.9371;
%     24.0 7.4825  4.3200  3.1644; 43.0 7.8115  4.5100  3.2697;
%     0.0 7.8115  4.5100  3.2697 ];
% model = [0 2.6000    1.5900    1.4600; 
%     2  6.1000    3.5300    2.7400;
%     12.0  6.9000    3.9300    2.9200;
%     22  8.1700    4.5300    3.3600];
% model = [0  6.1000    3.5300    2.7400; 20 8.1700    4.5300    3.3600];
%{
%% generate the synthetics
% convert depth to layer thickness, note the last layer needs to be 0
% thickness
thiki = [diff(model(:,1)); 0];
alpha = model(:,2);
beta = model(:,3);
rho = model(:,4);
% calcualte rho based on emperical relationship rho = 0.32*vp+0.77;
% rho = alpha*0.32+0.77;


% write out the model
dum1 = 0;
dum2 = 0;
dum3 = 0;
dum4 = 0;
% calculate the Poission's ratio
dum5 = vpovs_to_pr(alpha./beta,'forward');
% dum5 = 0.25;

% dir = '/home/yunfeng/new_rfun/RPB_stacking/';
dir = '/media/disk2/project/RF/synthetics/';
fid = fopen([dir,'crust1_model'],'w');
for n = 1:length(thiki)
    if n == 1
    fprintf(fid,'%2d %s\n',length(thiki),'crust1_model');
    end
    fprintf(fid,'%3d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n',...
             n,alpha(n),beta(n),rho(n),thiki(n),dum1,dum2,dum3,dum4,dum5(n));
end
fclose(fid);
% plot the velocity model
% plot_velocity_model(thiki,alpha,beta);
%% generate synthetics using Ammon's fortran code
modelpath = '/media/disk2/project/RF/synthetics/crust1_model';
gaussian  = 5;
waterlevel = 0.01;
rayp = 6.4/111.1775;
dt = 0.1;
sig_leng = 50;
ph = 5;
[RF] = generate_RF(modelpath,gaussian,waterlevel,rayp,dt,sig_leng,ph,'waterlevel');
figure;
subplot(311)
plot(RF.t,RF.rftn);hold on;
xlim([-10 50])
%% read in synthetic seismogram
% [t, synr, ~] = fget_sac([dir,'crust1_model_sp.r']);
% [t, synz, ~] = fget_sac([dir,'crust1_model_sp.z']);
% nt = length(t);
% VB = 0;
% [wlr,wlrms,nwl] = makeRFwater_ammon(synr,synz,ph,dt,nt,waterlevel,gaussian,VB);
% wltax = (dt*(0:1:length(t)-1)-ph)';
% plot(wltax,wlr,'r')
% xlim([-10 50])
%% generate synthetics using Jacobsen's matlab code
zs = model(2:end,1);
VPs = model(1:end,2);
VSs = model(1:end,3);
rhos = model(1:end,4);
p = 6.4;
ts = -5:0.1:50;
Gaussian_factor = 5;
rot_ang = NaN;
wlevel = 0.01;
addpath /media/disk2/project/RF/matlab_codes/forward_RF_Jacobsen
[Q_rf,L_rf,Q,L,ts_nice,lump]=...
    rf_forward2(zs,VPs,VSs,rhos,p,ts,Gaussian_factor,rot_ang);
[r_rf,up_rf,r,up,ts_nice,lump]=...
    rf_forward(zs,VPs,VSs,rhos,p,ts,Gaussian_factor,wlevel);
subplot(312)
plot(ts_nice,r_rf,'r');
% hold on;plot(ts_nice,Q_rf,'r')
%%
% seis = load('/usr/local/seis/IRFFMv1_linux/outmodel.0.0576.5.0.eqr.asc');
% subplot(313)
% plot(seis(:,1),seis(:,2))
% xlim([-10 50])
%}
