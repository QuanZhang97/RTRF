% Oct. 4, 2020, Yunfeng Chen, migrate the synthetic RFs from FD
% Oct. 6, 2020, migrate multiple shots
% Oct. 9, 2020, least-squares migration with CG
% Dec. 30, 2021, Quan Zhang, radon transform preprocessing

clear; clc; close all;
addpath ./data
addpath ./codes
%% define the grid
load('./model.mat');
datadir='./';
% interpolate the model
dx=1;
dy=dx/2;
x=x/1000; % convert from meters to kilometers
y=y/1000;
% pad the model
xpad=150;
xmax=max(x);
xi=0-xpad:dx:xmax+xpad;
yi = y(1):dy:y(end);

[X,Y]=meshgrid(x,y);
[XI,YI]=meshgrid(xi,yi);

keep1=xi>=0 & xi<=xmax;
Vp1=interp2(X,Y,Vp,XI(:,keep1),YI(:,keep1),'linear');
Vs1=interp2(X,Y,Vs,XI(:,keep1),YI(:,keep1),'linear');

Vp1_left=repmat(Vp1(:,1),1,sum(xi<0));
Vp1_right=repmat(Vp1(:,end),1,sum(xi>xmax));

Vs1_left=repmat(Vs1(:,1),1,sum(xi<0));
Vs1_right=repmat(Vs1(:,end),1,sum(xi>xmax));

Vp1=[Vp1_left Vp1 Vp1_right];
Vs1=[Vs1_left Vs1 Vs1_right];

% interp the Moho
dmoho=interp1(x,dmoho(1,:)/1000,xi);

% smooth the velocity model
N=31;
[Vp2,w]=moving_avg(Vp1,N,'constant',2);
[Vp2,w]=moving_avg(Vp2,N,'constant');
[Vs2,w]=moving_avg(Vs1,N,'constant',2);
[Vs2,w]=moving_avg(Vs2,N,'constant',2);

x=xi;
z=yi;
dz=dy;
vel=Vp2/1000;
vel_s=Vs2/1000;
[nz,nx]=size(vel);
img=zeros(nz,nx);
% dx=x(2)-x(1);
% dz=z(2)-z(1);
% x = [0:1:nx-1]*dx;
% z = [0:1:nz-1]*dz;
% pad
figure;
subplot(1,5,1)
plot(Vp1_left(:,1)/1000,z);hold on
plot(Vp1_right(:,1)/1000,z)
plot(Vs1_left(:,1)/1000,z)
plot(Vs1_right(:,1)/1000,z)
xlim([0 10])
%set(gcf,'Position',[0 0 250 500],'color','w')
title('Velocity model');
ylabel('Depth (km)')
xlabel('Velocity (km/s)')
ax=gca;ax.YDir='reverse';
set(gca,'fontsize',14)
legend("Left Vp","Right Vp","Left Vs","Right Vs",'Location','SouthEast')
subplot(1,5,2:5)
imagesc(x,z,Vp2/1000);
set(gcf,'Position',[0 0 900 500],'color','w')
title('Vp model profile');
xlim([0 400])
xlabel('Distance (km)')
ylabel('Depth (km)','rotation',-90)
set(gca,'fontsize',14)
set(gca,'YAxisLocation','right')
ch=colorbar('location','EastOutside');
set(get(ch,'title'),'string','Velocity (km/s)')
export_fig(['./figure/velocity_model.png']);

%% loop over all shots (events)
inc_angles=-20:4:20;
inc_angles(6)=[];
nshot=length(inc_angles);
dmig=zeros(nz,nx,nshot);
dmigls=zeros(nz,nx,nshot);
dmigls1=zeros(nz,nx,nshot);
resample_period=0.1;
flow=0.05;
fhigh=2.0;
% fhighs=1.0:0.5:5.0;
fhighs=1.2;
pwave=cell(1,nshot);
tshift=cell(1,nshot);
src=cell(1,nshot);
for ifreq=1:1
    fhigh=fhighs(ifreq);
    suffix=[];
for ishot=1:1
    disp(['Processing shot ', num2str(ishot)]);
    filename_x=['evt',num2str(ishot),'_vx.su'];
    filename_y=['evt',num2str(ishot),'_vy.su'];
    [sekr,hx] = readsegy(fullfile(datadir,filename_x));
    [sekv,hy] = readsegy(fullfile(datadir,filename_y));
    
    [nt, ntr]=size(sekr);
    rx=[hx.gx]/10^6; % convert to km
    rz=[hx.gy]/10^6; % convert to km
    dt=hx(1).dt/1000000;
    t=(0:nt-1)*dt;
    % deconvolution
    gauss=2.5;
    itmax = 100;
    minderr = 0.001;
    ph = 5; % phase delay
    VB = 0; % Allow verbose output
    % i1=8/dt+1;
    i1=1;
    sekr=sekr./max(sekr(:));
    % revere the sign if incidence angle>0
    if inc_angles(ishot)>0
        sekv=-sekv./max(sekv(:));
    else
        sekv=sekv./max(sekv(:));
    end
    sekv = bandpassSeis(sekv, dt, flow, fhigh, 4);
    sekr = bandpassSeis(sekr, dt, flow, fhigh, 4);
    
    [sekv, ~, ~] = resampleSeis( sekv, t, 0.1);
    [sekr, dt, t] = resampleSeis( sekr, t, 0.1);
    
    sekv=sekv(:,1:2:end);
    sekr=sekr(:,1:2:end);
    
    [nt, ntr]=size(sekr);
    itr=zeros(size(sekr));
       parfor n=1:ntr
        % bandpass filter
        seisr=sekr(i1:end,n);
        seisz=sekv(i1:end,n);
        seis=[seisr seisz];
        seis=taper(seis,0.05,0.05);
        [itr(:,n),itrms] = makeRFitdecon_la_norm(sekr(i1:end,n),sekv(i1:end,n),dt,nt-i1+1,ph,gauss,itmax,minderr);
    end
    ittax=(dt*(0:1:length(t)-i1)-ph);
    % save 'rf.mat' 'itr' 'ittax' 'rx' 'ry'
    % downsample the RF
    [itr_resample, dt, ittax] = resampleSeis( itr, ittax, resample_period );
    nt=length(ittax);
    
    fig1=figure;
    set(gcf,'Position',[0 0 450*4 300*2],'Color','w')
    subplot(241)
    imagesc(rx,t,sekv)
    title('Vertical (clean)');xlabel('Distance (km)');ylabel('Time (sec)')
    colorbar
    caxis([-0.2 0.2])
    ylim([0 40])
    text(-0.2,0.98,'a)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    colormap('gray');
 
    
    subplot(245)
    imagesc(rx,t,sekr)
    title('Horizontal (clean)');xlabel('Distance (km)');ylabel('Time (sec)')
    colorbar
    caxis([-0.2 0.2])
    ylim([0 40])
    text(-0.2,0.98,'e)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    colormap('gray');
end
end

%% add noise

k=0.03;

random_noise_1=randn(size(sekv))*k;
random_noise_2=randn(size(sekr))*k;

sekr_noise = sekr+random_noise_1;
sekv_noise = sekv+random_noise_2;

    % deconvolution
    gauss=2.5;
    itmax = 100;
    minderr = 0.001;
    ph = 5; % phase delay
    VB = 0; % Allow verbose output
    % i1=8/dt+1;
    i1=1;
    [nt, ntr]=size(sekr);
    sekr_noise=sekr_noise./max(sekr_noise(:));
    % revere the sign if incidence angle>0
    if inc_angles(ishot)>0
        sekv_noise=-sekv_noise./max(sekv_noise(:));
    else
        sekv_noise=sekv_noise./max(sekv_noise(:));
    end
    itr_noise=zeros(size(sekr_noise));
    for n=1:ntr
        % bandpass filter
        seisr=sekr_noise(i1:end,n);
        seisz=sekv_noise(i1:end,n);
        seis=[seisr seisz];
        seis=taper(seis,0.05,0.05);
        [itr_noise(:,n),itrms] = makeRFitdecon_la_norm(sekr_noise(i1:end,n),sekv_noise(i1:end,n),dt,nt-i1+1,ph,gauss,itmax,minderr);
    end
    nt=length(ittax);
    
%% plot noisy decon
    subplot(242)
    imagesc(rx,t,sekv_noise)
    title('Vertical (noisy)');xlabel('Distance (km)');ylabel('Time (sec)')
    colorbar
    caxis([-0.2 0.2])
    ylim([0 40])
    text(-0.2,0.98,'b)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    colormap('gray');
 
    
    subplot(246)
    imagesc(rx,t,sekr_noise)
    title('Horizontal (noisy)');xlabel('Distance (km)');ylabel('Time (sec)')
    colorbar
    caxis([-0.2 0.2])
    ylim([0 40])
    text(-0.2,0.98,'f)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    colormap('gray');

%% radon denoising


% Set radon parameters
    dp=0.0005;
    xx=rx(1:2:end);
    h = xx; nh = length(h);
    p = [-0.08:dp:0]; np = length(p);


    t=ittax;dt=t(2)-t(1);
    nt=length(t);
    m = zeros(nt,np);

    Param.h=h;
    Param.v=1./p;
    Param.nt=nt;
    Param.dt=dt;
    Param.type=1;
    ma=zeros(nt,np);
    N1 = 3;  % CG Iterations (Internal loop);
    N2 = 1;   % Update of weights for the sparse solution: N1 = 1 LS;  N2 > 3 for High Res (Sparse) solution
    [ma] = radon_op(sekr_noise,Param,-1);
    
    %denoise for sekr
    [mi,misfit] = yc_pcg(@radon_op,Param,sekr_noise,zeros(size(ma)),N1,N2,1);
    [dp] = radon_op(mi,Param,1);
    sekr_de = dp;
    
    subplot(247)
    imagesc(p,t,mi)
    ylim([0 40])
    title('Horizontal (high-resolution Radon domain)');xlabel('p');ylabel('\tau (sec)')
    colormap(seismic(1))
    colorbar
    text(-0.2,0.98,'g)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)

    
    %denoise for sekv
    [mi,misfit] = yc_pcg(@radon_op,Param,sekv_noise,zeros(size(ma)),N1,N2,1);
    [dp] = radon_op(mi,Param,1);
    sekv_de = dp;
    

    subplot(243)
    imagesc(p,t,mi)
    ylim([0 40])
    title('Vertical (high-resolution Radon domain)');xlabel('p');ylabel('\tau (sec)')
    colormap(seismic(1))
    colorbar
    text(-0.2,0.98,'c)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    
    [nt, ntr]=size(sekr);
    sekr_de=sekr_de./max(sekr_de(:));
    % revere the sign if incidence angle>0
    if inc_angles(ishot)>0
        sekv_de=-sekv_de./max(sekv_de(:));
    else
        sekv_de=sekv_de./max(sekv_de(:));
    end
    itr_de=zeros(size(sekr_de));
    for n=1:ntr
        % bandpass filter
        seisr=sekr_de(i1:end,n);
        seisz=sekv_de(i1:end,n);
        seisz = bandpassSeis(seisz, dt, flow, fhigh, 4);
        seisr = bandpassSeis(seisr, dt, flow, fhigh, 4);
        [itr_de(:,n),itrms] = makeRFitdecon_la_norm(seisr,seisz,dt,nt-i1+1,ph,gauss,itmax,minderr);
    end
    nt=length(ittax);
%% plot denoised decon

    subplot(244)
    imagesc(rx,t,sekv_de)
    title('Vertical (denoised)');xlabel('Distance (km)');ylabel('Time (sec)')
    colorbar
    caxis([-0.2 0.2])
    ylim([0 40])
    text(-0.2,0.98,'d)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    colormap('gray');
 
    
    subplot(248)
    imagesc(rx,t,sekr_de)
    title('Horizontal (denoised)');xlabel('Distance (km)');ylabel('Time (sec)')
    colorbar
    caxis([-0.2 0.2])
    ylim([0 40])
    text(-0.2,0.98,'h)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    colormap('gray');
    export_fig(['./figure/compare_synthetics.png']);
 

%% plot rf

    figure
    subplot(131)
    set(gcf,'Position',[0 0 450*3 300],'Color','w')
    imagesc(rx,ittax,itr_resample)
    xlim([120 280])
    ylim([-5 20])
    title('Receiver function (noisy)');xlabel('Distance (km)');ylabel('Time (sec)')
    colorbar
    caxis([-0.2 0.2])
    text(-0.2,0.98,'a)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    cm=colormap('gray');

    subplot(132)
    imagesc(rx,ittax,itr_noise)
    xlim([120 280])
    ylim([-5 20])
    title('Receiver function (noisy)');xlabel('Distance (km)');ylabel('Time (sec)')
    colorbar
    caxis([-0.2 0.2])
    text(-0.2,0.98,'b)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    cm=colormap('gray');

    subplot(133)
    imagesc(rx,ittax,itr_de)
    xlim([120 280])
    ylim([-5 20])
    
    title('Receiver function (noisy)');xlabel('Distance (km)');ylabel('Time (sec)')
    colorbar
    caxis([-0.2 0.2])
    text(-0.2,0.98,'c)','Units','normalized','FontSize',18)
    set(gca,'FontSize',12)
    cm=colormap('gray');
    export_fig(['./figure/compare_synthetics_rf.png']);
    
% %% compare a single trace
% figure;
% plot(ittax,itr(:,100)/max(itr(:,100)),ittax,itr_noise(:,100)/max(itr_noise(:,100)),ittax,itr_de(:,100)/max(itr_de(:,100)));
% legend({'Raw','Noise','Radon'})
