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
addpath ./codes

gauss=2.5;
load('newcolormap.mat')
load('event_1.mat')
for n=1:length(log)
       log(n).id = [num2str(log(n).year,'%04d'),num2str(log(n).jday,'%03d'),...
           num2str(log(n).H,'%02d'),num2str(log(n).m,'%02d')];
end
eventid=unique({log.id});
%% prepare common shot gather
lon1 = 85.80837236+0.; 
lat1 = 29.2673;
lon2 = 83.87268109+0.;
lat2 = 33.9673;
nlatlon = 100;
[latp,lonp] = gcwaypts(lat1,lon1,lat2,lon2,nlatlon);
[deg0,az0]= distance(lat1,lon1,latp,lonp);
% degree to distance
dist0 = deg0*2*pi*6371/360;
lognew=[];
% for ievt=1:length(eventid)
for ievt=1%6:16
    disp(['Processing event# ',num2str(ievt)])
    keep = strcmp({log.id},eventid{ievt});
    sub=log(keep);
    slat=[sub.slat];
    slon=[sub.slon];
    slatp=[];
    slonp=[];
    d=[];
    for i=1:length(slat)
        [slatp(i),slonp(i),d(i)] = proj_point_to_gcp(lat1,lon1,lat2,lon2,slat(i),slon(i));
    end
    % calculate the distance along the profile
    [deg0,az0]= distance(lat1,lon1,slatp,slonp);
    dist0 = deg0*2*pi*6371/360;
    
    % plot original RFs
    t=sub(1).ittax;
    figure;
    set(gcf,'Position',[100 100 1400 400],'Color','w');
    subplot(131)
    m_proj('lambert','long',[82 88],'lat',[28 35]); hold on;
    m_scatter(slon,slat,30,'^'); hold on;
    m_plot(lonp,latp,'k')
    m_scatter(slonp,slatp,30,'r^'); hold on;
    m_gshhs('l','line','color','k','linewidth',1)
    m_gshhs('lb2','line','color','k')
    m_grid('linewidth',2,'tickdir','out',...
        'xaxisloc','bottom','yaxisloc','right','fontsize',14);
    text(-0.1,0.98,'(a)','Units','normalized','FontSize',18)
    subplot(1,3,[2:3])
%     wigb([sub.itr],3,dist0,t);
    imagesc(dist0,t,[sub.itr]);
    colormap(seismic(3));
    caxis([-0.3 0.3])
    colorbar
    xlabel('Distance (km)');
    ylabel('Time (sec)');
    set(gca,'fontsize',14)
    title(['#RF=',num2str(length(sub))])
    ylim([-5 30])
    text(-0.1,0.98,'(b)','Units','normalized','FontSize',18)
    export_fig(['./figure/raw.png']);
    %% Set radon parameters
    dp=0.005;
    h = [sub.dist]; nh = length(h);
    p = [-2:dp:2]; np = length(p);

    t=sub(1).ittax;
    dt=t(2)-t(1); nt=length(t);
    m = zeros(nt,np);

    Param.h=h;
    Param.v=1./p;
    Param.nt=nt;
    Param.dt=dt;
    Param.type=1;
    ma=zeros(nt,np);

    %% Apply LS radon to Z component
    d = [sub.Z];
    N1 = 10;  % CG Iterations (Internal loop);
    N2 = 1;   % Update of weights for the sparse solution: N1 = 1 LS;  N2 > 3 for High Res (Sparse) solution
    [mi,misfit] = yc_pcg(@radon_op,Param,d,zeros(size(ma)),N1,N2,1);
    % Prediction (Use inverted velocity gather to estimate a preduction of the
    % data)
    [ma] = radon_op(d,Param,-1);
    [dp ] = radon_op(mi,Param,1);
    dp_z=dp;
    
    figure;
    set(gcf,'Position',[0 0 450*4 300*2],'color','w')
    subplot(241)
    wigbcyk(d,2,h,t)
    title('Vertical (raw)');
    xlabel('Distance (deg)')
    ylabel('Time (sec)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(a)','Units','normalized','FontSize',18)
    
    subplot(242)
    imagesc(p,t,ma/max(abs(ma(:))))
    title('Vertical (Radon)');
    xlabel('p')
    ylabel('\tau (sec)')
    set(gca,'fontsize',14)
    text(-0.2,0.98,'(b)','Units','normalized','FontSize',18)
    caxis([-1,1]);
    colormap(newcolormap)
    colorbar
    
    subplot(243)
    imagesc(p,t,mi/max(abs(mi(:))))
    title('Vertical (HRT)');
    xlabel('p')
    ylabel('\tau (sec)')
    set(gca,'fontsize',14)
    text(-0.2,0.98,'(c)','Units','normalized','FontSize',18)
    caxis([-1,1]);
    colormap(newcolormap)
    colorbar
    
    subplot(244)
    wigbcyk(dp,2,h,t)
    title('Vertical (denoised)');
    xlabel('Distance (km)')
    ylabel('Time (deg)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(d)','Units','normalized','FontSize',18)
    %% Apply LS radon to R component
    d=[sub.R];
    N1 = 10;  % CG Iterations (Internal loop);
    N2 = 1;   % Update of weights for the sparse solution: N1 = 1 LS;  N2 > 3 for High Res (Sparse) solution
    [mi,misfit] = yc_pcg(@radon_op,Param,d,zeros(size(ma)),N1,N2,1);
    % Prediction (Use inverted velocity gather to estimate a preduction of the
    % data)
    [ma] = radon_op(d,Param,-1);
    [dp ] = radon_op(mi,Param,1);
    dp_r=dp;
    

    subplot(245)
    wigbcyk(d,2,h,t)
    title('Horizontal (raw)');
    xlabel('Distance (deg)')
    ylabel('Time (sec)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(e)','Units','normalized','FontSize',18)
    
    
    subplot(246)
    imagesc(p,t,ma/max(abs(ma(:))))
    title('Horizontal (Radon)');
    xlabel('p')
    ylabel('\tau (sec)')
    set(gca,'fontsize',14)
    text(-0.2,0.98,'(f)','Units','normalized','FontSize',18)
    colormap(newcolormap)
    caxis([-1,1]);
    colorbar
    
    subplot(247)
    imagesc(p,t,mi/max(abs(mi(:))))
    title('Horizontal (HRT)');
    xlabel('p')
    ylabel('\tau (sec)')
    set(gca,'fontsize',14)
    text(-0.2,0.98,'(g)','Units','normalized','FontSize',18)
    colormap(newcolormap)
    caxis([-1,1]);
    colorbar
    
    subplot(248)
    wigbcyk(dp,2,h,t)
    title('Horizontal (denoised)');
    xlabel('Distance (deg)')
    ylabel('Time (sec)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(h)','Units','normalized','FontSize',18)
    export_fig(['./figure/compare_realdata.png']);
    %% deconvolution
    itmax = 400;
    minderr = 0.001;
    ph = 5; % phase delay
    VB = 0; % Allow verbose output
    for n=1:size(dp_r,2)
        R=dp_r(:,n);
        Z=dp_z(:,n);
        [itr,itrms] = makeRFitdecon_la_norm(R,Z,dt,nt,ph,gauss,itmax,minderr);
        sub(n).itr1=itr';
    end
    
    figure;
    set(gcf,'Position',[0 0 1200 700],'color','w')
    subplot(2,3,1:2)
    wigbcyk([sub.itr],4,h,t)
    title('Receiver function')
    xlabel('Distance (deg)')
    ylabel('Time (sec)')
    ylim([-5 100])
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(a)','Units','normalized','FontSize',18)


    subplot(2,3,4:5)
    wigbcyk([sub.itr1],4,h,t)
    title('Preprocessed receiver function')
    xlabel('Distance (deg)')
    ylabel('Time (sec)')
    ylim([-5 100])
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(b)','Units','normalized','FontSize',18)    
    
    
    subplot(2,3,3)
    wigbcyk([sub.itr],4,h,t)
    %title('Receiver function')
    xlabel('deg')
    ylabel('sec')
    xlim([29 31])
    ylim([-5 30])
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(c)','Units','normalized','FontSize',18)


    subplot(2,3,6)
    wigbcyk([sub.itr1],4,h,t)
    xlabel('deg')
    ylabel('sec')
    xlim([29 31])
    ylim([-5 30])
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(d)','Units','normalized','FontSize',18)    
    %% plot RFs
    t=sub(1).ittax;
    figure;
    set(gcf,'Position',[100 100 1400 400],'Color','w');
    subplot(131)
    m_proj('lambert','long',[82 88],'lat',[28 35]); hold on;
    m_scatter(slon,slat,30,'^'); hold on;
    m_plot(lonp,latp,'k')
    m_scatter(slonp,slatp,30,'r^'); hold on;
    m_gshhs('l','line','color','k','linewidth',1)
    m_gshhs('lb2','line','color','k')
    m_grid('linewidth',2,'tickdir','out',...
        'xaxisloc','bottom','yaxisloc','right','fontsize',14);
    text(-0.1,0.98,'(a)','Units','normalized','FontSize',18)

    subplot(1,3,[2:3])
%     wigb([sub.itr1],3,dist0,t);
    imagesc(dist0,t,[sub.itr1]);
    colormap(seismic(3));
    caxis([-0.3 0.3])
    colorbar
    xlabel('Distance (km)');
    ylabel('Time (sec)');
    set(gca,'fontsize',14)
    title(['#RF=',num2str(length(sub))])
    ylim([-5 30])
    text(-0.1,0.98,'(b)','Units','normalized','FontSize',18)
    export_fig(['./figure/rf_radon_decon.png']);
    
    %% save the new RFs
    lognew=[lognew sub];
   % close all;
end

%% compare post-processing
    figure
    set(gcf,'Position',[0 0 800 1050],'color','w')
    subplot(311)
    imagesc(h,t,[sub.itr])
    colormap(seismic(3));
    caxis([-0.3 0.3])
    title('Raw receiver function')
    xlabel('Distance (deg)')
    ylabel('Time (sec)')
    ylim([-5 30])
    colorbar
    set(gca,'fontsize',14)
    text(-0.1,0.98,'(a)','Units','normalized','FontSize',18)
    % deconvolution
    % plot original RFs
    subplot(312)
    d=[sub.itr];
    K = (1/10)*ones(2,5);
    d_smooth = conv2(d,K,'same');
    t=sub(1).ittax;
%     wigb([sub.itr],3,dist0,t);
    imagesc(dist0,t,d_smooth);
    colormap(seismic(3));
    caxis([-0.3 0.3])
    colorbar
    xlabel('Distance (deg)');
    ylabel('Time (sec)');
    set(gca,'fontsize',14)
    title('Moveing average')
    ylim([-5 30])
    text(-0.1,0.98,'(b)','Units','normalized','FontSize',18)
    % Apply LS radon to original RFs
    d=[sub.itr];
    N1 = 10;  % CG Iterations (Internal loop);
    N2 = 1;   % Update of weights for the sparse solution: N1 = 1 LS;  N2 > 3 for High Res (Sparse) solution
    [mi,misfit] = yc_pcg(@radon_op,Param,d,zeros(size(ma)),N1,N2,1);
    % Prediction (Use inverted velocity gather to estimate a preduction of the
    % data)
    [dp ] = radon_op(mi,Param,1);
    
    % plot RF after radon
    t=sub(1).ittax;
    subplot(313)
%     wigb(dp,3,dist0,t);
    imagesc(dist0,t,dp);
    colormap(seismic(3));
    caxis([-0.3 0.3])
    xlabel('Distance (deg)');
    ylabel('Time (sec)');
    text(-0.1,0.98,'(c)','Units','normalized','FontSize',18)
    set(gca,'fontsize',14)
    title('Radon transform')
    ylim([-5 30])
    colorbar
    export_fig(['./figure/post_processing_radon_decon.png']);
    
    
