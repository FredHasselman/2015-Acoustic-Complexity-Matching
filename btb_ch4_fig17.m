%%Supplementary Material to Hasselman (2014) "Classifying Complex Patterns into Speech Categories"
%
%%% Introduction
% This is a demonstration script accompanying the article "Classifying Complex Patterns into Speech Categories". Its purpose is to
% provide an example of how to use various freely available MATLAB sources on the web to extract variables from speech
% stimuli. This script should thus be used as an example to build your own scripts only. It is not a function or toolbox, not
% optimized for speed or functionality! Evaluate the code per Cell (using MATLAB's Cell Mode) and inspect the workspace to
% see what is going on.
%
% See OSF project for details: <https:osf.io/a8g32>
%
% Create Figures 1 - 7

%% Data / Toolboxes / Scripts used
%
% * CRP Toolbox by Norbert Marwan (<http://tocsy.pik-potsdam.de/CRPtoolbox/>)
% * Fred's toolbox: https://github.com/FredHasselman/toolboxML
% * Data is available in the GithHub repository: <https://github.com/FredHasselman/> and throuigh the OSF project:
% <https:osf.io/a8g32>

%% Author / Version / License
%
% Repository: <https://github.com/FredHasselman/https://github.com/FredHasselman/Hasselman2014-PeerJ-Classifying-Complex-Patterns-/ *Hasselman2014 on GitHub*>
%
% Created by: <http://www.fredhasselman.com/ *Fred Hasselman*> / January 2011
% Affiliations: <http://www.ru.nl/bsi _Behavioural Science Institute_> - <http://www.ru.nl _Radboud University Nijmegen_>

%% PREP
%
% Uncomment next line to clear and close everything... detergent grade!
% omo

%Change ... to the path on your machine where you stored the files
path=pwd;
cd(path)
load('Hasselman2014_stimfeatures_ORI.mat');


%% Figure 4.1: F2 Slope

h0=figure;
maximize(h0);

cm  = flipud(gray);
v   = [-200:10:0];
icm = 1:floor(length(cm)/length(v)):length(cm-20);
cmp = cm(icm,:);
colormap(cmp);
lvl   = .1;
 
for cnt=1:40
  
  % Get the Formant Tracks
  [t,m,n] = unique(Formants(1,cnt).tracks{1,1});
  
  IN = Formants(1,cnt).tracks{1,2};
  F1 = Formants(1,cnt).tracks{1,4};
  F2 = Formants(1,cnt).tracks{1,5};
  F3 = Formants(1,cnt).tracks{1,6};
  
  IN = IN(m); F1 = F1(m); F2 = F2(m); F3 = F3(m);
  
  Sweep(cnt).TI =  t(IN>=lvl);
  Sweep(cnt).F1 = F1(IN>=lvl);
  Sweep(cnt).F2 = F2(IN>=lvl);
  Sweep(cnt).F3 = F3(IN>=lvl);
  
  dsF1 = smooth(Sweep(cnt).F1,.6,'rloess');
  dsF2 = smooth(Sweep(cnt).F2,.6,'rloess');
  dsF3 = smooth(Sweep(cnt).F3,.6,'rloess');
  
  subplot(4,10,cnt);
  
  % Plot the spectrogram
  [~, h_spec] = contourf(STIM(cnt).T,STIM(cnt).F,20*log10(abs(STIM(cnt).S)+eps),v,'LineColor','none');hold on;
  %axis tight
  ax0 = gca;
  
  h_F1 = plot(ax0,Sweep(cnt).TI,dsF1,'-','Color',[.9 .9 .9],'LineWidth',3); hold on;
  h_F2 = plot(ax0,Sweep(cnt).TI,dsF2,'-','Color',[.9 .9 .9],'LineWidth',3); hold on;
  h_F3 = plot(ax0,Sweep(cnt).TI,dsF3,'-','Color',[.9 .9 .9],'LineWidth',3); hold on;
  
  set(ax0,'Ytick',[0:1000:3000],'YtickLabel',{'','1','2',''});
  
  ylim([1 3000]);
  
  grid on;
  
  [mxF2 In] = max(dsF2);
  Sweep(cnt).F2mx = mxF2; Sweep(cnt).tF2mx = Sweep(cnt).TI(In);
  [mnF2 mIn] = min(dsF2);
  Sweep(cnt).F2mn = mnF2; Sweep(cnt).F2mn = Sweep(cnt).TI(mIn); Sweep(cnt).swpF2 = (mxF2-mnF2)/(Sweep(cnt).TI(In)-Sweep(cnt).TI(mIn));
  
  % Note that the stimuli are not of equal length!
  % Obviously "Slowed Down" stimuli are longer
  if ismember(cnt,[1:10])
   title(num2str(cnt));
   xlim([0.01 .6]);
  end
  if ismember(cnt,[11:20])
   title(num2str(cnt-10));
   xlim([0.01 .9]);
  end
  if ismember(cnt,[21:30])
   title(num2str(cnt-20));
   xlim([0.01 .6]);
  end
  if ismember(cnt,[31:40])
   title(num2str(cnt-30));
   xlim([0.01 .9]);
  end
  
  if cnt==1
   Opos(1,:)=get(ax0,'Position');
   ylabel('Frequency (kHz)');
   xlabel('Time (s)');
  end
  if cnt==11
   Opos(2,:)=get(ax0,'Position');
   ylabel('Frequency (kHz)');
   xlabel('Time (s)');
  end
  if cnt==21
   Opos(3,:)=get(ax0,'Position');
   ylabel('Frequency (kHz)');
   xlabel('Time (s)');
  end
  if cnt==31
   Opos(4,:)=get(ax0,'Position');
   ylabel('Frequency (kHz)');
   xlabel('Time (s)');
  end
  
  h_t = text(.2,150,['\DeltaF2 = ',deblank(num2str(Sweep(cnt).swpF2/1000,'%#1.2f'))],'BackgroundColor',[.9 .9 .9],'Margin',0.01); %[.8 .8 .8] ,'EdgeColor','k'
 
  clear IN F1 F2 F3 t m n c h_spec h_s h_t h_F1 h_F2 h_F3 dsF1 dsF2 dsF3 In mIn mnF2 mxF2
  
end

Tpos = [-.1 -0.06 0 0];  

h_1= annotation('textbox',[Opos(1,:)+Tpos],'String','None','EdgeColor','none','FontSize',16);
set(h_1,'FitBoxToText','on');
h_2= annotation('textbox',[Opos(2,:)+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
set(h_2,'FitBoxToText','on');
h_3= annotation('textbox',[Opos(3,:)+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
set(h_3,'FitBoxToText','on');
h_4=annotation('textbox',[Opos(4,:)+Tpos],'String','Both','EdgeColor','none','FontSize',16);
set(h_4,'FitBoxToText','on');
h_5=annotation('textbox',[.146 .065 0 0],'String','/bAk/','EdgeColor','none','FontSize',16);
set(h_5,'FitBoxToText','on');
h_6=annotation('textbox',[.862 .065 0 0],'String','/dAk/','EdgeColor','none','FontSize',16);
set(h_6,'FitBoxToText','on');
h_7 = annotation('arrow',[.2 .84],[.05 .05],'HeadStyle','cback2','LineWidth',1.5);
 
keep Formants HNR RTent SPEC STIM rpSTATS rpTS stimuli rpMTRX

% Uncomment if you want to save a figure
%grab('Hasselman2014_Figure1',0)


%% Figure 4.2: maxENVELOPE Slope
%
% Slope till max formant amplitude from stimulus onset and formant onset

h0=figure;
maximize(h0);

cnt=0;
for cnt=1:40
  
  subplot(4,10,cnt);
  ax1 = gca;

  plot(STIM(cnt).IAT,  (stimuli(cnt).y./5),'LineWidth',.1,'Color',[.7 .7 .7]);hold on;
  axis tight;
  plot(STIM(cnt).IAT,[STIM(cnt).IAsm],'LineWidth',2,'Color',[.3 .3 .3]);hold on;
  axis tight;
  plot(STIM(cnt).IAT,-[STIM(cnt).IAsm],'LineWidth',2,'Color',[.3 .3 .3]);hold on;
  axis tight;
  
  %Plot min to max AMP line
  IAS = plot([STIM(cnt).IAT(1) STIM(cnt).IATmx],[STIM(cnt).IA(1) STIM(cnt).IAmx],'Color','k');hold on
  plot([STIM(cnt).IAT(1) STIM(cnt).IATmx],[STIM(cnt).IA(1) STIM(cnt).IAmx],'o','MarkerSize',4,'MarkerEdgeColor',[.3 .3 .3],'MarkerFaceColor',[.8 .8 .8]);
  
  axis tight;
  set(ax1,'Ytick',[-.5:.1:.5],'YtickLabel',{'','','','','','0','','','','',''});
  
  ylim([-.5 .5]);
  xlim([0 .9]);
  %grid on;
  
  %Print slope in figure
  text(.55,-.45,['\Delta = ',num2str(STIM(cnt).IASmxO,'%1.2f')]);
  
  %Garnish
  if ismember(cnt,[1:10])
   title(num2str(cnt));
  end
  if ismember(cnt,[11:20])
   title(num2str(cnt-10));
  end
  if ismember(cnt,[21:30])
   title(num2str(cnt-20));
  end
  if ismember(cnt,[31:40])
   title(num2str(cnt-30));
  end
  
  if cnt==1
   Opos(1,:)=get(ax1,'Position');
   ylabel('Amplitude (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==11
   Opos(2,:)=get(ax1,'Position');
   ylabel('Amplitude (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==21
   Opos(3,:)=get(ax1,'Position');
   ylabel('Amplitude (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==31
   Opos(4,:)=get(ax1,'Position');
   ylabel('Amplitude (a.u.)');
   xlabel('Time (s)');
  end
   
  clear IAS IASmx0
  
end

Tpos = [-.1 -0.06 0 0];
  
h_1= annotation('textbox',[Opos(1,:)+Tpos],'String','None','EdgeColor','none','FontSize',16);
set(h_1,'FitBoxToText','on');
h_2= annotation('textbox',[Opos(2,:)+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
set(h_2,'FitBoxToText','on');
h_3= annotation('textbox',[Opos(3,:)+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
set(h_3,'FitBoxToText','on');
h_4=annotation('textbox',[Opos(4,:)+Tpos],'String','Both','EdgeColor','none','FontSize',16);
set(h_4,'FitBoxToText','on');
h_5=annotation('textbox',[.146 .065 0 0],'String','/bAk/','EdgeColor','none','FontSize',16);
set(h_5,'FitBoxToText','on');
h_6=annotation('textbox',[.862 .065 0 0],'String','/dAk/','EdgeColor','none','FontSize',16);
set(h_6,'FitBoxToText','on');
h_7 = annotation('arrow',[.2 .84],[.05 .05],'HeadStyle','cback2','LineWidth',1.5);

keep Formants HNR RTent SPEC STIM rpSTATS rpTS stimuli rpMTRX

% Uncomment if you want to save a figure
%grab('Hasselman2014_Figure2',0);

%% Figure 4.3: Rise and Fall Time Entropy

h0=figure;
maximize(h0);

for cnt=1:40
 
  subplot(4,10,cnt);
  ax1 = gca;
  
  %Scale up the derivative of the smoothed envelope
  dsENV = derivative(STIM(cnt).IAsm).*600;
  
  plot(STIM(cnt).IAT, (stimuli(cnt).y./5),'LineWidth',.1,'Color',[.7 .7 .7]);hold on;
  axis tight;
 
  plot(STIM(cnt).IAT , dsENV,'LineWidth',2,'Color',[.3 .3 .3]);hold on;
  axis tight;
    
  % Plot the crossings
  [~,t0,s0] = crossing(dsENV,STIM(cnt).IAT);
  plot(t0,s0+.25,'x','MarkerSize',7,'MarkerEdgeColor','k');
  
  % Print entropy in figure
  text(.05,-.26,['RFTe = ',num2str(RTent(cnt),'%1.2f')]);
  
  
  set(ax1,'Ytick',[-.3:.1:.3],'YtickLabel',{'','','','0','','',''});
  
  ylim([-.3 .3]);
  xlim([0 .9]);
  %grid on;
  
  % Garnish
  if ismember(cnt,[1:10])
   title(num2str(cnt));
  end
  if ismember(cnt,[11:20])
   title(num2str(cnt-10));
  end
  if ismember(cnt,[21:30])
   title(num2str(cnt-20));
  end
  if ismember(cnt,[31:40])
   title(num2str(cnt-30));
  end
  
  if cnt==1
   Opos(1,:)=get(ax1,'Position');
   ylabel('Amplitude change (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==11
   Opos(2,:)=get(ax1,'Position');
   ylabel('Amplitude change (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==21
   Opos(3,:)=get(ax1,'Position');
   ylabel('Amplitude change (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==31
   Opos(4,:)=get(ax1,'Position');
   ylabel('Amplitude change (a.u.)');
   xlabel('Time (s)');
  end
   
  clear dsENV s0 t0 
  
end

Tpos = [-.1 -0.06 0 0];
  
h_1= annotation('textbox',[Opos(1,:)+Tpos],'String','None','EdgeColor','none','FontSize',16);
set(h_1,'FitBoxToText','on');
h_2= annotation('textbox',[Opos(2,:)+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
set(h_2,'FitBoxToText','on');
h_3= annotation('textbox',[Opos(3,:)+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
set(h_3,'FitBoxToText','on');
h_4=annotation('textbox',[Opos(4,:)+Tpos],'String','Both','EdgeColor','none','FontSize',16);
set(h_4,'FitBoxToText','on');
h_5=annotation('textbox',[.146 .065 0 0],'String','/bAk/','EdgeColor','none','FontSize',16);
set(h_5,'FitBoxToText','on');
h_6=annotation('textbox',[.862 .065 0 0],'String','/dAk/','EdgeColor','none','FontSize',16);
set(h_6,'FitBoxToText','on');
h_7 = annotation('arrow',[.2 .84],[.05 .05],'HeadStyle','cback2','LineWidth',1.5);

keep Formants HNR RTent SPEC STIM rpSTATS rpTS stimuli rpMTRX

% Uncomment if you want to save a figure
% grab('Hasselman2014_Figure3',0);


%% Figure 4.4: Phase Space Reconstruction Example

cnt=10; ds=1;
x = rpTS(cnt).ts(1:1024,2);
t = rpTS(cnt).ts(1:1024,1);
tau=6; m=3;

x=minplus1(x);

x1=x(1+(m-3)*tau:end-((m-1)*tau)+1); ts1=t(1+(m-3)*tau:end-((m-1)*tau)+1); %  t1=minplus1(t1);
x2=x(0+(m-2)*tau:end-((m-2)*tau));   ts2=t(0+(m-2)*tau:end-((m-2)*tau));   %  t2=minplus1(t2);
x3=x(0+(m-1)*tau:end-((m-3)*tau));   ts3=t(0+(m-1)*tau:end-((m-3)*tau));   %  t3=minplus1(t3);

t = -1:2/(length(x1)-1):1;

t1=t;
t2=t;
t3=t;

z1=-1+zeros(length(t1),1);
z2= 1+zeros(length(t2),1);
z3= 1+zeros(length(t3),1);

h0=figure;
maximize(h0);

ax_PS =axes('NextPlot','add');

h_ts1=plot3(x1(:),x2(:),z1,'-k');axis square;
xlim([-1 1]),ylim([-1 1]);zlim([-1 1]); view(3);
h_ts2=plot3(z2,x2(:),x3(:),'-k');axis square;
h_ts3=plot3(x1(:),z3,x3(:),'-k');axis square;

set(ax_PS,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
set([h_ts1 h_ts2 h_ts3],'Color',[.7 .7 .7]);

h_ps=plot3(x1(:), x2(:), x3(:),'-ko'); grid on; axis square;
set(h_ps,'MarkerFaceColor',[.5 .5 .5]);
xlabel('X_m_=_1'); ylabel('X_m_=_2'); zlabel('X_m_=_3');
title({'Reconstructed Phase Space of Stimulus 10 (first 1024 samples)', ['Delay Embedding with m = 3, \tau = 6, \epsilon = ',num2str(rpSTATS(cnt,1),1),' of maximum norm distance']});

s=.18;
point=525;
xc=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*s;xc=xc+x1(point)-(s/2);
yc=[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*s;yc=yc+x2(point)-(s/2);
zc=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*s;zc=zc+x3(point)-(s/2);

axes(ax_PS);
for i=1:6
 h=patch(xc(:,i),yc(:,i),zc(:,i),[.8 .8 .8]);
 set(h,'edgecolor','k','FaceAlpha',.5,'LineWidth',1)
end

xc1 = xc(:,1); yc1 = yc(:,1); zc1 = zc(:,1);
xc2 = xc(:,2); yc2 = yc(:,2); zc2 = zc(:,2);
xc3 = xc(:,3); yc3 = yc(:,3); zc3 = zc(:,3);
xc4 = xc(:,4); yc4 = yc(:,4); zc4 = zc(:,4);
xc5 = xc(:,5); yc5 = yc(:,5); zc5 = zc(:,5);
xc6 = xc(:,6); yc6 = yc(:,6); zc6 = zc(:,6);

axes(ax_PS);
h1=patch(xc5(:),yc5(:),[-1 -1 -1 -1]',[.8 .8 .8]);
h2=patch(xc6(:),[1 1 1 1]',zc1(:),[.8 .8 .8]);
h3=patch([1 1 1 1]',yc2(:),zc1(:),[.8 .8 .8]);

set([h1 h2 h3],'edgecolor','k','FaceAlpha',.5,'LineWidth',1)

text(xc(1,1)-(s/3),yc(1,1)+(s/3),zc(1,1)+.15,'\epsilon','FontSize',20);

axes(ax_PS);
hc1 = plot3(x1(point),x2(point),-1,'ok');
hc2 = plot3(1,x2(point),x3(point),'ok');
hc3 = plot3(x1(point),1,x3(point),'ok');
set([hc1 hc2 hc3],'MarkerFaceColor',[.2 .2 .2],'MarkerSize',8);

point3=795;
hc7 = plot3(x1(point3),x2(point3),-1,'sk');
hc8 = plot3(1,x2(point3),x3(point3),'sk');
hc9 = plot3(x1(point3),1,x3(point3),'sk');
set([hc7 hc8 hc9],'MarkerFaceColor',[.7 .7 .7],'MarkerSize',8);

keep Formants HNR RTent SPEC STIM rpSTATS rpTS stimuli rpMTRX

% Uncomment if you want to save a figure
% grab('Hasselman2014_Figure4',0);

%% Figure 4.5 (LEFT): RP example

cnt=10; tau=6; m=3; e=.01; thr= 'rr';
% Uncomment to calculate rahter than load the RP matrix
%[MTRX] = crp(rpTS(cnt).ts,m,tau,e,thr,'nonormalize','silent');

[MTRX] = rpMTRX(cnt).rp;
[r c]  = size(MTRX);

scrsz = get(0,'ScreenSize');
h0    = figure('Position',[scrsz(3)/4 0 scrsz(3)/2 scrsz(3)/2],'NextPlot','add');

RPdims  =[.35 .35 .6 .6];

TSHdims1=[.35 .06 .6 .06];
TSHdims2=[.35 .15 .6 .06];
TSHdims3=[.35 .24 .6 .06];

TSVdims1=[.06 .35 .06 .6];
TSVdims2=[.15 .35 .06 .6];
TSVdims3=[.24 .35 .06 .6];

%Recurrence Plot
ax_RP =axes('Position',RPdims);
spy(MTRX,'.k',1);
axis square; axis xy;
xlabel('Recurrent values in m-dimensional phase space');
title({'(Auto) Recurrence Plot of Stimulus 10:',['\tau = 6, m = 3, \epsilon = ',num2str(rpSTATS(cnt,1),1)]});
loi=line([0 r],[0 c],'Color','k','LineWidth',3);

x1=rpTS(cnt).ts(1+(m-3)*tau:end-((m-1)*tau),2); y1=rpTS(cnt).ts(1+(m-3)*tau:end-((m-1)*tau),1);
x2=rpTS(cnt).ts(0+(m-2)*tau:end-((m-2)*tau),2); y2=rpTS(cnt).ts(0+(m-2)*tau:end-((m-2)*tau),1);
x3=rpTS(cnt).ts(0+(m-1)*tau:end-((m-3)*tau),2); y3=rpTS(cnt).ts(0+(m-1)*tau:end-((m-3)*tau),1);

%Horizontal TS
ax_TSH1=axes('Position',TSHdims1);%drawnow
h_TSH1=line(y1,x1); axis tight
xlabel('Surrogate Dimensions: m Time Series Delayed by m*\tau');
ylabel('X_m_=_1');

ax_TSH2=axes('Position',TSHdims2);%drawnow
h_TSH2=line(y2,x2); axis tight
ylabel('X_m_=_2');

ax_TSH3=axes('Position',TSHdims3);%drawnow
h_TSH3=line(y3,x3); axis tight
ylabel('X_m_=_3');

%Vertical TS
ax_TSV1=axes('Position',TSVdims1);%drawnow
h_TSV1=line(x1,y1); axis tight
ylabel('Surrogate Dimensions: m Time Series Delayed by m*\tau');
xlabel('X_m_=_1');

ax_TSV2=axes('Position',TSVdims2);%drawnow
h_TSV2=line(x2,y2); axis tight
xlabel('X_m_=_2');

ax_TSV3=axes('Position',TSVdims3);%drawnow
h_TSV3=line(x3,y3); axis tight
xlabel('X_m_=_3');

set([h_TSH1,h_TSH2,h_TSH3,h_TSV1,h_TSV2,h_TSV3],'Color',[.5 .5 .5]);

set([ax_RP],...
 'FontSize',14,...
 'XTick',[1 round(r/2) r],...
 'YTick',[1 round(c/2) c],...
 'XTickLabel',{'1' num2str(round(r/2)) num2str(r)},...
 'YTickLabel',{'1' num2str(round(r/2)) num2str(r)});

set([ax_TSH1],...
 'FontSize',14,...
 'XTick',[y1(1) y1(round(length(y1)/2)) y1(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y1(1),4) num2str(y1(round(length(y1)/2)),4) num2str(y1(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSH2],...
 'FontSize',14,...
 'XTick',[y2(1) y2(round(length(y2)/2)) y2(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y2(1),4) num2str(y2(round(length(y2)/2)),4) num2str(y2(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSH3],...
 'FontSize',14,...
 'XTick',[y3(1) y3(round(length(y3)/2)) y3(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y3(1),4) num2str(y3(round(length(y1)/2)),4) num2str(y3(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSV1],...
 'FontSize',14,...
 'YTick',[y1(1) y1(round(length(y1)/2)) y1(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' '0' '1'},'Box','on');

set([ax_TSV2],...
 'FontSize',14,...
 'YTick',[y2(1) y2(round(length(y2)/2)) y2(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' 0'' '1'},'Box','on');

set([ax_TSV3],...
 'FontSize',14,...
 'YTick',[y3(1) y3(round(length(y3)/2)) y3(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' '0' '1'},'Box','on');

rpST = {['RQA measures:'],...
 [' '],...
 ['REC = ',num2str(rpSTATS(cnt,2),2)],...
 ['DET = ',num2str(rpSTATS(cnt,3),2)],...
 ['Lmn = ',num2str(rpSTATS(cnt,4),3)],...
 ['ENT = ',num2str(rpSTATS(cnt,6),3)],...
 ['LAM = ',num2str(rpSTATS(cnt,7),2)],...
 ['Vmn = ',num2str(rpSTATS(cnt,8),3)]};

h_s = annotation('textbox',[.06 .3 0 0],'String',rpST,'EdgeColor','none','FontName','Courier','FontSize',24);
set(h_s,'FitBoxToText','on');

keep Formants HNR RTent SPEC STIM rpSTATS rpTS stimuli rpMTRX

% Uncomment if you want to save a figure
% grab('Hasselman2014_Figure5',0);

%% Figure 4.5 (RIGHT): RP RANDOM example

cnt=10; tau=6; m=3; e=.01; thr= 'rr';

% Create a shuffled version
Y = shuffle(rpTS(cnt).ts(:,2));
[STATS]= crqa(Y, m,tau,e,thr,'nonormalize','silent');
[MTRX] = crp([rpTS(cnt).ts(:,1) Y],m,tau,e,thr,'nonormalize','silent');

[r c] = size(MTRX);

scrsz = get(0,'ScreenSize');
h0 = figure('Position',[scrsz(3)/4 0 scrsz(3)/2 scrsz(3)/2],'NextPlot','add');

RPdims  =[.35 .35 .6 .6];

TSHdims1=[.35 .06 .6 .06];
TSHdims2=[.35 .15 .6 .06];
TSHdims3=[.35 .24 .6 .06];

TSVdims1=[.06 .35 .06 .6];
TSVdims2=[.15 .35 .06 .6];
TSVdims3=[.24 .35 .06 .6];

%Recurrence Plot
ax_RP =axes('Position',RPdims);
spy(MTRX,'.k',1);
axis square; axis xy;
xlabel('Recurrent values in m-dimensional phase space');
title({'(Auto) Recurrence Plot of Stimulus 10 (RANDOMISED):',['\tau = 6, m = 3, \epsilon = ',num2str(STATS(1),1)]});
loi=line([0 r],[0 c],'Color','k','LineWidth',3);

x1=Y(1+(m-3)*tau:end-((m-1)*tau)); y1=rpTS(cnt).ts(1+(m-3)*tau:end-((m-1)*tau),1);
x2=Y(0+(m-2)*tau:end-((m-2)*tau)); y2=rpTS(cnt).ts(0+(m-2)*tau:end-((m-2)*tau),1);
x3=Y(0+(m-1)*tau:end-((m-3)*tau)); y3=rpTS(cnt).ts(0+(m-1)*tau:end-((m-3)*tau),1);

%Horizontal TS
ax_TSH1=axes('Position',TSHdims1);%drawnow
h_TSH1=line(y1,x1); axis tight
xlabel('Surrogate Dimensions: m Time Series Delayed by m*\tau');
ylabel('X_m_=_1');

ax_TSH2=axes('Position',TSHdims2);%drawnow
h_TSH2=line(y2,x2); axis tight
ylabel('X_m_=_2');

ax_TSH3=axes('Position',TSHdims3);%drawnow
h_TSH3=line(y3,x3); axis tight
ylabel('X_m_=_3');

%Vertical TS
ax_TSV1=axes('Position',TSVdims1);%drawnow
h_TSV1=line(x1,y1); axis tight
ylabel('Surrogate Dimensions: m Time Series Delayed by m*\tau');
xlabel('X_m_=_1');

ax_TSV2=axes('Position',TSVdims2);%drawnow
h_TSV2=line(x2,y2); axis tight
xlabel('X_m_=_2');

ax_TSV3=axes('Position',TSVdims3);%drawnow
h_TSV3=line(x3,y3); axis tight
xlabel('X_m_=_3');

set([h_TSH1,h_TSH2,h_TSH3,h_TSV1,h_TSV2,h_TSV3],'Color',[.5 .5 .5]);

set([ax_RP],...
 'FontSize',14,...
 'XTick',[1 round(r/2) r],...
 'YTick',[1 round(c/2) c],...
 'XTickLabel',{'1' num2str(round(r/2)) num2str(r)},...
 'YTickLabel',{'1' num2str(round(r/2)) num2str(r)});

set([ax_TSH1],...
 'FontSize',14,...
 'XTick',[y1(1) y1(round(length(y1)/2)) y1(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y1(1),4) num2str(y1(round(length(y1)/2)),4) num2str(y1(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSH2],...
 'FontSize',14,...
 'XTick',[y2(1) y2(round(length(y2)/2)) y2(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y2(1),4) num2str(y2(round(length(y2)/2)),4) num2str(y2(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSH3],...
 'FontSize',14,...
 'XTick',[y3(1) y3(round(length(y3)/2)) y3(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y3(1),4) num2str(y3(round(length(y1)/2)),4) num2str(y3(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSV1],...
 'FontSize',14,...
 'YTick',[y1(1) y1(round(length(y1)/2)) y1(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' '0' '1'},'Box','on');

set([ax_TSV2],...
 'FontSize',14,...
 'YTick',[y2(1) y2(round(length(y2)/2)) y2(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' 0'' '1'},'Box','on');

set([ax_TSV3],...
 'FontSize',14,...
 'YTick',[y3(1) y3(round(length(y3)/2)) y3(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' '0' '1'},'Box','on');

rpST = {['RQA measures:'],...
 [' '],...
 ['REC = ',num2str(STATS(2),2)],...
 ['DET = ',num2str(STATS(3),1)],...
 ['Lmn = ',num2str(STATS(4),3)],...
 ['ENT = ',num2str(STATS(6),1)],...
 ['LAM = ',num2str(STATS(7),1)],...
 ['Vmn = ',num2str(STATS(8),3)]};

h_s = annotation('textbox',[.06 .30 0 0],'String',rpST,'EdgeColor','none','FontName','Courier','FontSize',24);
set(h_s,'FitBoxToText','on');

keep Formants HNR RTent SPEC STIM rpSTATS rpTS stimuli rpMTRX

% Uncomment if you want to save a figure
% grab('Hasselman2014_Figure5r',0);


%% Figure 4.6: RP plots

ds = 2;
h0=figure;
maximize(h0);

for cnt=1:40
 
 subplot(4,10,cnt)
 tssz=length(rpTS(cnt).ts(:,1));
 
 rr = downsample(rpMTRX(cnt).rp,ds); cc = downsample(transpose(rr),ds);
 [r c] = size(cc);
 spy(cc,'.k',1);
 ax0 = gca;
 axis square; axis xy;
 xlabel(''); ylabel('');title(num2str(cnt));
 xlim([0 r(end)]);ylim([0 c(end)]);
 
 set(ax0,'XTick',[0 (round(r(end)/2)) r(end)],...
  'YTick',[0 (round(c(end)/2)) c(end)],...
  'XTickLabel',{'','',''},...
  'YTickLabel',{'','',''});
 
 if ismember(cnt,[1:10])
  title([num2str(cnt),'- \epsilon =',num2str(rpSTATS(cnt,1),1)]);
 end
 if ismember(cnt,[11:20])
  title([num2str(cnt-10),'- \epsilon =',num2str(rpSTATS(cnt,1),1)]);
 end
 if ismember(cnt,[21:30])
  title([num2str(cnt-20),'- \epsilon =',num2str(rpSTATS(cnt,1),1)]);
 end
 if ismember(cnt,[31:40])
  title([num2str(cnt-30),'- \epsilon =',num2str(rpSTATS(cnt,1),1)]);
 end
   
  if cnt==1
   Opos(1,:)=get(ax0,'Position');
  end
  if cnt==11
   Opos(2,:)=get(ax0,'Position');
  end
  if cnt==21
   Opos(3,:)=get(ax0,'Position');
  end
  if cnt==31
   Opos(4,:)=get(ax0,'Position');
  end
 
% Plot TS
 TSpos = get(ax0,'Position');
 ax_TS = axes('Position',[TSpos(1),TSpos(2)-.03,TSpos(3),TSpos(4)/4]);
 h_TSH = line(rpTS(cnt).ts(1:end,1),rpTS(cnt).ts(1:end,2),'Color',[.5 .5 .5]); axis tight
 set(ax_TS,'Visible','off');
 
end


Tpos = [-.1 -0.06 0 0];
  
h_1= annotation('textbox',[Opos(1,:)+Tpos],'String','None','EdgeColor','none','FontSize',16);
set(h_1,'FitBoxToText','on');
h_2= annotation('textbox',[Opos(2,:)+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
set(h_2,'FitBoxToText','on');
h_3= annotation('textbox',[Opos(3,:)+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
set(h_3,'FitBoxToText','on');
h_4=annotation('textbox',[Opos(4,:)+Tpos],'String','Both','EdgeColor','none','FontSize',16);
set(h_4,'FitBoxToText','on');
h_5=annotation('textbox',[.146 .065 0 0],'String','/bAk/','EdgeColor','none','FontSize',16);
set(h_5,'FitBoxToText','on');
h_6=annotation('textbox',[.862 .065 0 0],'String','/dAk/','EdgeColor','none','FontSize',16);
set(h_6,'FitBoxToText','on');
h_7 = annotation('arrow',[.2 .84],[.05 .05],'HeadStyle','cback2','LineWidth',1.5);
h_8 = annotation('textbox',[.39 .99 0 0],'String','(Auto) Recurrence Plots for all Stimuli','EdgeColor','none','FontSize',16);
set(h_8,'FitBoxToText','on'); 

keep Formants HNR RTent SPEC STIM rpSTATS rpTS stimuli rpMTRX

% Uncomment if you want to save a figure
% grab('Hasselman2014_Figure6',0);

%% Figure 4.7: LOGIT predictions

load('predictedLOGIT.mat');

f0=figure;
maximize(f0);
subplot(2,2,1)
ax0=gca; pos=get(ax0,'Position');

none  = errorbar([1:10]-0.1,data(:,1),[data(:,1)-data(:,2)],[data(:,3)-data(:,1)],'-');
hold on;

noned = errorbar([1:10]+0.1,data2(:,1),[data2(:,1)-data2(:,2)],[data2(:,3)-data2(:,1)],'-.');
hold on;
xlim([0.7 10.3]);
ylim([-0.05 1.05]);
threshold = line(xlim,[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--');

set([none, noned], ...
 'Marker'          , 'o'        , ...
 'MarkerSize'      , 7          , ...
 'MarkerEdgeColor' , 'k'     , ...
 'MarkerFaceColor' , [.5 .5 .5] , ...
 'Color'       , 'k');

set(none, ...
 'MarkerFaceColor' , [0 0 0]);

set(gca, ...
 'TickDir'     , 'out'     , ...
 'TickLength'  , [.01 .01] , ...
 'YTick'       , 0:.1:1, ...
 'XTickLabel'  , {'/bAk/','2','3','4','5','6','7','8','9','/dAk/'},...
 'Box','off', ...
 'LineWidth'   , 1         );

legend('Average reader','Dyslexic reader','Location','SouthEast');
title('A. None');
xlabel('Stimulus');
ylabel('\fontsize{14}\it\pi\fontsize{12}\rm with \itCI\fontsize{10}_{.95}\fontsize{10}\rm  for Perceiving /dAk/');


subplot(2,2,2)
none=errorbar([1:10]-0.1,data(:,4),[data(:,4)-data(:,5)],[data(:,6)-data(:,4)],'-');
hold on;
noned=errorbar([1:10]+0.1,data2(:,4),[data2(:,4)-data2(:,5)],[data2(:,6)-data2(:,4)],'-.');
hold on;
xlim([0.7 10.3]);
ylim([-0.05 1.05]);
threshold = line(xlim,[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--');

set([none, noned], ...
 'Marker'          , 'v'        , ...
 'MarkerSize'      , 7          , ...
 'MarkerEdgeColor' , 'k'     , ...
 'MarkerFaceColor' , [.5 .5 .5] , ...
 'Color'       , 'k');

set(none, ...
 'MarkerFaceColor' , [0 0 0]);

set(gca, ...
 'TickDir'     , 'out'     , ...
 'TickLength'  , [.01 .01] , ...
 'YTick'       , 0:.1:1, ...
 'XTickLabel'  , {'/bAk/','2','3','4','5','6','7','8','9','/dAk/'},...
 'Box','off', ...
 'LineWidth'   , 1         );


legend('Average reader','Dyslexic reader','Location','SouthEast');
title('B. Slowed down');
xlabel('Stimulus');
ylabel('\fontsize{14}\it\pi\fontsize{12}\rm with \itCI\fontsize{10}_{.95}\fontsize{10}\rm  for Perceiving /dAk/');


subplot(2,2,3)
none=errorbar([1:10]-0.1,data(:,7),[data(:,7)-data(:,8)],[data(:,9)-data(:,7)],'-');
hold on;
noned=errorbar([1:10]+0.1,data2(:,7),[data2(:,7)-data2(:,8)],[data2(:,9)-data2(:,7)],'-.');
hold on;
xlim([0.7 10.3]);
ylim([-0.05 1.05]);
threshold = line(xlim,[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--');

set([none, noned], ...
 'Marker'          , '^'        , ...
 'MarkerSize'      , 7          , ...
 'MarkerEdgeColor' , 'k'     , ...
 'MarkerFaceColor' , [.5 .5 .5] , ...
 'Color'       , 'k');

set(none, ...
 'MarkerFaceColor' , [0 0 0]);

set(gca, ...
 'TickDir'     , 'out'     , ...
 'TickLength'  , [.01 .01] , ...
 'YTick'       , 0:.1:1, ...
 'XTickLabel'  , {'/bAk/','2','3','4','5','6','7','8','9','/dAk/'},...
 'Box','off', ...
 'LineWidth'   , 1         );

legend('Average reader','Dyslexic reader','Location','SouthEast');
title('C. Amplified');
xlabel('Stimulus');
ylabel('\fontsize{14}\it\pi\fontsize{12}\rm with \itCI\fontsize{10}_{.95}\fontsize{10}\rm  for Perceiving /dAk/');


subplot(2,2,4)
none=errorbar([1:10]-0.1,data(:,10),[data(:,10)-data(:,11)],[data(:,12)-data(:,10)],'-');
hold on;
noned=errorbar([1:10]+0.1,data2(:,10),[data2(:,10)-data2(:,11)],[data2(:,12)-data2(:,10)],'-.');
hold on;
xlim([0.7 10.3]);
ylim([-0.05 1.05]);
threshold = line(xlim,[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--');

set([none, noned], ...
 'Marker'          , 'd'        , ...
 'MarkerSize'      , 7          , ...
 'MarkerEdgeColor' , 'k'     , ...
 'MarkerFaceColor' , [.5 .5 .5] , ...
 'Color'       , 'k');

set(none, ...
 'MarkerFaceColor' , [0 0 0]);

set(gca, ...
 'TickDir'     , 'out'     , ...
 'TickLength'  , [.01 .01] , ...
 'YTick'       , 0:.1:1, ...
 'XTickLabel'  , {'/bAk/','2','3','4','5','6','7','8','9','/dAk/'},...
 'Box','off', ...
 'LineWidth'   , 1         );

legend('Average reader','Dyslexic reader','Location','SouthEast');
title('D. Both');
xlabel('Stimulus');
ylabel('\fontsize{14}\it\pi\fontsize{12}\rm with \itCI\fontsize{10}_{.95}\fontsize{10}\rm  for Perceiving /dAk/');


keep Formants HNR RTent SPEC STIM rpSTATS rpTS stimuli rpMTRX data data2

% Uncomment if you want to save a figure
 grab('Hasselman2014_Figure7',0);
