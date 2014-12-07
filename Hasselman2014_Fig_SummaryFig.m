%% Supplementary Material to Hasselman (2014) "Classifying Complex Patterns into Speech Categories"
%
%%% Introduction
% This is a demonstration script accompanying the article "Classifying Complex Patterns into Speech Categories". Its purpose is to
% provide an example of how to use various freely available MATLAB sources on the web to extract variables from speech
% stimuli. This script should thus be used as an example to build your own scripts only. It is not a function or toolbox, not
% optimized for speed or functionality! Evaluate the code per Cell (using MATLAB's Cell Mode) and inspect the workspace to
% see what is going on.
%
% OSF project contains links to all the files: <https:osf.io/a8g32>

%% Data / Toolboxes / Scripts used
%
% * CRP Toolbox by Norbert Marwan (<http://tocsy.pik-potsdam.de/CRPtoolbox/>)
% * Fred's toolbox: https://github.com/FredHasselman/toolboxML
% *
% * Scripts and Wave File Stimuli are available in the GithHub repository: <https://github.com/FredHasselman/Hasselman2014-PeerJ-Classifying-Complex-Patterns>
% * Data Files are available at this OSF project page: <https:osf.io/a8g32>

%% Author / Version / License
%
% OSF Project Page: <https:osf.io/a8g32>
% Repository: <https://github.com/FredHasselman/Hasselman2014-PeerJ-Classifying-Complex-Patterns>// *Hasselman2014 on GitHub*>
%
% Created by: <http://www.fredhasselman.com/ *Fred Hasselman*> / 2011-2014
% Affiliations: <http://www.ru.nl/bsi _Behavioural Science Institute_> - <http://www.ru.nl _Radboud University Nijmegen_>

%% Create Summary Figure of Stimulus 1

%%  PREP
%
% Uncomment next line to clear and close everything... detergent grade!
% omo

% Change ... to the path on your machine where you stored the files
% If you copied the dropbox folder from the OSF, current path should be Ok (on machines that allow tilde expansion)
path='~/Dropbox/Hasselman2014-PeerJ-Classifying_Acoustic_Signals/';
% This datafile contains the values that were used in the article
load([path,'Hasselman2014_stimfeatures_ORI.mat']);

%% Figure 8 - Summary
h0=figure;
maximize(h0);

cm  = (gray);
v   = [-200:10:0];
icm = 1:floor(length(cm)/length(v)):length(cm-20);
cmp = cm(icm,:);
colormap(cmp);
lvl   = .1;

subplot(4,6,[2 11])
plot(stimuli(1).y,'-k')
axis off

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

subplot(4,6,13);

% Plot the spectrogram
[~, h_spec] = contourf(STIM(cnt).T,STIM(cnt).F,20*log10(abs(STIM(cnt).S)+eps),v,'LineColor','none');hold on;
axf = gca;

h_F1 = plot(axf,Sweep(cnt).TI,dsF1,'-','Color',[.5 .5 .5],'LineWidth',3); hold on;
h_F2 = plot(axf,Sweep(cnt).TI,dsF2,'-','Color',[.5 .5 .5],'LineWidth',3); hold on;
h_F3 = plot(axf,Sweep(cnt).TI,dsF3,'-','Color',[.5 .5 .5],'LineWidth',3); hold on;
 
set(axf,'Ylim',[1 3000],'Ytick',[],'YtickLabel','','Xlim',[0.1 .6],'Xtick',[],'XtickLabel','');

[mxF2 In] = max(dsF2);
Sweep(cnt).F2mx = mxF2; Sweep(cnt).tF2mx = Sweep(cnt).TI(In);
[mnF2 mIn] = min(dsF2);
Sweep(cnt).F2mn = mnF2; Sweep(cnt).F2mn = Sweep(cnt).TI(mIn); Sweep(cnt).swpF2 = (mxF2-mnF2)/(Sweep(cnt).TI(In)-Sweep(cnt).TI(mIn));

 ylabel('Frequency (kHz)');
 xlabel('Time (s)');

h_t = text(.45,180,['\Delta F2 = ',deblank(num2str(Sweep(cnt).swpF2/1000,'%#1.2f'))],'Margin',0.01); %[.8 .8 .8] ,'EdgeColor','k'

title('F2 Slope')

clear IN F1 F2 F3 t m n c h_spec h_s h_t h_F1 h_F2 h_F3 dsF1 dsF2 dsF3 In mIn mnF2 mxF2

subplot(4,6,14);
axh=gca;

[hnrp(cnt).S hnrp(cnt).F hnrp(cnt).T, hnrp(cnt).P] = spectrogram(rpTS(cnt).ts(:,2),SPEC.window,SPEC.noverlap,SPEC.f,SPEC.fs);
 [~, h_spec] = contourf(hnrp(cnt).T,hnrp(cnt).F,20*log10(abs(hnrp(cnt).S)+eps),v,'LineColor','none');hold on;

set(axh,'Ytick',[],'YtickLabel','','Xtick',[],'XtickLabel','');
hnr_t = text(.065,200,['HNR = ',deblank(num2str(HNR(cnt).HNR,'%#1.2f'))],'Margin',0.01); %[.8 .8 .8] ,'EdgeColor','k'

if cnt==1
 Opos(1,:)=get(axh,'Position');
 ylabel('Frequency (kHz)');
 xlabel('Time (s)');
end

title('HNR')


% Figure 2: maxENVELOPE Slope
%
% Slope till max formant amplitude from stimulus onset and formant onset

subplot(4,6,15);
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

set(ax1,'Ylim',[-5 .5],'Ytick',[],'YtickLabel','','Xlim',[0 .6],'Xtick',[],'XtickLabel','');

%Print slope in figure
text(.35,-.42,['\Delta maxENV = ',num2str(STIM(cnt).IASmxO,'%1.2f')]);


if cnt==1
 Opos(1,:)=get(ax1,'Position');
 ylabel('Amplitude (a.u.)');
 xlabel('Time (s)');
end

title('maxENV')
axis tight;

clear IAS IASmx0


% RFTe

subplot(4,6,16);
axr = gca;

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
text(.4,-.26,['RFTe = ',num2str(RTent(cnt),'%1.2f')]);

set(axr,'Ylim',[-.3 ,.3],'Ytick',[],'YtickLabel','','Xlim',[0 .6],'Xtick',[],'XtickLabel','');

ylabel('Amplitude change (a.u.)');
xlabel('Time (s)');

title('RFTe')

% Figure 6: RP plots

ds = 2;

subplot(4,6,17)
tssz=length(rpTS(cnt).ts(:,1));

rr = downsample(rpMTRX(cnt).rp,ds); cc = downsample(transpose(rr),ds);
[r c] = size(cc);
spy(cc,'.k',1);
ax6 = gca;
axis square; axis xy;title('DET / LAM');
xlabel(''); ylabel('');
xlim([0 r(end)]);ylim([0 c(end)]);

set(ax6,'XTick',[],...
 'YTick',[],...
 'ZTick',[],...
 'XTickLabel','',...
 'YTickLabel','',...
 'ZTickLabel','');

% Plot TS
TSpos = get(ax6,'Position');
ax_TS = axes('Position',[TSpos(1)+.01,TSpos(2)-.035,TSpos(3)-.02,TSpos(4)/6]);
h_TSH = line(rpTS(cnt).ts(1:end,1),rpTS(cnt).ts(1:end,2),'Color',[.5 .5 .5]); axis tight
set(ax_TS,'Visible','off');

ax_TSV = axes('Position',[TSpos(1)-.01,TSpos(2),TSpos(4)/9,TSpos(3)+.055]);
h_TSV = line(rpTS(cnt).ts(1:end,2),rpTS(cnt).ts(1:end,1),'Color',[.5 .5 .5]); axis tight
set(ax_TSV,'Visible','off');


% Figure 8 Multifractal Detrended Fluctuation Analysis

qmin=-10;
qmax=10;
qres=101;

qq = linspace(qmin,qmax,qres);

scmin=6;
scmax=12;
ressc=40;

scale=round(2.^[scmin:((scmax-scmin)/ressc):scmax]);

left = [1 4 7 10];
st   = [0 10 20 30];
stc  = [10 20 30 40 50 60 70 80 90 100];

cm  = gray(130);
i = 1;s=1;
subplot(4,6,18);
ax0 = gca;
ax2 = gca;

h(s)=plot(mf(cnt).hq,mf(cnt).Dq,'Color',cm(stc(s),:),'LineWidth',2); hold on;
title('Multifractal Spectrum');
xlabel('h(q)')
ylabel('D(q)')
ylim([-.05 1.05]);
xlim([.45 2.55]);

qzero=mf(i).Hq(qq==0);
qp2=mf(i).Hq(qq==2);
qm2=mf(i).Hq(qq==-2);
qp5=mf(i).Hq(qq==5);
qm5=mf(i).Hq(qq==-5);
plot(ax2,[qzero qzero],[-.05 2.55],':k');
plot(ax2,[qp2 qp2],[-.05 2.55],'--k');
plot(ax2,[qm2 qm2],[-.05 2.55],'--k');

set(ax2,'XTick',[qp2 qzero qm2],'XTickLabel',{'q=2','q=0','q=-2'});

title('CVhq+ / CVhq-')

subplot(4,6,[19 20]);
plot([0 0 1 1],[0 1 0 1],'.w');
text(0.5,0.5,'RTPDH','FontSize',30,'HorizontalAlignment','center')
axis off

subplot(4,6,[21 22]);
plot([0 0 1 1],[0 1 0 1],'.w');
text(0.5,0.5,'ATPDH','FontSize',30,'HorizontalAlignment','center')
axis off


subplot(4,6,[23 24]);
plot([0 0 1 1],[0 1 0 1],'.w');
text(0.5,0.5,'CMH','FontSize',30,'HorizontalAlignment','center')
axis off


h_4=annotation('textbox',[.52 .08 0 0],'String','Causal Ontology','EdgeColor','none','FontSize',18,'HorizontalAlignment','center');
set(h_4,'FitBoxToText','on');

h_5=annotation('textbox',[.08 .065 0 0],'String','Component Dominant','EdgeColor','none','FontSize',18);
set(h_5,'FitBoxToText','on');
h_6=annotation('textbox',[.862 .065 0 0],'String','Interaction Dominant','EdgeColor','none','FontSize',18);
set(h_6,'FitBoxToText','on');
h_7 = annotation('doublearrow',[.2 .84],[.05 .05],'HeadStyle','cback2','Head2Style','cback2','LineWidth',1.5);



grab('Hasselman2014_Figure8',0)