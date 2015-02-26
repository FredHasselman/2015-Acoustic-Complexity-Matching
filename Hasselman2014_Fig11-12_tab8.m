%% Supplementary Material to Hasselman (2014) "Classifying Acoustic Signals into Speech Categories"
%
%%%%%%%%%%%%%% MARKDOWN CODE %%%%%%%%%%%%%%
%
% ### Introduction
% This code performs the Quadratic Discriminant Analysis on extracted features of the speech signal.
% Use it to create Figure 11 and 12-13 and the data in Table 8
%
% This is not a proper function or toolbox, it is not optimized for speed or functionality or aestetics!
% Evaluate the code per Cell (using MATLAB's Cell Mode) and inspect the workspace to see what is going on.
%
% It is possible your monitor size will influence figure legibility.
%
% OSF project contains links to all the files: https://osf.io/a8g32
%
% ----------
%
% ### Data / Toolboxes / Scripts, etc. that need to be on the MATLAB PATH
%
% * [QDA code](http://www.maxlittle.net/home/index.php)  by Max Little used in: [M.A. Little, P.E. McSharry, S.J. Roberts, D.A.E. Costello, I.M. Moroz (2007). Exploiting nonlinear recurrence and fractal scaling properties for voice disorder detection, *BioMedical Engineering OnLine* 6:23](http://www.maxlittle.net/publications/bmeo.pdf)
% * [Fred's toolbox](https://github.com/FredHasselman/toolboxML) on GitHub
% * Scripts are available in a [GithHub repository](https://github.com/FredHasselman/Acoustic-Complexity-Matching)
% * Data Files are available at the [Open Science Framework](https://osf.io/a8g32/files) or directly from [Dropbox](https://www.dropbox.com/sh/i1vp0nlsj3mi3v7/AAD6V5WaQLdmEBTFRIzAqQm8a?dl=0)
%
%
% ### Author / Version / License
%
% Created by: [Fred Hasselman 2011-2014](http://www.fredhasselman.com)
% Affiliations: [School of Pedagogical and Educational Science](http://www.ru.nl/pwo) and [Behavioural Science Institute (Learning & Plasticity)](http://www.ru.nl/bsi) at the [Radboud University Nijmegen, the Netherlands](http://www.ru.nl)
%
%%%%%%%%%%%%%% MARKDOWN CODE %%%%%%%%%%%%%%

%% PREP

% Uncomment next line to clear and close everything... detergent grade!
% omo

% Method to use for Discriminant analysis
method = 'quadratic';

% Number of bootstrap repetitions for the classification
bootN = 1500;
% reproduce results
rng(12345)

% Change ... to the path on your machine where you stored the files
source='~/Dropbox/Hasselman2014-PeerJ-Classifying_Acoustic_Signals/';
datPath = [source,'DATA/'];

%% Load data used in the article
% If you extracted your own data the file will be saved as: Hasselman2014_stimfeatures.mat
load([datPath,'Hasselman2014_stimfeatures_ORI.mat']);
% This is the MLWiN prediction output for M4
load([datPath,'Hasselman2014_predictedLOGITM3M4.mat']);

%% Create the target labels
targetLabels = zeros(40,3);

% If lower CI.95> .5 the label changed from /bAk/ to /dAk/

% M3 Sample
N=observedLabels.M3(:,2)>.5;
S=observedLabels.M3(:,5)>.5;
A=observedLabels.M3(:,8)>.5;
B=observedLabels.M3(:,11)>.5;

targetLabels(:,1) = [N; S; A; B];

% M4 Average readers
N=observedLabels.M4.Average(:,2)>.5;
S=observedLabels.M4.Average(:,5)>.5;
A=observedLabels.M4.Average(:,8)>.5;
B=observedLabels.M4.Average(:,11)>.5;

targetLabels(:,2) = [N; S; A; B];

% M4 Dyslexic readers
N=observedLabels.M4.Dyslexic(:,2)>.5;
S=observedLabels.M4.Dyslexic(:,5)>.5;
A=observedLabels.M4.Dyslexic(:,8)>.5;
B=observedLabels.M4.Dyslexic(:,11)>.5;

targetLabels(:,3) = [N; S; A; B];

%save([datPath,'Hasselman2014_targetlabelsLOGITM3M4.mat']);

%% QDA based on sample labelling (M3) plot Figure 11
load([datPath,'Hasselman2014_targetlabelsLOGITM3M4.mat']);

h0 = figure; set(h0,'NextPlot','add') %,'Renderer','painters')
subplot(2,7,[1 2]); ax1 = gca; axis square;
subplot(2,7,[3 4]); ax2 = gca; axis square;
subplot(2,7,[5 6]); ax3 = gca; axis square;
subplot(2,7,[8 9]); ax4 = gca; axis square;
subplot(2,7,[10 11]); ax5 = gca; axis square;
subplot(2,7,[12 13]); ax6 = gca; axis square;
subplot(2,7,[7 14]); ax7 = gca; axis square;

set([ax1 ax2 ax3 ax4 ax5 ax6 ax7],'NextPlot','add');
maximize(h0);

warning off

markr= ['o','v','^','d'];
CLR  = [0 0 0;.5 .5 .5];
CLR2 = [.3 .3 .3;.7 .7 .7];

indM = [1:10;11:20;21:30;31:40];
indT = [1 10;11 20;21 30;31 40];

% Use column 1 of observedGR = targets based on M3
observed = targetLabels(:,1);

for s = 1:length(observed)
 if observed(s) == 1
  labels(s,1) = 1;
  labstr{s,1} = '/dAk/';
  clrT(s,:)   = CLR(2,:);
 else
  labels(s,1) = 0;
  labstr{s,1} = '/bAk/';
  clrT(s,:)   = CLR(1,:);
 end
end

% Get the measures
for i = 1:40
 MFMAX(i) = std(mf(i).hq(51:end),1)/mean(mf(i).hq(51:end));
 MFMIN(i) = std(mf(i).hq(1:50),1)/mean(mf(i).hq(51:end));
 DET(i)   = rpSTATS(i,3);
 LAM(i)   = rpSTATS(i,7);
 INH(i)   = HNR(i).HNR;
 NOISY(i) = RTent(i);
 FT       = Formants(i).tracks{5};
 TT       = Formants(i).tracks{1};
 bt       = find(TT>=rpTS(i).ts(1,1),1,'first');
 et       = find(TT>=.4,1,'first');
 F2       = polyfit(FT(bt:et),TT(bt:et),1);
 FS(i)    = F2(2);
 IA(i)    = STIM(i).IASmxO;
 
 LDr       = LAM(i)/DET(i);
 MFr       = MFMIN(i)/MFMAX(i);
 
end

% Convert to unit scale
MFMAX = unit(MFMAX);
MFMIN = unit(MFMIN);
DLD   = unit((DET-LAM));
DET   = unit(DET);
LAM   = unit(LAM);
INH   = unit(INH);
NOISY = unit(NOISY);
IA    = unit(IA);
FS    = unit(FS);

LDr  = unit(LDr);
MFr  = unit(MFr);

MEASURES  = {'CV Positive q (CVhq+)', 'CV Negative q (CVhq-)','Determinism (DET)','Laminarity (LAM)','Inharmonicity (HNR)',...
 'Rise Time Entropy (RFTe)','Max Envelope Slope (maxENV) ','Second Formant Sweep (\Delta F2)'};

for i = 1:6
 
 switch i
  case 1
   featvec = [MFMAX; MFMIN];
   rp1 = 1; rp2 = 2;
   ax = ax1;
  case 2
   featvec = [DET; LAM];
   rp1 = 3; rp2 = 4;
   ax = ax2;
  case 3
   featvec = [INH; NOISY];
   rp1 = 5; rp2 = 6;
   ax = ax3;
  case 4
   featvec = [IA; NOISY];
   rp1 = 7; rp2 = 6;
   ax = ax4;
  case 5
   featvec = [FS; INH];
   rp1 = 8; rp2 = 5;
   ax=ax5;
  case 6
   featvec = [FS; IA];
   rp1 = 8; rp2 = 7;
   ax=ax6;
 end
 
 gridsize=50;
 
 mn1=0; mx1=1; step1=(mx1-mn1)/gridsize;
 mn2=0; mx2=1; step2=(mx2-mn2)/gridsize;
 
 [x,y] = meshgrid(mn1:step1:mx1,mn2:step2:mx2);
 x1 = x(:);y1 = y(:);
 
 % NOTE: adapted gdalda_traintest to return a vector with mean class membership (Class)
 [TOTcl(i).Class, TOTcl(i).Perform, TOTcl(i).Confidence, TOTcl(i).mu0E,...
  TOTcl(i).mu1E, TOTcl(i).C0E, TOTcl(i).C1E] = qda_traintest(featvec, labels, bootN);
 
 axes(ax)
 
 [z, ld] = qdalda_classify([x1(:)'; y1(:)'], TOTcl(i).mu0E, TOTcl(i).mu1E,...
  TOTcl(i).C0E, TOTcl(i).C1E);
 
 z = reshape(z, length(x), length(y));
 [C,hz]=contour(x, y, z, [0 0], '-k');
 
 z=z(:);
 for zc = 1:length(z)
  if z(zc)<=0
   z1{zc,1} = '/dAk/';
  else
   z1{zc,1} = '/bAk/';
  end
 end
 
 if length(unique(z1))<2
  hi1 = 2; hi2=3;
 else
  hi1 = 4; hi2=6;
 end
 hg=gscatter(x1,y1,z1,CLR2,'..',4);
 
 for c = 1:length(TOTcl(i).Class)
  if round(TOTcl(i).Class(c))==0
   clrM(c,:) = CLR(1,:);
  else
   clrM(c,:) = CLR(2,:);
  end
 end
 
 for m = 1:4
  h1=scatter(ax,featvec(1,indM(m,:)),featvec(2,indM(m,:)),...
   32,clrM(indM(m,:),:),markr(m),'filled');
  h2=scatter(ax,featvec(1,indM(m,:)),featvec(2,indM(m,:)),...
   144,clrT(indM(m,:),:),markr(m)); set(h2,'LineWidth',2);
 end
 
 txt = {['Correct:'],['/bAk/: ',num2str(TOTcl(i).Perform(1)*100,'%3.1f%%')],...
  ['/dAk/: ',num2str(TOTcl(i).Perform(2)*100,'%3.1f%%')],...
  ['Total: ',num2str(TOTcl(i).Perform(3)*100,'%3.1f%%')]};
 
 axes(ax)
 nticks1=(mx1-mn1)/10;nticks2=(mx2-mn2)/10;
 set(ax,'XLim',[mn1-(nticks1/2) mx1+(nticks1/2)],'YLim',...
  [mn2-(nticks2/2) mx2+(nticks2/2)],'Box','on','XTick',...
  [mn1:nticks1:mx1],'YTick',[mn2:nticks2:mx2]);
 
 switch i
  case 1
   text(.1,.4,txt,'BackgroundColor','w','EdgeColor','k','Margin',5);
   title(ax,'Complex Temporal Patterns');
  case 2
   text(.1,.8,txt,'BackgroundColor','w','EdgeColor','k','Margin',5);
   title(ax,'Complex Temporal Patterns');
  case 3
   text(.8,.1,txt,'BackgroundColor','w','EdgeColor','k','Margin',5);
   title(ax,'Periodicity Measures');
  case 4
   text(.8,.1,txt,'BackgroundColor','w','EdgeColor','k','Margin',5);
   title(ax,'Periodicity & Component Process Measures');
  case 5
   text(.7,.8,txt,'BackgroundColor','w','EdgeColor','k','Margin',5);
   title(ax,'Periodicity & Component Process Measures');
  case 6
   text(.7,.8,txt,'BackgroundColor','w','EdgeColor','k','Margin',5);
   title(ax,'Component Process Measures');
 end
 
 xlabel(MEASURES{rp1});ylabel(MEASURES{rp2});
 legend off
 
 clear featvec x x1 y y1 z z1 zc m c clrM ldmn1 mx1 step1 mn2 mx2 step2 h1 h2
 clear nticks1 nticks2 rp1 rp2
 
end

% Create "legend" plot

subplot(2,7,[7 14]); ax7 = gca; axis square;
ax=ax7;

s1=plot(ax,.41,.8,'o','Color',CLR(1,:),'MarkerSize',12,'LineWidth',2); hold on;
s2=plot(ax,.41,.7,'v','Color',CLR(1,:),'MarkerSize',12,'LineWidth',2); hold on;
s3=plot(ax,.41,.6,'^','Color',CLR(1,:),'MarkerSize',12,'LineWidth',2); hold on;
s4=plot(ax,.41,.5,'d','Color',CLR(1,:),'MarkerSize',12,'LineWidth',2); hold on;


s9 =plot(ax,.41,.35,'o','Color',CLR(1,:),'MarkerSize',6,'MarkerFaceColor',CLR(1,:)); hold on;
s10=plot(ax,.41,.25,'v','Color',CLR(1,:),'MarkerSize',6,'MarkerFaceColor',CLR(1,:)); hold on;
s11=plot(ax,.41,.15,'^','Color',CLR(1,:),'MarkerSize',6,'MarkerFaceColor',CLR(1,:)); hold on;
s12=plot(ax,.41,.05,'d','Color',CLR(1,:),'MarkerSize',6,'MarkerFaceColor',CLR(1,:)); hold on;

s5=plot(ax,.56,.8,'o','Color',CLR(2,:),'MarkerSize',12,'LineWidth',2); hold on;
s6=plot(ax,.56,.7,'v','Color',CLR(2,:),'MarkerSize',12,'LineWidth',2); hold on;
s7=plot(ax,.56,.6,'^','Color',CLR(2,:),'MarkerSize',12,'LineWidth',2); hold on;
s8=plot(ax,.56,.5,'d','Color',CLR(2,:),'MarkerSize',12,'LineWidth',2); hold on;


s13=plot(ax,.56,.35,'o','Color',[.5 .5 .5],'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]); hold on;
s14=plot(ax,.56,.25,'v','Color',[.5 .5 .5],'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]); hold on;
s15=plot(ax,.56,.15,'^','Color',[.5 .5 .5],'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]); hold on;
s16=plot(ax,.56,.05,'d','Color',[.5 .5 .5],'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]); hold on;

text(.35,.9,'/dAk/');
text(.50,.9,'/bAk/');
text(.7,.9,'Manipulation');

text(0,.65,'Observed as');
text(0,.2,'Classified as');

text(.7,.8,'None');
text(.7,.7,'Slowed Down');
text(.7,.6,'Amplified');
text(.7,.5,'Both');

text(.7,.35,'None');
text(.7,.25,'Slowed Down');
text(.7,.15,'Amplified');
text(.7,.05,'Both');

title('Legend');
ylim([0 1]); xlim([0 1]);
set(gca,'XTickLabel',[],'YTickLabel',[]);
axis off

% Uncomment both lines to save the Figure
%patch2ind(gcf) % This is necessary in order to use 'Renderer','painters' with patch object.
%grab('Hasselman2014_Figure11',0)

keep source datPath bootN method observedLabels targetLabels MFMAX MFMIN DET LAM INH NOISY IA FS TOTcl

% Use these results for table 8 and Figure 12
%save([datPath,'Hasselman2014_QDAresults.mat']);

%% QDA based on labeling of average readers and dyslexic readers (M4)
% Uncomment to load variables
load([datPath,'Hasselman2014_QDAresults.mat']);

% reproduce results
rng(12345)

for i = 1:6
 
 switch i
  case 1
   featvec = [MFMAX; MFMIN];
  case 2
   featvec = [DET; LAM];
  case 3
   featvec = [INH; NOISY];
  case 4
   featvec = [IA; NOISY];
  case 5
   featvec = [FS; INH];
  case 6
   featvec = [FS; IA];
 end
 
 % AVERAGE labeling results
 observed = targetLabels(:,2);
 
 [AVEcl(i).Class, AVEcl(i).Perform, AVEcl(i).Confidence, AVEcl(i).mu0E,...
  AVEcl(i).mu1E, AVEcl(i).C0E, AVEcl(i).C1E] = qda_traintest(featvec, observed, bootN);
 
 % DYSLECTIC labeling results
 observed = targetLabels(:,3);
 
 [DYScl(i).Class, DYScl(i).Perform, DYScl(i).Confidence, DYScl(i).mu0E,...
  DYScl(i).mu1E, DYScl(i).C0E, DYScl(i).C1E] = qda_traintest(featvec, observed, bootN);
 
end

keep tab tabt observedLabels targetLabels TOTcl AVEcl DYScl Group Feature datPath

% Create Table 8 with Classifier results to paste into wordprocessor

f=length(TOTcl);
for i = 1:f
 tab(i,:) = [TOTcl(i).Perform(1) TOTcl(i).Confidence(1),...
  TOTcl(i).Perform(2) TOTcl(i).Confidence(2) TOTcl(i).Perform(3) TOTcl(i).Confidence(3)];
 tab(i+f,:) = [AVEcl(i).Perform(1) AVEcl(i).Confidence(1)....
  AVEcl(i).Perform(2) AVEcl(i).Confidence(2) AVEcl(i).Perform(3) AVEcl(i).Confidence(3)];
 tab(i+2*f,:) = [DYScl(i).Perform(1) DYScl(i).Confidence(1),...
  DYScl(i).Perform(2) DYScl(i).Confidence(2) DYScl(i).Perform(3) DYScl(i).Confidence(3)];
end

tab = tab.*100;
Group = {'Sample','Average Readers','Dyslexic Readers'};
Feature={'CVhq+ / CVhq-','LAM / DET','HNR / RFTe','mxENV / RFTe','F2 / HNR','F2 / mxENV'};


[r c]=size(tab);
cnt=1;
cnt2=1;
% Copy and paste tabt into Excel or Numbers
for i=1:r
 if ((i==c)||(i==2*c))
  cnt=cnt+1;
 end
 tabt{i,1} = Group{cnt};
 tabt{i,2} = Feature{cnt2};
 for j=3:8
  tabt{i,j} = num2str(tab(i,j-2),'%3.1f%%');
 end
 if cnt2==c
  cnt2=1;
 else
  cnt2=cnt2+1;
 end
end

clear i j

keep tab tabt observedLabels targetLabels TOTcl AVEcl DYScl Group Feature datPath
% save([datPath,'Hasselman2014_Table8.mat']);

%% Create Figure 12
% Uncomment to load results
%load([datPath,'Hasselman2014_Table8.mat']);

markr= ['o','v','^','d'];
CLR  = [0 0 0;.5 .5 .5];
CLR2 = [.7 .7 .7;.85 .85 .85];

h0 = figure;
maximize(h0);

stims = [1:10; 11:20; 21:30; 31:40];
cols = [2:3; 5:6; 8:9; 11:12];


x = 1:10;

averageM4 = [];
for i = 1:length(cols)
 averageM4  = [averageM4; observedLabels.M4.Average(:,cols(i))];
end

cnt=0;
for i=1:4
 for j=1:6
  cnt=cnt+1;
  subplot(4,6,cnt)
  ciplot(observedLabels.M4.Average(:,cols(i,1)),observedLabels.M4.Average(:,cols(i,2)),x',CLR2(2,:),'none'); hold on;
  plot([1 10],[.5 .5],':k'); hold on;
  hl(1)=plot(1:10,targetLabels(stims(i,:),2),'-','Color',CLR(2,:),'LineWidth',2,'Marker',markr(i),'MarkerSize',8,'MarkerFaceColor',CLR2(1,:)); hold on;
  hl(2)=plot(1:10,AVEcl(j).Class(stims(i,:))','-','Color',CLR(1,:),'LineWidth',1,'Marker',markr(i),'MarkerSize',8,'MarkerFaceColor',CLR(2,:));hold on;
  ylim([-0.1 1.1]);
  xlim([0.9 10.1]);
  set(gca,'YTick',[0 .5 1],'YTickLabel',{'0','.5','1'},...
   'XTick',[1:10],'XTickLabel',{'/bAk/','2','3','4','5','6','7','8','9','/dAk/'});
  
  id1=round(averageM4(stims(i,:)))' == round(AVEcl(j).Class(stims(i,:)));
  if any(id1==0)
   plot(x(id1==0)',AVEcl(j).Class(stims(i,id1==0))','xr','MarkerSize',10,'LineWidth',2); hold on;
  end
  
  ax0=gca;
  
  if i==1
   title(Feature{j},'FontSize',16)
   
   switch j
    case 1
     title(ax0,'Complex Temporal Patterns');
    case 2
     title(ax0,'Complex Temporal Patterns');
    case 3
     title(ax0,'Periodicity Measures');
    case 4
     title(ax0,'Periodicity & Component Process Measures');
    case 5
     title(ax0,'Periodicity & Component Process Measures');
    case 6
     title(ax0,'Component Process Measures');
   end
   
  end
  
  if ismember(cnt,[1,7,13,19])
   Opos(i,:)=get(ax0,'Position');
   ylabel('\pi /dAk/')
   legend(hl,{'Target','Estimate'},'Location','East');
   %    legend(hf,{'Predicted CI'},'Location','South');
   legend('boxoff')
  end
  
  if cnt >= 19
   xlabel('Stimulus');
   Hpos(j,:) = get(ax0,'Position');
  end
  
 end
end
colormap(gray(2));
Tpos = [-.1 -0.07 0 0];

h_1= annotation('textbox',[Opos(1,:)+Tpos],'String','None','EdgeColor','none','FontSize',16);
set(h_1,'FitBoxToText','on');
h_2= annotation('textbox',[Opos(2,:)+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
set(h_2,'FitBoxToText','on');
h_3= annotation('textbox',[Opos(3,:)+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
set(h_3,'FitBoxToText','on');
h_4=annotation('textbox',[Opos(4,:)+Tpos],'String','Both','EdgeColor','none','FontSize',16);
set(h_4,'FitBoxToText','on');

for f=1:6
 h_5(f)= annotation('textbox',[(Hpos(f,1)+.05) 0.07 0 0],'String',Feature{f},'EdgeColor','none','FontSize',16,'HorizontalAlignment','center');
 set(h_5(f),'FitBoxToText','on');
 h_6(f)= annotation('textbox',[(Hpos(f,1)+.05) .98 0 0],'String',Feature{f},'EdgeColor','none','FontSize',16,'HorizontalAlignment','center');
 set(h_6(f),'FitBoxToText','on');
end

h_7=annotation('textbox',[Hpos(1,1)-.1 0.98 0 0],'String','Average Readers','EdgeColor','none','FontSize',16,'FontWeight','bold');
set(h_7,'FitBoxToText','on');

keep tab tabt observedLabels targetLabels TOTcl AVEcl DYScl Group Feature datPath

% Uncomment both lines to save Figure
% patch2ind(gcf)
% grab('Hasselman2014_Figure12',0);

%% Create Figure 13
% Uncomment to load results
%load([datPath,'Hasselman2014_Table8.mat']);

markr= ['o','v','^','d'];
CLR  = [0 0 0;.5 .5 .5];
CLR2 = [.7 .7 .7;.85 .85 .85];

h0 = figure;
maximize(h0);

stims = [1:10; 11:20; 21:30; 31:40];
cols  = [2:3; 5:6; 8:9; 11:12];

x = 1:10;

dyslexicM4 = [];
for i = 1:length(cols)
 dyslexicM4 = [dyslexicM4; observedLabels.M4.Dyslexic(:,cols(i))];
end

cnt=0;
for i=1:4
 for j=1:6
  cnt=cnt+1;
  subplot(4,6,cnt)
  ciplot(observedLabels.M4.Dyslexic(:,cols(i,1)),observedLabels.M4.Dyslexic(:,cols(i,2)),x',CLR2(2,:),'none'); hold on;
  plot([1 10],[.5 .5],':k'); hold on;
  hl(1)=plot(1:10,targetLabels(stims(i,:),3),'-','Color',CLR(2,:),'LineWidth',2,'Marker',markr(i),'MarkerSize',8,'MarkerFaceColor',CLR2(1,:)); hold on;
  hl(2)=plot(1:10,DYScl(j).Class(stims(i,:))','-','Color',CLR(1,:),'LineWidth',1,'Marker',markr(i),'MarkerSize',8,'MarkerFaceColor',CLR(2,:));hold on;
  ylim([-0.1 1.1]);
  xlim([0.9 10.1]);
  set(gca,'YTick',[0 .5 1],'YTickLabel',{'0','.5','1'},...
   'XTick',[1:10],'XTickLabel',{'/bAk/','2','3','4','5','6','7','8','9','/dAk/'});
  
  id1=round(dyslexicM4(stims(i,:)))' == round(DYScl(j).Class(stims(i,:)));
  if any(id1==0)
   plot(x(id1==0)',DYScl(j).Class(stims(i,id1==0))','xr','MarkerSize',10,'LineWidth',2); hold on;
  end
  
  ax0=gca;
  
  if i==1
   title(Feature{j},'FontSize',16)
   
   switch j
    case 1
     title(ax0,'Complex Temporal Patterns');
    case 2
     title(ax0,'Complex Temporal Patterns');
    case 3
     title(ax0,'Periodicity Measures');
    case 4
     title(ax0,'Periodicity & Component Process Measures');
    case 5
     title(ax0,'Periodicity & Component Process Measures');
    case 6
     title(ax0,'Component Process Measures');
   end
   
  end
  
  if ismember(cnt,[1,7,13,19])
   Opos(i,:)=get(ax0,'Position');
   ylabel('\pi /dAk/')
   legend(hl,{'Target','Estimate'},'Location','East');
   legend('boxoff')
  end
  
  if cnt >= 19
   xlabel('Stimulus');
   Hpos(j,:) = get(ax0,'Position');
  end
  
 end
end

colormap(gray(2));
Tpos = [-.1 -0.07 0 0];

h_1= annotation('textbox',[Opos(1,:)+Tpos],'String','None','EdgeColor','none','FontSize',16);
set(h_1,'FitBoxToText','on');
h_2= annotation('textbox',[Opos(2,:)+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
set(h_2,'FitBoxToText','on');
h_3= annotation('textbox',[Opos(3,:)+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
set(h_3,'FitBoxToText','on');
h_4=annotation('textbox',[Opos(4,:)+Tpos],'String','Both','EdgeColor','none','FontSize',16);
set(h_4,'FitBoxToText','on');

for f=1:6
 h_5(f)= annotation('textbox',[(Hpos(f,1)+.05) 0.07 0 0],'String',Feature{f},'EdgeColor','none','FontSize',16,'HorizontalAlignment','center');
 set(h_5(f),'FitBoxToText','on');
 h_6(f)= annotation('textbox',[(Hpos(f,1)+.05) .98 0 0],'String',Feature{f},'EdgeColor','none','FontSize',16,'HorizontalAlignment','center');
 set(h_6(f),'FitBoxToText','on');
end

h_7=annotation('textbox',[Hpos(1,1)-.1 0.98 0 0],'String','Dyslexic Readers','EdgeColor','none','FontSize',16,'FontWeight','bold');
set(h_7,'FitBoxToText','on');

keep tab tabt observedLabels targetLabels TOTcl AVEcl DYScl Group Feature datPath

% Uncomment both lines to save Figure
% patch2ind(gcf)
% grab('Hasselman2014_Figure13',0);
