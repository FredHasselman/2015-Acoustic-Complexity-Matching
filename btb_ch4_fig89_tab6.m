%% Supplementary Material to "Beyond the Boundary - Chapter 4"
%
%%% Introduction
% This code performs the Quadratic Discriminant analysis on extracted features of the speech signal.
% Use it to create Figure 4.8 and 4.9 and the data in Table 4.6
%
% This is not a proper function or toolbox, it is not optimized for speed or functionality or aestetics! 
% Evaluate the code per Cell (using MATLAB's Cell Mode) and inspect the workspace to see what is going on.
%
% It is possible your monitor size will influence figure legibility.

%% Data / Toolboxes / Scripts, etc. that need to be on the MATLAB PATH
%
% * QDA code by Max Little (<http://www.maxlittle.net/home/index.php Homepage>) used in:
% M.A. Little, P.E. McSharry, I.M. Moroz, S.J. Roberts (2006). Nonlinear,
% Biophysically-Informed Speech Pathology Detection in Proceedings of IEEE ICASSP
% 2006, IEEE Publishers: Toulouse, France.
% * Several helper scripts and data files created and authored by Fred Hasselman included in the GitHub repository
% * Fred's toolbox: https://github.com/FredHasselman/toolboxML

%% Author / Version / License
%
% Repository: <https://github.com/FredHasselman/BTB-Supplemental-Material-CHAPTER4/ *BTB Chapter 4 on GitHub*>
%
% Created by: <http://www.fredhasselman.com/ *Fred Hasselman*> / January 2011
% Affiliations: <http://www.ru.nl/bsi _Behavioural Science Institute_> - <http://www.ru.nl _Radboud University Nijmegen_>

%% PREP

% Uncomment next line to clear and close everything... detergent grade!
% omo

% Change ... to the path on your machine where you stored the files
% cd('...');
load('btb_ch4_stimfeatures_ORI.mat');
load('observedLOGIT.mat');

% Method to use for Discriminant analysis
method = 'quadratic';

% Number of bootstrap repetitions for the classification
bootN = 1500;

%% QDA based on sample labelling (M3) plot Figure 4.8

h0 = figure; set(h0,'NextPlot','add')
subplot(2,6,[1 2]); ax1 = gca; axis square;
subplot(2,6,[3 4]); ax2 = gca; axis square;
subplot(2,6,[5 6]); ax3 = gca; axis square;
subplot(2,6,[7 8]); ax4 = gca; axis square;
subplot(2,6,[9 10]); ax5 = gca; axis square;
set([ax1 ax2 ax3 ax4 ax5],'NextPlot','add');
maximize(h0);

markr= ['o','v','^','d'];
CLR  = [0 0 0;.5 .5 .5];
CLR2 = [.3 .3 .3;.7 .7 .7];

indM = [1:10;11:20;21:30;31:40];
indT = [1 10;11 20;21 30;31 40];

% Use column 1 of observedGR = targets based on M3
observed = observedGR(:,1);

for i = 1:length(observed)
 if round(observed(i)) == 1
  labels(i,1) = 1;
  labstr{i,1} = '/dAk/';
  clrT(i,:)   = CLR(2,:);
 else
  labels(i,1) = 0;
  labstr{i,1} = '/bAk/';
  clrT(i,:)   = CLR(1,:);
 end
end


% Get the measures
for i = 1:40
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
end

% Convert to unit scale
DET   = unit(DET);
LAM   = unit(LAM);
INH   = unit(INH);
NOISY = unit(NOISY);
IA    = unit(IA);
FS    = unit(FS);

MEASURES  = {'DETERMINISM','LAMINARITY','Inharmonicity (HNR)',...
 'Rise Time Entropy','Max Envelope Slope','Formant Sweep (F2 Slope)'};

for i = 1:5
 
 switch i
  case 1
   featvec = [DET; LAM];
   rp1 = 1; rp2 = 2;
   ax = ax1;
  case 2
   featvec = [INH; NOISY];
   rp1 = 3; rp2 = 4;
   ax = ax2;
  case 3
   featvec = [IA; NOISY];
   rp1 = 5; rp2 = 4;
   ax = ax3;
  case 4
   featvec = [FS; INH];
   rp1 = 6; rp2 = 3;
   ax=ax4;
  case 5
   featvec = [FS; IA];
   rp1 = 6; rp2 = 5;
   ax=ax5;
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
 [C,hz]=contour(x, y, z, [0 0], 'k-');
 
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
 
 if i == 1
  text(.01,.8,txt,'BackgroundColor','w','EdgeColor','k','Margin',5);
 elseif i==5
  text(.7,.8,txt,'BackgroundColor','w','EdgeColor','k','Margin',5);
 else
  text(.7,.2,txt,'BackgroundColor','w','EdgeColor','k','Margin',5);
 end
 xlabel(MEASURES{rp1});ylabel(MEASURES{rp2});
 legend off

 clear featvec x x1 y y1 z z1 zc m c clrM ld h1 h2 mn1 mx1 step1 mn2 mx2 step2
 clear nticks1 nticks2 rp1 rp2
  
end

% Create "legend" plot
subplot(2,6,[11 12]); ax6 = gca; axis square;
ax=ax6;

s1=plot(ax,.4,.8,'o','Color',CLR(1,:),'MarkerSize',12,'LineWidth',2); hold on;
s2=plot(ax,.4,.7,'v','Color',CLR(1,:),'MarkerSize',12,'LineWidth',2); hold on;
s3=plot(ax,.4,.6,'^','Color',CLR(1,:),'MarkerSize',12,'LineWidth',2); hold on;
s4=plot(ax,.4,.5,'d','Color',CLR(1,:),'MarkerSize',12,'LineWidth',2); hold on;
s5=plot(ax,.55,.8,'o','Color',CLR(2,:),'MarkerSize',12,'LineWidth',2); hold on;
s6=plot(ax,.55,.7,'v','Color',CLR(2,:),'MarkerSize',12,'LineWidth',2); hold on;
s7=plot(ax,.55,.6,'^','Color',CLR(2,:),'MarkerSize',12,'LineWidth',2); hold on;
s8=plot(ax,.55,.5,'d','Color',CLR(2,:),'MarkerSize',12,'LineWidth',2); hold on;

s9 =plot(ax,.4,.35,'o','Color',CLR(1,:),'MarkerSize',6,'MarkerFaceColor',CLR(1,:)); hold on;
s10=plot(ax,.4,.25,'v','Color',CLR(1,:),'MarkerSize',6,'MarkerFaceColor',CLR(1,:)); hold on;
s11=plot(ax,.4,.15,'^','Color',CLR(1,:),'MarkerSize',6,'MarkerFaceColor',CLR(1,:)); hold on;
s12=plot(ax,.4,.05,'d','Color',CLR(1,:),'MarkerSize',6,'MarkerFaceColor',CLR(1,:)); hold on;
s13=plot(ax,.55,.35,'o','Color',CLR(2,:),'MarkerSize',6,'MarkerFaceColor',CLR(2,:)); hold on;
s14=plot(ax,.55,.25,'v','Color',CLR(2,:),'MarkerSize',6,'MarkerFaceColor',CLR(2,:)); hold on;
s15=plot(ax,.55,.15,'^','Color',CLR(2,:),'MarkerSize',6,'MarkerFaceColor',CLR(2,:)); hold on;
s16=plot(ax,.55,.05,'d','Color',CLR(2,:),'MarkerSize',6,'MarkerFaceColor',CLR(2,:)); hold on;
 
text(.36,.9,'/dAk/');
text(.51,.9,'/bAk/');
text(.7,.9,'Manipulation');

text(.05,.65,'Observed as');
text(.05,.2,'Classified as');

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

keep bootN observedGR DET LAM INH NOISY IA FS TOTcl

%Print to file (for some reason cannot use 'painters' renderer)
% shg;
% set(gcf, 'PaperPositionMode', 'auto')
% print('-deps2','-r900','btb_ch4_Figure8')

%save('btb_ch4_QDAresults.mat')

%% QDA based on labeling of average readers and dyslexic readers (M4)

for i = 1:5
 
 switch i
  case 1
   featvec = [DET; LAM];
   rp1 = 1; rp2 = 2;
  case 2
   featvec = [INH; NOISY];
   rp1 = 3; rp2 = 4;
  case 3
   featvec = [IA; NOISY];
   rp1 = 5; rp2 = 4;
  case 4
   featvec = [FS; INH];
   rp1 = 6; rp2 = 3;
  case 5
   featvec = [FS; IA];
   rp1 = 6; rp2 = 5;
 end
 
 % AVERAGE labeling results
 observed = observedGR(:,2);
 
 [AVEcl(i).Class AVEcl(i).Perform, AVEcl(i).Confidence, AVEcl(i).mu0E,...
  AVEcl(i).mu1E, AVEcl(i).C0E, AVEcl(i).C1E] = qda_traintest(featvec, observed, bootN);

 % DYSLECTIC labeling results
 observed = observedGR(:,3);

 [DYScl(i).Class DYScl(i).Perform, DYScl(i).Confidence, DYScl(i).mu0E,...
  DYScl(i).mu1E, DYScl(i).C0E, DYScl(i).C1E] = qda_traintest(featvec, observed, bootN);
 
end

keep observedGR TOTcl AVEcl DYScl

%% Create Table 4.6 with Classifier results to paste into wordprocessor
 
 for i = 1:5
  tab(i,:) = [TOTcl(i).Perform(1) TOTcl(i).Confidence(1),...
   TOTcl(i).Perform(2) TOTcl(i).Confidence(2) TOTcl(i).Perform(3) TOTcl(i).Confidence(3)];
  tab(i+5,:) = [AVEcl(i).Perform(1) AVEcl(i).Confidence(1)....
   AVEcl(i).Perform(2) AVEcl(i).Confidence(2) AVEcl(i).Perform(3) AVEcl(i).Confidence(3)];
  tab(i+10,:) = [DYScl(i).Perform(1) DYScl(i).Confidence(1),...
   DYScl(i).Perform(2) DYScl(i).Confidence(2) DYScl(i).Perform(3) DYScl(i).Confidence(3)];
 end
 
 tab = tab.*100;
 
 % Copy and paste tabt into Excel or Numbers
 for i=1:15
  for j=1:6
   tabt{i,j} = num2str(tab(i,j),'%3.1f%%');
  end
 end
 
 clear i j 
 
 save('btb_ch4_Table6.mat');
 
 %% Create Figure 4.9

 h1 = figure;
 maximize(h1);

 for i = 1:5
  tabl(i,:) = [TOTcl(i).Class(1:10) TOTcl(i).Class(11:20),...
   TOTcl(i).Class(21:30) TOTcl(i).Class(31:40)];
  tabl(i+5,:) = [AVEcl(i).Class(1:10) AVEcl(i).Class(11:20),....
   AVEcl(i).Class(21:30) AVEcl(i).Class(31:40)];
  tabl(i+10,:) = [DYScl(i).Class(1:10) DYScl(i).Class(11:20),...
   DYScl(i).Class(21:30) DYScl(i).Class(31:40)];
 end

stims = [1:10; 11:20; 21:30; 31:40];
cnt=0;
for i=1:5
 for j=1:4
  cnt=cnt+1;
  subplot(5,4,cnt)
  plot(1:10,tabl(i,stims(j,:)),'-ko','LineWidth',3);hold on;
  plot(1:10,observedGR(stims(j,:)',1),'-k^','Color',[.5 .5 .5],'LineWidth',3); hold on;
  plot([1 10],[.5 .5],':k'); xlabel('Stimulus');
  ylim([-0.1 1.1]); xlim([0.9 10.1]);
  set(gca,'XTick',[1:10]);
  if i==1
   switch j
    case 1
     title('None');
    case 2
     title('Slowed Down');
    case 3
     title('Amplified');
    case 4
     title('Both');
   end
  end
  
  switch cnt
   case 1
    ylabel('LAM / DET');
   case 5
    ylabel('HNR / RFTe');
   case 9
    ylabel('mxENV / RFTe');
   case 13
    ylabel('F2 / HNR');
   case 17
    ylabel('F2 / mxENV');
  end
  
 end
end

%grab('btb_ch4_Figure9',0);
