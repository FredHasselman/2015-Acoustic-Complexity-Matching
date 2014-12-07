%%

scmin=6;
scmax=12;
ressc=60;

scale=round(2.^[scmin:((scmax-scmin)/ressc):scmax]);

y=cumsum(stimuliMF(1).IA-mean(stimuliMF(1).IA))';
m=1;

for ns=1:length(scale),
    segments(ns)=floor(length(y)/scale(ns));
    for v=1:segments(ns),
        Index=((((v-1)*scale(ns))+1):(v*scale(ns)));
        C=polyfit(Index,y(Index),m);
        fit1=polyval(C,Index);     
        RMS_scale{ns}(v)=sqrt(mean((y(Index)-fit1).^2));
    end
    F(ns)=sqrt(mean(RMS_scale{ns}.^2));
end
Ch = polyfit(log2(scale),log2(F),1);
H = Ch(1);
RegLine = polyval(Ch,log2(scale));

%%
h0=figure;
maximize(h0);

ns  = [21 31 51 61];
col = [9 10 11 12] ;
ids = [53 27 7 4];

subplot(8,4,[1 8])
plot(stimuliMF(1).IA,'-k'); hold on
plot(1:length(stimuliMF(1).IA),cumsum(stimuliMF(1).IA-mean(stimuliMF(1).IA))./1000,'-','Color',[.5 .5 .5],'LineWidth',2); hold on

axis off
grid off

for s = 1:length(ns)

subplot(8,4,col(s))

sc = [1:scale(ns(s)):length(y)]';
x  = (sc(1:end-1)+sc(2:end))./2;
yv = unit(RMS_scale{ns(s)});

plot(y,'-','Color',[.5 .5 .5]); hold on
plot([sc(ids(s)):sc(ids(s)+1)],y([sc(ids(s)):sc(ids(s)+1)]),'-k','LineWidth',2); hold on
set(gca,'XTick',[1:scale(ns(s)):length(y)],'YTick',[],'XTickLabel','','XGrid','on','YGrid','off');
axis tight

ylabel('Profile');
xlabel('');
title(['s = ',num2str(scale(ns(s))),' (scale)  |  Ns = ',num2str(length(sc)),' (segments v)']);

subplot(8,4,col(s)+4)
plot([sc(ids(s)):sc(ids(s)+1)],detrend(y(sc(ids(s)):sc(ids(s)+1))),'-k'); hold on
set(gca,'XTick',[sc(ids(s)) sc(ids(s)+1)],'XTickLabel',{num2str(sc(ids(s)).*(1/stimuli(1).fs),2),num2str(sc(ids(s)+1).*(1/stimuli(1).fs),2)},'YTick',0,'YGrid','on');

ylabel('');
title(['Detrended Segment (v = ',num2str(ids(s)),')'])
axis tight


subplot(8,4,col(s)+8)
plot(x,yv,'ok'); hold on
plot(x(ids(s)),yv(ids(s)),'ok','MarkerFaceColor',[.4 .4 .4],'MarkerSize',8); hold on
set(gca,'XTick',[1:scale(ns(s)):length(y)],'YTick',[],'XTickLabel','','XGrid','on','YGrid','off');
ylabel('F^2(s,v)');
xlabel(['RMS variation of F^2 (s = ',num2str(scale(ns(s))),', N_s = ',num2str(length(sc)),') = ',num2str(F(ns(s)),2)]);

axis tight

end

subplot(8,4,[21 32])
plot(log2(scale),log2(F),'sk'); hold on;
plot(log2(scale(ns)),log2(F(ns)),'sk','MarkerFaceColor',[.4 .4 .4],'MarkerSize',8); hold on;
plot(log2(scale),RegLine,'-k','LineWidth',2); hold on;
xlim([scmin-.5 scmax+.5]);
set(gca,'XTick',[scmin:scmax],'XTickLabel',[2.^[scmin:scmax]],'YTick',floor(min(log2(F))):ceil(max(log2(F))),'YTickLabel',2.^[floor(min(log2(F))):ceil(max(log2(F)))]);
ylabel('RMS variation of [F^2(s,v)]')
xlabel('Scale (s)')

% Uncomment to save for further processing in Vector Graphics Software (e.g., Adobe Illustrator to add arrows and lines)
% grab('Hasselman2014_Figure9',0)
