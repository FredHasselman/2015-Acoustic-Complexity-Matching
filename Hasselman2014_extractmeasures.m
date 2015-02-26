% Supplementary Material to Hasselman (2014) "Classifying Acoustic Signals into Speech Categories"

%%%%%%%%%%%%%% MARKDOWN CODE %%%%%%%%%%%%%%
%
% ### Introduction
% This is a demonstration script accompanying the fourth chapter of my dissertation (Beyond the Boundary). Its purpose is to
% provide an example of how to use various freely available MATLAB sources on the web to extract variables from speech
% stimuli. This script should thus be used as an example to build your own scripts only. It is not a function or toolbox, not
% optimized for speed or functionality! Evaluate the code per Cell (using MATLAB's Cell Mode) and inspect the workspace to
% see what is going on.
%
% Extracts spectral and pitch data, formant sweeps and amplitude envelopes and various measures of temporal evolution of the
% signal.
%
% Results should be comparable to those obtained by [Praat](http://www.praat.org) (Boersma & Weeninck).
%
% ### Data / Toolboxes / Scripts, etc. that need to be on the MATLAB PATH
%
% * The [Signal processing toolbox](http://www.mathworks.com) by the Mathworks
% * [`hilbert2.m` and `derivative.m`](http://www.mathworks.com/matlabcentral/fileexchange/authors/110216) by Scott McKinney
% * The [Sound Processing Toolbox](http://note.sonots.com/SciSoftware/Pitch.html) by Naotoshi Seo
% * [MIR toolbox](https://www.jyu.fi/hum/laitokset/musiikki/en/research/coe/materials/mirtoolbox/mirtoolbox) created by departement of Music at the University of Jyvskyl Finland
% * [`crossing.m`](http://www.mathworks.nl/matlabcentral/fileexchange/2432-crossing) by Steffen Brueckner
% * [CRP Toolbox](http://tocsy.pik-potsdam.de/CRPtoolbox/) by Norbert Marwan. Note that an edited scriptfile `crp_edit.m` in the [GithHub repository](https://github.com/FredHasselman/Acoustic-Complexity-Matching) is needed to extract the RQA threshold values.
% * `MFDFA1.m` available as part of the [Multifractal Toolbox](http://www.ntnu.edu/inm/geri/software). Download the [Introduction to MFDFA (zip-file)](http://www.ntnu.edu/documents/170234/1315232/Introduction_to_MFDFA4.zip). Detailed instruction are available in [Ihlen, E.A.F. (2012). FrontiersIn Fractal Physiology, 3(141), 1-18](http://www.ntnu.edu/documents/170234/1315232/Introduction_to_MFDFA.pdf)
% * [Fred's toolbox](https://github.com/FredHasselman/toolboxML) on GitHub
% * **DATA** available at [Open Science Framework](https://osf.io/hpjse/files/) or directly from [Dropbox](https://www.dropbox.com/sh/i1vp0nlsj3mi3v7/AAD6V5WaQLdmEBTFRIzAqQm8a?dl=0)
% * **STIMULI** available at [Open Science Framework](https://osf.io/hpjse/files/) or directly from [Dropbox](https://www.dropbox.com/sh/i1vp0nlsj3mi3v7/AAD6V5WaQLdmEBTFRIzAqQm8a?dl=0)
% * Scripts are available in a [GithHub repository](https://github.com/FredHasselman/Acoustic-Complexity-Matching)
%
% ### Author / Version / License
%
% Created by: [Fred Hasselman 2011-2014](http://www.fredhasselman.com)
% Affiliations: [School of Pedagogical and Educational Science](http://www.ru.nl/pwo) and [Behavioural Science Institute (Learning & Plasticity)](http://www.ru.nl/bsi) at the [Radboud University Nijmegen, the Netherlands](http://www.ru.nl)
%
%%%%%%%%%%%%%% MARKDOWN CODE %%%%%%%%%%%%%%

%% PREP: SETTINGS
% Reads .WAV data and performs spectral, envelope and pitch analyses.

% Change these variables to reflect your datafile names and locations.
% Here there are four different manipulations (manis) of a
% 10-step continuum (stims), with files named accordingly.

manis     = {'BAKDAK', 'BAKDAKV', 'BAKDAKA', 'BAKDAKB'};
stims     = 10;

% Change 'source' to the path on your machine where you stored the files
% If you copied the dropbox folder from the OSF to your own dropbox, the current path should be Ok (on machines that allow tilde expansion)
source='~/Dropbox/Hasselman2014-PeerJ-Classifying_Acoustic_Signals/';
wavPath   = [source '/STIMULI/'];
txtPath   = [source '/FORMANT TABLES/'];
datPath   = [source '/DATA/'];
ext       = '.WAV';

DELIMITER = '\t';
SKIP      = '--undefined--';

%Spectrogram settings put into structure SPEC
SPEC.fs        = 44100;     % Or read from file, see for loop below
SPEC.ts        = .002;      % Praat settings: Time step         (s)
SPEC.wl        = .005;      % Praat settings: Window length     (s)
SPEC.mf        = 5500;      % Praat settings: Maximum frequency (Hz)
SPEC.sf        = 20;        % Praat settings: Frequency step    (Hz)

%%% Praat default behaviour
%
% Praat changes Frequencey step (|sf|) and Time step (|ts|) to the most
% sensible minimal values given current Window length (wl) in the case they
% were set too low by the user.

if SPEC.ts < SPEC.wl/(8*sqrt(pi))
 SPEC.ts = SPEC.wl/(8*sqrt(pi));
end
if SPEC.sf < (sqrt(pi)/8/SPEC.wl)
 SPEC.sf = (sqrt(pi)/8/SPEC.wl);
end

SPEC.f       = 1:SPEC.sf:SPEC.mf; % Spectrogram will be est. at freqs in SPEC.f
SPEC.wintype = 'gausswin';        % Praat settings: Gaussian window is used

SPEC.nfft    = round(SPEC.wl*SPEC.fs); % Convert s to samples
SPEC.noverlap= round(SPEC.ts*SPEC.fs); % Convert s to samples
SPEC.window  = eval(sprintf('%s(SPEC.nfft)', SPEC.wintype)); % Create window

%Formant tracking with LPC, female speaker of 5 formants
SPEC.fl = 25;                     % Frame length  (ms)
SPEC.fo = 5;                      % Frame Overlap (ms)
SPEC.po = 10;                     % Prediction order (LPC poles / coefs)
SPEC.d  = gcd(2*SPEC.mf,SPEC.fs); % Get smallest integer conversion ratio
SPEC.p  = (2*SPEC.mf)/SPEC.d;     % Resample to   factor for LPC
SPEC.q  = SPEC.fs/SPEC.d;         % Resample from factor for LPC

% CHECK SETTINGS
%
% * *Check frequency vector:* [f(1) f(end)]
% * *Check the window properties:* wvtool(SPEC.window)

%% Load the .WAV files into MATLAB

cnt=0;
cnt2=0;
for mani = 1:length(manis)
 for stim = 1:stims
  
  wavFile = [wavPath char(manis(mani)) num2str(stim) ext];
  txtFile = [txtPath char(manis(mani)) num2str(stim) '.txt'];
  
  if exist(wavFile,'file')
   cnt=cnt+1;
   
   % Get wavefile and put into structure [stimuli(cnt)._]
   [stimuli(cnt).y stimuli(cnt).fs stimuli(cnt).nbits,stimuli(cnt).opts]=wavread(wavFile);
   stimuli(cnt).name = [char(manis(mani)) num2str(stim)];
   
  else
   disp(wavFile)
  end
  
  % Read Praat Formant track files
  if exist(txtFile,'file')
   cnt2=cnt2+1;
   fid=fopen(txtFile,'r');
   Formants(cnt2).tracks = textscan(fid,'%*s %f %f %f %f %f %f %f %f','TreatAsEmpty',SKIP,'HeaderLines',1);
   Formants(cnt2).name = txtFile;
   fclose(fid);
  else
   disp(txtFile);
  end
  
 end
end

clear cnt cnt2 wavFile wavPath ext txtPath txtFile DELIMITER SKIP fid

save([datPath,'Hasselman2104_stimfeatures.mat'],'stimuli','Formants');

%% Get Envelope and Formant measures

cnt=0;
for mani = 1:length(manis)
 for stim = 1:stims
  
  cnt=cnt+1;
  
  % Perform spectrogram analysis (SP toolbox) and put results into structure
  % NOTE: If you added the MIRtoolbox to your path it is possible there will be a naming conflict!
  %       Just move the MIRtoolbox - AuditoryToolbox to the bottom of the search path
  [STIM(cnt).S STIM(cnt).F STIM(cnt).T, STIM(cnt).P] = spectrogram(stimuli(cnt).y,SPEC.window,SPEC.noverlap,SPEC.f,SPEC.fs);
  STIM(cnt).name = [char(manis(mani)) num2str(stim)];
  
  % Get immediate amplitude envelope and frequency by Hilbert transform
  % (HILBERT2 by Scott McKinney) and put results into structure STIM(cnt).IA .IF
  [STIM(cnt).IA, STIM(cnt).IF] = hilbert2(stimuli(cnt).y,SPEC.fs);
  
  %Smooth the Envelope
  STIM(cnt).IAsm = smooth(STIM(cnt).IA,.1,'loess');
  STIM(cnt).IAT  = [1:length(STIM(cnt).IAsm)]./stimuli(cnt).fs;
  STIM(cnt).IAdf = diff(STIM(cnt).IAT);
  
  %Get MAX envelope amplitude at formant
  [mxIA mxIn] = max(STIM(cnt).IAsm);
  STIM(cnt).IAmx = mxIA; STIM(cnt).IATmx = STIM(cnt).IAT(mxIn);
  
  %Slope till MAX amplitude envelope (from stimulus onset)
  STIM(cnt).IASmxO = (STIM(cnt).IAmx-STIM(cnt).IAsm(1))/(STIM(cnt).IATmx-STIM(cnt).IAT(1));
  
  % Get formant tracks using LPC (Sound Processing Toolbox by Naotoshi Seo).
  % Resample so Nyquist = max frequency before LPC
  % (this means sample rate should be 2*SPEC.mf)
  y_rs = resample(stimuli(cnt).y,SPEC.p,SPEC.q);
  [STIM(cnt).FT, STIM(cnt).FTt] = spFormantsTrackLpc(y_rs,2*SPEC.mf,SPEC.po,...
   SPEC.fl,SPEC.fo,SPEC.wintype,0);
  
  % Get pitch track by using autocorrelation method
  % (Sound Processing Toolbox by Naotoshi Seo)
  [STIM(cnt).F0, STIM(cnt).F0t, STIM(cnt).F0r] = spPitchTrackCorr(stimuli(cnt).y,...
   SPEC.fs,SPEC.fl,SPEC.fo,[],0);
  STIM(cnt).name = [char(manis(mani)) num2str(stim)];
  
  clear mxIA mxIn y_rs
  
 end % for stims
end % for mani

clear cnt mani manis stim stims

save([datPath,'Hasselman2104_stimfeatures.mat'],'SPEC','STIM','-append');

%% Get Rise and Fall Time Entropy using MIR toolbox (MIR will smooth and resample Hilbert transform)

for cnt = 1:40
 
 % Create miraudio
 [mirAFile] = miraudio(stimuli(cnt).y);
 
 % Get MIRenvelope (Hilbert transform)
 [mirENV] = mirenvelope(mirAFile);
 sENV     = mirgetdata(mirENV);
 t        = (1:length(stimuli(cnt).y))./stimuli(cnt).fs;
 
 dsENV = derivative(sENV);
 lENV  = length(dsENV);
 
 % Resample time vector
 step = floor(length(stimuli(cnt).y)/lENV);
 tx   = prep(decimate(t,step),lENV);
 
 % In figure 3 the envelope is exaggerated (scaled)
 dsENV = dsENV.*150;
 
 % Get the zero crossings
 [~,t0,~] = crossing(dsENV,tx);
 t0 = [0 t0];
 RTent(cnt) = entropy(nonzeros(sort(round(diff(t0)*1000),'ascend')));
 
 clear dsENV sENV lENV step t0 t tx mirAFile mirENV
 
end

clear cnt

save([datPath,'Hasselman2104_stimfeatures.mat'],'RTent','-append');

%% Create TS with equal samples, normalized to [-1 1] for HNR and RQA

for i = 1:40
 
 fs = stimuli(i).fs;
 x  = stimuli(i).y;
 
 % Find transition part
 indb = find(x.^2>=.2,1,'first');
 xlr  = flipud(x);
 inde = find(xlr.^2>=.2,1,'first');
 xc   = x(indb:(length(x)-inde)); ts = indb/fs;
 t    = ([1:length(xc)]'./fs)+ts;
 
 % Create TS of size 4096
 step = floor(length(xc)/4096);
 xx   = prep(decimate(xc,step),4096);
 tx   = prep(decimate(t,step),4096);
 
 %Normalize waveform to [-1 1]
 rng=(max(xx)-min(xx));
 mid=(max(xx)+min(xx))/2;
 xu = (xx-(mid))/(rng/2);
 
 %Store into structure rpTS
 rpTS(i).ts = [tx xu];
 
 clear t ts tx fs x xc xx xn xu xlr en b step indb inde rng mid j
 
end

clear i

save([datPath,'Hasselman2104_stimfeatures.mat'],'rpTS','-append');

%% Calculate HNR based on resampled waveforms

for i= 1:1
 
 %Create miraudio (MIR toolbox University of Jyv?skyl?)
 [mirAFile] = miraudio(rpTS(i).ts(:,2));
 [mirSFile] = mirspectrum(mirAFile,'Frame',SPEC.ts,'Min',...
  SPEC.f(1),'Max',max(SPEC.f),'Window',SPEC.wintype,'Res',SPEC.sf);
 
 % Get Pitch track by using MIR toolbox (also ac). Pitch floor
 % and ceiling frequencies are Praat defaults.
 [mirF0] = mirpitch(mirAFile,'Mono','Min',75,'Max',600);
 HNR(i).F0mir = mirgetdata(mirF0);
 
 % Get inharmonicity based on pitch estimates
 % mirSTATS(i).F0acm = median(mirSTATS(i).F0);
 [hnr] = mirinharmonicity(mirAFile,'f0',mirF0);
 HNR(i).HNR = mirgetdata(hnr);
 
 clear mirF0 mirAFile mirSFile hnr;
 
end

clear i

save([datPath,'Hasselman2104_stimfeatures.mat'],'HNR','-append');


%% Calculate RECURRENCE MEASURES
% !!!!!WARNING!!!!!!
% This will take a long time to compute!

rpSTATS = zeros(40,14);

%RQA settings (based on Mutual Information and Nearest Neigbour analysis)
tau=6;m=3; e=.01;W=[];WS=[];LMIN=2;VMIN=2;TW=0;thr= 'rr';

for i=1:40
 
 % Get RP for plotting and threshold of Fixed RR using an adaptation of the
 % crp.m code in Marwan's toolbox. The edited version is included in the zip
 % file available at GitHub as crp_edit.m. Copy crp_edit.m to the CRP toolbox folder
 % Uncomment if you want the RP matrix and threshold
 % [rpMTRX(i).rp rpSTATS(i,1)] = crp_edit(rpTS(i).ts(:,2),m,tau,e,thr,'nonormalize','silent');
 
 %RQA measures
 rpSTATS(i,2:14) = crqa(rpTS(i).ts(:,2),m,tau,e,W,WS,LMIN,VMIN,TW,thr,'nonormalize','nogui');
 
end

clear tau thr e i m W WS LMIN VMIN TW

save([datPath,'Hasselman2104_stimfeatures.mat'],'rpSTATS','-append');

%%  Get MULTIFRACTAL SPECTRUM
% Makes use of the Multifractal Detrended Fluctuation Analysis code MFDFA1.m by E. Ihlen.
% The code is availale here: http://www.ntnu.edu/inm/geri/software
%
% An article descibing the MFDFA technique and how to use the code:
% Ihlen(2012). Introduction to MFDFA, FrontiersIn Fractal Physiology, 3(141), 1-18.
% http://www.ntnu.edu/documents/170234/1315232/Introduction_to_MFDFA.pdf

scmin=6;
scmax=12;
ressc=30;

scale=round(2.^[scmin:((scmax-scmin)/ressc):scmax]);

qmin=-10;
qmax=10;
qres=101;

qq = linspace(qmin,qmax,qres);

m=1;

for i = 1:40
 
 stimuliMF(i).y = stimuli(i).y;
 stimuliMF(i).t = [0:(1/stimuli(i).fs):(length(stimuliMF(i).y)-1)/stimuli(i).fs]';
 
 [stimuliMF(i).IA, stimuliMF(i).IF] = hilbert2(stimuliMF(i).y,stimuli(i).fs);
 
 [mf(i).HqE,mf(i).Hq,mf(i).tq,mf(i).hq,mf(i).Dq,mf(i).Fq]=MFDFA1(stimuliMF(i).IA,scale,qq,m,0);
 %plot(log2(scale),log2(mf(i).Fq(qq==1,:)./scale))
 
end

save([datPat,'Hasselman2104_stimfeatures.mat'],'mf','stimuliMF','-append');
%% Get the residual error Table S1

for i = 1:40
 for j=1:101
  nr(i,j) = mf(1,i).HqE(1,j).normr;
 end
end

outnr=[mean(nr(1:10,:),2), std(nr(1:10,:),1,2),mean(nr(11:20,:),2), std(nr(11:20,:),1,2),mean(nr(21:30,:),2), std(nr(21:30,:),1,2), mean(nr(31:40,:),2), std(nr(31:40,:),1,2)];
