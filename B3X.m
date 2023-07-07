function [SBP,DBP] = B3X(cuff,osc,ks,fs)
% This function implements the B³X algorithm for blood pressure estimation
% from oscillometric wave and Korotkoff sounds recordings.
% 
% The B3X algorithm uses the beat-by-beat cross-correlation (B³X) of
% adjacent Korotkoff sounds and the first derivative of B³X (B³XD) to
% locate the time instants at which the pressure of the BP brachial cuff
% reaches the systolic (SBP) and diastolic (DBP) blood pressure values.
% The SBP and DBP timings are eventually used to recover the actual SBP and
% DBP values from the cuff pressure curve.
%
% INPUTS:
% - cuff: 1D array that contains the brachial cuff pressure readings
% - osc: 1D array that contains the oscillometric wave recording
% - ks: 1D array that contains the Korotkoff sounds recording
% - fs: sampling frequency of recorded signals
%
% OUTPUTS:
% - SBP: estimated systolic blood pressure
% - DBP: estimated diastolic blood pressure


normabs = @(x) x/max(abs(x(:)));


%% PARAMETERS OF SEARCH ALGORITHM
SBP_threshold = 0.02;  % perc. threshold for B³XD
DBP_threshold = 0.02;  % perc. threshold for B³XD
DBP_threshold2 = 0.02;  % perc. threshold for B³X



%% FILTERS

% LPF to filter out spurious oscillations from OSC signal
fcut_osc = 3;
[b,a] = butter(4,fcut_osc/(fs/2));

% BPF for Korotkoff sounds
fcut_KS = [30 80];
[bs,as] = butter(2,fcut_KS/(fs/2));



%% %%%%%%%%%%%%%%%%%%%%%%%%% OSC SIGNAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%% %

% OSC LP filtering
lp_osc = filtfilt(b,a,osc);                      % zero-phase filtering

% Find peaks in the LPF OSC signal
[~,wa] = findpeaks(lp_osc,'MinPeakProminence',0.04);

% Exclude peaks in the first and last 1-second intervals of the signal
wa(wa<fs*1) = [];
wa(wa>length(osc)-fs*1) = [];


%% %%%%%%%%%%%%%%%%%%%%%% KOROTKOFF SOUNDS PROCESSING %%%%%%%%%%%%%%%%%%%%%% %

% KS BP filtering
ks_filt = filtfilt(bs,as,ks);              % zero-phase filtering


% LOCALIZATION OF KOROTKOFF SOUNDS PEAKS
% For each OSC peak location, the algorithm searches for the absolute
% maximum in a window of Npre samples before and Npost samples after the
% OSC peak. 

Tpre = 500; % ms
Tpost = 20; % ms
Npre = round(Tpre*fs/1000);
Npost = round(Tpost*fs/1000);
ps = zeros(size(wa));
for j = 1:length(wa)
    win = wa(j)-Npre:wa(j)+Npost;
    [~,im] = max(ks_filt(win));
    ps(j) = win(1)+im-1;
end


% SEGMENTATION OF KOROTKOFF SOUNDS (KS)
% Segments of Nwin samples before and after each KS peak are extracted
Twin = 120; % ms
Nwin = Twin*fs/1000;

for j = 1:length(ps)
    win = ps(j)-Nwin:ps(j)+Nwin;
    ks_win(:,j) = win;
    ks_seg(:,j) = ks_filt(win);
end


% COMPUTATION OF BEAT-BY-BEAT CROSS-CORRELATION (B³X) AND ITS FIRST DERIVATIVE (B³XD)
nomean = @(x) x-mean(x);
ks_xc = zeros(length(ps)-1,1);
for k = 1:length(ps)-1
    xc = xcorr(nomean(ks_seg(:,k)),nomean(ks_seg(:,k+1)));
    ks_xc(k) = max(xc);
end
ksxd = normabs([0; diff(ks_xc)]);



 %% %%%%%%%%%%%%%%%%%%%%%% BLOOD PRESSURE ESTIMATION %%%%%%%%%%%%%%%%%%%%%% %%

%---------------------- SYSTOLIC BLOOD PRESSURE ----------------------%
% SBP corresponds to the first point of B³XD that is above the positive
% threshold XTp and is preceded by at least two points below XTp. 

% Absolute maximum of B³XD
[max_ksxd,imax1] = max(ksxd);


% Exclude unreasonable data
while((cuff(wa(imax1))<60) || (cuff(wa(imax1))> 300) || imax1<3)
    ksxd(imax1) = 0;
    [max_ksxd,imax1] = max(ksxd);
end

% Positive threshold for B³XD
XTp = SBP_threshold * max_ksxd;


% BACKWARD SEARCH FROM MAXIMUM OF B³XD
% The index of the heartbeat corresponding to SBP is found by searching
% for the last sample of B³XD (preceding the absolute maximum) that:
%    1) is below the positive threshold XTp
%    2) has at least two preceding samples below XTp


ksxd1=[0;0;ksxd];
[~,imax] = max(ksxd1);
SBPcond = (abs(ksxd1(3:imax)) < XTp) & (abs(ksxd1(2:imax-1)) < XTp) & (abs(ksxd1(1:imax-2)) < XTp);
iSP = find(SBPcond,1,'last')+1;   % location of SBP

if isempty(iSP)
    iSP = 1;
end

% Estimated SBP
SBP = cuff(wa(iSP));               



%---------------------- DIASTOLIC BLOOD PRESSURE ---------------------%
% DBP corresponds to the last point at which B³XD is below the negative
% threshold XTn and B³X is above the positive threshold Tp.  

% Absolute minimum of B³XD
[min_ksxd,imin] = min(ksxd);

% Exclude unreasonable data
while((cuff(ps(imin))<30) || (cuff(ps(imin))> 120) || ((length(ksxd)-imin)<3))
    ksxd(imin) = 0;
    [min_ksxd,imin] = min(ksxd);
end
    

% Negative threshold for B³XD
XTn = DBP_threshold * min_ksxd;  % absolute threshold

% Threshold for B³X
Tp = DBP_threshold2 * max(ks_xc); 


% MIXED FORWARD SEARCH FROM MINIMUM OF B³XD
% Starting from the absolute minimum, the index of the heartbeat
% corresponding to DBP is found by searching for the index that:
%    1) corresponds to a value of the B³XD above the negative threshold XTn 
%    2) corresponds to a value of B³X below the positive threshold Tp

DBPcond = (ksxd(imin:length(ksxd)) > XTn) & (ks_xc(imin:length(ksxd)) < Tp);
iDP = find(DBPcond, 1) + imin - 1;
if (isempty(iDP))
    iDP = length(ksxd);
end
if (iDP > length(ksxd))
    iDP = length(ksxd);
end

% Estimated DBP
DBP = cuff(wa(iDP));     

end