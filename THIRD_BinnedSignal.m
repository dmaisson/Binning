%% Bin the signal

%% load fileclear; clc;
clear; clc;

dir='X:\08. Lab personnel\Current\David\Projects\Ephys\HC Modulation - Re_Vs_dcMEC LFP\2. Output\Ephys\Wasnt\7';
load(strcat(dir,'\HC10_TaskPhases7_Binned'));

% user adjusted input
% SLFP=MEC_SLFP;
% CLFP=MEC_CLFP;
% STime=MEC_STime;
% CTime=MEC_CTime;
%% Mark start and end of each bin per trial

% Column 2 within each level of the bins struct ure array will be the first
% timestamp for each bin. Column 3 will be the last.

%Stem
for j=1:length(Stembins);
for i=1:length(Stembins(j).Sample);
    Stembins(j).Sample{i,2}=Stembins(j).Sample{i,1}(1,1);
    Stembins(j).Choice{i,2}=Stembins(j).Choice{i,1}(1,1);
end
end

for j=1:length(Stembins);
for i=1:length(Stembins(j).Sample);
    Stembins(j).Sample{i,3}=Stembins(j).Sample{i,1}(1,end);
    Stembins(j).Choice{i,3}=Stembins(j).Choice{i,1}(1,end);
end
end

% Right
for j=1:length(RightBins);
for i=1:length(RightBins(j).Sample);
    if isempty(RightBins(j).Sample)==1
        RightBins(j).Sample=[];
    else
        RightBins(j).Sample{i,2}=RightBins(j).Sample{i,1}(1,1);
    end
end
end

for j=1:length(RightBins);
for i=1:length(RightBins(j).Sample);
    if isempty(RightBins(j).Sample)==1
        RightBins(j).Sample=[];
    else
        RightBins(j).Sample{i,3}=RightBins(j).Sample{i,1}(1,end);
    end
end
end

for j=1:length(RightBins);
for i=1:length(RightBins(j).Choice);
    if isempty(RightBins(j).Choice)==1
        RightBins(j).Choice=[];
    else
        RightBins(j).Choice{i,2}=RightBins(j).Choice{i,1}(1,1);
    end
end
end

for j=1:length(RightBins);
for i=1:length(RightBins(j).Choice);
    if isempty(RightBins(j).Choice)==1
        RightBins(j).Choice=[];
    else
        RightBins(j).Choice{i,3}=RightBins(j).Choice{i,1}(1,end);
    end
end
end

% Left
for j=1:length(LeftBins);
for i=1:length(LeftBins(j).Sample);
    if isempty(LeftBins(j).Sample)==1
        LeftBins(j).Sample=[];
    else
        LeftBins(j).Sample{i,2}=LeftBins(j).Sample{i,1}(1,1);
    end
end
end

for j=1:length(LeftBins);
for i=1:length(LeftBins(j).Sample);
    if isempty(LeftBins(j).Sample)==1
        LeftBins(j).Sample=[];
    else
        LeftBins(j).Sample{i,3}=LeftBins(j).Sample{i,1}(1,end);
    end
end
end

for j=1:length(LeftBins);
for i=1:length(LeftBins(j).Choice);
    if isempty(LeftBins(j).Choice)==1
        LeftBins(j).Choice=[];
    else
        LeftBins(j).Choice{i,2}=LeftBins(j).Choice{i,1}(1,1);
    end
end
end

for j=1:length(LeftBins);
for i=1:length(LeftBins(j).Choice);
    if isempty(LeftBins(j).Choice)==1
        LeftBins(j).Choice=[];
    else
        LeftBins(j).Choice{i,3}=LeftBins(j).Choice{i,1}(1,end);
    end
end
end

%% Index Signal for timestamps between start and end of each bin in each trial
%Stem - Sample
for i=1:length(Stembins) %for each trial (24)
    for j=1:size(Stembins(i).Sample,1) %for each bin 
        SBinnedLFP(i).Stem{j,:}=SLFP(i).trav(STime(i).trav>Stembins(i).Sample{j,2} & STime(i).trav<Stembins(i).Sample{j,3});
    end
end

%Stem - Choice
for i=1:length(Stembins) %for each trial (24)
    for j=1:size(Stembins(i).Choice,1)
        CBinnedLFP(i).Stem{j,:}=CLFP(i).trav(CTime(i).trav>Stembins(i).Choice{j,2} & CTime(i).trav<Stembins(i).Choice{j,3});
    end
end

%Right - Sample
for i=1:length(RightBins) %for each trial (24)
    for j=1:size(RightBins(i).Sample,1)
        if isempty(RightBins(i).Sample)==1
            SBinnedLFP(i).Right{j,:}=[];
        else
            SBinnedLFP(i).Right{j,:}=SLFP(i).trav(STime(i).trav>RightBins(i).Sample{j,2} & STime(i).trav<RightBins(i).Sample{j,3});
        end
    end
end

%Left - Sample
for i=1:length(LeftBins)
    for j=1:size(LeftBins(i).Sample,1)
        if isempty(LeftBins(i).Sample)==1
            SBinnedLFP(i).Left{j,:}=[];
        else
            SBinnedLFP(i).Left{j,:}=SLFP(i).trav(STime(i).trav>LeftBins(i).Sample{j,2} & STime(i).trav<LeftBins(i).Sample{j,3});
        end
    end
end

%Right - Choice
for i=1:length(RightBins)
    for j=1:size(RightBins(i).Choice,1)
        if isempty(RightBins(i).Choice)==1
            CBinnedLFP(i).Right{j,:}=[];
        else
            CBinnedLFP(i).Right{j,:}=CLFP(i).trav(CTime(i).trav>RightBins(i).Choice{j,2} & CTime(i).trav<RightBins(i).Choice{j,3});
        end
    end
end

%Left - Choice
for i=1:length(LeftBins)
    for j=1:size(LeftBins(i).Choice,1)
        if isempty(LeftBins(i).Choice)==1
            CBinnedLFP(i).Left{j,:}=[];
        else
            CBinnedLFP(i).Left{j,:}=CLFP(i).trav(CTime(i).trav>LeftBins(i).Choice{j,2} & CTime(i).trav<LeftBins(i).Choice{j,3});
        end
    end
end

clear CLFP SLFP;
%% Reorganize
% Render cell arrays of Trial(y) X Bin(x) structure for raw signal

% Sample
% Stem Bins
for i=1:length(SBinnedLFP)
    for j=1:length(SBinnedLFP(i).Stem)
        SLFP{i,j}=SBinnedLFP(i).Stem{j,1};
    end
end

% Right Bins
for i=1:length(SBinnedLFP)
    for j=1:length(SBinnedLFP(i).Right)
        SLFP{i,(length(SBinnedLFP(i).Stem)+j)}=SBinnedLFP(i).Right{j,1};
    end
end

% Left Bins
for i=1:length(SBinnedLFP)
    for j=1:length(SBinnedLFP(i).Left)
        SLFP{i,(length(SBinnedLFP(i).Stem)+j)}=SBinnedLFP(i).Left{j,1};
    end
end

% Choice
% Stem Bins
for i=1:length(CBinnedLFP)
    for j=1:length(CBinnedLFP(i).Stem)
        CLFP{i,j}=CBinnedLFP(i).Stem{j,1};
    end
end

% Right Bins
for i=1:length(CBinnedLFP)
    for j=1:length(CBinnedLFP(i).Right)
        CLFP{i,(length(CBinnedLFP(i).Stem)+j)}=CBinnedLFP(i).Right{j,1};
    end
end

% Left Bins
for i=1:length(CBinnedLFP)
    for j=1:length(CBinnedLFP(i).Left)
        CLFP{i,(length(CBinnedLFP(i).Stem)+j)}=CBinnedLFP(i).Left{j,1};
    end
end

%% Save matrices
% use if you are running this on multiple signals.
% Adjust all for each signal it's run on before moving to Clean-up

SLFP_HC=SLFP;
CLFP_HC=CLFP;
STime_HC=STime;
CTime_HC=CTime;

%% Clean-up
clearvars -except CLFP_HC CTime_HC DLFP DTime GBbins params SLFP_HC STime_HC

%% Save
cd 'X:\08. Lab personnel\Current\David\Projects\Ephys\HC Modulation - Re_Vs_dcMEC LFP\2. Output\Ephys\Wasnt\7';
save ('HC10_Signal_BinPhase.mat','-v7.3');