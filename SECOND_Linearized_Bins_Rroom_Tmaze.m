%% Run this script section-by-section

dir='X:\08. Lab personnel\Current\David\Projects\Ephys\HC Modulation - Re_Vs_dcMEC LFP\2. Output\Ephys\Wasnt\7';
load(strcat(dir,'\HC10_TaskPhases7_Binned.mat'));

%user defined input needed
SLFP=HC_SLFP; %change based on name of signal matrix
%% Linearize the bins for left
%ordinalize the stembins
for i=1:length(stembins);
    LBins_Lin(i)=stembins(i);
end
%add the CP bin
for i=length(stembins)+1:length(stembins)+length(CPbin);
    LBins_Lin(i)=CPbin;
end

%save the length of array before starting to cycle in reward arm
j=length(LBins_Lin); 

%add ordinalized reward arm bins
for i=(length(LBins_Lin)+1):(length(LBins_Lin)+length(Lrbins))
        LBins_Lin(i)=Lrbins(i-j);
end

% add Left Rew Zone
for i=length(LBins_Lin)+1:length(LBins_Lin)+length(RZbin(1))
    LBins_Lin(i)=RZbin(1);
end

%save the length of array before starting to cycle in reward arm
j=length(LBins_Lin); 

%add left return arm
for i=length(LBins_Lin)+1:length(LBins_Lin)+length(Lretbins)
    LBins_Lin(i)=Lretbins(i-j);
end

n_Lbins=length(LBins_Lin);

%% Linearize the bins for Right
%ordinalize the stembins
for i=1:length(stembins);
    RBins_Lin(i)=stembins(i);
end
%add the CP bin
for i=length(stembins)+1:length(stembins)+length(CPbin);
    RBins_Lin(i)=CPbin;
end

%save the length of array before starting to cycle in reward arm
j=length(RBins_Lin); 

%add ordinalized reward arm bins
for i=(length(RBins_Lin)+1):(length(RBins_Lin)+length(Rrbins));
        RBins_Lin(i)=Rrbins(i-j);
end

% add Left Rew Zone
for i=length(RBins_Lin)+1:length(RBins_Lin)+length(RZbin(2));
    RBins_Lin(i)=RZbin(2);
end

%save the length of array before starting to cycle in reward arm
j=length(RBins_Lin); 

%add left return arm
for i=length(RBins_Lin)+1:length(RBins_Lin)+length(Rretbins);
    RBins_Lin(i)=Rretbins(i-j);
end

n_Rbins=length(RBins_Lin);

%% index the timestamps of bin occupancies

%%%%% user adjusted input required %%%%%

%generates collection of all timestamps associated with X,Y coordinates for
%each bin. The product is a collection of times for the occupancy of each
%bin, without regard for traversal

for i=1:length(LBins_Lin); %for the number of bins in the stem and CP
LBinX(i,:)=(pos_x>LBins_Lin(i).Xmin & pos_x<LBins_Lin(i).Xmax); %return boolean for case in which pos_x is within bin's X-range
LBinY(i,:)=(pos_y>LBins_Lin(i).Ymin & pos_y<LBins_Lin(i).Ymax); %return boolean for cases in which pos_y is within bin's Y-range
LBinFlag(i,:)=(LBinX(i,:)==1 & LBinY(i,:)==1); %return boolean for cases in which both above booleans are both true "TRUE"

RBinX(i,:)=(pos_x>RBins_Lin(i).Xmin & pos_x<RBins_Lin(i).Xmax); %return boolean for case in which pos_x is within bin's X-range
RBinY(i,:)=(pos_y>RBins_Lin(i).Ymin & pos_y<RBins_Lin(i).Ymax); %return boolean for cases in which pos_y is within bin's Y-range
RBinFlag(i,:)=(RBinX(i,:)==1 & RBinY(i,:)==1); %return boolean for cases in which both above booleans are both true "TRUE"

Bin(i).Left=strfind(LBinFlag(i,:),1); %index cell numbers for above cases
BinStamps(i).Left=pos_t((Bin(i).Left));%index pos_t values for those cases
Bin(i).Right=strfind(RBinFlag(i,:),1); %index cell numbers for above cases
BinStamps(i).Right=pos_t(Bin(i).Right);%index pos_t values for those cases

end

% Stem/CP bins between right and light traversals overlap. This pulls out the
% stem/CP bins into a separate variable and groups the rest together by
% Left/Right arms

for i=1:2 %adjust if needed to meet different quantity of stemp/CP bins
    StemCPBins(i).stem=BinStamps(i).Left;
end

for i=1:3 %adjust if needed to meet different quantity of reward/return bins
    RestBins(i).Left=BinStamps(i+length(StemCPBins)).Left;
    RestBins(i).Right=BinStamps(i+length(StemCPBins)).Right;
end

%% Segment Binned Timestamps to Sample and Choice Traversals

% This section goes back through each bin and separates the bin occupancy
% by the sample and choice traversals. It does this by pulling out the
% binned timstamps that fit within the Int_Sample and Int_Choice start/end
% parameters. 
for i=1:length(StemCPBins);
    for j=1:length(SLFP);
        SCPBins(i).Sample{j,:}=StemCPBins(i).stem(StemCPBins(i).stem>Int_Sample(j,1) & StemCPBins(i).stem<Int_Sample(j,8));
        SCPBins(i).Choice{j,:}=StemCPBins(i).stem(StemCPBins(i).stem>Int_Choice(j,1) & StemCPBins(i).stem<Int_Choice(j,8));
    end
end

for i=1:length(RestBins);
    for j=1:length(SLFP);
        LRestBins(i).Sample{j,:}=RestBins(i).Left(RestBins(i).Left>Int_Sample(j,1) & RestBins(i).Left<Int_Sample(j,8));
        LRestBins(i).Choice{j,:}=RestBins(i).Left(RestBins(i).Left>Int_Choice(j,1) & RestBins(i).Left<Int_Choice(j,8));
        RRestBins(i).Sample{j,:}=RestBins(i).Right(RestBins(i).Right>Int_Sample(j,1) & RestBins(i).Right<Int_Sample(j,8));
        RRestBins(i).Choice{j,:}=RestBins(i).Right(RestBins(i).Right>Int_Choice(j,1) & RestBins(i).Right<Int_Choice(j,8));
    end
end

%% Invert formatting to have trials in struc array house the bins in cell array
%STEM
for j=1:length(SLFP); %adjust for different number of trials
for i=1:length(SCPBins);
        Stembins(j).Sample{i,:}=SCPBins(i).Sample{j,:};
        Stembins(j).Choice{i,:}=SCPBins(i).Choice{j,:};
end
end

%Right
for j=1:length(SLFP); %adjust for different number of trials
for i=1:length(RRestBins);
        RightBins(j).Sample{i,:}=RRestBins(i).Sample{j,:};
        RightBins(j).Choice{i,:}=RRestBins(i).Choice{j,:};
end
end

%Left
for j=1:length(SLFP); %adjust for different number of trials
for i=1:length(LRestBins);
        LeftBins(j).Sample{i,:}=LRestBins(i).Sample{j,:};
        LeftBins(j).Choice{i,:}=LRestBins(i).Choice{j,:};
end
end

%% cleanup
clearvars -except Stembins RightBins LeftBins CTime CLFP SLFP DLFP DTime GBbins Int_Choice Int_Sample params pos_t pos_x pos_y STime Times;

%% Manually empty lines from Bins structure array that have empty cells

%% save
cd 'X:\08. Lab personnel\Current\David\Projects\Ephys\HC Modulation - Re_Vs_dcMEC LFP\2. Output\Ephys\Wasnt\7';
save ('HC10_TaskPhases7_Binned.mat','-v7.3');
