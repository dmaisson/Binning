%% load file after running processing script
clear; clc;

dir='X:\08. Lab personnel\Current\David\Projects\Ephys\HC Modulation - Re_Vs_dcMEC LFP\2. Output\Ephys\Wasnt\7';
load(strcat(dir,'\HC10_TaskPhases7.mat'));
%% 
% function [Rbins, Lbins] = Binmaker_f()

%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stem
numbins = 1;

%User input 
stem_xmax = 610; %Enter x coordinate for the border between the stem and choice point
stem_center =  250; %Enter y coordinate in the center of the stem

%Calculate the rest of the stem coordinates
stem_xmin = stem_xmax - 372; %218;
stem_ymin = stem_center - 20; %200; 
stem_ymax = stem_center + 20; %250;

%Check to make sure that the box is centered on the maze stem
plot(pos_y,pos_x','.'), hold on
plot([stem_ymin stem_ymax stem_ymax stem_ymin stem_ymin],[stem_xmax stem_xmax stem_xmin stem_xmin stem_xmax],'r')


binsize = (stem_xmax - stem_xmin)/numbins;

x1 = stem_xmin:binsize:stem_xmax-binsize;
x2 = stem_xmin+binsize:binsize:stem_xmax;
x3 = x2;
x4 = x1;
x5 = x1;

y1 = zeros(1,numbins) + stem_ymin;
y2 = y1;
y3 = zeros(1,numbins) + stem_ymax;
y4 = y3;
y5 = y1;

for k=1:numbins;
    stembins(k).Xmin=x1(k);
    stembins(k).Xmax=x2(k);
    stembins(k).Ymin=y1(k);
    stembins(k).Ymax=y3(k);
end

%Plotting. View the bins that you created.
for (k=1:numbins)
   plot([y1(k),y2(k),y3(k),y4(k),y1(k)],[x1(k),x2(k),x3(k),x4(k),x1(k)],'r');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Right Goal Arm

numbinRr = 1;

Rreward_xminM = stem_xmax;  % 590; Nearest to choice point
Rreward_yminM = stem_ymax; %250;
Rreward_ymaxM = Rreward_yminM+10; %260;
Rreward_xminR = 575; %556.8; Nearest to reward %%User input
Rreward_yminR = stem_ymax+155; %405;
Rreward_ymaxR = Rreward_ymaxM+155; %415;
% 
Rreward_xmaxR = Rreward_xminR+50; %606.8;
Rreward_xmaxM = stem_xmax+ 65; %640;
% 
Rreward_x = [Rreward_xminM Rreward_xmaxM Rreward_xmaxR Rreward_xminR Rreward_xminM];
Rreward_y = [Rreward_yminM Rreward_ymaxM Rreward_ymaxR Rreward_yminR Rreward_yminM];
plot(Rreward_y,Rreward_x,'g')

Rr_binsizeX = (Rreward_xminM - Rreward_xminR)/numbinRr;

x1Rr = Rreward_xminM:-Rr_binsizeX:Rreward_xminR+Rr_binsizeX;
x4Rr = Rreward_xminM-Rr_binsizeX:-Rr_binsizeX:Rreward_xminR;
x3Rr = x4Rr + 50;
x2Rr = x1Rr + 50;
x5Rr = x1Rr;

Rr_binsizeY = (Rreward_yminR - Rreward_yminM)/numbinRr;

y1Rr = Rreward_yminM:Rr_binsizeY:Rreward_yminR-Rr_binsizeY;
y2Rr = Rreward_ymaxM:Rr_binsizeY:Rreward_ymaxR-Rr_binsizeY;
y3Rr = y2Rr+Rr_binsizeY;
y4Rr = y1Rr+Rr_binsizeY;
y5Rr = y1Rr;

for k=1:numbinRr;
    Rrbins(k).Xmin=x4Rr(k);
    Rrbins(k).Xmax=x2Rr(k);
    Rrbins(k).Ymin=y1Rr(k);
    Rrbins(k).Ymax=y3Rr(k);
end

%Plot bins
for (k=1:numbinRr)
   plot([y1Rr(k),y2Rr(k),y3Rr(k),y4Rr(k),y5Rr(k)],[x1Rr(k),x2Rr(k),x3Rr(k),x4Rr(k),x5Rr(k)],'r');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left Goal Arm

numbinLr = 1;

Lreward_xminM = stem_xmax; %590;
Lreward_ymaxM = stem_ymin; %200;
Lreward_yminM = Lreward_ymaxM-10; %190;
Lreward_xminL = 585; %556.8; User input
Lreward_ymaxL = Lreward_ymaxM-155; %45;
Lreward_yminL = Lreward_yminM-155; %35;

Lreward_xmaxL = Lreward_xminL+65; %606.8;
Lreward_xmaxM = stem_xmax+50; %640;

Lreward_x = [Lreward_xminM Lreward_xmaxM Lreward_xmaxL Lreward_xminL Lreward_xminM];
Lreward_y = [Lreward_ymaxM Lreward_yminM Lreward_yminL Lreward_ymaxL Lreward_ymaxM];
plot(Lreward_y,Lreward_x,'g')

Lr_binsizeX = (Lreward_xminM - Lreward_xminL)/numbinLr;

x1Lr = Lreward_xminM:-Lr_binsizeX:Lreward_xminL+Lr_binsizeX;
x2Lr = x1Lr + 50;
x4Lr = Lreward_xminM-Lr_binsizeX:-Lr_binsizeX:Lreward_xminL;
x3Lr = x4Lr + 50;
x5Lr = x1Lr;

Lr_binsizeY = (Lreward_ymaxM - Lreward_ymaxL)/numbinLr;

y2Lr = Lreward_yminM:-Lr_binsizeY:Lreward_yminL+Lr_binsizeY;
y1Lr = Lreward_ymaxM:-Lr_binsizeY:Lreward_ymaxL+Lr_binsizeY;
y4Lr = y1Lr-Lr_binsizeY;
y3Lr = y2Lr-Lr_binsizeY;
y5Lr = y1Lr;

for k=1:numbinLr;
    Lrbins(k).Xmin=x4Lr(k);
    Lrbins(k).Xmax=x3Lr(k);
    Lrbins(k).Ymin=y4Lr(k);
    Lrbins(k).Ymax=y1Lr(k);
end

for k=1:numbinLr
   plot([y1Lr(k),y2Lr(k),y3Lr(k),y4Lr(k),y5Lr(k)],[x1Lr(k),x2Lr(k),x3Lr(k),x4Lr(k),x5Lr(k)],'r');
%pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choice point%%%
% Triangles (Turn from stem to reward arm

midTri_x = stem_xmax; %590;
midLTri_y = Lreward_ymaxM+30; %230;
midRTri_y = Rreward_ymaxM-40; %220;

L_Tri_y = Lreward_ymaxM; %200;
R_Tri_y =  stem_ymax; %250;

outMidTri_x = midTri_x+52; %642;
outMidLTri_y = L_Tri_y+20; %220;
outMidRTri_y = R_Tri_y-20; %230;

outTri_x = 640;
outLTri_y = Lreward_yminM; %190;
outRTri_y = Rreward_ymaxM; 260;


trix1L = [midTri_x midTri_x];
trix2L = [midTri_x outMidTri_x];
trix3L = [outMidTri_x outTri_x];
trix4L = trix1L;
trix5L = trix2L;

triy1L = [L_Tri_y L_Tri_y];
triy2L = [midLTri_y outMidLTri_y];
triy3L = [outMidLTri_y outLTri_y];
triy4L = triy1L;
triy5L = triy2L;

trix1R = [midTri_x midTri_x];
trix2R = [midTri_x outMidTri_x];
trix3R = [outMidTri_x outTri_x];
trix4R = trix1R;
trix5R = trix2R;

triy1R = [R_Tri_y R_Tri_y];
triy2R = [midRTri_y outMidRTri_y];
triy3R = [outMidRTri_y outRTri_y];
triy4R = triy1R;
triy5R = triy2R;

CPbin(1).Xmin = trix1L(1,1);
CPbin(1).Xmax = trix2R(1,2);
CPbin(1).Ymin = triy1L(1,1);
CPbin(1).Ymax = triy1R(1,1);

plot([L_Tri_y outMidLTri_y outLTri_y L_Tri_y],[midTri_x outMidTri_x outTri_x midTri_x],'g') %Left toward reward
plot([L_Tri_y midLTri_y outMidLTri_y L_Tri_y],[midTri_x midTri_x outMidTri_x midTri_x],'g')

plot([R_Tri_y outMidRTri_y outRTri_y R_Tri_y],[midTri_x outMidTri_x outTri_x midTri_x],'g')
plot([R_Tri_y midRTri_y outMidRTri_y R_Tri_y],[midTri_x midTri_x outMidTri_x midTri_x],'g')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left reward zone
Lrz_xA=530;
Lrz_yA=50;

Lrz_xB=630;
Lrz_yB=50;

Lrz_xC=630;
Lrz_yC=0;

Lrz_xD=530;
Lrz_yD=0;

RZbin(1).Xmin = Lrz_xA;
RZbin(1).Xmax = Lrz_xB;
RZbin(1).Ymin = Lrz_yC;
RZbin(1).Ymax = Lrz_yB;

% Right reward zone
Rrz_xA=530;
Rrz_yA=500;

Rrz_xB=630;
Rrz_yB=500;

Rrz_xC=630;
Rrz_yC=440;

Rrz_xD=530;
Rrz_yD=440;

RZbin(2).Xmin = Rrz_xA;
RZbin(2).Xmax = Rrz_xB;
RZbin(2).Ymin = Rrz_yC;
RZbin(2).Ymax = Rrz_yB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left return arm

numbinLb = 1;

Lreturn_Ax = stem_xmin-27.5; Lreturn_Ay = stem_ymin-41; %180.5;158;
Lreturn_Bx = stem_xmin; Lreturn_By = stem_ymin;%208 199

Lreturn_Cx = Lreturn_Ax+300; Lreturn_Cy = Lreturn_Ay-108; %50; %493.5;
Lreturn_Dx = 475; Lreturn_Dy = 0;

Lreturn_x = [Lreturn_Ax Lreturn_Bx Lreturn_Cx Lreturn_Dx Lreturn_Ax];
Lreturn_y = [Lreturn_Ay Lreturn_By Lreturn_Cy Lreturn_Dy Lreturn_Ay];
plot(Lreturn_y,Lreturn_x,'g')

Lb_binsizeX = (Lreturn_Dx - Lreturn_Ax)/numbinLb;

x1Lb = Lreturn_Dx-Lb_binsizeX:-Lb_binsizeX:Lreturn_Ax;
x2Lb = x1Lb + Lb_binsizeX;
x4Lb = Lreturn_Cx-Lb_binsizeX:-Lb_binsizeX:Lreturn_Bx;
x3Lb = x4Lb + Lb_binsizeX;
x5Lb = x1Lb;

Lb_binsizeY = (Lreturn_Ay - Lreturn_Dy)/numbinLb;

y1Lb = Lreturn_Dy+Lb_binsizeY:Lb_binsizeY:Lreturn_Ay;
y2Lb = y1Lb-Lb_binsizeY;
y4Lb = Lreturn_Cy+Lb_binsizeY:Lb_binsizeY:Lreturn_By;
y3Lb = y4Lb-Lb_binsizeY;
y5Lb = y1Lb;

for k=1:numbinLb;
    Lretbins(k).Xmin=x4Lb(k);
    Lretbins(k).Xmax=x3Lb(k);
    Lretbins(k).Ymin=y2Lb(k);
    Lretbins(k).Ymax=y1Lb(k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Right return arm

numbinRb = 1;


Rreturn_Ax = stem_xmin-27.5; Rreturn_Ay = stem_ymax+41; %180.5;158;
Rreturn_Bx = 475; Rreturn_By = 500;

Rreturn_Cx = Rreturn_Ax+300; Rreturn_Cy = Rreturn_Ay+108; %50; %493.5;
Rreturn_Dx = stem_xmin; Rreturn_Dy = stem_ymax;%208 199

Rreturn_x = [Rreturn_Ax Rreturn_Bx Rreturn_Cx Rreturn_Dx Rreturn_Ax];
Rreturn_y = [Rreturn_Ay Rreturn_By Rreturn_Cy Rreturn_Dy Rreturn_Ay];
plot(Rreturn_y,Rreturn_x,'g')

Rb_binsizeX = (Rreturn_Dx - Rreturn_Ax)/numbinRb;

x1Rb = Rreturn_Dx-Rb_binsizeX:-Rb_binsizeX:Rreturn_Ax;
x2Rb = x1Rb + Rb_binsizeX;
x4Rb = Rreturn_Cx-Rb_binsizeX:-Rb_binsizeX:Rreturn_Bx;
x3Rb = x4Rb + Rb_binsizeX;
x5Rb = x1Rb;

Rb_binsizeY = (Rreturn_Ay - Rreturn_Dy)/numbinRb;

y1Rb = Rreturn_Dy+Rb_binsizeY:Rb_binsizeY:Rreturn_Ay;
y2Rb = y1Rb-Rb_binsizeY; y2Rb=flipud(y2Rb');y2Rb=y2Rb';
y4Rb = Rreturn_Cy+Rb_binsizeY:Rb_binsizeY:Rreturn_By;
y3Rb = y4Rb-Rb_binsizeY; y3Rb=flipud(y3Rb');y3Rb=y3Rb';
y5Rb = y1Rb;

for k=1:numbinRb;
    Rretbins(k).Xmin=x4Lb(k);
    Rretbins(k).Xmax=x3Lb(k);
    Rretbins(k).Ymin=y2Rb(k);
    Rretbins(k).Ymax=y3Rb(k);
end
% stembins = [x1;x2;x3;x4;x5;y1;y2;y3;y4;y5];
% Lretbins = [x1Lb;x2Lb;x3Lb;x4Lb;x5Lb;y1Lb;y2Lb;y3Lb;y4Lb;y5Lb];
% Rretbins = [x1Rb;x2Rb;x3Rb;x4Rb;x5Rb;y1Rb;y2Rb;y3Rb;y4Rb;y5Rb];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GB quadrants
GB_Ax=180; GB_Ay=130;
GB_Bx=180; GB_By=350;
GB_Cx=25; GB_Cy=350;
GB_Dx=25; GB_Dy=130;

GBbins(1).Xmin=GB_Dx; GBbins(1).Xmax=((GB_Ax-GB_Cx)/2)+GB_Dx; GBbins(1).Ymin=GB_Dy; GBbins(1).Ymax=((GB_By-GB_Dy)/2)+GB_Dy;

GBbins(2).Xmin=GB_Dx; GBbins(2).Xmax=((GB_Ax-GB_Cx)/2)+GB_Dx; GBbins(2).Ymin=((GB_By-GB_Dy)/2)+GB_Dy; GBbins(2).Ymax=GB_By;

GBbins(3).Xmin=((GB_Ax-GB_Cx)/2)+GB_Dx;GBbins(3).Xmax=GB_Ax;GBbins(3).Ymin=((GB_By-GB_Dy)/2)+GB_Dy;GBbins(3).Ymax=GB_By;

GBbins(4).Xmin=((GB_Ax-GB_Cx)/2)+GB_Dx;GBbins(4).Xmax=GB_Ax;GBbins(4).Ymin=GB_Dy;GBbins(4).Ymax=((GB_By-GB_Dy)/2)+GB_Dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1R = [x1 trix1R x1Rr x1Rb];
x2R = [x2 trix2R x2Rr x2Rb];
x3R = [x3 trix3R x3Rr x3Rb];
x4R = [x4 trix4R x4Rr x4Rb];
x5R = [x5 trix5R x5Rr x5Rb];

y1R = [y1 triy1R y1Rr y1Rb];
y2R = [y2 triy2R y2Rr y2Rb];
y3R = [y3 triy3R y3Rr y3Rb];
y4R = [y4 triy4R y4Rr y4Rb];
y5R = [y5 triy5R y5Rr y5Rb];

x1L = [x1 trix1L x1Lr x1Lb];
x2L = [x2 trix2L x2Lr x2Lb];
x3L = [x3 trix3L x3Lr x3Lb];
x4L = [x4 trix4L x4Lr x4Lb];
x5L = [x5 trix5L x5Lr x5Lb];

y1L = [y1 triy1L y1Lr y1Lb];
y2L = [y2 triy2L y2Lr y2Lb];
y3L = [y3 triy3L y3Lr y3Lb];
y4L = [y4 triy4L y4Lr y4Lb];
y5L = [y5 triy5L y5Lr y5Lb];

% Rbins = [x1R; x2R; x3R; x4R; x5R; y1R; y2R; y3R; y4R; y5R];
% Lbins = [x1L; x2L; x3L; x4L; x5L; y1L; y2L; y3L; y4L; y5L];



% Plotting
% 
stem_x = [stem_xmin stem_xmax stem_xmax stem_xmin stem_xmin];
stem_y = [stem_ymax stem_ymax stem_ymin stem_ymin stem_ymax];

figure;
plot(pos_x,pos_y,'.'), set(gca,'YDir','reverse')
hold on

plot(stem_x,stem_y,'r')
hold on
plot([Lrz_xA Lrz_xB Lrz_xC Lrz_xD Lrz_xA],[Lrz_yA Lrz_yB Lrz_yC Lrz_yD Lrz_yA],'r');
plot([Rrz_xA Rrz_xB Rrz_xC Rrz_xD Rrz_xA],[Rrz_yA Rrz_yB Rrz_yC Rrz_yD Rrz_yA],'r');

for (k=1:numbins)
   plot([x1(k),x2(k),x3(k),x4(k),x1(k)],[y1(k),y2(k),y3(k),y4(k),y1(k)],'r');
end

for (k=1:numbinRr)
   plot([x1Rr(k),x2Rr(k),x3Rr(k),x4Rr(k),x5Rr(k)],[y1Rr(k),y2Rr(k),y3Rr(k),y4Rr(k),y5Rr(k)],'r');
end

for k=1:numbinLr
   plot([x1Lr(k),x2Lr(k),x3Lr(k),x4Lr(k),x5Lr(k)],[y1Lr(k),y2Lr(k),y3Lr(k),y4Lr(k),y5Lr(k)],'r');
%pause
end

plot([midTri_x outMidTri_x outTri_x midTri_x],[L_Tri_y outMidLTri_y outLTri_y L_Tri_y],'r') %Left toward reward
plot([midTri_x midTri_x outMidTri_x midTri_x],[L_Tri_y midLTri_y outMidLTri_y L_Tri_y],'r')

plot([midTri_x outMidTri_x outTri_x midTri_x],[R_Tri_y outMidRTri_y outRTri_y R_Tri_y],'r')
plot([midTri_x midTri_x outMidTri_x midTri_x],[R_Tri_y midRTri_y outMidRTri_y R_Tri_y],'r')

% for (k=1:numbinLb)
%    plot([x1Lb(k),x2Lb(k),x3Lb(k),x4Lb(k),x5Lb(k)],[y1Lb(k),y2Lb(k),y3Lb(k),y4Lb(k),y5Lb(k)],'r');
% end
% 
% for (k=1:numbinRb)
%    plot([x1Rb(k),x2Rb(k),x3Rb(k),x4Rb(k),x5Rb(k)],[y1Rb(k),y2Rb(k),y3Rb(k),y4Rb(k),y5Rb(k)],'r');
% end

%% clean-up
clearvars -except pos_x pos_y pos_t GBbins stembins RZbin Rretbins Rrbins Lretbins Lrbins CPbin CLFP CTime SLFP STime Int_Choice Int_Sample params;
