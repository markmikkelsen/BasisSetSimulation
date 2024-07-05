function [ sysETH,sysAcetate,sysAcac,sysSucc,sysGlyc,sysVal,sysAceton,sysbHBHM] = define_spin_systems( )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Ethanol
sysETH.shifts=[1.19;1.19;1.19;3.67;3.67];
sysETH.J=[0,0,0,7.1,7.1; ...
    0,0,0,7.1,7.1; ...
    0,0,0,7.1,7.1; ...
    0,0,0,0,0; ...
    0,0,0,0,0];
sysETH.scaleFactor=1;
sysETH.name='Ethanol';

% Acetate
sysAcetate.shifts=1.904;
sysAcetate.J=0;
sysAcetate.scaleFactor=1;
sysAcetate.name='Acetate';

% Acetoacetate
sysAcac(1).shifts=2.27;
sysAcac(1).J=0;
sysAcac(1).scaleFactor=3;
sysAcac(1).name='AcAc_1';
%
sysAcac(2).shifts=3.43;
sysAcac(2).J=0;
sysAcac(2).scaleFactor=2;
sysAcac(2).name='AcAc_2';
%
% Succinate
sysSucc.shifts=2.3920;
sysSucc.J=0;
sysSucc.scaleFactor=4;
sysSucc.name='Succinate';

% Glycerol
sysGlyc.shifts=[3.5522; 3.6402; 3.7704; 3.5522; 3.6402];

sysGlyc.J=[0,-11.715,4.427,0,0; ...
           0, 0,6.485, 0,0; ...
           0, 0, 0,4.427,6.485; ...
           0, 0, 0, 0, -11.715; ...
           0, 0, 0, 0, 0];
sysGlyc.scaleFactor=1;
sysGlyc.name='Glycerol';


% Valine
sysVal.shifts=[3.5953; 2.2577; 1.0271; 1.0271; 1.0271;0.9764;0.9764;0.9764];

sysVal.J=[0,4.405,0,0,0,0,0,0; ...
          0, 0, 6.971, 6.971, 6.971, 7.071, 7.071, 7.071; ...
          0, 0, 0, 0, 0, 0, 0, 0; ...
          0, 0, 0, 0, 0, 0, 0, 0; ...
          0, 0, 0, 0, 0, 0, 0, 0; ...
          0, 0, 0, 0, 0, 0, 0, 0; ...
          0, 0, 0, 0, 0, 0, 0, 0; ...
          0, 0, 0, 0, 0, 0, 0, 0];
sysVal.scaleFactor=1;
sysVal.name='Valine';

% Aceton
sysAceton.shifts=2.22;
sysAceton.J=0;
sysAceton.scaleFactor=1;
sysAceton.name='Aceton';

% bHBHM
% Based on the Human Metabolome Database
% http://www.hmdb.ca/spectra/nmr_one_d/1027
sysbHBHM.shifts=[2.402; 2.326; 4.16; 1.204; 1.204; 1.204];
sysbHBHM.J=[0,14.4,7.3,0,0,0; ...
            0,0,6.3,0,0,0; ...
            0,0,0,6.26,6.26,6.26; ...
            0,0,0,0,0,0; ...
            0,0,0,0,0,0; ...
            0,0,0,0,0,0];
sysbHBHM.scaleFactor=1;
sysbHBHM.name='bHBHM';

end

