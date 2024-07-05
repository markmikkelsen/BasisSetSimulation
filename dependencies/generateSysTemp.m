clear;
close all;
% Angepasst fuer FID-A 2
%-------------------------------------------------------------------------
% Edit here
%--------------------------------------------------------------------------
act_temp=10;
%-------------------------------------------------------------------------
zielpfad_base='C:\Users\admin\CloudStation\Programme\FID-A-1.2\simulationTools';
zielpfad=[zielpfad_base,'\metabolites_T_',sprintf('%d',act_temp)];
if (exist(zielpfad,'dir')==0)
            mkdir(zielpfad);
else
    rmdir(zielpfad,'s')
    mkdir(zielpfad);
end
    
%--------------------------------------------------------------------------
smaller_than_ref=-1;
greater_than_ref=1;

ref_temp=37;
diff_temp=act_temp-ref_temp;
%--------------------------------------------------------------------------
load('C:\Users\admin\CloudStation\Programme\FID-A-1.2\simulationTools\metabolites\spinSystems')
%--------------------------------------------------------------------------
slopeNAA_1=2.5070e-4; %ppm/K
%slopeNAA_2=-0.0234e-4; %ppm/K
slopeNAA_3=-5.999e-4; %ppm/K
%
sysNAA(2).shifts=[7.8205;4.3817+greater_than_ref*slopeNAA_1*diff_temp;2.6727+greater_than_ref*slopeNAA_3*diff_temp;2.4863];
%
%--------------------------------------------------------------------------
slopeAla=1.5402e-4; %ppm/K
sysAla.shifts=[1.4667+smaller_than_ref*slopeAla*diff_temp; ...
    1.4667+smaller_than_ref*slopeAla*diff_temp; ...
    1.4667+smaller_than_ref*slopeAla*diff_temp; ...
    3.7746];
%--------------------------------------------------------------------------
%sys GABA
%--------------------------------------------------------------------------
% nach neustem Govin
%https://onlinelibrary.wiley.com/doi/epdf/10.1002/nbm.3336
sysGABA.shifts=[2.2840;2.2840;1.889;1.889;3.0128;3.0128];
sysGABA.J=[0,-10.744,7.775,6.173,0,0; ...
    0,0,7.432,7.933,0,0; ...
    0,0,0,-13.121,5.372,10.578; ...
    0,0,0,0,7.127,6.982; ...
    0,0,0,0,0,-12.021; ...
    0,0,0,0,0,0];
slopeGABA=7.2800e-4; %ppm/K
sysGABA.shifts=[2.2840;2.2840;1.889;1.889;3.0128;3.0128+greater_than_ref*slopeGABA*diff_temp];
%--------------------------------------------------------------------------
%Asp
%neglegted%slopeAsp_1=-0.0691e-4; %ppm/K 
slopeAsp_2=4.1970e-4; %ppm/K
sysAsp.shifts=[3.8914;2.8011; 2.6533+smaller_than_ref*slopeAsp_2*diff_temp];
%--------------------------------------------------------------------------
% Cr
slopeCr=-6.2166e-4; %ppm/K
sysCr.shifts=[6.6490;3.9130+greater_than_ref*slopeCr*diff_temp;3.9130+greater_than_ref*slopeCr*diff_temp;3.0270;3.0270;3.0270];
%--------------------------------------------------------------------------
slopePCr=-6.6944e-4; %ppm/K
sysPCr(1).shifts=[3.9300+greater_than_ref*slopePCr*diff_temp; ...
    3.9300+greater_than_ref*slopePCr*diff_temp; ...
    6.5810;7.2960];
%--------------------------------------------------------------------------
% Glu
%--------------------------------------------------------------------------
slopeGlu_1=-3.4740e-4; %ppm/K
slopeGlu_2=4.8710e-4; %ppm/K
slopeGlu_3=-0.3854e-4; %ppm/K
slopeGlu_4=-3.4450e-4; %ppm/K
%
sysGlu.shifts=[3.7433+greater_than_ref*slopeGlu_1*diff_temp; ...
               2.0375+smaller_than_ref*slopeGlu_2*diff_temp; ...
               2.1200+smaller_than_ref*slopeGlu_3*diff_temp; ...
               2.3378+smaller_than_ref*slopeGlu_4*diff_temp; ...
               2.3520];
%--------------------------------------------------------------------------
slopeGln_1=-0.5810e-4; %ppm/K
slopeGln_2=3.4090e-4; %ppm/K
slopeGln_3=7.6150e-4; %ppm/K
slopeGln_4=2.9690e-4; %ppm/K
sysGln(1).shifts=[3.7530+greater_than_ref*slopeGln_1*diff_temp; ...
    2.1290+greater_than_ref*slopeGln_2*diff_temp; ...
    2.1090; ...
    2.4320+greater_than_ref*slopeGln_3*diff_temp; ... 
    2.4540+greater_than_ref*slopeGln_4*diff_temp];
%--------------------------------------------------------------------------
slopeIns_1=-3.1710e-4; %ppm/K
slopeIns_2=3.1430e-4; %ppm/K
slopeIns_3=3.5590e-4; %ppm/K^ 
sysIns.shifts=[3.5217; ...
               4.0538+greater_than_ref*slopeIns_3*diff_temp; ...
               3.5217; ...
               3.6144+greater_than_ref*slopeIns_2*diff_temp; ...
               3.2690+smaller_than_ref*slopeIns_1*diff_temp; ...
               3.6144+greater_than_ref*slopeIns_2*diff_temp];
%--------------------------------------------------------------------------
slopeLac=-3.1866e-4; %ppm/K
sysLac.shifts=[4.0974+greater_than_ref*slopeLac*diff_temp; ...
    1.3142;1.3142;1.3142];
%--------------------------------------------------------------------------
slopeTau=3.2098e-4; %ppm/K
sysTau.shifts=[ 3.4206;3.4206;3.2459+smaller_than_ref*slopeTau*diff_temp;3.2459+smaller_than_ref*slopeTau*diff_temp];
%--------------------------------------------------------------------------

save([zielpfad,'\spinSystems'],'sys*') 
