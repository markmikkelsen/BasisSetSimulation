clear all;
close all;
% run_simSemiLASERShaped_fast.m
% Jamie Near, McGill University 2015.
% Dana Goerzen and Jamie Near, McGill University 2021.
% Fast version by Muhammad G Saleh (Johns Hopkins University School of Medicine, 2019)

% DESCRIPTION & MODIFICATIONS:
% This script was modified by Niklaus Zölch & Jessica Archibald to include: 

%     a.	A 0 reference peak-  useful depending on the fiting software
%     b.	A loop to run through the selected metabolites
%     c.	Saving .raw, .png, and .mat files
%     d.	Output a .pdf and .basis file -> using modified functions from Osprey- Georg Oeltzschner, Johns Hopkins University 2019.

% USAGE:
% This script simulates a semi-LASER  experiment with fully shaped refocusing 
% pulses. Coherence order filtering is employed to only simulate desired signals
% This results in a 4x speed up compared to phase cycling 
% (see deprecated run_simSemiLASERShaped_fast_phCyc.m)
% Furthermore, simulations are run at various locations in space to account for the 
% within-voxel spatial variation of the metabolite signal.  Summation 
% across spatial positions is performed. The MATLAB parallel computing toolbox 
% (parfor loop) was used to accelerate the simulations.  Acceleration 
% is currently performed in the direction of the slice selective pulse along
% the x-direction, but this can be changed.  Up to a factor of 12 acceleration
% can be achieved using this approach. To achieve 
% faster perfomance compared to the original 'run_simSemiLASER_shaped.m' function,
% this code uses the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.  Some additional acceleration is currently performed using parfor 
% loops in both x and y directions.  To enable the use of the MATLAB
% parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes.  If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.
%
% 
% INPUTS: 
% To run this script, there is technically only one input argument:
% spinSys           = spin system to simulate 
% However, the user should also edit the following parameters as 
% desired before running the function:
% refocWaveform     = name of refocusing pulse waveform.
% refTp             = duration of refocusing pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% Bfield            = magnetic field strength [Tesla]
% lw                = linewidth of the output spectrum [Hz]
% thkX              = slice thickness of x refocusing pulse [cm]
% thkY              = slice thickness of y refocusing pulse [cm]
% fovX              = full simulation FOV in the x direction [cm]
% fovY              = full simulation FOV in the y direction [cm]
% nX                = number of spatial grid points to simulate in x-direction
% nY                = number of spatial grid points to simulate in y-direction
% taus              = vector of pulse sequence timings  [ms]
%
% OUTPUTS:
% out               = Simulation results, summed over all space.

ToolboxCheck
% ************INPUT PARAMETERS**********************************
% 
% Define the variable Basis_name at the beginning of your script
basis_name = 'lcm_gamma_new.basis' ; %keep "_gamma_" 
pathtofida='/Users/jessicaarchibald/Desktop/GLUATLAS/FID-A-master2022';
addpath(genpath(pathtofida));
folder_to_save=[pathtofida,'/GeneratedRawFiles/sLASER_JESS_2024/'];
save_result=true;
%refocWaveform='AM_REFOMAN6.pta'; %name of refocusing pulse waveform.
Waveform='oit_800_6500.pta';
flip_angle=180;
%refocWaveform='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
refTp=4.4496; %duration of refocusing pulses[ms]
Npts=4096; %number of spectral points
sw=4000; %spectral width [Hz]
lw=2; %linewidth of the output spectrum [Hz]
Bfield=3; %Magnetic field strength in [T]
thkX=2.4; %slice thickness of x refocusing pulse [cm]
thkY=2.2; %slice thickness of y refocusing pulse [cm]
%
fovX=3; %size of the full simulation Field of View in the x-direction [cm]
fovY=3; %size of the full simulation Field of View in the y-direction [cm]
%
nX=40; %Number of grid points to simulate in the x-direction
nY=40; %Number of grid points to simulate in the y-direction
%
% full voxel
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
% eff voxel
%x=linspace(-0.31,0.35,nX); %X positions to simulate [cm]
% full voxel
y=linspace(-fovY/2,fovY/2,nY);
% effvoxel
%y=linspace(-1.06,1.20,nY); %y positions to simulate [cm]%
%x=0;
%y=0;
te=32;%timing of the pulse sequence [ms]
centreFreq=2.02; %Centre frequency of MR spectrum [ppm]
%
%PhCyc1=[0 0 0 0];  %phase cycling scheme of first refocusing pulse
%PhCyc2=[0 0 90 90]; %phase cycling scheme of second refocusing pulse
%PhCyc3=[0 0 0 0]; %phase cycling scheme of third refocusing pulse
%PhCyc4=[0 90 0 90]; %phase cycling scheme of fourth refocusing pulse
%
fovX=-x(1)+x(end);
fovY=-y(1)+y(end);
% spin systems 
spinSysList={'PE', 'Asc', 'Scyllo','Glu','Cr','NAA','NAAG','PCr','GSH','Gly','Glc','GPC',...
    'PCh','Ala','Asp','GABA', 'Gln', 'Ins', 'Lac', 'Tau'};

% shift
shift_in_ppm=(4.65-centreFreq);
% ************ END OF INPUT PARAMETERS BY USER **********************************
%--------------------------------------------------------------------------
%Load RF waveform
%--------------------------------------------------------------------------
%Niklaus : set inv / macht keinen unterschied oder?
rfPulse=io_loadRFwaveform(Waveform,'inv',0);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
sysRef.J=0;
sysRef.shifts=0;
sysRef.scaleFactor=1;
sysRef.name='Ref_0ppm';
%
sysRef.centreFreq=centreFreq;
%
[ ref] = run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sysRef,flip_angle);
%
tau1=15; tau2=13;%fake timing
refjustforppmrange=sim_press(Npts,sw,Bfield,lw,sysRef,tau1,tau2);
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% shift 
% https://ch.mathworks.com/matlabcentral/newsreader/view_thread/243061 
ppm_range=ref.ppm(1)-ref.ppm(end);
ppm_per_point=ppm_range/size(ref.ppm,2);
shift_in_points=round(shift_in_ppm/ppm_per_point);
ref.fids=ref.fids.*exp(-1i*2*pi*shift_in_points*(0:1:(size(ref.fids,1)-1)).'/(size(ref.fids,1)));
%--------------------------------------------------------------------------
% Additional Metabolites
[sysETH,sysAcetate,sysAcac,sysSucc,sysGlyc,sysVal,sysAceton,sysbHBHM]=define_spin_systems();

%-------------------------------------------------------------------------
%Load spin systems (for the rest)
load([pathtofida,'/simulationTools/metabolites/spinSystems'])
%-------------------------------------------------------------------------

for met_nr=1:size(spinSysList,2);
    %
    spinSys=spinSysList{met_nr}; %spin system to simulate
    sys=eval(['sys' spinSys]);
    % Schreibe die einfach im ersten rein
    sys(1).centreFreq=centreFreq;
    
    %-------------------------------------------------------------------------
    % Simulation
    %-------------------------------------------------------------------------
    [ out] = run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sys,flip_angle);
    % save out as
    % Save before the shift -
    save_out_mat=[folder_to_save,'matfiles_pre'];
    if (exist(save_out_mat,'dir')==0)
                 mkdir(save_out_mat);
    end
    save([save_out_mat,'/',spinSys],'out')

    %----------------------------------------------------------------------
    % effective voxel size
    %----------------------------------------------------------------------
    %-------------------------------------------------------------------------
    % Add shift here and later for every simulated 
    %
    %-------------------------------------------------------------------------
    out.fids=out.fids.*exp(-1i*2*pi*shift_in_points*(0:1:(size(out.fids,1)-1)).'/(size(out.fids,1)));
    %
    %effvox.fids=effvox.fids.*exp(-1i*2*pi*shift_in_points*(0:1:(size(effvox.fids,1)-1)).'/(size(effvox.fids,1)));
    %
    %-------------------------------------------------------------------------
    % add tms ref
    %-------------------------------------------------------------------------
    out=op_addScans(out,ref);
    %effvox=op_addScans(effvox,ref);
    %
    save_figure=[folder_to_save,'figures'];
    if (exist(save_figure,'dir')==0)
                 mkdir(save_figure);
    end
    % figure
    figure;plot(refjustforppmrange.ppm,real(ifftshift(ifft(out.fids))),'b')
    set(gca,'xdir','reverse')
    colormap;set(gcf,'color','w');
    xlim([-1 5])
    xlabel('ppm');
    title(['figure with ref',spinSys])
    print('-dpng','-r300',[save_figure,'/',spinSys])
    % effvoxel
    %figure;plot(ref.ppm,real(ifftshift(ifft(effvox.fids))),'b')
    %set(gca,'xdir','reverse')
    %title(['figure effvox with ref ppm ',spinSys])
    %print('-dpng','-r300',[save_figure,'\effvox_',spinSys])
    out.name=spinSys;
    out.centreFreq=centreFreq; % This is needed for the check within fit_LCMmakeBasis.
    save_raw=[folder_to_save,'raw'];

    if save_result
        if (exist(save_raw,'dir')==0)
                 mkdir(save_raw);
        end

        RF=io_writelcmraw(out,[save_raw, '/', spinSys '.raw'],spinSys);
       
    end
 % Saving after shift   
    save_out_mat_end=[folder_to_save,'matfiles_post'];
    if (exist(save_out_mat_end,'dir')==0)
                 mkdir(save_out_mat_end);
    end
    save([save_out_mat_end,'/',spinSys],'out')
    
end
disp('Running fit_makeLCMBasis...');

BASIS=fit_makeLCMBasis_2Jess(save_out_mat_end, false, [folder_to_save,'/', basis_name],'Philips','sLASER');
% %Vizualize your created basis set 
% figure;plot(BASIS.ppm,real(BASIS.specs));legend(BASIS.name)
% set(gca,'xdir','reverse','XGrid','on')
% so schon 

rmpath(genpath(pathtofida));
disp('Done');

       