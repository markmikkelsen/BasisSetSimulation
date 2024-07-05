% adopted from
% run_simSemiLASERShaped_fast.m
% Dana Goerzen and Jamie Near, McGill University 2021.
% Fast version by Muhammad G Saleh (Johns Hopkins University School of Medicine, 2019)

% USAGE:
% out=run_simPressShaped_fast(spinSys);
% 
% DESCRIPTION:
% This script simulates a PRESS experiment with fully shaped refocusing 
% pulses. Coherence order filtering is employed to only simulate desired signals
% This results in a 4x speed up compared to phase cycling (see deprecated run_simSemiLASERShaped_fast_phCyc.m)
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
function [ out] = run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sys,flipAngle)
%function out=run_simPressShaped_fast(spinSys)

% INPUTS:
% To run this script, there is technically only one input argument:
% spinSys           = spin system to simulate 
%
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
% refPhCyc1         = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2         = vector of phase cycling steps for 2nd refocusing pulse [degrees]
%
% OUTPUTS:
% out               = Simulation results, summed over all space.
%
nX=size(x,2);
nY=size(y,2);
fovX=-x(1)+x(end);
fovY=-y(1)+y(end);

%

gamma=42577000; %gyromagnetic ratio
%Niklaus
centreFreq=sys(1).centreFreq;
%

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
rfPulse=rf_resample(rfPulse,100);
%
if ~rfPulse.isGM
    %Non-gradient modulated pulse - Calculating the x and y gradient 
    %strengths for the desired slice thickness
    Gx=(rfPulse.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
    Gy=(rfPulse.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
else
    %Gradient modulated pulse
    %1.  Calculating the unitless scaling factor for the GM waveform.
    Gx=(rfPulse.tthk/(refTp/1000))/thkX;
    Gy=(rfPulse.tthk/(refTp/1000))/thkY;
end

%Initialize structures:
% out_posxy_rpc=cell(length(x),length(y),length(ph1));
%out_posx_rpc =cell(length(x),length(PhCyc1));
out_posx_rpc =cell(length(x),1);
%d=cell(length((PhCyc1)));
d=cell(1,1);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%if you do not have the parallel computing toolbox, uncomment the first
%for loop and delete "parfor X=1:length(x)"


% parfor X=1:length(x)
% %  for X=1:length(x)
% 
%        % disp(['Executing X-position ' num2str(X) '!!' ]);
%         out_posx_rpc{X}=sim_sLASER_shaped_Ref1(Bfield,sys,te,rfPulse,refTp,x(X),Gx,flipAngle,centreFreq);
% %                            sim_sLASER_shaped_Ref1(Bfield,sys,te,RF,       tp,  dx, Gx,ph1,   ph2,  centreFreq)
% end

% Define the number of iterations
numIterations = length(x);

% Loop through the iterations
parfor X = 1:numIterations
    % Calculate the progress percentage
    progressPercentage = X / numIterations * 100;
    
    % Create a progress bar string
    progressBar = ['[', repmat('=', 1, round(progressPercentage)), repmat(' ', 1, 100 - round(progressPercentage)), ']'];
    
    % Display the progress bar in the command window
    fprintf('Progress X-position: %s %.1f%%\r', progressBar, progressPercentage);

    % Perform your computation here
    out_posx_rpc{X} = sim_sLASER_shaped_Ref1(Bfield, sys, te, rfPulse, refTp, x(X), Gx, flipAngle, centreFreq);
    % sim_sLASER_shaped_Ref1(Bfield,sys,te,RF,tp,dx,Gx,ph1,ph2,centreFreq)
end

% Print a new line to clear the progress bar
fprintf('\n');

%calculate the average density matrix (Doing this inside a separate for
%loop because I couldn't figure out how to do this inside the parfor loop):
for X=1:length(x)
        d{1}=sim_dAdd(d{1},out_posx_rpc{X});
end

% %Initialize structures:
out_posy_rpc =cell(length(x),1);
out=struct([]);

%Now loop through y direction (second refoc pulse only);
parfor Y=1:length(y) %Use this if you do have the MATLAB parallel processing toolbox
%for Y=1:length(y) %Use this if you don't have the MATLAB parallel processing toolbox
 % Calculate the progress percentage
    progressPercentage = Y / numIterations * 100;
  % Create a progress bar string
    progressBar = ['[', repmat('=', 1, round(progressPercentage)), repmat(' ', 1, 100 - round(progressPercentage)), ']'];
    
    % Display the progress bar in the command window
    fprintf('Progress Y-position: %s %.1f%%\r', progressBar, progressPercentage);
     %   disp(['Executing Y-position ' num2str(Y) '!!' ]);
     % Perform your computation here
        out_posy_rpc{Y}=sim_sLASER_shaped_Ref2(d{1},Npts,sw,Bfield,lw,sys,te,rfPulse,refTp,y(Y),Gy,flipAngle,centreFreq);
%                            sim_sLASER_shaped_Ref2(d,   n,sw,Bfield,linewidth,sys,te,RF,       tp, dy,  Gy,ph3,    ph4,  centreFreq)
end

% Print a new line to clear the progress bar
fprintf('\n');


%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
for Y=1:length(y)
        out=op_addScans(out,out_posy_rpc{Y});
end


%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims=(nX*nY);
out=op_ampScale(out,1/numSims);

%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
voxRatio=(thkX*thkY)/(fovX*fovY);
out=op_ampScale(out,1/voxRatio);


end















       