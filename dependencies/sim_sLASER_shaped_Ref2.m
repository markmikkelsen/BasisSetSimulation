%% Simulate in Y-direction only
function out = sim_sLASER_shaped_Ref2(d,n,sw,Bfield,linewidth,sys,te,RF,tp,dy,Gy,flipAngle,centreFreq)

if nargin<13%15 falsch bei Jamie
    centreFreq=2.3;
    if nargin<12%14 falsch bei Jamie
        flipAngle=180;
    end
end

%Check if this is a gradient modulated pulse.  If so, set Gy equal to zero:
if RF.isGM
    %Scale the GM waveform by the factor Gy and then set Gy equal to zero:
    RF=rf_scaleGrad(RF,Gy);
    Gy=0;
end

if (te/4)<(tp/1000)
    error('ERROR: the duration of the refocusing pulse cannot be longer than a quarter of the echo time! ABORTING!!');
end

%initialize evolution times
tau1=(te/4-tp)/2;
tau2=te/4-tp;

%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H]=sim_Hamiltonian(sys,Bfield);

%BEGIN sLASER PULSE SEQUENCE************
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dy,Gy);          %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_COF(H,d,1);
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dy,Gy);          %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1

[out,~]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along +y axis (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='semi-LASER';
out.te=te;
out.sim='shaped';

%Additional fields for compatibility with FID-A processing tools.
out.sz=size(out.specs);
out.date=date;
out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;
out.dims.extras=0;
out.averages=1;
out.rawAverages=1;
out.subspecs=1;
out.rawSubspecs=1;
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=1;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.isISIS=0;
end


