%% Simulate in X-direction only
function d = sim_sLASER_shaped_Ref1(Bfield,sys,te,RF,tp,dx,Gx,flipAngle,centreFreq)

if nargin<9 %11 Falsch bei Jamie
    centreFreq=2.3;
    if nargin<8 %10 Falsch bei Jamie
        flipAngle=180;
    end
end

%Check if this is a gradient modulated pulse.  If so, set Gx equal to zero:
if RF.isGM
    %Scale the GM waveform by the factor Gx and then set Gx equal to zero:
    RF=rf_scaleGrad(RF,Gx);
    Gx=0;
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
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN sLASER PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                  %EXCITE instantaneously
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dx,Gx);          %1st shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_COF(H,d,1);
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dx,Gx);          %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2

end