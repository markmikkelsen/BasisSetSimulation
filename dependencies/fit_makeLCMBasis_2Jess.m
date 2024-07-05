% fit_makeBasis.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [BASIS] = fit_makeBasis(folder, addMMFlag, sequence, editTarget)
% 
% DESCRIPTION:
% Generates a basis set in FID-A structure. The code will search all *.mat
% files in the input folder for FID-A structures with simulated spectra. It
% also performs sanity checks on the simulation parameters, and returns
% warnings if parameters are not identical for all parameters.
% 
% INPUTS:
% folder    = folder containing *.mat files representing FID-A structures
% addMMFlag = Flag to decide whether MM and lipid basis functions should be
%               added to the basis set.
%             OPTIONS:  1 = Add MM+lip (Default)
%                       0 = Don't add MM+lip
% sequence  = sequence type
%             OPTIONS:  'unedited' (default)
%                       'MEGA'
%                       'HERMES'
%                       'HERCULES'
% editTarget= Target molecule of edited data.
%             OPTIONS:  'GABA'
%                       'GSH'
%                       '

%
% OUTPUTS:
% BASIS     = Simulated basis set in FID-A structure format. 

% Kann nicht soviel wie fit_makeBasis von Georg
% Wenn man das einbauen will einfach kopieren, aber Achtung, wegen
% Unterschiede in ifft und fft


function BASIS = fit_makeLCMBasis_2Jess(folder, addMMFlag,fullpath_to_save_basis,vendor,sequence)

% folder that contains all matfiles form simulation
% addMMflag
% sequence

%
% DC offset correction?
%
do_offset_correction=true;
if do_offset_correction
    disp('DC offset correction : PPMOFF')
end
%
% Collect *.mat filenames from input folder
mat_files       = dir([folder filesep '*.mat']);
mat_filenames   = strcat(folder, filesep, {mat_files.name});
idx = contains(mat_filenames, 'Ref');
mat_filenames(idx) = [];
nMets           = length(mat_filenames);
%
fprintf('Number of Metabolites : %d\n',nMets)
%
% Loop over all *.mat filenames, load their data, store in a buffer
%
for kk = 1:nMets

    temp = load(mat_filenames{kk});
    %
    % Multiplexed experiments (e.g. MEGA/HERMES) have more than one sub-basis 
    % simulated per metabolite. Find out how many:
    basisFct = fieldnames(temp);
    % Load the signals, DC-correct and store them in separate dimensions
    for ll = 1:length(basisFct)
        if isfield(temp.(basisFct{1}), 'centerFreq')
            buffer.centerFreq = temp.(basisFct{1}).centerFreq;
        else
            temp.(basisFct{1}).centerFreq = 3;
            buffer.centerFreq = 3;
        end
        % orig
        %temp.(basisFct{ll}).specs      = fftshift(fft(temp.(basisFct{ll}).fids, [], 1), 1);
        % aber in fida eigentlich so
        %fftshift(ifft(fids,[],dims.t),dims.t);
        temp.(basisFct{ll}).specs      = fftshift(ifft(temp.(basisFct{ll}).fids, [], 1), 1);
        
        spectralwidth = temp.(basisFct{ll}).spectralwidth;
        sz = temp.(basisFct{ll}).sz;
        Bo = temp.(basisFct{ll}).Bo;
        f=(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)));
        ppm=f/(Bo*42.577);
        ppm=-(ppm-4.65); % achtung 4.68 before

        %temp.(basisFct{ll}).ppm = ppm - (4.68 - temp.(basisFct{1}).centerFreq);
        temp.(basisFct{ll}).ppm=ppm;
        % Niklaus
        % now we can use it with ifft above
        % otherwise there was an error
        % but this is the thing that causes the difference between
        % makebasis and this script in the first point of the fid (which
        % turned out to be not that relevant)
        % if in makebasis PPMOFF is set, then something similiar is done 
        if do_offset_correction
            temp.(basisFct{ll})        = op_dccorr(temp.(basisFct{ll}),'p'); 
        end
        buffer.fids(:,kk,ll)           = temp.(basisFct{ll}).fids;
        buffer.specs(:,kk,ll)          = temp.(basisFct{ll}).specs;
    end
    % The following field should always be the same for all sub-bases
    % (unless something is seriously flawed with the simulation code)
    buffer.t(:,kk)              = temp.(basisFct{1}).t;
    buffer.ppm(:,kk)            = temp.(basisFct{1}).ppm;
    buffer.spectralwidth(kk)    = temp.(basisFct{1}).spectralwidth;
    buffer.dwelltime(kk)        = temp.(basisFct{1}).dwelltime;
    buffer.n(kk)                = temp.(basisFct{1}).sz(1);
    buffer.linewidth(kk)        = temp.(basisFct{1}).linewidth;
    buffer.Bo(kk)               = temp.(basisFct{1}).Bo;
    if iscell(temp.(basisFct{1}).seq)
        buffer.seq{kk}              = temp.(basisFct{1}).seq{1};
    else
        buffer.seq{kk}              = temp.(basisFct{1}).seq;
    end
    if isfield(temp.(basisFct{1}),'name')
        buffer.name{kk}             = temp.(basisFct{1}).name;
    else
        C = strsplit(mat_files(kk).name,'_');
        C = C{end};
        buffer.name{kk} = strrep(C,'.mat','');
    end
    buffer.te(kk)               = temp.(basisFct{1}).te;
    buffer.dims                 = temp.(basisFct{1}).dims;
    buffer.flags                = temp.(basisFct{1}).flags;
end

% Test whether parameters are the same across all basis functions; flag
% warning if they are not; write into basis set struct if they are.
seq_params = {'spectralwidth','dwelltime','n','linewidth','Bo','seq','te', 'centerFreq'};
for pp = 1:length(seq_params)
    unique_params = unique(buffer.(seq_params{pp}));
    if length(unique_params) > 1
        error('WARNING! One or more sequence parameters are not the same across all input basis functions.');
    else
        BASIS.(seq_params{pp}) = unique_params;
    end
end

% Test whether ppm and t aves are the same across all basis functions; flag
% error if they are not; write into basis set struct if they are.
seq_params = {'ppm','t'};
for pp = 1:length(seq_params)
    unique_params = unique(buffer.(seq_params{pp}),'stable');
    if length(unique_params) ~= BASIS.n
        error('WARNING! One or more sequence parameters are not the same across all input basis functions.');
    else
        BASIS.(seq_params{pp}) = unique_params';
    end
end

% If chosen, add MM
% this is not tested 
% could be useful
if addMMFlag
    n = BASIS.n;
    sw = BASIS.spectralwidth;
    Bo = BASIS.Bo;
    centerFreq = BASIS.centerFreq;
    % The amplitude and FWHM values are determined as for the LCModel and
    % TARQUIN algorithms (see Wilson et al., MRM 2011).
    hzppm = Bo*42.577;
    
    % To scale the amplitudes correctly, we first need to determine the
    % area of the 3.027 ppm CH3 signal of creatine
    CrArea = detCrArea(buffer);
    oneProtonArea = CrArea/3;
    
    % Next, we determine the area of a Gaussian singlet with nominal area 1
    testGaussian    = op_gaussianPeak(n,sw,Bo,centerFreq,0.1*hzppm,0,1);
    testGaussian    = op_dccorr(testGaussian,'p');
    gaussianArea    = sum(real(testGaussian.specs));
    
    % Now we know the scaling factor to generate MM/lipid signals with the
    % correct relative scaling with respect to the CH3 signal
    MM09            = op_gaussianPeak(n,sw,Bo,centerFreq,0.14*hzppm,0.91,3*oneProtonArea/gaussianArea);
    MMBase.MM09     = op_dccorr(MM09,'p');
    MM12            = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.21,2*oneProtonArea/gaussianArea);
    MMBase.MM12     = op_dccorr(MM12,'p');
    MM14            = op_gaussianPeak(n,sw,Bo,centerFreq,0.17*hzppm,1.43,2*oneProtonArea/gaussianArea);
    MMBase.MM14     = op_dccorr(MM14,'p');
    MM17            = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.67,2*oneProtonArea/gaussianArea);
    MMBase.MM17     = op_dccorr(MM17,'p');
    MM20a           = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,2.08,1.33*oneProtonArea/gaussianArea);
    MM20b           = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm,2.25,0.33*oneProtonArea/gaussianArea);
    MM20c           = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.95,0.33*oneProtonArea/gaussianArea);
    MM20d           = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm,3.0,0.4*oneProtonArea/gaussianArea);
    MM20            = op_addScans(MM20a,MM20b); MM20 = op_addScans(MM20,MM20c); MM20 = op_addScans(MM20,MM20d);
    MMBase.MM20     = op_dccorr(MM20,'p');
    Lip09           = op_gaussianPeak(n,sw,Bo,centerFreq,0.14*hzppm,0.89,3*oneProtonArea/gaussianArea);
    MMBase.Lip09    = op_dccorr(Lip09,'p');
    Lip13a          = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.28,2*oneProtonArea/gaussianArea);
    Lip13b          = op_gaussianPeak(n,sw,Bo,centerFreq,0.89*hzppm,1.28,2*oneProtonArea/gaussianArea);
    Lip13           = op_addScans(Lip13a,Lip13b);
    MMBase.Lip13    = op_dccorr(Lip13,'p');
    Lip20a          = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,2.04,1.33*oneProtonArea/gaussianArea);
    Lip20b          = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,2.25,0.67*oneProtonArea/gaussianArea);
    Lip20c          = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm,2.8,0.87*oneProtonArea/gaussianArea);
    Lip20           = op_addScans(Lip20a,Lip20b); Lip20 = op_addScans(Lip20,Lip20c);
    MMBase.Lip20    = op_dccorr(Lip20,'p');
    MMLips = {'MM09','MM12','MM14','MM17','MM20','Lip09','Lip13','Lip20'};
    
    % Now copy over the names, fids, and specs into the basis set structure
    for rr = 1:length(MMLips)
        buffer.name{nMets+rr} = MMLips{rr};
        for qq = 1:length(basisFct)
            buffer.fids(:,nMets+rr,qq)  = MMBase.(MMLips{rr}).fids;
            buffer.specs(:,nMets+rr,qq) = MMBase.(MMLips{rr}).specs;
        end
    end
    
    BASIS.flags.addedMM     = 1;
    BASIS.nMM               = length(MMLips);
    save_str = '_MM';
else
    BASIS.flags.addedMM     = 0;
    BASIS.nMM               = 0;
    save_str = '_noMM';
end



% Copy over the FID, specs, dims, and the metabolite names
BASIS.fids              = buffer.fids;
BASIS.specs             = buffer.specs;
BASIS.name              = buffer.name;
BASIS.dims              = buffer.dims;
BASIS.flags             = buffer.flags;
BASIS.nMets             = nMets;
BASIS.sz                = size(BASIS.fids);

% Normalize basis set
%BASIS.scale = max(max(max(real(buffer.specs))));
%BASIS.fids  = BASIS.fids ./ BASIS.scale;
%BASIS.specs = BASIS.specs ./ BASIS.scale;

%
% save
%

print_basis(BASIS,strrep(fullpath_to_save_basis,'.basis','.pdf'))
% Vorsicht noch alles genau anschauen
io_writelcmBASIS(BASIS,fullpath_to_save_basis,vendor,sequence);

% Save as *.mat file
% Nizo rausgenommen vorerst
%save(['BASIS' save_str '.mat'], 'BASIS');

end


%% io_writelcmBASIS
%   This function creates a LCM usable .BASIS file from an Osprey Basisset.
%   MMs are not included in the output
%
%
%   USAGE:
%       RF=io_writelcmBASIS(in,outfile,vendor,resample);
%
%   INPUT:      in      = Osprey BASIS file
%               outfile = path and name of the LCModel .BASIS file
%               vendor  = String with Vendor name for consistent naming 
%               SEQ     = name of the sequence
%
%   OUTPUT:     RF is unused, but .BASIS file is created
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-02-11: First version of the code.


function RF=io_writelcmBASIS(in,outfile,vendor,SEQ)

% metabList = fit_createMetabList({'full'});
% 
%  % Add basis spectra (if they were removed to reduce thhe file size)
% if ~isfield(in,'specs')
%     [in]=osp_recalculate_basis_specs(in);
% end
% 
% % Create the modified basis set without macro molecules 
% basisSet = fit_selectMetabs(in, metabList, 0);
basisSet=in;
Bo=basisSet.Bo;
HZPPPM=42.577*Bo;
FWHMBA = basisSet.linewidth/HZPPPM;
ECHOT = basisSet.te;
  
BADELT=basisSet.dwelltime;
NDATAB= basisSet.sz(1)*2; % Achtung wegen zerofilling

XTRASH = 0;
%write to txt file
fid=fopen(outfile,'w+');
fprintf(fid,' $SEQPAR');
fprintf(fid,'\n FWHMBA = %5.6f,',FWHMBA);
fprintf(fid,'\n HZPPPM = %5.6f,',HZPPPM);
fprintf(fid,'\n ECHOT = %2.2f,',ECHOT);
fprintf(fid,'\n SEQ = ''%s''',SEQ);
fprintf(fid,'\n $END');
fprintf(fid,'\n $BASIS1');
%
fprintf(fid,'\n IDBASI = ''%s'',',[vendor ' ' SEQ ' ' num2str(ECHOT) '/ (c) Jess and Niklaus']);
%
fprintf(fid,'\n FMTBAS = ''(2E15.6)'',');
% 6E13.5
%fprintf(fid,'\n FMTBAS = ''(6E13.5)'',');
fprintf(fid,'\n BADELT = %5.6f,',BADELT);
fprintf(fid,'\n NDATAB = %i', NDATAB);
fprintf(fid,'\n $END\n');
for i = 1 : basisSet.nMets
    if ~strcmp(basisSet.name{i}, 'CrCH2') && ~strcmp(basisSet.name{i}, 'H2O')
        RF = shift_centerFreq(basisSet,i);
        fprintf(fid,' $NMUSED');
        % This doesnt change anything when changed 
        % but they need to be set equally for the makebasis to get the same
        % results 
        fprintf(fid,'\n AUTOPH = %s','F');
        fprintf(fid,'\n AUTOSC = %s','F');
        fprintf(fid,'\n NOSHIF = %s','T');
        %
        fprintf(fid,'\n XTRASH = %2.2f',XTRASH); % as now effect as noting is implemented with this here 
        fprintf(fid,'\n $END');
        fprintf(fid,'\n $BASIS');
        fprintf(fid,'\n ID = ''%s'',',basisSet.name{i});
        fprintf(fid,'\n METABO = ''%s'',',basisSet.name{i});
        fprintf(fid,'\n CONC = 1.,');%
        fprintf(fid,'\n TRAMP = 1.,');
        fprintf(fid,'\n VOLUME = 1.,');
        fprintf(fid,'\n ISHIFT = 0'); % 0 with this shift the basis shift can be shifted during the fit
        fprintf(fid,'\n $END\n');
        %2E15.6
        fprintf(fid,' %7.6e  %7.6e\n',RF');
        % if  6E13.5
		%RFprint=RF';
        %for n =1:(size(RFprint,2))
        %    if mod(n,3)==0
        %        fprintf(fid,'%13.5E%13.5E\n',RFprint(1,n),RFprint(2,n));
        %    else
        %        fprintf(fid,'%13.5E%13.5E',RFprint(1,n),RFprint(2,n)); 
        %    end
        %end 
        % last row
        %fprintf(fid,'\n');
   end
end

fclose(fid);
end

% nizo Be carefull
function [RF] = shift_centerFreq(data_struct,idx)

    t=repmat(data_struct.t',[1 data_struct.sz(2:end,1)]);
    hzpppm = data_struct.Bo*42.577;
    %f = (4.68-data_struct.centerFreq)*hzpppm;%nizo removed
    fids = data_struct.fids(:,idx);
    %fids=fids.*exp(-1i*t*f*2*pi);% nizo removed
    %Take the complex conjugate becuase the sense of rotation in LCModel seems to
    %be opposite to that used in FID-A.
    fids = conj(fids);
    %
    %scaling the fids % not relevant
    % zerofilling this is done in the makebasis as well
    samples=data_struct.sz(1);
    fids_zf=zeros(2*samples,1);
    fids_zf(1:samples,1)=fids;
    %
    %
    specs=(fft(fids_zf,[],data_struct.dims.t));
    % vorher
    %specs=(ifft(fids,[],data_struct.dims.t));
    RF=zeros(length(specs(:)),2);
    RF(:,1)=real(specs(:));
    RF(:,2)=imag(specs(:));

end

function print_basis(BASIS,fullpath_to_basis)
%
xlim_range=[-1 10.0];
%
a=figure;
plot(BASIS.ppm,real(BASIS.specs));legend(BASIS.name,'Location','eastoutside')
set(gca,'xdir','reverse','XGrid','on')
%
text(0.1,0.8,['Echo Time: ',sprintf('%d',BASIS.te)],'Units','normalized')
text(0.1,0.75,[BASIS.seq{1}],'Units','normalized')
text(0.1,0.70,['No. Mets: ',sprintf('%d',BASIS.nMets)],'Units','normalized')
text(0.1,0.65,['LW: ',sprintf('%0.2f',BASIS.linewidth)],'Units','normalized')
text(0.1,0.60,['SpectralW: ',sprintf('%d',BASIS.spectralwidth)],'Units','normalized')
%
ax=gca;
ax.XAxis.MinorTick       = 'on';
ax.XAxis.MinorTickValues = xlim_range(1):0.5:xlim_range(2);
ax.XMinorGrid = 'on';
xlim(xlim_range)
%
exportgraphics(a,fullpath_to_basis, 'Append', false);

for jj=1:BASIS.nMets
    a=figure;
    subplot(2,1,1);
    plot(BASIS.ppm,real(BASIS.specs(:,jj)));
    set(gca,'xdir','reverse','XGrid','on')
    ax=gca;
    ax.XAxis.MinorTick       = 'on';
    ax.XAxis.MinorTickValues = xlim_range(1):0.5:xlim_range(2);
    ax.XMinorGrid = 'on';
    xlim(xlim_range)
    %
    subplot(2,1,2);
    plot(BASIS.ppm,imag(BASIS.specs(:,jj)));
    set(gca,'xdir','reverse','XGrid','on')
    ax=gca;
    ax.XAxis.MinorTick       = 'on';
    ax.XAxis.MinorTickValues = xlim_range(1):0.5:xlim_range(2);
    ax.XMinorGrid = 'on';
    xlim(xlim_range)
    xlim(xlim_range)
    %
    sgtitle(BASIS.name{jj});
    %
    exportgraphics(a,fullpath_to_basis, 'Append', true);
end

end


% detCrArea.m
% Georg Oeltzschner, Johns Hopkins University 2020
% 
% USAGE:
% [CrArea] = detCrArea(buffer);
% 
% DESCRIPTION:
% Finds the creatine spectrum in the temporary basis set buffer, then fits
% a Lorentzian to the 3.027 ppm CH3 creatine singlet to determine its area.
% Subsequently, macromolecule and lipid basis functions are scaled
% accordingly.
% 
% INPUTS:
% in        = a temporary buffer containing simulated basis functions
%
% OUTPUTS:
% CrArea    = Estimated area under the 3.027 ppm CH3 Cr singlet.

function [CrArea] = detCrArea(in)

% Find the creatine basis function
idx_Cr          = find(strcmp(in.name,'Cr'));
if isempty(idx_Cr)
    error('No basis function with nametag ''Cr'' found! Abort!');
end

%[~, idx_3027]   = min(abs(buffer.ppm(:,1)-3.027));

% Determine the window where we are going to look for the peak.
ppm = in.ppm(:,1);
ppmmin = 3.027 - 0.4;
ppmmax = 3.027 + 0.4;
refWindow = in.specs(ppm>ppmmin & ppm<ppmmax, idx_Cr);
ppmWindow = in.ppm(ppm>ppmmin & ppm<ppmmax);

% Find the maximum and its index
maxRef_index    = find(abs(real(refWindow)) == max(abs(real((refWindow)))));
maxRef          = real(refWindow(maxRef_index));

% Determine an initial estimate for the FWHM
% Peak lines can be super narrow, so overestimate it slightly
gtHalfMax   = find(abs(real(refWindow)) >= 0.4*abs(maxRef));
FWHM1       = abs(ppmWindow(gtHalfMax(1)) - ppmWindow(gtHalfMax(end)));
FWHM1       = FWHM1*(42.577*in.Bo(1));  %Assumes proton.

% Determine an initial estimate for the center frequency of the Cr peak
crFreq = ppmWindow(maxRef_index);

% Set up the fit
parsGuess=zeros(1,5);
parsGuess(1) = maxRef;  % amplitude
parsGuess(2) = (5*in.Bo/3)/(42.577*in.Bo); %FWHM.  Assumes Proton.  LW = 5/3 Hz/T.   % FWHM. Assumes Proton.
parsGuess(3) = crFreq;  % center frequency
parsGuess(4) = 0;       % baseline offset
parsGuess(5) = 0;       % phase
    
% Run first guess
yGuess  = op_lorentz(parsGuess, ppmWindow);
parsFit = nlinfit(ppmWindow, real(refWindow'), @op_lorentz, parsGuess);
yFit    = op_lorentz(parsFit, ppmWindow);
    
% figure;
% plot(ppmWindow,refWindow,'.',ppmWindow,yGuess,':',ppmWindow,yFit);
% legend('data','guess','fit');

CrArea = sum(yFit);

end
