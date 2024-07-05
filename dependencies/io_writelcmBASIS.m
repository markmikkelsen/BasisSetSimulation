%% io_writelcmBASIS
%   This function creates a LCM usable .BASIS file from an Osprey in.
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

%%%%% August 2022 Modified by Niklaus Zölch to run in
%%%%% 'run-simLASERShapred_fast_JNZ' 

function io_writelcmBASIS(in,outfile,vendor,SEQ)

% This we need to discuss how we want to organize that
% But so far we set the list of metabolites that should be simulated and then all of them % are in the basis.
% metabList = fit_createMetabList({'full'});

%  Add basis spectra (if they were removed to reduce thhe file size)
% if ~isfield(in,'specs')
%     [in]=osp_recalculate_basis_specs(in);
% end

% Create the modified basis set without macro molecules 
% in = fit_selectMetabs(in, metabList, 0);

Bo=in.Bo;
HZPPPM=42.577*Bo;
FWHMBA = in.linewidth/HZPPPM;
ECHOT = in.te;
  
BADELT=in.dwelltime;
NDATAB= in.sz(1);

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
fprintf(fid,'\n IDBASI = ''%s'',',[vendor ' ' SEQ ' ' num2str(ECHOT) ' ms Osprey']);
fprintf(fid,'\n FMTBAS = ''(2E15.6)'',');
fprintf(fid,'\n BADELT = %5.6f,',BADELT);
fprintf(fid,'\n NDATAB = %i', NDATAB);
fprintf(fid,'\n $END\n');
for i = 1 : in.nMets
    if ~strcmp(in.name{i}, 'CrCH2') && ~strcmp(in.name{i}, 'H2O')
        RF = shift_centerFreq(in,i);
        fprintf(fid,' $NMUSED');
        fprintf(fid,'\n XTRASH = %2.2f',XTRASH);
        fprintf(fid,'\n $END');
        fprintf(fid,'\n $BASIS');
        fprintf(fid,'\n ID = ''%s'',',in.name{i});
        fprintf(fid,'\n METABO = ''%s'',',in.name{i});
        fprintf(fid,'\n CONC = 1.,');
        fprintf(fid,'\n TRAMP = 1.,');
        fprintf(fid,'\n VOLUME = 1.,');
        fprintf(fid,'\n ISHIFT = 0');
        %J
        fprintf(fid,'\n AUTOPH = %s','F');
        fprintf(fid,'\n AUTOSC = %s','F');
        fprintf(fid,'\n NOSHIF = %s','T');
        fprintf(fid,'\n $END\n');
        fprintf(fid,' %7.6e  %7.6e\n',RF');
    end
end

fclose(fid);
end

function [RF] = shift_centerFreq(data_struct,idx)

    t=repmat(data_struct.t',[1 data_struct.sz(2:end,1)]);
    hzpppm = data_struct.Bo*42.577;
    f = (4.68-data_struct.centerFreq)*hzpppm;
    fids = data_struct.fids(:,idx);
    % fids=fids.*exp(-1i*t*f*2*pi);- %% Remove the shift since we do it in the simulation % 
    %Take the complex conjugate becuase the sense of rotation in LCModel seems to
    %be opposite to that used in FID-A.
    fids = conj(fids);
    % vorher
    %specs=(fft(fids,[],data_struct.dims.t));
    specs=(ifft(fids,[],data_struct.dims.t));
    RF=zeros(length(specs(:)),2);
    RF(:,1)=real(specs(:));
    RF(:,2)=imag(specs(:));

end

function [in]=osp_recalculate_basis_specs(in)
    % This function recalculates the basis spectra and ppm-axis of the
    % basis set

    in.specs = fftshift(fft(in.fids,[],1),1);

    % Calcualte ppm-axis
    f = [(-in.spectralwidth/2)+(in.spectralwidth/(2*in.sz(1))):in.spectralwidth/(in.sz(1)):(in.spectralwidth/2)-(in.spectralwidth/(2*in.sz(1)))];
    in.ppm = f/(in.Bo*42.577);
    in.ppm=in.ppm + in.centerFreq;
end