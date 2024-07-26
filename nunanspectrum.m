


function [varargout] = nunanspectrum(TS, T, TIME_UNITS, varargin)
%% Eliminate missing data (nufft does not need them)

if isreal(TS) % real (scalar) time series
    REAL = true;
else % complex (vector, e.g. current) time series u + i*v
    REAL = false;
end

TS = TS(:);
T = T(:);
T = T(isfinite(TS));
TS = TS(isfinite(TS));

%% Sort out optional inputs

if ~isempty(varargin)
    % Turn option/value pairs into a structure:
    Struct = struct(varargin{:});
    Names = fieldnames(Struct);
    % Fields should be drawn from:
    AllowedVars = {'Segments', 'Plot_option', 'Plot', 'Window', 'Freq'};
    % Note: Freq should be defined by df:df:f_Ny, not the 0:df:[2*f_Ny - df] that nufft assumes
    % (those adjustments are made automatically).
    for ii=1:length(Names)
        if ismember(Names{ii},AllowedVars)
        else
            error(['''',Names{ii},''' is not a possible option for this function. Please see documentation.'])
        end
    end
    % If that test is passed, reassign to variables without the structure
    % prefix (this should not affect memory in MATLAB because assigning two
    % variables to the same array does not copy the array until one of them is
    % changed). Nevertheless, clear "Str" after.
    for ii=1:length(Names)
        eval([Names{ii},' = Struct.',Names{ii},';'])
    end
    % Note: dynamically assigning variables like above is generally bad
    % practice; however, in this case the variables are constrained to the few
    % named in "AllowedVars", so there's no chance of this spiraling out of
    % control.
end
clear Struct
% Set any variables that were undefined:
if ~exist('Segments','var')
    Segments = 1;
end
if ~exist('Plot_option','var')
    Plot_option = '.-';
end
if ~exist('Plot','var')
    Plot = false;
end
if ~exist('Window','var')
    Window = 'rectwin';
end
if ~exist('Freq','var')
    df = Segments/[T(end) - T(1)];
    % dt_avg = (T(end) - T(1))/(length(TS) - 1); % bad for large gaps
    dt_avg = median(diff(T),'omitnan');
    f_Ny = 1/[2*dt_avg];
    Freq = [df:df:f_Ny]';
    Freq_ = [0; Freq(1:[end-1])];
    FreqFreq = [df:df:(2*f_Ny)]';
    FreqFreq_ = [0; FreqFreq(1:[end-1])];
else
    Freq = Freq(:);
    FreqFreq = [Freq ; Freq+Freq(end)];
    Freq_ = [0; Freq(1:[end-1])];
    FreqFreq_ = [0; FreqFreq(1:[end-1])];
end

%% Segment

TS_reshape = cell(2*Segments - 1,1);
T_reshape = cell(2*Segments - 1,1);
N_TS = length(TS);
seg_len = round(N_TS/Segments);
if Segments == 1
    TS_reshape{1} = TS;
    T_reshape{1} = T;
else
    if mod(seg_len,2) % odd
        seg_i = -1;
        for ii = 1:2:[length(TS_reshape) - 2]
            seg_i = seg_i + 1;
            % edge segments
            TS_reshape{ii} = TS( (1 + seg_i*seg_len):((seg_i + 1)*seg_len) );
            T_reshape{ii} =   T( (1 + seg_i*seg_len):((seg_i + 1)*seg_len) );
            % overlap segments
            TS_reshape{ii+1} = TS( (1 + seg_i*seg_len + ([seg_len - 1]/2)):((seg_i + 1)*seg_len + ([seg_len - 1]/2)) );
            T_reshape{ii+1} =   T( (1 + seg_i*seg_len + ([seg_len - 1]/2)):((seg_i + 1)*seg_len + ([seg_len - 1]/2)) );
        end
        TS_reshape{ii+2} = TS(((seg_i + 1)*seg_len + 1):end);
        T_reshape{ii+2} =   T(((seg_i + 1)*seg_len + 1):end);
    else % even
        seg_i = -1;
        for ii = 1:2:[length(TS_reshape) - 2]
            seg_i = seg_i + 1;
            % edge segments
            TS_reshape{ii} = TS( (1 + seg_i*seg_len):((seg_i + 1)*seg_len) );
            T_reshape{ii} =   T( (1 + seg_i*seg_len):((seg_i + 1)*seg_len) );
            % overlap segments
            TS_reshape{ii+1} = TS( (1 + seg_i*seg_len + ([seg_len]/2)):((seg_i + 1)*seg_len + ([seg_len]/2)) );
            T_reshape{ii+1} =   T( (1 + seg_i*seg_len + ([seg_len]/2)):((seg_i + 1)*seg_len + ([seg_len]/2)) );
        end
        TS_reshape{ii+2} = TS(((seg_i + 1)*seg_len + 1):end);
        T_reshape{ii+2} =   T(((seg_i + 1)*seg_len + 1):end);
    end
end

%% Detrend data segments

% Remove the linear trend in each segment via a least squares fit
for ii = 1:length(TS_reshape)
    HH = [ones(size(T_reshape{ii})) [T_reshape{ii} - T_reshape{ii}(1)]];
    TREND_LINE_COEF = HH\TS_reshape{ii};
    TS_reshape{ii} = TS_reshape{ii} - HH*TREND_LINE_COEF;
end

%% Window
% For future window options, define "NormFactor" as the following:
% NormFactor = 1/mean([0; Window(T_reshape{1})]);
% where T_reshape{1} could be any other element of T_reshape

if REAL
    mean_var = 0;
    for ii = 1:length(TS_reshape)
        TSii_var = var(TS_reshape{ii});
        mean_var = mean_var + TSii_var/length(TS_reshape);
    end
else
    mean_var = [0 0];
    for ii = 1:length(TS_reshape)
        % % % TSii_var = var(real(TS_reshape{ii}));
        % % % mean_var(1) = mean_var(1) + TSii_var/length(TS_reshape);
        % % % TSii_var = var(imag(TS_reshape{ii}));
        % % % mean_var(2) = mean_var(2) + TSii_var/length(TS_reshape);

        % We can use the regular FFT here because we are only looking for
        % the CCW vs. CW variance, which is encoded in the FFT even if the
        % frequencies are wrong:
        FFT_TSii = fft(TS_reshape{ii});
        CCW_FFT_TSii = FFT_TSii(2:[floor(length(FFT_TSii)/2 + 1)]);
        CW_FFT_TSii  = FFT_TSii([floor(length(FFT_TSii)/2 + 1) + 1]:end);
        CCW_var = sum(CCW_FFT_TSii.*conj(CCW_FFT_TSii))/[length(TS_reshape{ii})^2];
        CW_var  = sum(CW_FFT_TSii.*conj(CW_FFT_TSii))/[length(TS_reshape{ii})^2];
        mean_var(1) = mean_var(1) + CCW_var/length(TS_reshape);
        mean_var(2) = mean_var(2) + CW_var/length(TS_reshape);
    end
end

if strcmp(Window,'rectwin')
    % no change to the data
    NormFactor = 1;
elseif strcmp(Window,'hanning')
    for ii = 1:length(TS_reshape)
        TS_reshape{ii} = TS_reshape{ii}.*[sin([T_reshape{ii} - T_reshape{ii}(1)]*pi/[T_reshape{ii}(end) - T_reshape{ii}(1)]).^2];
    end
    NormFactor = 1/mean([0; sin([T_reshape{1} - T_reshape{1}(1)]*pi/[T_reshape{1}(end) - T_reshape{1}(1)]).^2].^2);
else
    error('At this time, only von Hann (''hanning'') and rectangular (''rectwin'') windows are implemented in this function.')
end
% for ii = 1:length(TS_reshape); plot(T_reshape{ii},TS_reshape{ii},'.-');hold on; end % for testing purposes

%% Calculate power spectrum with nufft

NUFFT_TS = cell(2*Segments - 1,1);
mean_NUFFT_TS_squared = zeros(size(FreqFreq_));
for ii = 1:length(TS_reshape)
    NUFFT_TS{ii} = nufft(TS_reshape{ii}, T_reshape{ii}, FreqFreq_)/[length(T_reshape{ii})];
    mean_NUFFT_TS_squared = mean_NUFFT_TS_squared + [abs(NUFFT_TS{ii}).^2]/length(TS_reshape);
end
% figure;for ii = 1:length(NUFFT_TS); loglog(FreqFreq_,abs(NUFFT_TS{ii}).^2,'.-');hold on; end % for testing purposes
% semilogy(FreqFreq_,mean_NUFFT_TS_squared,'k.-','MarkerSize',15)% for testing purposes
% figure; semilogy(mean_NUFFT_TS_squared,'k.-','MarkerSize',15)% for testing purposes

if REAL
    if mod(length(mean_NUFFT_TS_squared),2) % odd
        N_spec = length(mean_NUFFT_TS_squared);
        SPEC =   [mean_NUFFT_TS_squared(2:[[N_spec + 1]/2]) + flip(mean_NUFFT_TS_squared([[N_spec + 1]/2 + 1]:end)) ]/2;
    else % even
        N_spec = length(mean_NUFFT_TS_squared);
        SPEC = [ [mean_NUFFT_TS_squared(2:[N_spec/2]) + flip(mean_NUFFT_TS_squared([N_spec/2 + 2]:end)) ]/2 ; 2*mean_NUFFT_TS_squared([N_spec/2 + 1])];
    end
    SPEC = 2*NormFactor*(N_spec*FreqFreq(end))*SPEC;
    % The above is the correct normalization step, but because the effect of
    % gaps in the data are difficult to account for, simply force the spectrum
    % to fulfill Parseval's Theorem:
    SPEC = mean_var*SPEC/sum(SPEC*Freq(1));
else
    if mod(length(mean_NUFFT_TS_squared),2) % odd
        N_spec = length(mean_NUFFT_TS_squared);
        SPEC(:,1) =        mean_NUFFT_TS_squared(2:[[N_spec + 1]/2]);
        SPEC(:,2) = + flip(mean_NUFFT_TS_squared([[N_spec + 1]/2 + 1]:end));
    else % even
        N_spec = length(mean_NUFFT_TS_squared);
        SPEC(:,1) = [      mean_NUFFT_TS_squared(2:[N_spec/2]) ;        mean_NUFFT_TS_squared([N_spec/2 + 1])];
        SPEC(:,2) = [ flip(mean_NUFFT_TS_squared([N_spec/2 + 2]:end)) ; mean_NUFFT_TS_squared([N_spec/2 + 1])];
    end
    SPEC = NormFactor*(N_spec*FreqFreq(end))*SPEC;
    % The above is the correct normalization step, but because the effect of
    % gaps in the data are difficult to account for, simply force the spectrum
    % to fulfill Parseval's Theorem:
    SPEC(:,1) = mean_var(1)*SPEC(:,1)/sum(SPEC(:,1)*Freq(1));
    SPEC(:,2) = mean_var(2)*SPEC(:,2)/sum(SPEC(:,2)*Freq(1));
end

%% Error

% In short, do this calculation for 95% confidence ratio:
EffectiveSegments = 2*Segments - 1;% degrees of freedom; for overlapping
                     % Hann windowed segments, this is the total number of
                     % segments, i.e.:
                     % 2*SEGMENTS - 1
                     % This may only work perfectly with the Hann window.
                     % With other windows, it's not entirely clear, but my
                     % guess is, if the window is "Hann-like", it will not
                     % make a big difference.

err_high = 2*EffectiveSegments/chi2inv(.05/2,2*EffectiveSegments);
err_low = 2*EffectiveSegments/chi2inv(1-.05/2,2*EffectiveSegments);
Err = [err_low err_high];

%% Plot if requested

if Plot
    if REAL
        % Do not set as a new figure in case the user wants to combine several
        % applications of nunanspectrum.
        % sum(SPEC*Freq(1)) % for testing purposes
        loglog(Freq,SPEC,Plot_option)
        xlabel(['Frequency (cycles per ' TIME_UNITS ')'])
        ylabel(['Spectral Power Density [(time\_series\_units)^2 (cycles per ',...
            TIME_UNITS,')^{-1}]']); hold on
        % error bar plotted
        loglog([1.1 1.1]*Freq(end),[1 (err_high/err_low)]*min(SPEC),'k','LineWidth',3);
        legend('Spectrum','95% error')
    else
        % Do not set as a new figure in case the user wants to combine several
        % applications of nunanspectrum.
        % sum(SPEC(:,1)*Freq(1)) % for testing purposes
        % sum(SPEC(:,2)*Freq(1)) % for testing purposes
        loglog(Freq,SPEC(:,1),Plot_option);hold on
        loglog(Freq,SPEC(:,2),Plot_option)
        xlabel(['Frequency (cycles per ' TIME_UNITS ')'])
        ylabel(['Spectral Power Density [(time\_series\_units)^2 (cycles per ',...
            TIME_UNITS,')^{-1}]']); hold on
        % error bar plotted
        loglog([1.1 1.1]*Freq(end),[1 (err_high/err_low)]*min(SPEC(:)),'k','LineWidth',3);
        legend('CCW Spectrum','CW Spectrum','95% error')
    end
else
end

%% Outputs

if nargout == 0
elseif nargout == 1
    varargout{1} = SPEC;
elseif nargout == 2
    varargout{1} = SPEC;
    varargout{2} = Freq;
elseif nargout == 3
    varargout{1} = SPEC;
    varargout{2} = Freq;
    varargout{3} = Err;
else
    varargout{1} = SPEC;
    varargout{2} = Freq;
    varargout{3} = Err;
    warning('Only 1, 2, or 3 outputs are expected')
end

end

%% Examples

% nunanspectrum(AR_make(0.95,10000), [1:10000]', 'X', ...
%     'Segments',11,'Window','hanning','Plot',true)
% nunanspectrum(sin([1:10001]'*2*pi/10) + AR_make(0.95,10001), [1:10001]' - 1, 'X', ...
%     'Segments',10,'Window','hanning','Plot',true)
% nunanspectrum(sin([1:5000 6000:10000]'*2*pi/10) + AR_make(0.95,length([1:5000 6000:10000]')), [1:5000 6000:10000]', 'X', ...
%     'Segments',11,'Window','hanning','Plot',true)
