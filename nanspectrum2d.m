% nanspectrum2d(DATA, DT, DX, TSEGMENTS, XSEGMENTS, TIME_UNITS, SPACE_UNITS, PLOT_BOOLEAN, INTERPMETHOD)
%
% Makes a 2D spectrum given a data set (e.g. one which could be made into a
% Hovmöller diagram).
% 
% Note: This documentation talks about the two dimensions of DATA as if one
% needs to be time and the other needs to be space, but they could both be
% space.
% 
% % IN:  DATA = Data set wherein time progresses along the x-axis and space
% %             progresses along the y-axis, i.e. a single row is a time
% %             series and a single column is a spatial record at fixed time.
% % IN:  DT = Time between measurements
% % IN:  DX = Distance between measurements
% % IN:  TSEGMENTS = Number of segments into which to divide the temporal
% %                  aspect of the data (must be >1)
% % IN:  XSEGMENTS = Number of segments into which to divide the spatial
% %                  aspect of the data (can be = 1)
% % IN:  TIME_UNITS = The units of "DT". Should be in singular form, e.g.
% %                   'second'
% % IN:  SPACE_UNITS = The units of "DX". Should be in singular form, e.g.
% %                   'meter'
% % IN:  PLOT_BOOLEAN = Choose "1" to plot, "0" to not plot
% % IN:  INTERPMETHOD = A one- or two-element cell array, of which the
% %                     first element is the choice that the called
% %                     "intertp1" function will use for the interpolation
% %                     (only matters for data with gaps). See "help
% %                     interp1" for more information. OR: Set to 0 for
% %                     zero-padding (i.e. NaN -> 0), set to 1 to replaces
% %                     NaN's with randn noise with the same  variance as
% %                     the data series (an issue with this is that it will
% %                     give a different answer each time it's used).
% %                     Otherwise use the appropriate string from "interp1"
% %                     options.
% %                    
% %                     The second element in "INTERPMETHOD" is one of two
% %                     strings:
% %                         'time' if you want to interpolate each time
% %                                series separately, or:
% %                         'space' if you want to interpolate each spatial
% %                                record separately
% %                         'both' if you want a 2D interpolation (WARNING:
% %                                this is not yet implemented, but is
% %                                included here as a placeholder)
% % 
% % OUT: Spectrum2D = A 2D spectrum the same size as one of the segments,
% %                   except if one or both of size(segment) is even, in
% %                   which case a redundant(?) row and/or column is
% %                   removed so as to match f_vec and k_vec.
% % OUT: f_vec = the frequencies corresponding to the columns of Spectrum2D
% % OUT: k_vec = the wavenumbers corresponding to the rows of Spectrum2D
% % OUT: DATA_treated = what the function used to get the spectrum, (i.e.
%                       after interpolating) which in an ideal world would
%                       be the same as "DATA".
% 
% [Spectrum2D, f_vec, k_vec, DATA_treated] = nanspectrum2d(DATA, DT, DX, TWINDOWS, XWINDOWS, TIME_UNITS, SPACE_UNITS, PLOT_BOOLEAN, INTERPMETHOD)
% 
% 

function [Spectrum2D, f_vec, k_vec, DATA_treated] = nanspectrum2d(...
    DATA, DT, DX, TSEGMENTS, XSEGMENTS, TIME_UNITS, SPACE_UNITS, ...
    PLOT_BOOLEAN, INTERPMETHOD)
%% Nested function (for plane fit)
function OUT = funct3(tt,xx)
	OUT = ones(length(tt),3);
	OUT(:,2) = tt;
	OUT(:,3) = xx;
end
%% Fast instructions printed to terminal if zero inputs are given
if nargin == 0
    disp('nanspectrum2d(DATA, DT, DX, TSEGMENTS, XSEGMENTS, TIME_UNITS, SPACE_UNITS, PLOT_BOOLEAN, INTERPMETHOD)')
else
end
%% Interpolate

N_t = size(DATA,2); % number of points in time
N_x = size(DATA,1); % number of points in space ("x" is generic space coordinate, not specifically "zonal")
T = 1:N_t;
X = 1:N_x;

if length(INTERPMETHOD) ~= 2
    if isempty(INTERPMETHOD)
        % Nothing, because this option means no interpolation is done
    else
        eval('help nanspectrum2d')
        error('Please format the variable "INTERPMETHOD" correctly.')
    end
elseif strcmp(INTERPMETHOD{2},'time')
    % i.e. in time
    if ischar(INTERPMETHOD{1})
        for i=1:N_x
            TS = DATA(i,:);
            DATA(i,:) = interp1(T(isfinite(TS)),TS(isfinite(TS)),T,INTERPMETHOD{1});
        end
    elseif INTERPMETHOD{1} == 0
        for i=1:N_x
            % Do a quick linear detrend and demean (not possible in "detrend",
            % and no "nandetrend" available). Solve TS' ? A*m -> m = A\TS'
            TS = DATA(i,:);
            TS_nonan = TS(isfinite(TS)); T_nonan = T(isfinite(TS));
%                 if i==1 % define the fitting matrix only once
                    A = ones(length(TS_nonan),2); A(:,2) = T_nonan';
%                 else
%                 end
            m = A\(TS_nonan');
            nancoords = find(~isfinite(TS));
            for j=1:length(nancoords)
                TS(nancoords(j)) = m(1) + m(2)*T(nancoords(j));
            end
            DATA(i,:) = TS;
        end
        
    elseif INTERPMETHOD{1} == 1
        for i=1:N_x
            % Do a quick linear detrend and demean (not possible in "detrend",
            % and no "nandetrend" available). Solve TS' ? A*m -> m = A\TS'
            TS = DATA(i,:);
            TS_nonan = TS(isfinite(TS)); T_nonan = T(isfinite(TS));
%                 if i==1 % define the fitting matrix only once
                    A = ones(length(TS_nonan),2); A(:,2) = T_nonan';
%                 else
%                 end
            m = A\(TS_nonan');
            TS_nonandetrended = TS_nonan - (A*m)'; % not kept around because we eventually
            % detrend the padded segments
            TS_nandetrended = TS; TS_nandetrended(isfinite(TS)) = TS_nonandetrended;
            TS_var = nanvar(TS_nandetrended);
            nancoords = find(~isfinite(TS));
            for j=1:length(nancoords)
                TS(nancoords(j)) = sqrt(TS_var)*randn + m(1) + m(2)*T(nancoords(j));
            end
            DATA(i,:) = TS;
        end
    else
        eval('help NAN_SPECTRUM_MAKE_LK')
        error(['The first element of the variable "INTERPMETHOD" must be 0, 1, or one of several',...
            ' specific strings. See the documentation above.'])
    end
    %%%
    
elseif strcmp(INTERPMETHOD{2},'space')
    % i.e. in space
    if ischar(INTERPMETHOD{1})
        for n=1:N_t
            SS = DATA(:,n)'; % Space series
            DATA(:,n) = interp1(X(isfinite(SS)),SS(isfinite(SS)),X,INTERPMETHOD{1})';
        end
    elseif INTERPMETHOD{1} == 0
        for n=1:N_t
            % Do a quick linear detrend and demean (not possible in "detrend",
            % and no "nandetrend" available). Solve TS' \approx A*m -> m = A\TS'
            SS = DATA(:,n)';
            SS_nonan = SS(isfinite(SS)); X_nonan = X(isfinite(SS));
            A = ones(length(SS_nonan),2); A(:,2) = X_nonan';
            m = A\(SS_nonan');
            nancoords = find(~isfinite(SS));
            for j=1:length(nancoords)
                SS(nancoords(j)) = m(1) + m(2)*X(nancoords(j));
            end
            DATA(:,n) = SS';
        end
        
    elseif INTERPMETHOD{1} == 1
        for n=1:N_t
            % Do a quick linear detrend and demean (not possible in "detrend",
            % and no "nandetrend" available). Solve TS' \approx A*m -> m = A\TS'
            SS = DATA(:,n)';
            SS_nonan = SS(isfinite(SS)); X_nonan = X(isfinite(SS));
            A = ones(length(SS_nonan),2); A(:,2) = X_nonan';
            m = A\(SS_nonan');
            SS_nonandetrended = SS_nonan - (A*m)';
            % ^ Not kept around because we eventually detrend the padded
            % segments.
            SS_nandetrended = SS; SS_nandetrended(isfinite(SS)) = SS_nonandetrended;
            SS_var = mean(nanvar(SS_nandetrended));
            nancoords = find(~isfinite(SS));
            for j=1:length(nancoords)
                SS(nancoords(j)) = sqrt(SS_var)*randn + m(1) + m(2)*X(nancoords(j));
            end
            DATA(:,n) = SS';
        end
    else
        eval('help NAN_SPECTRUM_MAKE_LK')
        error(['The first element of the variable "INTERPMETHOD" must be 0, 1, or one of several',...
            ' specific strings. See the documentation above.'])
    end
elseif strcmp(INTERPMETHOD{2},'both')
    % interp2 (this probably should only be used if you are calculating
    % a 2D wavenumber spectrum, i.e. x/y or x/z, instead of x/t or z/t)
    error(['This option is not yet implemented. See',...
           ' https://github.com/lkach/Spectral_Analysis_Package for',...
           ' updates or feel free to implement your own solution.',...
           ' Otherwise, one of the other otions will probably be',...
           ' sufficient.'])
else
    eval('help nanspectrum2d')
    error('Please format the variable "INTERPMETHOD" correctly.')
end

%% Segmenting
% This allows 1 segment in space if necessary (but not in time)

% N_t = size(DATA,2); % number of points in time
% N_x = size(DATA,1); % number of points in space ("x" is generic space coordinate, not specifically "zonal")

if TSEGMENTS == 1
    error('TWINDOWS cannot = 1')
elseif XSEGMENTS == 1
    if ~mod(N_t,TSEGMENTS)
        % build segments segmented only in time
        
        TSegLeng = N_t/TSEGMENTS; % length of a time segment
        XSegLeng = N_x/XSEGMENTS; % length of a space segment
        Hann = repmat(hanning(TSegLeng)', N_x, 1); % size N_x by TWinLeng
        % chop up the DATA set
        seg_DATA = zeros(N_x,TSegLeng,2*TSEGMENTS - 1);
        for s=1:TSEGMENTS % just the non-overlapping ones (overlapping ones come next)
            seg_DATA(:,:,s) = DATA(:,((s-1)*TSegLeng + 1):(s*TSegLeng));
        end
        for s=1:(TSEGMENTS-1)% overlapping
            seg_DATA(:,:,s+TSEGMENTS) = DATA(:,((s-1)*TSegLeng + 1 + floor(TSegLeng/2)):(s*TSegLeng + floor(TSegLeng/2)));
        end
        
        % Detrend each segment, pre-windowing (probably unimportant, Hanning has no trend).
        xxx = repmat((1:size(seg_DATA,1))', size(seg_DATA,2), 1); ttt = [];
        for i=1:size(seg_DATA,2)
            ttt = [ttt,linspace(i,i,size(seg_DATA,1))];
        end
        A = funct3(ttt,xxx);
        for s=1:size(seg_DATA,3)
            m = A\reshape(seg_DATA(:,:,s), XSegLeng*TSegLeng, 1);
            Deplaned = A*m;
            seg_DATA(:,:,s) = seg_DATA(:,:,s) - reshape(Deplaned,XSegLeng,TSegLeng);
        end
        % Detrend done
        
        hann_seg_DATA = seg_DATA.*Hann; clear seg_DATA % just to free up memory (hopefully this isn't too expensive)
        
    else
        error(['You must be able to divide the length of the time series into ',...
            num2str(TSEGMENTS),' windows.'])
    end
else % i.e. XWINDOWS ~= 1
    if ~mod(N_t,TSEGMENTS) && ~mod(N_x,XSEGMENTS)
        % build segments windowed in time and space, e.g. repmat(hanning(20)',30,1).*hanning(30)
        
        TSegLeng = N_t/TSEGMENTS; % length of a time segment
        XSegLeng = N_x/XSEGMENTS; % length of a space segment
        Hann = repmat(hanning(TSegLeng)', XSegLeng, 1).*hanning(XSegLeng); % size N_x by TWinLeng
        % chop up the DATA set
        NumSegs = (XSEGMENTS*TSEGMENTS) + ((XSEGMENTS-1)*(TSEGMENTS-1)); % Total number of segments
        seg_DATA = zeros(XSegLeng,TSegLeng,NumSegs);
        s = 1;
        for n=1:TSEGMENTS % just the non-overlapping ones (overlapping ones come next)
            for i=1:XSEGMENTS
                seg_DATA(:,:,s) = DATA(((i-1)*XSegLeng + 1):(i*XSegLeng),((n-1)*TSegLeng + 1):(n*TSegLeng));
                s = s + 1;
            end
        end
        s = 1;
        for n=1:(TSEGMENTS-1) % just the non-overlapping ones (overlapping ones come next)
            for i=1:(XSEGMENTS-1)
                seg_DATA(:,:,s+(XSEGMENTS*TSEGMENTS)) = DATA(((i-1)*XSegLeng + 1 + floor(XSegLeng/2)):(i*XSegLeng + floor(XSegLeng/2)),((n-1)*TSegLeng + 1 + floor(TSegLeng/2)):(n*TSegLeng + floor(TSegLeng/2)));
                s = s + 1;
            end
        end

        % Detrend each segment, pre-windowing (probably unimportant, Hanning has no trend).
        xxx = repmat((1:size(seg_DATA,1))', size(seg_DATA,2), 1); ttt = [];
        for i=1:size(seg_DATA,2)
            ttt = [ttt,linspace(i,i,size(seg_DATA,1))];
        end
        A = funct3(ttt,xxx);
        for s=1:size(seg_DATA,3)
            m = A\reshape(seg_DATA(:,:,s), XSegLeng*TSegLeng, 1);
            Deplaned = A*m;
            seg_DATA(:,:,s) = seg_DATA(:,:,s) - reshape(Deplaned,XSegLeng,TSegLeng);
        end
        % Detrend done
        
        hann_seg_DATA = seg_DATA.*Hann; clear seg_DATA % just to free up memory (hopefully this isn't too expensive)
        
    else
        error(['You must be able to divide the length of the time series into ',...
            num2str(TSEGMENTS), ' windows AND the length of the space series into ',...
            num2str(XSEGMENTS),' windows.'])
    end
end

%% Apply FFT

Sheets = size(hann_seg_DATA,3);

fft2_hann_seg_DATA = zeros(size(hann_seg_DATA));
for s = 1:Sheets
    fft2_hann_seg_DATA(:,:,s) = fftshift(fft2(hann_seg_DATA(:,:,s)));
end
clear hann_seg_DATA % just to free up memory (hopefully this isn't too expensive)

%% Abs -> Average -> Normalize

amp = zeros(size(fft2_hann_seg_DATA));
for s = 1:Sheets
    amp(:,:,s) = abs( fft2_hann_seg_DATA(:,:,s) ).^2;
end
clear fft2_hann_seg_DATA % just to free up memory (hopefully this isn't too expensive)

% amp = amp*sqrt(8/3)/(size(amp,1)*size(amp,2));
amp = (8/3)*amp/(numel(amp(:,:,1)));
    warning('Normalization should be correct to satisfy Parselval''s Theorem, but this should be verified for serious use.')

Spectrum2D = mean(amp,3);

%% f_vec and k_vec

% The frequency space in which you're working:
f_s = 1/DT; % sample rate
timelength_record = ((N_t-1)*DT)/TSEGMENTS; % needed to find the lowest frequency
                                            % OF A SEGMENT, not the whole time series
f_Ny = 0.5*f_s;
df = 1/timelength_record;
f_vec = df:df:f_Ny;
    f_vec = [-flip(f_vec), 0, f_vec];

% The wavenumber space in which you're working:
k_s = 1/DX; % sample rate
spacelength_record = ((N_x-1)*DX)/XSEGMENTS; % needed to find the lowest frequency
                                            % OF A SEGMENT, not the whole time series
k_Ny = 0.5*k_s;
dk = 1/spacelength_record;
k_vec = dk:dk:k_Ny;
    k_vec = [-flip(k_vec), 0, k_vec];

%% Trim the left-most column (the redundant one) if TSegLeng is even
% % ^ it's already ok if TSegLeng is odd

if ~mod(TSegLeng,2) % i.e. "IF TSegLeng IS EVEN"
    Spectrum2D = Spectrum2D(:,2:end);
else
end

%% Trim the top row (the redundant one) if XSegLeng is even
% % ^ it's already ok if XSegLeng is odd

if ~mod(XSegLeng,2) % i.e. "IF XSegLeng IS EVEN"
    Spectrum2D = Spectrum2D(2:end,:);
else
end

%% Plotting

if PLOT_BOOLEAN
    % No "figure", that's up to the user to put in when they use this script
    imagesc(f_vec,k_vec,log10(Spectrum2D))
    title('log$_{10}\textbf{S}$','interpreter','latex','fontsize',16)
    xlabel(['frequency (cycles per ',TIME_UNITS,')'],'interpreter','latex','fontsize',16)
    ylabel(['wavenumber (cycles per ',SPACE_UNITS,')'],'interpreter','latex','fontsize',16)
    colormap('jet')
    c = colorbar;
    c.YLabel.String = 'log_1_0S';
    ax = gca; ax.YDir = 'normal';
else
end

%% Error
% These will certainly be different depending on whether or not we segment
% in space, but for now I'll keep like this.
warning(['Error estimate is an educated guess based on 1D spectral principles,',...
        ' inspect code and determine if this is correct and/or sufficient for your needs.'])

% In short, do this calculation for 95% confidence ratio:
NumSegs = (XSEGMENTS*TSEGMENTS) + ((XSEGMENTS-1)*(TSEGMENTS-1)); % Total number of segments
% degrees of freedom; for Hanning windowed overlapping
% This only works because of the Hanning filter

err_high = 2*NumSegs/chi2inv(.05/2,2*NumSegs);
err_low = 2*NumSegs/chi2inv(1-.05/2,2*NumSegs);
err = [err_low err_high];

disp(err)

%% Trouble shooting:
DATA_treated = DATA;
end

%% Compare to a simple fft2:

% figure;subplot(211);imagesc(fftshift(log10(abs(real(fft2(DATA))))));subplot(212);imagesc(fftshift(log10(abs(imag(fft2(DATA))))))

%% Suggested labels for k-l-spectra

% xlabel('$2\pi k = \frac{1}{\lambda_x}$ (cycles per meter)','interpreter','latex','fontsize',14)
% ylabel('$2\pi l = \frac{1}{\lambda_y}$ (cycles per meter)','interpreter','latex','fontsize',14)
