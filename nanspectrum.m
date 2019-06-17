% [Spectrum, f_vec, err] = nanspectrum(TS, DT, TIME_UNITS, SEGMENTS, PLOT_OPTION, PLOT_BOOLEAN, INTERPMETHOD, WINDOWMETHOD)
% 
% Spectrum estimator, able to handle time series with NaN's, Inf's, etc.
% Based on material taught in Sarah Gille's SIOC 221a, notes accessible at
% <http://pordlabs.ucsd.edu/sgille/sioc221a/>.
% This appears to be Welch's method.
% Author: Luke Kachelein
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % 
% % IN:  TS =           time series (row or column vector)
% % IN:  DT =           time step (scalar)
% % IN:  TIME_UNITS =   the units of DT (string, use singular form, e.g. 'day')
% % IN:  SEGMENTS =     The number of segments into which the time series is
% %                     broken for averaging purposes (scalar). If you
% %                     choose 1, then no windowing will apply, and you
% %                     will just be taking the absolute value of the fft
% %                     of the detrended time series.
% % IN:  PLOT_OPTION =  format for the plot (the string argument one would
% %                     use in "plot", e.g. '--r' for a red dashed line.
% %                     enter 0 for '.-' and the default color.
% % IN:  PLOT_BOOLEAN = 1 for an automatically generated plot, 0 for none
% % IN:  INTERPMETHOD = The choice that the called MATLAB function "interp1"
% %                     will use for the interpolation (only matters for data
% %                     with gaps). See "help interp1" for more information.
% %                     Set to 0 for zero-padding (i.e. NaN -> 0, after
% %                     considering trend and offset), set to 1 to replaces
% %                     NaN's with randn noise with the same variance as the
% %                     time series (an issue with this is that it will give
% %                     a different answer each time it's used). Otherwise
% %                     use the appropriate string from "interp1" options.
% % IN:  WINDOWMETHOD = (Optional argument: string or vector)
% %                      - If string: The choice of window function applied
% %                     to each segment (string). If left blank, a "hann"
% %                     window is used, unless SEGMENTS == 1, in which case
% %                     "rectwin" is used. Options include:
% % {bartlett, blackman, boxcar, rectwin, chebwin, hamming, hann, hanning, kaiser, triang}
% %                     It would be possible to use other filters,
% %                     including user-made ones, as long as they are
% %                     defined functions.
% %                      - If vector: The user can enter an arbitrary
% %                     window without needing to define a filter function.
% %                     This vector needs to be exactly as long as the
% %                     segment length = length(TS)/SEGMENTS.
% %                                 
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % 
% % OUT: Spectrum =     The power spectrum (column vector), units of
% %                     [TS_units^2 f_vec_units^-1]
% % OUT: f_vec =        A column vector of the frequencies corresponding to
% %                     the Fourier coefficients in "Spectrum". Units of:
% %                     cycles per TIME_UNITS (not an angular frequency;
% %                     for that, define omega_vec = 2*pi*f_vec)
% % OUT: err =          A two element column vector, where 
% %                     err(1) = "err_low" and err(2) = "err_high".
% %                     See section "%% Error". 95% confidence interval.
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 

function [Spectrum, f_vec, err] = nanspectrum(TS,DT,TIME_UNITS,...
    SEGMENTS,PLOT_OPTION,PLOT_BOOLEAN,INTERPMETHOD,varargin)

%% For troubleshooting (part 1):
% raw_TS = TS;

%% Alert the user if the time series cannot be properly segmented

leng_TS = length(TS);
error_Msg = ['You need to make sure that your time series can',...
        ' actually be divided into ',num2str(SEGMENTS),' windows.',...
        ' As it is given, the time series has ',num2str(leng_TS),...
        ' data points. This can be factored into [',...
        num2str(factor(leng_TS)),'].'];
    
% Give yourself an error if you can't split TS into SEGMENTS segments:
if ~mod(leng_TS,SEGMENTS)
% Using one segment is not advised, because it will just result in the
% absolute value of the fft, which will have a lot of noise, but if you
% do need to use the full time series without segmenting, it would
% probably be best to not dampen the edges and just use 'rectwin' for
% the variable WINDOWMETHOD ("rectwin" gives all ones). A warning will be
% given when WINDOWMETHOD is assigned.
else
    warning(error_Msg)
    error('The time series is of the wrong length.')
end

%% Handle WINDOWMETHOD variable:

N = leng_TS/SEGMENTS;

if nargin == 7
    if SEGMENTS == 1
        WINDOWMETHOD = 'rectwin';
    else
        WINDOWMETHOD = 'hann';
    end
elseif nargin == 8
    WINDOWMETHOD = varargin{1};
    if ischar(WINDOWMETHOD)
        if SEGMENTS == 1
            if strcmp(WINDOWMETHOD,'rectwin') || strcmp(WINDOWMETHOD,'boxcar')
            else
                warning(['You have chosen not to segment your time'...
                    'series. In this case, it is advised (but not'...
                    'required) that you choose "rectwin" or "boxcar" '...
                    'for WINDOWMETHOD, because they do not deemphasize'...
                    'the edges of the time series.'])
            end
        else
        end
    elseif isrow(WINDOWMETHOD) && length(WINDOWMETHOD) == N
        WINDOWMETHOD = WINDOWMETHOD';
    elseif iscolumn(WINDOWMETHOD) && length(WINDOWMETHOD) == N
    else
        error('The optional eighth argument "WINDOWMETHOD" needs to be a string for a window function or vector of length length(TS)/SEGMENTS.')
    end
else
    error('"nanspectrum" only takes 7 or 8 arguments.')
end

%% Make a column into a row

if isrow(TS)
else % reformat if given as a column
    TS = TS';
end

%% Interpolate over NaN's using the chosen method "INTERPMETHOD"

T = 1:leng_TS;
if ischar(INTERPMETHOD)
    TS = interp1(T(isfinite(TS)),TS(isfinite(TS)),T,INTERPMETHOD);
elseif INTERPMETHOD == 0
    % Do a quick linear detrend and demean (not possible in "detrend",
    % and no "nandetrend" available). Solve TS' \approx A*m -> m = A\TS'
    TS_nonan = TS(isfinite(TS)); T_nonan = T(isfinite(TS));
    
    A = ones(length(TS_nonan),2); A(:,2) = T_nonan';
    m = A\(TS_nonan');
    nancoords = find(~isfinite(TS));
    for i=1:length(nancoords)
        TS(nancoords(i)) = m(1) + m(2)*T(nancoords(i));
    end
    
elseif INTERPMETHOD == 1
    % Do a quick linear detrend and demean (not possible in "detrend",
    % and no "nandetrend" available). Solve TS' ? A*m -> m = A\TS'
    TS_nonan = TS(isfinite(TS)); T_nonan = T(isfinite(TS));
    
    A = ones(length(TS_nonan),2); A(:,2) = T_nonan';
    m = A\(TS_nonan');
    TS_nonandetrended = TS_nonan - (A*m)'; % not kept around because we eventually
                                           % detrend the padded segments
    TS_nandetrended = TS; TS_nandetrended(isfinite(TS)) = TS_nonandetrended;
    TS_var = nanvar(TS_nandetrended);
    nancoords = find(~isfinite(TS));
    for i=1:length(nancoords)
        TS(nancoords(i)) = sqrt(TS_var)*randn + m(1) + m(2)*T(nancoords(i));
    end
else
    eval('help nanspectrum')
    error(['The variable "INTERPMETHOD" must be 0, 1, or one of several',...
        ' specific strings. See the documentation above.'])
end

%% Segment (reduce noise), window (account for edge effects), and calculated "Spectrum"

% Build window matrix (The total number of windows including
% overlapping ones is 2*SEGMENTS-1)
if ischar(WINDOWMETHOD)
    eval(['Window = repmat(',WINDOWMETHOD,'(N), 1, 2*SEGMENTS - 1);'])
else % There is no need to ask if it's a vector of the right length; this
     % was already checked before.
    Window = repmat(WINDOWMETHOD, 1, 2*SEGMENTS - 1);
end

% The normalization factor for the spectrum (to ensure parseval's theorem
% is satisfied) is 1/(the mean of the square of the filter):
norm_factor = 1/mean([0;Window(:,1)].^2);
% ^ I am not sure if this should possibly be "1/mean(Window(:,1).^2)"
% instead; however, I chose the extra 0 because it results in 8/3 exactly
% if the Window is built using "hanning" (this 8/3 agrees with the
% literature, Bendat & Piersol). I am guessing that this has to do with the
% window being periodic and only hitting zero once per cycle.

% We do not want to discount the very beginning and end of the time series,
% which if we use "Window" as is, will surpress the contribution of the
% first segLen/2 points and last segLen/2 points compared to all of the
% others. However, if we decided to not surpress them (e.g. by defining the
% first half of the first window as being all ones then the second half
% being the window function, and the last window being the first half of
% the window function, followed by all ones), we get significant peak
% spreading, so we will just have to live with the disadvantage of the
% first few and last few data being deemphasized. In the case that there
% are a lot of windows, this effect will be pretty small anyway.

% Divide up the time series into rows in a matrix:
% The old way was to use "buffer", which was imperfect:
% seg_TS = buffer(TS, leng_TS/SEGMENTS,floor(0.5*leng_TS/SEGMENTS));
% Therefore, we build the matrix "seg_TS" by fundamental means:
seg_TS = zeros(N, 2*SEGMENTS - 1); % Initialize
for i = 1:2:(2*SEGMENTS - 1) % Odd numbered columns
    seg_TS(:,i) = TS( (1 + ((i-1)/2)*N):(((i+1)/2)*N) );
end
Offset = floor(0.5*N); % Midpoint at which the even-numbered columns
                       % will start; trivial for even leng_TS, and for
                       % odd leng_TS, it starts in the middle.
for i = 1:2:(2*SEGMENTS - 2) % Even numbered columns
    seg_TS(:,i+1) = TS( (Offset + 1 + ((i-1)/2)*N):(Offset + ((i+1)/2)*N) );
end
% ^ This works for odd and even time series and segments

seg_TS = detrend(seg_TS); % Column-wise by default
% From Sarah's lecture 9:
% "First you must demean your data, otherwise, the window will shift energy
% from the mean into other frequencies. If you're working in segments, you
% should demean (and detrend) each segment before you do anything further."
% Element-wise mulyiply:
window_seg_TS = Window.*seg_TS;

%%%%%%%%%%% fft:
fft_window_seg_TS = fft(window_seg_TS);
%%%%%%%%%%% mod squared:
if mod(N,2) == 0 % even N
    amp = abs(fft_window_seg_TS(1:(N/2+1),:)).^2;
elseif mod(N,2) == 1 % odd N
    amp = abs(fft_window_seg_TS(1:((N+1)/2),:)).^2;
else
    error('Fatal error, figure it out yourself')
end
%%%%%%%%%%% add back lost energy
if mod(N,2) == 0 % even N
    amp(2:(end-1),:) = 2*amp(2:(end-1),:); % for even N
elseif mod(N,2) == 1 % odd N
    amp(2:end,:) = 2*amp(2:end,:); % for odd N
else
    error('Fatal error, figure it out yourself')
end

%%%%%%%%%%% The frequency space in which you're working:
f_s = 1/DT; % sample rate
timelength_record = ((leng_TS)*DT)/SEGMENTS;
% ^ Lowest frequency OF A WINDOW, not the whole time series

f_Ny = 0.5*f_s;
df = 1/timelength_record;
%%%%%%%%%%% Build frequency vector:
if mod(N,2) == 0 % even N
    f_vec = (df:df:f_Ny)';
elseif mod(N,2) == 1 % odd N
    f_vec = (df:df:(f_Ny - df/2))';
    % Because f_Ny is not a multiple of df if N is odd, but (f_Ny - df/2) is
else
    error('Fatal error, figure it out yourself')
end

amp = norm_factor*amp*DT/N;
% Divide by N because of the normalization which MATLAB uses in "fft". The
% DT ensures that any unit for time will work for normalizing.

% The spectrum:
Spectrum = mean(amp(2:end,:),2);

%% Error

% In short, do this calculation for 95% confidence ratio:
EffectiveSegments = 2*SEGMENTS - 1;% degrees of freedom; for overlapping
                     % Hann windowed segments, this is the total number of
                     % segments, i.e.:
                     % 2*SEGMENTS - 1
                     % This may only work perfectly with the Hann window.
                     % With other windows, it's not entirely clear, but my
                     % guess is, if the window is "Hann-like", it will not
                     % make a big difference.

err_high = 2*EffectiveSegments/chi2inv(.05/2,2*EffectiveSegments);
err_low = 2*EffectiveSegments/chi2inv(1-.05/2,2*EffectiveSegments);
err = [err_low err_high];

%% Plot
if PLOT_BOOLEAN
    % sort out blank plot option
    if PLOT_OPTION == 0
        PLOT_OPTION = '.-';
    else
    end

    % You need to provide "figure" already; this is to facilitate plotting
    % several spectra together called only this function.
    loglog(f_vec,Spectrum,PLOT_OPTION); hold on
    xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
    ylabel(['Spectral density [(time\_series\_units)$^2$ (cycles per ',...
            TIME_UNITS,')$^{-1}$]'],'Interpreter','LaTeX')
    
    % error bar plotted
    loglog([1.1 1.1]*f_vec(end),[1 (err_high/err_low)]*min(Spectrum),'k');
    
else
end
% figure;imagesc(Window)
%% Troubleshooting (part 2)
% % Good way to check that your time series are behaving properly:
% figure
% subplot(211);plot(raw_TS);title('raw')
% subplot(212);plot(TS);title('post-processing')

end

