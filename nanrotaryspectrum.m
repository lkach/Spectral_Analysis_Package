% [S_cw, S_ccw, f_vec, err] =
% nanrotaryspectrum(u, v, DT, TIME_UNITS, SEGMENTS, PLOT_OPTION, PLOT_BOOLEAN, INTERPMETHOD, varargin)
% 
% Rotary spectrum estimator, able to handle time series with NaN's, Inf's,
% etc. Based on material taught in Sarah Gille's SIOC 221a, notes
% accessible at:
% <http://pordlabs.ucsd.edu/sgille/sioc221a/>.
% This appears to be Welch's method.
% Author: Luke Kachelein
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % 
% % IN:  u =            x-velocity time series (row or column vector)
% % IN:  v =            y-velocity time series (row or column vector)
% %                     length(u) must be equal to length(v)
% % IN:  DT =           time step (scalar)
% % IN:  TIME_UNITS =   the units of DT (string, use singular form, e.g. 'day')
% % IN:  SEGMENTS =     The number of segments into which the time series is
% %                     broken for averaging purposes (scalar). If you
% %                     choose 1, then no windowing will apply, and you
% %                     will just be taking the absolute value of the fft
% %                     of the detrended time series.
% % IN:  PLOT_OPTION =  format for the plot (the string arguments one would
% %                     give for CW and CCW spectra, in the format of a cell,
% %                     e.g. {'b.-','g.-'} for blue connected dots (CW) and
% %                     green connected dots (CCW). Enter 0 if you want it
% %                     to be {'.-','.-'}
% % IN:  PLOT_BOOLEAN = 1 for an automatically generated plot, 0 for none
% % IN:  INTERPMETHOD = The choice that the called MATLAB function "interp1"
% %                     will use for the interpolation (only matters for
% %                     data with gaps). See "help interp1" for more information.
% %                     Set to 0 for zero-padding (i.e. NaN -> 0), set to
% %                     1 to replaces NaN's with randn noise with the same
% %                     variance as the time series (an issue with this is
% %                     that it will give a different answer each time it's
% %                     used). Otherwise use the appropriate string from
% %                     "interp1" options.
% % IN:  WINDOWMETHOD = (Optional) The choice of window function applied to
% %                     each segment (strong). If left blank, a "hann"
% %                     window is used, unless SEGMENTS == 1, in which case
% %                     "rectwin" is used. Options include:
% % {bartlett, blackman, boxcar, rectwin, chebwin, hamming, hann, hanning, kaiser, triang}
% %                     It would be possible to use other filters,
% %                     including user-made ones, as long as they are
% %                     defined functions.
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % 
% % OUT: S_cw =         The clockwise power spectrum (column vector).
% %                     Units of:
% %                     [TS_units^2 f_vec_units^-1]
% % OUT: S_ccw =        The counter-clockwise power spectrum (column vector).
% %                     This has the same units as S_cw, of course.
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

function [S_cw, S_ccw, f_vec, err] = nanrotaryspectrum(u,v,DT,TIME_UNITS,...
    SEGMENTS,PLOT_OPTION,PLOT_BOOLEAN,INTERPMETHOD,varargin)

%% For troubleshooting (part 1):
% raw_TS = TS;

%% Alert the user if the time series cannot be properly segmented

leng_TS = length(u);

if length(v) ~= leng_TS
    error('u and v must be the same length')
else
end

error_Msg = ['You need to make sure that your time series can',...
        ' actually be divided into ',num2str(SEGMENTS),' windows.',...
        ' As it is given, the time series each have ',num2str(leng_TS),...
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
    error('The time series are of the wrong length.')
end

%% Before wasting computation, warn user if PLOT_OPTION is incorrectly formatted
if isnumeric(PLOT_OPTION) && PLOT_OPTION == 0
    PLOT_OPTION = {'.-','.-'};
elseif iscell(PLOT_OPTION) && length(PLOT_OPTION) == 2
else
    eval('help nanrotaryspectrum')
    error('The sixth argument "PLOT_OPTION" is incorrectly formatted. See documentation above.')
end

%% Handle WINDOWMETHOD variable:

N = leng_TS/SEGMENTS;

if nargin == 8
    if SEGMENTS == 1
        WINDOWMETHOD = 'rectwin';
    else
        WINDOWMETHOD = 'hann';
    end
elseif nargin == 9
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
    error('"nanrotaryspectrum" only takes 8 or 9 arguments.')
end

%% Make a column into a row

if isrow(u)
else % reformat if given as a column
    u = u';
end
if isrow(v)
else % reformat if given as a column
    v = v';
end

%% Interpolate over NaN's using the chosen method "INTERPMETHOD"

T = 1:leng_TS;
if ischar(INTERPMETHOD)
    u = interp1(T(isfinite(u)),u(isfinite(u)),T,INTERPMETHOD);
    v = interp1(T(isfinite(v)),v(isfinite(v)),T,INTERPMETHOD);
elseif INTERPMETHOD == 0
    % Do a quick linear detrend and demean (not possible in "detrend",
    % and no "nandetrend" available). Solve TS' \approx A*m -> m = A\TS'
% % % % % % % % % % % % % % % %     U     % % % % % % % % % % % % % % % % %
    u_nonan = u(isfinite(u)); Tu_nonan = T(isfinite(u));
    
    Au = ones(length(u_nonan),2); Au(:,2) = Tu_nonan';
    mu = Au\(u_nonan');
    nancoords_u = find(~isfinite(u));
    for i=1:length(nancoords_u)
        u(nancoords_u(i)) = mu(1) + mu(2)*T(nancoords_u(i));
    end
% % % % % % % % % % % % % % % %     V     % % % % % % % % % % % % % % % % %
    v_nonan = v(isfinite(v)); Tv_nonan = T(isfinite(v));
    
    Av = ones(length(v_nonan),2); Av(:,2) = Tv_nonan';
    mv = Av\(v_nonan');
    nancoords_v = find(~isfinite(v));
    for i=1:length(nancoords_v)
        v(nancoords_v(i)) = mv(1) + mv(2)*T(nancoords_v(i));
    end
    
elseif INTERPMETHOD == 1
    % Do a quick linear detrend and demean (not possible in "detrend",
    % and no "nandetrend" available). Solve TS' ? A*m -> m = A\TS'
    
% % % % % % % % % % % % % % % %     U     % % % % % % % % % % % % % % % % %
    u_nonan = u(isfinite(u)); Tu_nonan = T(isfinite(u));
    
    Au = ones(length(u_nonan),2); Au(:,2) = Tu_nonan';
    mu = Au\(u_nonan');
    u_nonandetrended = u_nonan - (Au*mu)'; % not kept around because we eventually
                                           % detrend the padded segments
    u_nandetrended = u; u_nandetrended(isfinite(u)) = u_nonandetrended;
    u_var = nanvar(u_nandetrended);
    nancoords_u = find(~isfinite(u));
    for i=1:length(nancoords_u)
        u(nancoords_u(i)) = sqrt(u_var)*randn + mu(1) + mu(2)*T(nancoords_u(i));
    end
% % % % % % % % % % % % % % % %     V     % % % % % % % % % % % % % % % % %
    v_nonan = v(isfinite(v)); Tv_nonan = T(isfinite(v));
    
    Av = ones(length(v_nonan),2); Av(:,2) = Tv_nonan';
    mv = Av\(v_nonan');
    v_nonandetrended = v_nonan - (Av*mv)'; % not kept around because we eventually
                                           % detrend the padded segments
    v_nandetrended = v; v_nandetrended(isfinite(v)) = v_nonandetrended;
    v_var = nanvar(v_nandetrended);
    nancoords_v = find(~isfinite(v));
    for i=1:length(nancoords_v)
        v(nancoords_v(i)) = sqrt(v_var)*randn + mv(1) + mv(2)*T(nancoords_v(i));
    end
else
    eval('help nanrotaryspectrum')
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
seg_u = zeros(N, 2*SEGMENTS - 1); % Initialize
seg_v = zeros(N, 2*SEGMENTS - 1);
for i = 1:2:(2*SEGMENTS - 1) % Odd numbered columns
    seg_u(:,i) = u( (1 + ((i-1)/2)*N):(((i+1)/2)*N) );
    seg_v(:,i) = v( (1 + ((i-1)/2)*N):(((i+1)/2)*N) );
end
Offset = floor(0.5*N); % Midpoint at which the even-numbered columns
                       % will start; trivial for even leng_TS, and for
                       % odd leng_TS, it starts in the middle.
for i = 1:2:(2*SEGMENTS - 2) % Even numbered columns
    seg_u(:,i+1) = u( (Offset + 1 + ((i-1)/2)*N):(Offset + ((i+1)/2)*N) );
    seg_v(:,i+1) = v( (Offset + 1 + ((i-1)/2)*N):(Offset + ((i+1)/2)*N) );
end
% ^ This works for odd and even time series and segments

seg_u = detrend(seg_u); % Column-wise by default
seg_v = detrend(seg_v);
% From Sarah's lecture 9:
% "First you must demean your data, otherwise, the window will shift energy
% from the mean into other frequencies. If you're working in segments, you
% should demean (and detrend) each segment before you do anything further."
% Element-wise mulyiply:
window_seg_u = Window.*seg_u;
window_seg_v = Window.*seg_v;

%%%%%%%%%%% fft:
fft_window_seg_u = fft(window_seg_u);
fft_window_seg_v = fft(window_seg_v);

%%%%%%%%%%% Isolate imaginary and real components:
if mod(N,2) == 0 % even N
    A = real(fft_window_seg_u(1:(N/2+1),:)); % for even N
    B = imag(fft_window_seg_u(1:(N/2+1),:));
    C = real(fft_window_seg_v(1:(N/2+1),:));
    D = imag(fft_window_seg_v(1:(N/2+1),:));
elseif mod(N,2) == 1 % odd N
    A = real(fft_window_seg_u(1:(N+1)/2,:)); % for even N
    B = imag(fft_window_seg_u(1:(N+1)/2,:));
    C = real(fft_window_seg_v(1:(N+1)/2,:));
    D = imag(fft_window_seg_v(1:(N+1)/2,:));
else
    error('Fatal error, figure it out yourself')
end

%%%%%%%%%%% Combine into the CW and CCW parts
fft_window_seg_CW = 0.5*(A + D + 1i*(C - B));
fft_window_seg_CCW  = 0.5*(A - D + 1i*(C + B));


%% CW

%%%%%%%%%%% mod squared:
if mod(N,2) == 0 % even N
    ampCW = abs(fft_window_seg_CW(1:(N/2+1),:)).^2; % for even N
elseif mod(N,2) == 1 % odd N
    ampCW = abs(fft_window_seg_CW(1:((N+1)/2),:)).^2; % for odd N
else
    error('Fatal error, figure it out yourself')
end
%%%%%%%%%%% add back lost energy
if mod(N,2) == 0 % even N
    ampCW(2:end-1,:) = 2*ampCW(2:(end-1),:); % for even N
elseif mod(N,2) == 1 % odd N
    ampCW(2:end,:) = 2*ampCW(2:end,:); % for odd N
else
    error('Fatal error, figure it out yourself')
end
ampCW = norm_factor*ampCW*DT/N;
% Divide by N because of the normalization which MATLAB uses in "fft". The
% DT ensures that any unit for time will work for normalizing.

% The spectrum:
S_cw = mean(ampCW(2:end,:),2);

%% CCW

%%%%%%%%%%% mod squared:
if mod(N,2) == 0 % even N
    ampCCW = abs(fft_window_seg_CCW(1:(N/2+1),:)).^2; % for even N
elseif mod(N,2) == 1 % odd N
    ampCCW = abs(fft_window_seg_CCW(1:((N+1)/2),:)).^2; % for odd N
else
    error('Fatal error, figure it out yourself')
end
%%%%%%%%%%% add back lost energy
if mod(N,2) == 0 % even N
    ampCCW(2:end-1,:) = 2*ampCCW(2:(end-1),:); % for even N
elseif mod(N,2) == 1 % odd N
    ampCCW(2:end,:) = 2*ampCCW(2:end,:); % for odd N
else
    error('Fatal error, figure it out yourself')
end
ampCCW = norm_factor*ampCCW*DT/N;
% Divide by N because of the normalization which MATLAB uses in "fft". The
% DT ensures that any unit for time will work for normalizing.

% The spectrum:
S_ccw = mean(ampCCW(2:end,:),2);

%% Frequency

%%%%%%%%%%% The frequency space in which you're working:
f_s = 1/DT; % sample rate
timelength_record = ((leng_TS)*DT)/SEGMENTS;
% ^ Lowest frequency OF A WINDOW, not the whole time series

f_Ny = 0.5*f_s; %
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

%% Error

% In short, do this calculation for 95% confidence ratio:
EffectiveSegments = 2*SEGMENTS - 1;% degrees of freedom; for Hanning windowed overlapping
                     % windows, this is the total number of windows, i.e.
                     % 2*WINDOWS - 1
                     % This only works because of the Hanning filter

err_high = 2*EffectiveSegments/chi2inv(.05/2,2*EffectiveSegments);
err_low = 2*EffectiveSegments/chi2inv(1-.05/2,2*EffectiveSegments);
err = [err_low err_high];

%% Plot
if PLOT_BOOLEAN
    
    % You need to provide "figure" already; this is to facilitate plotting
    % several spectra together called only this function.
    loglog(f_vec,S_cw,PLOT_OPTION{1}); hold on
    loglog(f_vec,S_ccw,PLOT_OPTION{2})
    xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
    ylabel(['Spectral density [(time\_series\_units)$^2$ (cycles per ',...
            TIME_UNITS,')$^{-1}$]'],'Interpreter','LaTeX')
    
    % error bar plotted
    loglog([1.1 1.1]*f_vec(length(S_cw)),[1 (err_high/err_low)]*min([min(S_cw) min(S_ccw)]),'k');
    legend('CW','CCW','95% confidence')
    
else
end

%% Troubleshooting (part 2)
% % Good way to check that your time series are behaving properly:
% figure
% subplot(211);plot(raw_TS);title('raw')
% subplot(212);plot(TS);title('post-processing')

end
%% Inertial currents (CW in NH, CCW in SH)

% see:
% Chapter 7: Dynamical Processes for Descriptive Ocean Circulation - 2011 Descriptive Physical Oceanography Sixth Edition
% Lynne Talley
