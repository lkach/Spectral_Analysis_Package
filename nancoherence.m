% [f_vec,Coherence_linear, Coherence_rotary, Phase_linear, Phase_rotary,Beta, S_u_v1,S_u_v2,S_cw_ccw1,S_cw_ccw2,err] =
% NANCOHERENCE(w1, w2, DT, TIME_UNITS, SEGMENTS, PLOT_OPTION, PLOT_BOOLEAN, INTERPMETHOD, WINDOWMETHOD, Alpha)
%                                                                                         ^ optional ^ ^ ^ ^
% 
% Linear and rotary spectrum and coherence estimator, able to handle time
% series with NaN's, Inf's, etc. Based on material taught in Sarah Gille's
% SIOC 221a, notes accessible at:
% <http://pordlabs.ucsd.edu/sgille/sioc221a/>.
% Welch's method
% Author: Luke Kachelein
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % 
% % IN:  w1 =           first complex time series  = u1 + i*v1 (row or column vector)
% % IN:  w2 =           second complex time series = u2 + i*v2 (row or column vector)
% %                     length(w1) must be equal to length(w2)
% %                     NOTE: if both w1 and w2 are real, then any figures
% %                     requested will only show results for the "u"
% %                     component, and all figures that would have shown
% %                     rotary quantities are not generated.
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
% % IN:  PLOT_BOOLEAN = 1 for automatically generated plots for spectra and
% %                     coherence, 0 for none
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
% % IN:  Alpha =        (Optional, but if Alpha is included you must also
% %                     include WINDOWMETHOD, i.e. Alpha is 10th argument)
% %                     The threshold fraction at
% %                     which we are willing to tolerate random chance for
% %                     claiming that coherence exists. E.G. Alpha = 0.05
% %                     means that for Coherence == Beta (see output
% %                     variables), there is a 5% chance that the coherence
% %                     is due to random noise, and a 95% chance that it's
% %                     meaningful. Default Alpha = 0.05.
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % 
% % OUT: f_vec =        A column vector of the frequencies corresponding to
% %                     the Fourier coefficients in "Spectrum". Units of:
% %                     cycles per TIME_UNITS (not an angular frequency;
% %                     for that, define omega_vec = 2*pi*f_vec)
% % Coherence_linear =  [Coh_u1u2,   Coh_v1v2,     Coh_u1v2,    Coh_v1u2   ]
% % Coherence_rotary =  [Coh_CW1CW2, Coh_CCW1CCW2, Coh_CW1CCW2, Coh_CCW1CW2]
% % Phase_linear =      [Pha_u1u2,   Pha_v1v2,     Pha_u1v2,    Pha_v1u2   ]
% % Phase_rotary =      [Pha_CW1CW2, Pha_CCW1CCW2, Pha_CW1CCW2, Pha_CCW1CW2]
% % Beta =              Threshold value for coherence corresponding to Alpha
% % OUT: S_u_v1 =       Spectrum of w1, = [S_u1 S_v1]
% % OUT: S_u_v2 =       Spectrum of w2, = [S_u2 S_v2]
% % OUT: S_cw_ccw1 =    Rotary spectrum of w1, = [S_cw1 S_ccw1]
% %                     [TS_units^2 f_vec_units^-1]
% % OUT: S_cw_ccw2 =    Rotary spectrum of w2, = [S_cw2 S_ccw2]
% %                     [TS_units^2 f_vec_units^-1]
% % OUT: err =          A two element column vector, where 
% %                     err(1) = "err_low" and err(2) = "err_high".
% %                     See section "%% Error". 95% confidence interval.
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% Coherence is in the range [0 1]
% Phase is in the range [-pi pi]
% 

function [f_vec,Coherence_linear, Coherence_rotary, Phase_linear, Phase_rotary,Beta, S_u_v1,S_u_v2,S_cw_ccw1,S_cw_ccw2,err] = ...
         nancoherence(w1,w2,DT,TIME_UNITS,SEGMENTS,PLOT_OPTION,PLOT_BOOLEAN,INTERPMETHOD,varargin)

if nargin == 0
    disp('  [f_vec,Coherence_linear, Coherence_rotary, Phase_linear, Phase_rotary, S_u_v1,S_u_v2,S_cw_ccw1,S_cw_ccw2,err] =')
    disp('  nancoherence(w1, w2, DT, TIME_UNITS, SEGMENTS, PLOT_OPTION, PLOT_BOOLEAN, INTERPMETHOD, WINDOWMETHOD)')
else


%% Take complex inputs and convert to u and v

u1 = real(w1);
v1 = imag(w1);

u2 = real(w2);
v2 = imag(w2);

%% For troubleshooting (part 1):
% raw_TS = TS;

%% Alert the user if the time series cannot be properly segmented

leng_TS = length(u1);

if length(v1) ~= leng_TS && length(v2) ~= leng_TS
    error('first and second u+iv must be the same length')
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
elseif nargin > 8
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
end

Alpha = 0.05; % default
if nargin == 10
    Alpha = varargin{2};
elseif nargin > 10
    error('"nancoherence" only takes 8, 9, or 10 arguments.')
else
end

%% Make a column into a row

if isrow(u1)
else % reformat if given as a column
    u1 = u1';
end
if isrow(v1)
else % reformat if given as a column
    v1 = v1';
end

if isrow(u2)
else % reformat if given as a column
    u2 = u2';
end
if isrow(v2)
else % reformat if given as a column
    v2 = v2';
end

%% W1 START
%% Interpolate over NaN's using the chosen method "INTERPMETHOD"

T = 1:leng_TS;
if ischar(INTERPMETHOD)
    u1 = interp1(T(isfinite(u1)),u1(isfinite(u1)),T,INTERPMETHOD);
    v1 = interp1(T(isfinite(v1)),v1(isfinite(v1)),T,INTERPMETHOD);
elseif INTERPMETHOD == 0
    % Do a quick linear detrend and demean (not possible in "detrend",
    % and no "nandetrend" available). Solve TS' \approx A*m -> m = A\TS'
% % % % % % % % % % % % % % % %     U     % % % % % % % % % % % % % % % % %
    u1_nonan = u1(isfinite(u1)); Tu1_nonan = T(isfinite(u1));
    
    Au1 = ones(length(u1_nonan),2); Au1(:,2) = Tu1_nonan';
    mu1 = Au1\(u1_nonan');
    nancoords_u1 = find(~isfinite(u1));
    for i=1:length(nancoords_u1)
        u1(nancoords_u1(i)) = mu1(1) + mu1(2)*T(nancoords_u1(i));
    end
% % % % % % % % % % % % % % % %     V     % % % % % % % % % % % % % % % % %
    v1_nonan = v1(isfinite(v1)); Tv1_nonan = T(isfinite(v1));
    
    Av1 = ones(length(v1_nonan),2); Av1(:,2) = Tv1_nonan';
    mv1 = Av1\(v1_nonan');
    nancoords_v1 = find(~isfinite(v1));
    for i=1:length(nancoords_v1)
        v1(nancoords_v1(i)) = mv1(1) + mv1(2)*T(nancoords_v1(i));
    end
    
elseif INTERPMETHOD == 1
    % Do a quick linear detrend and demean (not possible in "detrend",
    % and no "nandetrend" available). Solve TS' ? A*m -> m = A\TS'
    
% % % % % % % % % % % % % % % %     U     % % % % % % % % % % % % % % % % %
    u1_nonan = u1(isfinite(u1)); Tu1_nonan = T(isfinite(u1));
    
    Au1 = ones(length(u1_nonan),2); Au1(:,2) = Tu1_nonan';
    mu1 = Au1\(u1_nonan');
    u1_nonandetrended = u1_nonan - (Au1*mu1)'; % not kept around because we eventually
                                           % detrend the padded segments
    u1_nandetrended = u1; u1_nandetrended(isfinite(u1)) = u1_nonandetrended;
    u1_var = nanvar(u1_nandetrended);
    nancoords_u1 = find(~isfinite(u1));
    for i=1:length(nancoords_u1)
        u1(nancoords_u1(i)) = sqrt(u1_var)*randn + mu1(1) + mu1(2)*T(nancoords_u1(i));
    end
% % % % % % % % % % % % % % % %     V     % % % % % % % % % % % % % % % % %
    v1_nonan = v1(isfinite(v1)); Tv1_nonan = T(isfinite(v1));
    
    Av1 = ones(length(v1_nonan),2); Av1(:,2) = Tv1_nonan';
    mv1 = Av1\(v1_nonan');
    v1_nonandetrended = v1_nonan - (Av1*mv1)'; % not kept around because we eventually
                                           % detrend the padded segments
    v1_nandetrended = v1; v1_nandetrended(isfinite(v1)) = v1_nonandetrended;
    v1_var = nanvar(v1_nandetrended);
    nancoords_v1 = find(~isfinite(v1));
    for i=1:length(nancoords_v1)
        v1(nancoords_v1(i)) = sqrt(v1_var)*randn + mv1(1) + mv1(2)*T(nancoords_v1(i));
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
seg_u1 = zeros(N, 2*SEGMENTS - 1); % Initialize
seg_v1 = zeros(N, 2*SEGMENTS - 1);
for i = 1:2:(2*SEGMENTS - 1) % Odd numbered columns
    seg_u1(:,i) = u1( (1 + ((i-1)/2)*N):(((i+1)/2)*N) );
    seg_v1(:,i) = v1( (1 + ((i-1)/2)*N):(((i+1)/2)*N) );
end
Offset = floor(0.5*N); % Midpoint at which the even-numbered columns
                       % will start; trivial for even leng_TS, and for
                       % odd leng_TS, it starts in the middle.
for i = 1:2:(2*SEGMENTS - 2) % Even numbered columns
    seg_u1(:,i+1) = u1( (Offset + 1 + ((i-1)/2)*N):(Offset + ((i+1)/2)*N) );
    seg_v1(:,i+1) = v1( (Offset + 1 + ((i-1)/2)*N):(Offset + ((i+1)/2)*N) );
end
% ^ This works for odd and even time series and segments

seg_u1 = detrend(seg_u1); % Column-wise by default
seg_v1 = detrend(seg_v1);
% From Sarah's lecture 9:
% "First you must demean your data, otherwise, the window will shift energy
% from the mean into other frequencies. If you're working in segments, you
% should demean (and detrend) each segment before you do anything further."
% Element-wise mulyiply:
window_seg_u1 = Window.*seg_u1;
window_seg_v1 = Window.*seg_v1;

%%%%%%%%%%% fft:
fft_window_seg_u1 = fft(window_seg_u1);
fft_window_seg_v1 = fft(window_seg_v1);

%%%%%%%%%%% Isolate imaginary and real components:
if mod(N,2) == 0 % even N
    A1 = real(fft_window_seg_u1(1:(N/2+1),:)); % for even N
    B1 = imag(fft_window_seg_u1(1:(N/2+1),:));
    C1 = real(fft_window_seg_v1(1:(N/2+1),:));
    D1 = imag(fft_window_seg_v1(1:(N/2+1),:));
elseif mod(N,2) == 1 % odd N
    A1 = real(fft_window_seg_u1(1:(N+1)/2,:)); % for even N
    B1 = imag(fft_window_seg_u1(1:(N+1)/2,:));
    C1 = real(fft_window_seg_v1(1:(N+1)/2,:));
    D1 = imag(fft_window_seg_v1(1:(N+1)/2,:));
else
    
end

%%%%%%%%%%% Combine into the CW and CCW parts
fft_window_seg_CW1 = 0.5*(A1 + D1 + 1i*(C1 - B1));
fft_window_seg_CCW1  = 0.5*(A1 - D1 + 1i*(C1 + B1));
% ^ backwards from the sin-cos formulation

%% CW

%%%%%%%%%%% mod squared:
if mod(N,2) == 0 % even N
    ampCW1 = abs(fft_window_seg_CW1(1:(N/2+1),:)).^2; % for even N
elseif mod(N,2) == 1 % odd N
    ampCW1 = abs(fft_window_seg_CW1(1:((N+1)/2),:)).^2; % for odd N
else
    
end
%%%%%%%%%%% add back lost energy
if mod(N,2) == 0 % even N
    ampCW1(2:end-1,:) = 2*ampCW1(2:(end-1),:); % for even N
elseif mod(N,2) == 1 % odd N
    ampCW1(2:end,:) = 2*ampCW1(2:end,:); % for odd N
else
    
end
ampCW1 = norm_factor*ampCW1*DT/N;
% Divide by N because of the normalization which MATLAB uses in "fft". The
% DT ensures that any unit for time will work for normalizing.

% The spectrum:
S_cw1 = mean(ampCW1(2:end,:),2);

%% CCW

%%%%%%%%%%% mod squared:
if mod(N,2) == 0 % even N
    ampCCW1 = abs(fft_window_seg_CCW1(1:(N/2+1),:)).^2; % for even N
elseif mod(N,2) == 1 % odd N
    ampCCW1 = abs(fft_window_seg_CCW1(1:((N+1)/2),:)).^2; % for odd N
else
    
end
%%%%%%%%%%% add back lost energy
if mod(N,2) == 0 % even N
    ampCCW1(2:end-1,:) = 2*ampCCW1(2:(end-1),:); % for even N
elseif mod(N,2) == 1 % odd N
    ampCCW1(2:end,:) = 2*ampCCW1(2:end,:); % for odd N
else
    
end
ampCCW1 = norm_factor*ampCCW1*DT/N;
% Divide by N because of the normalization which MATLAB uses in "fft". The
% DT ensures that any unit for time will work for normalizing.

% The spectrum:
S_ccw1 = mean(ampCCW1(2:end,:),2);

%% W1 END
%% W2 START
%% Interpolate over NaN's using the chosen method "INTERPMETHOD"

if ischar(INTERPMETHOD)
    u2 = interp1(T(isfinite(u2)),u2(isfinite(u2)),T,INTERPMETHOD);
    v2 = interp1(T(isfinite(v2)),v2(isfinite(v2)),T,INTERPMETHOD);
elseif INTERPMETHOD == 0
    % Do a quick linear detrend and demean (not possible in "detrend",
    % and no "nandetrend" available). Solve TS' \approx A*m -> m = A\TS'
% % % % % % % % % % % % % % % %     U     % % % % % % % % % % % % % % % % %
    u2_nonan = u2(isfinite(u2)); Tu2_nonan = T(isfinite(u2));
    
    Au2 = ones(length(u2_nonan),2); Au2(:,2) = Tu2_nonan';
    mu2 = Au2\(u2_nonan');
    nancoords_u2 = find(~isfinite(u2));
    for i=1:length(nancoords_u2)
        u2(nancoords_u2(i)) = mu2(1) + mu2(2)*T(nancoords_u2(i));
    end
% % % % % % % % % % % % % % % %     V     % % % % % % % % % % % % % % % % %
    v2_nonan = v2(isfinite(v2)); Tv2_nonan = T(isfinite(v2));
    
    Av2 = ones(length(v2_nonan),2); Av2(:,2) = Tv2_nonan';
    mv2 = Av2\(v2_nonan');
    nancoords_v2 = find(~isfinite(v2));
    for i=1:length(nancoords_v2)
        v2(nancoords_v2(i)) = mv2(1) + mv2(2)*T(nancoords_v2(i));
    end
    
elseif INTERPMETHOD == 1
    % Do a quick linear detrend and demean (not possible in "detrend",
    % and no "nandetrend" available). Solve TS' ? A*m -> m = A\TS'
    
% % % % % % % % % % % % % % % %     U     % % % % % % % % % % % % % % % % %
    u2_nonan = u2(isfinite(u2)); Tu2_nonan = T(isfinite(u2));
    
    Au2 = ones(length(u2_nonan),2); Au2(:,2) = Tu2_nonan';
    mu2 = Au2\(u2_nonan');
    u2_nonandetrended = u2_nonan - (Au2*mu2)'; % not kept around because we eventually
                                           % detrend the padded segments
    u2_nandetrended = u2; u2_nandetrended(isfinite(u2)) = u2_nonandetrended;
    u2_var = nanvar(u2_nandetrended);
    nancoords_u2 = find(~isfinite(u2));
    for i=1:length(nancoords_u2)
        u2(nancoords_u2(i)) = sqrt(u2_var)*randn + mu2(1) + mu2(2)*T(nancoords_u2(i));
    end
% % % % % % % % % % % % % % % %     V     % % % % % % % % % % % % % % % % %
    v2_nonan = v2(isfinite(v2)); Tv2_nonan = T(isfinite(v2));
    
    Av2 = ones(length(v2_nonan),2); Av2(:,2) = Tv2_nonan';
    mv2 = Av2\(v2_nonan');
    v2_nonandetrended = v2_nonan - (Av2*mv2)'; % not kept around because we eventually
                                           % detrend the padded segments
    v2_nandetrended = v2; v2_nandetrended(isfinite(v2)) = v2_nonandetrended;
    v2_var = nanvar(v2_nandetrended);
    nancoords_v2 = find(~isfinite(v2));
    for i=1:length(nancoords_v2)
        v2(nancoords_v2(i)) = sqrt(v2_var)*randn + mv2(1) + mv2(2)*T(nancoords_v2(i));
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
seg_u2 = zeros(N, 2*SEGMENTS - 1); % Initialize
seg_v2 = zeros(N, 2*SEGMENTS - 1);
for i = 1:2:(2*SEGMENTS - 1) % Odd numbered columns
    seg_u2(:,i) = u2( (1 + ((i-1)/2)*N):(((i+1)/2)*N) );
    seg_v2(:,i) = v2( (1 + ((i-1)/2)*N):(((i+1)/2)*N) );
end
Offset = floor(0.5*N); % Midpoint at which the even-numbered columns
                       % will start; trivial for even leng_TS, and for
                       % odd leng_TS, it starts in the middle.
for i = 1:2:(2*SEGMENTS - 2) % Even numbered columns
    seg_u2(:,i+1) = u2( (Offset + 1 + ((i-1)/2)*N):(Offset + ((i+1)/2)*N) );
    seg_v2(:,i+1) = v2( (Offset + 1 + ((i-1)/2)*N):(Offset + ((i+1)/2)*N) );
end
% ^ This works for odd and even time series and segments

seg_u2 = detrend(seg_u2); % Column-wise by default
seg_v2 = detrend(seg_v2);
% From Sarah's lecture 9:
% "First you must demean your data, otherwise, the window will shift energy
% from the mean into other frequencies. If you're working in segments, you
% should demean (and detrend) each segment before you do anything further."
% Element-wise mulyiply:
window_seg_u2 = Window.*seg_u2;
window_seg_v2 = Window.*seg_v2;

%%%%%%%%%%% fft:
fft_window_seg_u2 = fft(window_seg_u2);
fft_window_seg_v2 = fft(window_seg_v2);

%%%%%%%%%%% Isolate imaginary and real components:
if mod(N,2) == 0 % even N
    A2 = real(fft_window_seg_u2(1:(N/2+1),:)); % for even N
    B2 = imag(fft_window_seg_u2(1:(N/2+1),:));
    C2 = real(fft_window_seg_v2(1:(N/2+1),:));
    D2 = imag(fft_window_seg_v2(1:(N/2+1),:));
elseif mod(N,2) == 1 % odd N
    A2 = real(fft_window_seg_u2(1:(N+1)/2,:)); % for even N
    B2 = imag(fft_window_seg_u2(1:(N+1)/2,:));
    C2 = real(fft_window_seg_v2(1:(N+1)/2,:));
    D2 = imag(fft_window_seg_v2(1:(N+1)/2,:));
else
    
end

%%%%%%%%%%% Combine into the CW and CCW parts
fft_window_seg_CW2 = 0.5*(A2 + D2 + 1i*(C2 - B2));
fft_window_seg_CCW2  = 0.5*(A2 - D2 + 1i*(C2 + B2));
% ^ backwards from the sin-cos formulation

%% CW

%%%%%%%%%%% mod squared:
if mod(N,2) == 0 % even N
    ampCW2 = abs(fft_window_seg_CW2(1:(N/2+1),:)).^2; % for even N
elseif mod(N,2) == 1 % odd N
    ampCW2 = abs(fft_window_seg_CW2(1:((N+1)/2),:)).^2; % for odd N
else
    
end
%%%%%%%%%%% add back lost energy
if mod(N,2) == 0 % even N
    ampCW2(2:end-1,:) = 2*ampCW2(2:(end-1),:); % for even N
elseif mod(N,2) == 1 % odd N
    ampCW2(2:end,:) = 2*ampCW2(2:end,:); % for odd N
else
    
end
ampCW2 = norm_factor*ampCW2*DT/N;
% Divide by N because of the normalization which MATLAB uses in "fft". The
% DT ensures that any unit for time will work for normalizing.

% The spectrum:
S_cw2 = mean(ampCW2(2:end,:),2);

%% CCW

%%%%%%%%%%% mod squared:
if mod(N,2) == 0 % even N
    ampCCW2 = abs(fft_window_seg_CCW2(1:(N/2+1),:)).^2; % for even N
elseif mod(N,2) == 1 % odd N
    ampCCW2 = abs(fft_window_seg_CCW2(1:((N+1)/2),:)).^2; % for odd N
else
    
end
%%%%%%%%%%% add back lost energy
if mod(N,2) == 0 % even N
    ampCCW2(2:end-1,:) = 2*ampCCW2(2:(end-1),:); % for even N
elseif mod(N,2) == 1 % odd N
    ampCCW2(2:end,:) = 2*ampCCW2(2:end,:); % for odd N
else
    
end
ampCCW2 = norm_factor*ampCCW2*DT/N;
% Divide by N because of the normalization which MATLAB uses in "fft". The
% DT ensures that any unit for time will work for normalizing.

% The spectrum:
S_ccw2 = mean(ampCCW2(2:end,:),2);
%% W2 END

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
    
end

%% Error

% In short, do this calculation for 95% confidence ratio:
EffectiveSegments = 2*SEGMENTS - 1;% degrees of freedom; for Hanning windowed overlapping
                     % windows, this is the total number of windows, i.e.
                     % 2*WINDOWS - 1
                     % This only works because of the Hanning filter

err_high = 2*EffectiveSegments/chi2inv(Alpha/2, 2*EffectiveSegments);
err_low = 2*EffectiveSegments/chi2inv(1 - Alpha/2, 2*EffectiveSegments);
err = [err_low err_high];

Beta = 1 - Alpha^(1/(EffectiveSegments-1));

%% Linear Spectra and Coherence



%%%%%%%%%%% mod squared:
if mod(N,2) == 0 % even N
    S_u1 = fft_window_seg_u1(1:(N/2+1),:).*conj(fft_window_seg_u1(1:(N/2+1),:)); % for even N
    S_v1 = fft_window_seg_v1(1:(N/2+1),:).*conj(fft_window_seg_v1(1:(N/2+1),:));
    S_u2 = fft_window_seg_u2(1:(N/2+1),:).*conj(fft_window_seg_u2(1:(N/2+1),:));
    S_v2 = fft_window_seg_v2(1:(N/2+1),:).*conj(fft_window_seg_v2(1:(N/2+1),:));
    
    S_u1u2 = fft_window_seg_u1(1:(N/2+1),:).*conj(fft_window_seg_u2(1:(N/2+1),:));
    S_v1v2 = fft_window_seg_v1(1:(N/2+1),:).*conj(fft_window_seg_v2(1:(N/2+1),:));
    S_u1v2 = fft_window_seg_u1(1:(N/2+1),:).*conj(fft_window_seg_v2(1:(N/2+1),:));
    S_v1u2 = fft_window_seg_v1(1:(N/2+1),:).*conj(fft_window_seg_u2(1:(N/2+1),:));
elseif mod(N,2) == 1 % odd N
    S_u1 = fft_window_seg_u1(1:((N+1)/2),:).*conj(fft_window_seg_u1(1:((N+1)/2),:)); % for odd N
    S_v1 = fft_window_seg_v1(1:((N+1)/2),:).*conj(fft_window_seg_v1(1:((N+1)/2),:));
    S_u2 = fft_window_seg_u2(1:((N+1)/2),:).*conj(fft_window_seg_u2(1:((N+1)/2),:));
    S_v2 = fft_window_seg_v2(1:((N+1)/2),:).*conj(fft_window_seg_v2(1:((N+1)/2),:));
    
    S_u1u2 = fft_window_seg_u1(1:((N+1)/2),:).*conj(fft_window_seg_u2(1:((N+1)/2),:));
    S_v1v2 = fft_window_seg_v1(1:((N+1)/2),:).*conj(fft_window_seg_v2(1:((N+1)/2),:));
    S_u1v2 = fft_window_seg_u1(1:((N+1)/2),:).*conj(fft_window_seg_v2(1:((N+1)/2),:));
    S_v1u2 = fft_window_seg_v1(1:((N+1)/2),:).*conj(fft_window_seg_u2(1:((N+1)/2),:));
else
end
%%%%%%%%%%% add back lost energy
if mod(N,2) == 0 % even N
    ampCW2(2:end-1,:) = 2*ampCW2(2:(end-1),:); % for even N
    
    S_u1(2:end-1,:) = 2*S_u1(2:(end-1),:); % for even N
    S_v1(2:end-1,:) = 2*S_v1(2:(end-1),:);
    S_u2(2:end-1,:) = 2*S_u2(2:(end-1),:);
    S_v2(2:end-1,:) = 2*S_v2(2:(end-1),:);
    
    S_u1u2(2:end-1,:) = 2*S_u1u2(2:(end-1),:);
    S_v1v2(2:end-1,:) = 2*S_v1v2(2:(end-1),:);
    S_u1v2(2:end-1,:) = 2*S_u1v2(2:(end-1),:);
    S_v1u2(2:end-1,:) = 2*S_v1u2(2:(end-1),:);
    % ^ this *2 may be unnecessary for cross spectra, because of the
    % asymmetry of the complex fft*conj(fft), but this is necessary for the
    % ratios to be correct for coherence and phase
elseif mod(N,2) == 1 % odd N
    ampCW2(2:end,:) = 2*ampCW2(2:end,:); % for odd N
    
    S_u1(2:end,:) = 2*S_u1(2:end,:); % for odd N
    S_v1(2:end,:) = 2*S_v1(2:end,:);
    S_u2(2:end,:) = 2*S_u2(2:end,:);
    S_v2(2:end,:) = 2*S_v2(2:end,:);
    
    S_u1u2(2:end,:) = 2*S_u1u2(2:end,:);
    S_v1v2(2:end,:) = 2*S_v1v2(2:end,:);
    S_u1v2(2:end,:) = 2*S_u1v2(2:end,:);
    S_v1u2(2:end,:) = 2*S_v1u2(2:end,:);
    % ^ this *2 may be unnecessary for cross spectra, because of the
    % asymmetry of the complex fft*conj(fft), but this is necessary for the
    % ratios to be correct for coherence and phase
else
end

% Coherence and phase:
Coh_u1u2 = [real(mean( S_u1u2 , 2)).^2 + imag(mean( S_u1u2 , 2)).^2]./[mean( S_u1 , 2).*mean( S_u2 , 2)]; Coh_u1u2 = Coh_u1u2(2:end);
Coh_v1v2 = [real(mean( S_v1v2 , 2)).^2 + imag(mean( S_v1v2 , 2)).^2]./[mean( S_v1 , 2).*mean( S_v2 , 2)]; Coh_v1v2 = Coh_v1v2(2:end);
Coh_u1v2 = [real(mean( S_u1v2 , 2)).^2 + imag(mean( S_u1v2 , 2)).^2]./[mean( S_u1 , 2).*mean( S_v2 , 2)]; Coh_u1v2 = Coh_u1v2(2:end);
Coh_v1u2 = [real(mean( S_v1u2 , 2)).^2 + imag(mean( S_v1u2 , 2)).^2]./[mean( S_v1 , 2).*mean( S_u2 , 2)]; Coh_v1u2 = Coh_v1u2(2:end);
% 
Pha_u1u2 = atan2( -imag(mean( S_u1u2 , 2)), real(mean( S_u1u2 , 2))); Pha_u1u2 = Pha_u1u2(2:end);
Pha_v1v2 = atan2( -imag(mean( S_v1v2 , 2)), real(mean( S_v1v2 , 2))); Pha_v1v2 = Pha_v1v2(2:end);
Pha_u1v2 = atan2( -imag(mean( S_u1v2 , 2)), real(mean( S_u1v2 , 2))); Pha_u1v2 = Pha_u1v2(2:end);
Pha_v1u2 = atan2( -imag(mean( S_v1u2 , 2)), real(mean( S_v1u2 , 2))); Pha_v1u2 = Pha_v1u2(2:end);

% Collapse all the spectra into vectors now that their matrix forms are no
% longer needed:
S_u1 = mean(S_u1(2:end,:),2);
S_v1 = mean(S_v1(2:end,:),2);
S_u2 = mean(S_u2(2:end,:),2);
S_v2 = mean(S_v2(2:end,:),2);

S_u_v1 = [S_u1, S_v1];
S_u_v2 = [S_u2, S_v2];

Coherence_linear = [Coh_u1u2, Coh_v1v2, Coh_u1v2, Coh_v1u2];
Phase_linear =     [Pha_u1u2, Pha_v1v2, Pha_u1v2, Pha_v1u2];

%% Rotary Spectra (already defined) and Coherence

%%%%%%%%%%% mod squared:
if mod(N,2) == 0 % even N
    S_CW1 = fft_window_seg_CW1(1:(N/2+1),:).*conj(fft_window_seg_CW1(1:(N/2+1),:)); % for even N
    S_CCW1 = fft_window_seg_CCW1(1:(N/2+1),:).*conj(fft_window_seg_CCW1(1:(N/2+1),:));
    S_CW2 = fft_window_seg_CW2(1:(N/2+1),:).*conj(fft_window_seg_CW2(1:(N/2+1),:));
    S_CCW2 = fft_window_seg_CCW2(1:(N/2+1),:).*conj(fft_window_seg_CCW2(1:(N/2+1),:));
    
    S_CW1CW2 = fft_window_seg_CW1(1:(N/2+1),:).*conj(fft_window_seg_CW2(1:(N/2+1),:));
    S_CCW1CCW2 = fft_window_seg_CCW1(1:(N/2+1),:).*conj(fft_window_seg_CCW2(1:(N/2+1),:));
    S_CW1CCW2 = fft_window_seg_CW1(1:(N/2+1),:).*conj(fft_window_seg_CCW2(1:(N/2+1),:));
    S_CCW1CW2 = fft_window_seg_CCW1(1:(N/2+1),:).*conj(fft_window_seg_CW2(1:(N/2+1),:));
elseif mod(N,2) == 1 % odd N
    S_CW1 = fft_window_seg_CW1(1:((N+1)/2),:).*conj(fft_window_seg_CW1(1:((N+1)/2),:)); % for odd N
    S_CCW1 = fft_window_seg_CCW1(1:((N+1)/2),:).*conj(fft_window_seg_CCW1(1:((N+1)/2),:));
    S_CW2 = fft_window_seg_CW2(1:((N+1)/2),:).*conj(fft_window_seg_CW2(1:((N+1)/2),:));
    S_CCW2 = fft_window_seg_CCW2(1:((N+1)/2),:).*conj(fft_window_seg_CCW2(1:((N+1)/2),:));
    
    S_CW1CW2 = fft_window_seg_CW1(1:((N+1)/2),:).*conj(fft_window_seg_CW2(1:((N+1)/2),:));
    S_CCW1CCW2 = fft_window_seg_CCW1(1:((N+1)/2),:).*conj(fft_window_seg_CCW2(1:((N+1)/2),:));
    S_CW1CCW2 = fft_window_seg_CW1(1:((N+1)/2),:).*conj(fft_window_seg_CCW2(1:((N+1)/2),:));
    S_CCW1CW2 = fft_window_seg_CCW1(1:((N+1)/2),:).*conj(fft_window_seg_CW2(1:((N+1)/2),:));
else
end
%%%%%%%%%%% add back lost energy
if mod(N,2) == 0 % even N
    ampCW2(2:end-1,:) = 2*ampCW2(2:(end-1),:); % for even N
    
    S_CW1(2:end-1,:) = 2*S_CW1(2:(end-1),:); % for even N
    S_CCW1(2:end-1,:) = 2*S_CCW1(2:(end-1),:);
    S_CW2(2:end-1,:) = 2*S_CW2(2:(end-1),:);
    S_CCW2(2:end-1,:) = 2*S_CCW2(2:(end-1),:);
    
    S_CW1CW2(2:end-1,:) = 2*S_CW1CW2(2:(end-1),:);
    S_CCW1CCW2(2:end-1,:) = 2*S_CCW1CCW2(2:(end-1),:);
    S_CW1CCW2(2:end-1,:) = 2*S_CW1CCW2(2:(end-1),:);
    S_CCW1CW2(2:end-1,:) = 2*S_CCW1CW2(2:(end-1),:);
    % ^ this *2 may be unnecessary for cross spectra, because of the
    % asymmetry of the complex fft*conj(fft), but this is necessary for the
    % ratios to be correct for coherence and phase
elseif mod(N,2) == 1 % odd N
    ampCW2(2:end,:) = 2*ampCW2(2:end,:); % for odd N
    
    S_CW1(2:end,:) = 2*S_CW1(2:end,:); % for odd N
    S_CCW1(2:end,:) = 2*S_CCW1(2:end,:);
    S_CW2(2:end,:) = 2*S_CW2(2:end,:);
    S_CCW2(2:end,:) = 2*S_CCW2(2:end,:);
    
    S_CW1CW2(2:end,:) = 2*S_CW1CW2(2:end,:);
    S_CCW1CCW2(2:end,:) = 2*S_CCW1CCW2(2:end,:);
    S_CW1CCW2(2:end,:) = 2*S_CW1CCW2(2:end,:);
    S_CCW1CW2(2:end,:) = 2*S_CCW1CW2(2:end,:);
    % ^ this *2 may be unnecessary for cross spectra, because of the
    % asymmetry of the complex fft*conj(fft), but this is necessary for the
    % ratios to be correct for coherence and phase
else
end

% Coherence and phase:
Coh_CW1CW2 = [real(mean( S_CW1CW2 , 2)).^2 + imag(mean( S_CW1CW2 , 2)).^2]./[mean( S_CW1 , 2).*mean( S_CW2 , 2)];           Coh_CW1CW2 = Coh_CW1CW2(2:end);
Coh_CCW1CCW2 = [real(mean( S_CCW1CCW2 , 2)).^2 + imag(mean( S_CCW1CCW2 , 2)).^2]./[mean( S_CCW1 , 2).*mean( S_CCW2 , 2)];   Coh_CCW1CCW2 = Coh_CCW1CCW2(2:end);
Coh_CW1CCW2 = [real(mean( S_CW1CCW2 , 2)).^2 + imag(mean( S_CW1CCW2 , 2)).^2]./[mean( S_CW1 , 2).*mean( S_CCW2 , 2)];       Coh_CW1CCW2 = Coh_CW1CCW2(2:end);
Coh_CCW1CW2 = [real(mean( S_CCW1CW2 , 2)).^2 + imag(mean( S_CCW1CW2 , 2)).^2]./[mean( S_CCW1 , 2).*mean( S_CW2 , 2)];       Coh_CCW1CW2 = Coh_CCW1CW2(2:end);
% Rotary phase is multiplied by -1, as it seems that this ensures negative
% lag corresponds to w1 leading w2 (which is the result for linear phase
% above for u's and v's)
Pha_CW1CW2 =   -atan2( -imag(mean( S_CW1CW2 , 2)),   real(mean( S_CW1CW2 , 2)));   Pha_CW1CW2 = Pha_CW1CW2(2:end);
Pha_CCW1CCW2 = -atan2( -imag(mean( S_CCW1CCW2 , 2)), real(mean( S_CCW1CCW2 , 2))); Pha_CCW1CCW2 = Pha_CCW1CCW2(2:end);
Pha_CW1CCW2 =  -atan2( -imag(mean( S_CW1CCW2 , 2)),  real(mean( S_CW1CCW2 , 2)));  Pha_CW1CCW2 = Pha_CW1CCW2(2:end);
Pha_CCW1CW2 =  -atan2( -imag(mean( S_CCW1CW2 , 2)),  real(mean( S_CCW1CW2 , 2)));  Pha_CCW1CW2 = Pha_CCW1CW2(2:end);

% Collapse all the spectra into vectors now that their matrix forms are no
% longer needed:
S_CW1 =  mean(S_CW1(2:end,:),2);
S_CCW1 = mean(S_CCW1(2:end,:),2);
S_CW2 =  mean(S_CW2(2:end,:),2);
S_CCW2 = mean(S_CCW2(2:end,:),2);

% Spectra:
S_cw_ccw1 = [S_CW1, S_CCW1];
S_cw_ccw2 = [S_CW2, S_CCW2];

Coherence_rotary = [Coh_CW1CW2, Coh_CCW1CCW2, Coh_CW1CCW2, Coh_CCW1CW2];
Phase_rotary =     [Pha_CW1CW2, Pha_CCW1CCW2, Pha_CW1CCW2, Pha_CCW1CW2];

%% Plot
if PLOT_BOOLEAN
    
    LEG_LOC = 'northoutside';
    
    % LINEAR SPECTRA
    figure('Color',[1 1 1])
    subplot(2,1,1)
    loglog(f_vec,S_u1,PLOT_OPTION{1}); hold on
    loglog(f_vec,S_v1,PLOT_OPTION{2})
    xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
    ylabel(['Spectral density [(time\_series\_units)$^2$ (cycles per ',...
        TIME_UNITS,')$^{-1}$]'],'Interpreter','LaTeX')
    % error bar plotted
    loglog([1.1 1.1]*f_vec(length(S_u1)),[1 (err_high/err_low)]*min([min(S_u1) min(S_v1)]),'k');
    legend({'u','v','95% confidence'},'Location',LEG_LOC,'Orientation','horizontal'); title('w_1 Spectra: u & v')
    subplot(2,1,2) % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    loglog(f_vec,S_u2,PLOT_OPTION{1}); hold on
    loglog(f_vec,S_v2,PLOT_OPTION{2})
    xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
    ylabel(['Spectral density [(time\_series\_units)$^2$ (cycles per ',...
        TIME_UNITS,')$^{-1}$]'],'Interpreter','LaTeX')
    % error bar plotted
    loglog([1.1 1.1]*f_vec(length(S_u2)),[1 (err_high/err_low)]*min([min(S_u2) min(S_v2)]),'k');
    legend({'u','v','95% confidence'},'Location',LEG_LOC,'Orientation','horizontal'); title('w_2 Spectra: u & v')

    % ROTARY SPECTRA
    if ~isreal(w1) || ~isreal(w2) % i.e. only plot this if either input time series is complex
        figure('Color',[1 1 1])
        subplot(2,1,1)
        loglog(f_vec,S_CW1,PLOT_OPTION{1}); hold on
        loglog(f_vec,S_CCW1,PLOT_OPTION{2})
        xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
        ylabel(['Spectral density [(time\_series\_units)$^2$ (cycles per ',...
            TIME_UNITS,')$^{-1}$]'],'Interpreter','LaTeX')
        % error bar plotted
        loglog([1.1 1.1]*f_vec(length(S_CW1)),[1 (err_high/err_low)]*min([min(S_CW1) min(S_CCW1)]),'k');
        legend({'CW','CCW','95% confidence'},'Location',LEG_LOC,'Orientation','horizontal'); title('w_1 Rotary Spectra: cw & ccw')
        subplot(2,1,2) % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        loglog(f_vec,S_CW2,PLOT_OPTION{1}); hold on
        loglog(f_vec,S_CCW2,PLOT_OPTION{2})
        xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
        ylabel(['Spectral density [(time\_series\_units)$^2$ (cycles per ',...
            TIME_UNITS,')$^{-1}$]'],'Interpreter','LaTeX')
        % error bar plotted
        loglog([1.1 1.1]*f_vec(length(S_CW2)),[1 (err_high/err_low)]*min([min(S_CW2) min(S_CCW2)]),'k');
        legend({'CW','CCW','95% confidence'},'Location',LEG_LOC,'Orientation','horizontal'); title('w_2 Rotary Spectra: cw & ccw')
    end



    % COHERENCE AND PHASE

    % LINEAR, PARALLEL
    figure('Color',[1 1 1])
    subplot(2,1,1)
    semilogx(f_vec,Coh_u1u2,PLOT_OPTION{1}); hold on
    semilogx(f_vec,Coh_v1v2,PLOT_OPTION{2})
    semilogx([f_vec(1), f_vec(end)],Beta*[1 1],'-k')
    set(gca,'ylim',[0 1])
    xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
    ylabel('Coherence [unitless]','Interpreter','LaTeX')
    legend({'\gamma^2(u1,u2) (f)','\gamma^2(v1,v2) (f)',[num2str(100 - 100*Alpha),'% confidence threshold']},'Location',LEG_LOC,'Orientation','horizontal')
    subplot(2,1,2)
    semilogx(f_vec,Pha_u1u2,PLOT_OPTION{1}); hold on
    semilogx(f_vec,Pha_v1v2,PLOT_OPTION{2})
    set(gca,'ylim',[-pi pi])
    xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
    ylabel('Phase difference [radians]','Interpreter','LaTeX')
    legend({'\phi(u1,u2) (f)','\phi(v1,v2) (f)'},'Location',LEG_LOC,'Orientation','horizontal')

    % LINEAR, PERPENDICULAR
    if ~isreal(w1) || ~isreal(w2) % i.e. only plot this if either input time series is complex
        figure('Color',[1 1 1])
        subplot(2,1,1)
        semilogx(f_vec,Coh_u1v2,PLOT_OPTION{1}); hold on
        semilogx(f_vec,Coh_v1u2,PLOT_OPTION{2})
        semilogx([f_vec(1), f_vec(end)],Beta*[1 1],'-k')
        set(gca,'ylim',[0 1])
        xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
        ylabel('Coherence [unitless]','Interpreter','LaTeX')
        legend({'\gamma^2(u1,v2) (f)','\gamma^2(v1,u2) (f)',[num2str(100 - 100*Alpha),'% confidence threshold']},'Location',LEG_LOC,'Orientation','horizontal')
        subplot(2,1,2)
        semilogx(f_vec,Pha_u1v2,PLOT_OPTION{1}); hold on
        semilogx(f_vec,Pha_v1u2,PLOT_OPTION{2})
        set(gca,'ylim',[-pi pi])
        xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
        ylabel('Phase difference [radians]','Interpreter','LaTeX')
        legend({'\phi(u1,v2)(f)','\phi(v1,u2)(f)'},'Location',LEG_LOC,'Orientation','horizontal')

        % ROTARY, PARALLEL
        figure('Color',[1 1 1])
        subplot(2,1,1)
        semilogx(f_vec,Coh_CW1CW2,PLOT_OPTION{1}); hold on
        semilogx(f_vec,Coh_CCW1CCW2,PLOT_OPTION{2})
        semilogx([f_vec(1), f_vec(end)],Beta*[1 1],'-k')
        set(gca,'ylim',[0 1])
        xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
        ylabel('Coherence [unitless]','Interpreter','LaTeX')
        legend({'\gamma^2(CW1,CW2) (f)','\gamma^2(CCW1,CCW2) (f)',[num2str(100 - 100*Alpha),'% confidence threshold']},'Location',LEG_LOC,'Orientation','horizontal')
        subplot(2,1,2)
        semilogx(f_vec,Pha_CW1CW2,PLOT_OPTION{1}); hold on
        semilogx(f_vec,Pha_CCW1CCW2,PLOT_OPTION{2})
        set(gca,'ylim',[-pi pi])
        xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
        ylabel('Phase difference [radians]','Interpreter','LaTeX')
        legend({'\phi(CW1,CW2) (f)','\phi(CCW1,CCW2) (f)'},'Location',LEG_LOC,'Orientation','horizontal')

        % ROTARY, PERPENDICULAR
        figure('Color',[1 1 1])
        subplot(2,1,1)
        semilogx(f_vec,Coh_CW1CCW2,PLOT_OPTION{1}); hold on
        semilogx(f_vec,Coh_CCW1CW2,PLOT_OPTION{2})
        semilogx([f_vec(1), f_vec(end)],Beta*[1 1],'-k')
        set(gca,'ylim',[0 1])
        xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
        ylabel('Coherence [unitless]','Interpreter','LaTeX')
        legend({'\gamma^2(CW1,CCW2) (f)','\gamma^2(CCW1,CW2) (f)',[num2str(100 - 100*Alpha),'% confidence threshold']},'Location',LEG_LOC,'Orientation','horizontal')
        subplot(2,1,2)
        semilogx(f_vec,Pha_CW1CCW2,PLOT_OPTION{1}); hold on
        semilogx(f_vec,Pha_CCW1CW2,PLOT_OPTION{2})
        set(gca,'ylim',[-pi pi])
        xlabel(['Cycles per ',TIME_UNITS],'Interpreter','LaTeX')
        ylabel('Phase difference [radians]','Interpreter','LaTeX')
        legend({'\phi(CW1,CCW2)(f)','\phi(CCW1,CW2)(f)'},'Location',LEG_LOC,'Orientation','horizontal')
    end
else
end


%% Troubleshooting (part 2)
% % Good way to check that your time series are behaving properly:
% figure
% subplot(211);plot(raw_TS);title('raw')
% subplot(212);plot(TS);title('post-processing')

end
end
%% Inertial currents (CW in NH, CCW in SH)

% see:
% Chapter 7: Dynamical Processes for Descriptive Ocean Circulation - 2011 Descriptive Physical Oceanography Sixth Edition
% Lynne Talley
