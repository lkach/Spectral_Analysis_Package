% [Spec, f_spec] = WKT_spectrum(X,maxlag)
%
% Estimates the power spectrum using the Wiener–Khinchin theorem (Fourier
% transform of the autocovariance), as opposed to Welch's method (averaged
% periodigrams).
%
% IN:   X = Nx1 vector time series. May contain NaN's and/or Inf's, but
%           data must be evenly spaced.
% IN:   maxlags = Scalar, the maximum lag (multiple of the time step)
%           considered by the internally calculated autocovariance.
% IN:   WindowMethod = (Optional) If it's a string, the function used to
%           build the window by which the autocovariance is multiplied. The
%           default is "hann". For no adjustment to the autocovariance, use
%           "rectwin".
%           If it's a vector, it will be used as window itself, e.g.
%           fftshift(hann(2*maxlag,'periodic'));
%
% OUT:  Spec = The output power spectrum, which has units of X^2 / freq.
%           Multiply by f_spec(1) to get variance.
% OUT:  f_spec = The frequencies corresponding to the elements of "Spec".
%           This has units of 1/[units of time step], where
%           1 time step = 1 time unit (e.g. if dt = 1 hr, f_spec is in
%           1/hr; likewise, if dt = 3 seconds, f_spec is in 1/(3 sec);
%           convert accordingly)
%

function [Spec, f_spec] = WKT_spectrum(X,maxlag,varargin)

X(isinf(X)) = nan;

N_x = length(X);
X = X - nanmean(X);

autocov_x = zeros(1 + maxlag,1); % the autocovariance
% autocov_std_x = zeros(1 + maxlag,1); % the std of the points that went into each element in pre-shaped "autocov_x" (starting with zero-lag)
% n_autocov_x = zeros(1 + maxlag,2); % how many pairs were finite and how many pairs could have been finite if perfect
for i=1:(1 + maxlag) % lag = (i-1)*dt
    autocov_x_i = zeros(N_x - i + 1,1);
    for j=1:(N_x - i + 1)
        autocov_x_i(j) = X(j)*X(j+i-1);
    end
    autocov_x(i) = nanmean(autocov_x_i); % autocov
    %     autocov_std_x(i) = nanstd(autocov_x_i); % std of points that made each element of autocov
    %     n_autocov_x(i,1) = sum(isfinite(autocov_x_i)); % how many pairs were finite
    %     n_autocov_x(i,2) = N_x - i + 1; % how many pairs could have been finite if perfect
end
% NOTE: "autocov_std_x" and "n_autocov_x" are unused, and are
% retained here for troubleshooting and possible future inclusion.
autocov_x = [autocov_x;flip(autocov_x(2:(end-1)))];

if nargin == 2
    Window = fftshift(hann(length(autocov_x),'periodic'));
elseif nargin == 3
    if isstr(varargin{1})
        if strcmp(varargin{1},'blackman') || ...
                strcmp(varargin{1},'flattopwin') || ...
                strcmp(varargin{1},'hamming') || ...
                strcmp(varargin{1},'hann') || ...
                strcmp(varargin{1},'hanning') || ...
                strcmp(varargin{1},'blackmanharris')
            Window = eval(['fftshift(',varargin{1},'(length(autocov_x),''periodic''));']);
        else
            Window = eval(['fftshift(',varargin{1},'(length(autocov_x)));']);
        end
    else
        Window = varargin{1};
    end
else
    error('Incorrect number of inputs; please read documentation.')
end

Spec = fft(autocov_x.*Window); % Window applied for smoothing

% figure; subplot(211); plot(Window,'.-'); subplot(212); plot(autocov_x,'.-') % uncomment for debugging

Spec = Spec/length(Spec); % Satisfy Parseval's Theorem
if var(real(Spec)) < 10000*var(imag(Spec)) % make sure that the autocovariance is periodic and thus that its fft is real
    warning(['Imaginary part [removed] of fft(autocov_x) was ',num2str(var(imag(Spec))/var(real(Spec)))],' times that of the real part.')
end
Spec = real(Spec);
if sum(Spec < 0) % if the WKT spectrum is negative at any frequency
    warning(['The internally calculated spectrum (using the WKT) results',...
        ' in a partly negative spectrum.']);
    % error('See warning.')
    % % ^^^ Uncomment if the possibility of negative spectrum is strongly undesired.
else
end
% Fold the spectrum in on itself to give it symmetry:
if mod(length(Spec),2) % odd length
    Spec = Spec(1:ceil(length(Spec)/2));
elseif ~mod(length(Spec),2) % even length
    Spec = Spec(1:ceil(length(Spec)/2 + 0.5));
end
f_Ny = 1/2; % "T" will be evenly-spaced if the user follows the directions in the documentation
Spec(2:end) = 2*Spec(2:end);
Spec = Spec(2:end);
df_cov = f_Ny/length(Spec);
f_spec = (df_cov:df_cov:f_Ny)';
Spec = Spec/df_cov; % Make this a proper power spectrum
end

%% Example

% T = 0:999;
% X = 3*sin(2*pi*T/50) + randn(size(T)) + 10;
% 
% [WelchS, WelchF, ~] = nanspectrum_LK(X, 1, 'unit', 5, '.-', 0, 0, 'hanning');
% [wktS, wktF] = WKT_spectrum(X,length(T)/10);
% 
% figure
% plot(wktF,wktS,'.-'); hold on
% plot(WelchF,WelchS,'.-')
% set(gca,'yscale','log')
