% X, maxlag, 
function [spec, f_spec] = WKT_spectrum(X,maxlag)
% Make spectrum, or use the one given if one is given:
N_x = length(X);
X = X - nanmean(X);
% if ...
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
autocov_x = [autocov_x;flip(autocov_x(2:(end-1)))]; % autocov_x = [autocov_x;0;flip(autocov_x(2:end))];
Window = fftshift(hann(length(autocov_x),'periodic'));
spec = fft(autocov_x.*Window); % Window applied for smoothing
figure; subplot(211); plot(Window,'.-'); subplot(212); plot(autocov_x,'.-')%$debugging
spec = spec/length(spec); % Satisfy Parseval's Theorem
if var(real(spec)) < 10000*var(imag(spec)) % make sure that the autocovariance is periodic and thus that its fft is real
    warning(['Imaginary part [removed] of fft(autocov_x) was ',num2str(var(imag(spec))/var(real(spec)))],' times that of the real part.')
end
spec = real(spec);
if sum(spec < 0) % if the WKT spectrum is negative at any frequency
    warning(['The internally calculated spectrum (using the WKT) results',...
        ' in a partly negative spectrum.']);
    % error('See warning.') % Only if we don't want this as a possibility.
else
end
% Fold the spectrum in on itself du to symmetry:
if mod(length(spec),2) % odd length
    spec = spec(1:ceil(length(spec)/2));%%%
elseif ~mod(length(spec),2) % even length
    spec = spec(1:ceil(length(spec)/2 + 0.5));%%%
end
f_Ny = 1/2; % "T" will be evenly-spaced if the user follows the directions in the documentation
    % Eliminated because it shifted and aligned frequencies incorrectly:
    % df_cov = f_Ny/length(spec);
    % f_spec = [df_cov:df_cov:f_Ny]';
    % spec(2:end) = 2*spec(2:end);
spec(2:end) = 2*spec(2:end);
spec = spec(2:end);
df_cov = f_Ny/length(spec);
f_spec = [df_cov:df_cov:f_Ny]';

end

%% Test Case

% T = 0:999;
% X = 3*sin(2*pi*T/50) + randn(size(T)) + 10;
% 
% [WelchS, WelchF, ~] = nanspectrum_LK(X, 1, 'unit', 5, '.-', 0, 0, 'hanning');
% [wktS, wktF] = WKT_spectrum(X,length(T)/10);
% %%
% figure
% plot(wktF,wktS,'.-'); hold on
% plot(WelchF,WelchS*(5/(T(end)-T(1))),'.-')
% % set(gca,'yscale','log')
