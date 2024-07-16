% [SPEC, F12, err] = NANSPECTRUM2(DATA, [D1  D2], [SEG1 SEG2], {UNITS1, UNITS2}, PLOT_BOOLEAN, WINDOWMETHOD)
% 
% Two-dimensional power spectrum estimation. Arbitrary dimensions for the
% two inputs (e.g. space-time or space-space).
% 
% -------------------------------------------------------------------------
% REQUIRED INPUTS:
%   DATA = Data matrix, size nd by md. May be complex for corrent = u + i*v
%   DD   = 2-element vector where:
%          DD(1) = step length in the horizontal dimension of DATA
%          DD(2) = step length in the vertical dimension of DATA
%   SEGS = 2-element vector where:
%          SEGS(1) = number of segments in the horizontal dimension of DATA
%          SEGS(2) = number of segments in the vertical dimension of DATA
% 
% OPTIONAL INPUTS:
%   UNITS = 2-element cell where:
%          UNITS{1} = physical units of the horizontal dimension of DATA
%          UNITS{2} = physical units of the vertical dimension of DATA
%   PLOT_BOOLEAN = true or false, decides whether or not to plot (default false)
%   WINDOWMETHOD = string, name of the windowing function applied to
%          segments in both dimensions (default 'hanning')
% 
% OPTIONAL INPUT NAME-VALUE PAIRS:
%   'Detrend'   default true, decides if each segment should have a planar
%               fit removed (true) or not (false)
%   'Overlap'   default true, decides if overlapping segments should be
%               used
%   'Colormap'  default 'parula', the colormap used for plots (only used if
%               PLOT_BOOLEAN == true)
% -------------------------------------------------------------------------
% OUTPUTS:
%   SPEC = Power spectrum, size nd/SEGS(2) by md/SEGS(1), units of
%          [data_units]^2 * [UNITS{1}] * [UNITS{2}]
%   F12  = 3-dimensional array, size nd/SEGS(2) by md/SEGS(1) by 2, where
%          F12(:,:,1) = frequencies corresponding to the horizontal
%                       dimension of DATA
%          F12(:,:,2) = frequencies corresponding to the vertical
%                       dimension of DATA
%   err  = (OPTIONAL) error ratio, where err(1) = 1 and err(2) > 1.
%          Calculated using chi2inv (see function for more information).
% -------------------------------------------------------------------------

function [SPEC, F12, varargout] = nanspectrum2(DATA, DD, SEGMENTS, varargin)
%% Variables defined
Detrend = true;
Overlap = true;
Colormap = 'parula';
if nargin == 3
    UNITS        = {'',''};
    PLOT_BOOLEAN = false;
    WINDOWMETHOD = 'hanning';
elseif nargin == 4
    UNITS        = varargin{1};
    PLOT_BOOLEAN = false;
    WINDOWMETHOD = 'hanning';
elseif nargin == 5
    UNITS        = varargin{1};
    PLOT_BOOLEAN = varargin{2};
    WINDOWMETHOD = 'hanning';
elseif nargin == 6
    UNITS        = varargin{1};
    PLOT_BOOLEAN = varargin{2};
    WINDOWMETHOD = varargin{3};
elseif nargin > 6
    % Other inputs given as name-value pairs:
    UNITS        = varargin{1};
    PLOT_BOOLEAN = varargin{2};
    WINDOWMETHOD = varargin{3};
    OTHER_ARGS_  = varargin;
    OTHER_ARGS = cell(1,nargin - 6);
    for ii=4:length(OTHER_ARGS_)
        OTHER_ARGS{ii-3} = OTHER_ARGS_{ii};
    end
    Str = struct(OTHER_ARGS{:});
    Names = fieldnames(Str);
    % Verify that variables are only the ones that are allowed:
    AllowedVars = {'Detrend','Overlap','Colormap'}';
    for i=1:length(Names)
        if ismember(Names{i},AllowedVars)
        else
            error(['''',Names{i},''' is not a possible option for nanspectrum2.' char(13) ...
                   'Valid options are: ''Detrend'', ''Overlap'', and ''Colormap''' char(13) ...
                   'Please see documentation.'])
        end
    end
    for i=1:length(Names)
        eval([Names{i},' = Str.',Names{i},';'])
    end
else
    error('Incorrect number of inputs')
end

N1 = size(DATA,1); % Y's
N2 = size(DATA,2); % X's

% Confusingly, this flips the expected order of the inputs:
D1 = DD(2); % delta in the Y direction
D2 = DD(1); % delta in the X direction
SEG1 = SEGMENTS(2); % segments in the Y direction
SEG2 = SEGMENTS(1); % segments in the X direction

error_Msg = ['Make sure that the data can actually be divided into the given number of windows.' char(13) ...
        'As it is given, the data has ' num2str(N1) 'data points in the vertical direction and ' char(13) ...
        num2str(N2) ' in the horizontal direction. These can be factored as follows: ' char(13) ...
        'VERTICAL   [' num2str(factor(N1)) ']' char(13) ...
        'HORIZONTAL [' num2str(factor(N2)) ']'];

if mod(N1,SEG1) || mod(N2,SEG2) % either or both fail
    warning(error_Msg)
    error('Incompatible number of segments given.')
else
end

%% Segment

% % % Non-overlapping part:
if Overlap
    % for both non-overlapping and overlapping segments
    DATA_reshape = nan([ N1/SEG1, N2/SEG2, SEG1*SEG2 + [SEG1-1]*[SEG2-1] ]);
else
    % for only non-overlapping segments
    DATA_reshape = nan([ N1/SEG1, N2/SEG2, SEG1*SEG2 ]);
end
Iys = reshape(mod(0:[[SEG1*SEG2]-1],SEG1)+1,[SEG1,SEG2]); % Index y segments
Ixs = reshape(mod(0:[[SEG1*SEG2]-1],SEG2)+1,[SEG2,SEG1])'; % Index x segments
Iys = Iys(:); Ixs = Ixs(:);
for ii = 1:[SEG1*SEG2]
    DATA_reshape(:,:,ii) = DATA([1:(N1/SEG1)] + (Iys(ii)-1)*N1/SEG1 , ...
                                [1:(N2/SEG2)] + (Ixs(ii)-1)*N2/SEG2);
    % % % For index verification:
    % [min([1:(N1/SEG1)] + (Iys(ii)-1)*N1/SEG1), max([1:(N1/SEG1)] + (Iys(ii)-1)*N1/SEG1);...
    %  min([1:(N2/SEG2)] + (Ixs(ii)-1)*N2/SEG2), max([1:(N2/SEG2)] + (Ixs(ii)-1)*N2/SEG2)]
end
% % % Overlapping part:
if Overlap
    Iys = reshape(mod(0:[[[SEG1-1]*[SEG2-1]]-1],SEG1-1)+1,[SEG1-1,SEG2-1]); % Index y segments
    Ixs = reshape(mod(0:[[[SEG1-1]*[SEG2-1]]-1],SEG2-1)+1,[SEG2-1,SEG1-1])'; % Index x segments
    Iys = Iys(:); Ixs = Ixs(:);
    for ii = 1:[[SEG1-1]*[SEG2-1]]
        DATA_reshape(:,:,ii+[SEG1*SEG2]) = DATA([1:(N1/SEG1)] + (Iys(ii)-1)*N1/SEG1 + floor(0.5*[N1/SEG1]+1)-1, ...
                                                [1:(N2/SEG2)] + (Ixs(ii)-1)*N2/SEG2 + floor(0.5*[N2/SEG2]+1)-1);
        % % % For index verification:
        % [min([1:(N1/SEG1)] + (Iys(ii)-1)*N1/SEG1 + floor(0.5*[N1/SEG1]+1)-1), max([1:(N1/SEG1)] + (Iys(ii)-1)*N1/SEG1 + floor(0.5*[N1/SEG1]+1)-1);...
        %  min([1:(N2/SEG2)] + (Ixs(ii)-1)*N2/SEG2 + floor(0.5*[N2/SEG2]+1)-1), max([1:(N2/SEG2)] + (Ixs(ii)-1)*N2/SEG2 + floor(0.5*[N2/SEG2]+1)-1)]
    end
else
end

%% Segment-by-segment interpolate over NaN's using a plane fit (real and imaginary separately)

[XX,YY] = meshgrid(1:[N2/SEG2],1:[N1/SEG1]);
% data = H*coef + residual ---> H\data = coef_estimated
HH = [XX(:), YY(:), ones(size(XX(:)))];

for seg_i = 1:size(DATA_reshape,3)
    % Fill gaps if there are any:
    DATA_slice = DATA_reshape(:,:,seg_i);
    IsFiniteData = isfinite(DATA_slice); IsFiniteData = IsFiniteData(:);
    if sum(~IsFiniteData) || Detrend
        DATA_col = DATA_slice; DATA_col = DATA_col(:);
        Coef_R = HH(IsFiniteData,:)\real(DATA_col(IsFiniteData));
        Coef_I = HH(IsFiniteData,:)\imag(DATA_col(IsFiniteData));
        DATA_plane = reshape(HH*Coef_R, N1/SEG1, N2/SEG2) + 1i*reshape(HH*Coef_I, N1/SEG1, N2/SEG2);
        DATA_slice(~isfinite(DATA_slice)) = DATA_plane(~isfinite(DATA_slice));
        DATA_reshape(:,:,seg_i) = DATA_slice;
    else
    end
    % De-plane
    if Detrend
        DATA_reshape(:,:,seg_i) = DATA_reshape(:,:,seg_i) - DATA_plane;
    else
    end
end

%% Window

WINDOW = Window2D(WINDOWMETHOD,squeeze(DATA_reshape(:,:,1)));

for seg_i = 1:size(DATA_reshape,3)
    DATA_reshape(:,:,seg_i) = WINDOW.*DATA_reshape(:,:,seg_i);
end

% The normalization factor for the spectrum (to ensure parseval's theorem
% is fulfilled) is 1/(the mean of the square of the filter):
norm_factor = 1/mean([0;WINDOW(:)].^2);

%% FFT2 and calculate spectrum

fftDATA = nan(size(DATA_reshape));
for seg_i = 1:size(DATA_reshape,3)
    fftDATA(:,:,seg_i) = fft2(DATA_reshape(:,:,seg_i));
end

SPEC = sum(abs(fftDATA).^2, 3)/size(fftDATA,3);
SPEC = fftshift(SPEC,1);
SPEC = fftshift(SPEC,2);

%% Frequencies

F1 = [-1/2]:[1/[N1/SEG1]]:[1/2]; F1 = F1/D1;
F2 = [-1/2]:[1/[N2/SEG2]]:[1/2]; F2 = F2/D2;

% [F2,F1] = meshgrid(F2(1:[end-1]),F1(1:[end-1]));
[F2,F1] = meshgrid(F2(1:[end-1]),F1(2:end));

F1 = flip(F1);

F12 = cat(3,F2,F1);

% figure;subplot(121);imagesc(F12(:,:,1));subplot(122);imagesc(F12(:,:,2))

%% Normalize:

SPEC = norm_factor*SPEC*D1*D2/[[N1/SEG1]*[N2/SEG2]];

%% Plot (if requested)

if PLOT_BOOLEAN
    % whos F1 F2 SPEC
    figure
    pcolor_centered(F2,F1,log10(SPEC)); shading flat
    clim([prctile(log10(SPEC(:)),1) max(log10(SPEC(:)))])
    xlabel(['Cycles per ' UNITS{1}])
    ylabel(['Cycles per ' UNITS{2}])
    colormap(Colormap)
    c = colorbar;
    c.YLabel.String = 'log_{10}S';

    % % % Data variance:
    disp('Data variance:')
    disp(var(DATA(:),'omitnan'))
    % % % Sum of the spectrum/df (should be the same as variance according to Parseval's theorem):
    disp('Sum of the spectrum/df:')
    disp(sum(SPEC(:) * [1/[N1/SEG1]] * [1/[N2/SEG2]] * 1/[D1*D2]))

end

%% Error
% These will certainly be different depending on whether or not we segment
% in space, but for now I'll keep like this.

% In short, do this calculation for 95% confidence ratio:
if Overlap
    NumSegs = (SEG1*SEG2) + ((SEG1-1)*(SEG2-1)); % Total number of segments
else
    NumSegs = (SEG1*SEG2); % Total number of segments
end
% degrees of freedom; for Hanning windowed overlapping
% This only works because of the Hanning filter

err_high = 2*NumSegs/chi2inv(   0.05/2, 2*NumSegs);
err_low = 2*NumSegs/chi2inv(1 - 0.05/2, 2*NumSegs);
err = [err_low err_high]/err_low;
% Divided by "err_low" so that err(1) = 1. The quantity "err" is only
% useful as a ratio.

% disp(err)

%% Outputs

if nargout < 3
elseif nargout == 3
    warning(['Error estimate is an educated guess based on 1D spectral principles,',...
             ' inspect code and determine if this is correct and/or sufficient for your needs.'])
    varargout{1} = err;
end

%% Auxiliary function
function WIN2D = Window2D(WIN,INPUT)
W1 = eval([WIN '(size(INPUT,1));']); % vertical dim.
W2 = eval([WIN '(size(INPUT,2));']); % horiz. dim.
[WW2,WW1] = meshgrid(W2,W1);
WIN2D = WW2.*WW1;
end

function PC = pcolor_centered(X,Y,D)
dx = X(1,2) - X(1,1);
X_ = [X , X(:,end) + dx ; ...
      X(end,:) , X(end,end) + dx];
dy = Y(2,1) - Y(1,1);
Y_ = [Y , Y(:,end) ; ...
      Y(end,:) + dy , Y(end,end) + dy];
D_ = [D , nan(size(D,1),1) ; ...
      nan(1,size(D,2)) , nan ];
PC = pcolor(X_ - dx/2, Y_ - dy/2, D_);
end
end
%% Example

% close all
% 
% dT = 1;
% TT = 0:dT:999;
% T_wave = 50; Omega_wave = 2*pi/T_wave;
% 
% dX = 0.1;
% XX = [0:dX:99.9]';
% L_wave = 5; K_wave = 2*pi/L_wave;
% PHASE_SPEED = Omega_wave/K_wave;
% 
% WAVES = real( exp(1i*[K_wave*XX - Omega_wave*TT]) );
% WAVES = WAVES + randn(size(WAVES))*std(WAVES(:),'omitnan')*0.05;
% 
% % % % Animate:
% % figure
% % for ti = 1:length(TT)
% %     plot(XX,WAVES(:,ti),'.-')
% %     pause(0.05)
% % end
% 
% figure
% imagesc(WAVES); colorbar
% 
% tic
% [SPEC,FF] = nanspectrum2(WAVES, [dT dX], [5 2], {'second','meter'}, true, 'hanning', ...
%     'Colormap','parula' , 'Overlap',true , 'Detrend',false);
% toc
