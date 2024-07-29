%%

% Simple demonstrations of nunanspectrum

%%

TIME = [1:5000 7000:10000]';
% TIME = [1:3000]';
DATA =     cos(TIME*2*pi/10) + AR_make(0.95,length(TIME)) + ...
       1i*[-sin(TIME*2*pi/10) + AR_make(0.95,length(TIME))];
% DATA =     cos(TIME*2*pi/10) + AR_make(0.95,length(TIME));
close all
figure
plot(TIME,real(DATA)); hold on
plot(TIME,imag(DATA))
figure
[SPEC,FF] = nunanspectrum(DATA, ...
TIME, 'X', ...
'Segments',4,'Window','hanning','Plot',true,'PlotSegments',true);

%%

TIME = [1:5000 9000:10000]';
% TIME = [1:3000]';
DATA =     cos(TIME*2*pi/10) + AR_make(0.95,length(TIME)) + ...
       1i*[sin(TIME*2*pi/10) + AR_make(0.95,length(TIME))];
DATA =     0*cos(TIME*2*pi/10) + AR_make(0.95,length(TIME));
close all
figure
plot(TIME,real(DATA)); hold on
plot(TIME,imag(DATA))
figure
[SPEC,FF] = nunanspectrum(DATA, ...
TIME, 'X', ...
'Segments',1,'Window','rectwin','Plot',true,'PlotSegments',true);

(var(real(DATA)) + var(imag(DATA))) ./ sum(SPEC(:)*FF(1))
