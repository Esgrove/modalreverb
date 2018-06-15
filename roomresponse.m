%% ELEC-E5630 - Acoustics and Audio Technology Seminar
%  Juri Lukkarila
%  2016

%% Frequency response

% Import measurement data from text file: 
filename = 'D:\Dropbox\KOULU\AKU SEMINAR\Latex\figures\PC_stereo_uus.csv';
delimiter = ',';
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
    'MultipleDelimsAsOne', true);
fclose(fileID);
freq = dataArray{:, 1};
mag = dataArray{:, 2};

% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

% apply 1/96 octave smoothing
mag = smooth_spectrum(mag,freq,96);

%% Plot frequency response
figure()
semilogx(freq, mag); grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
axis([20 22000 30 110])
set(gca,'XTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
set(gcf,'PaperUnits','centimeters','PaperPosition',[-1.6 0 24 8],...
        'PaperSize',[20.6 7.8])
print(gcf, 'freqresponse', '-dpdf', '-painters'); 

%% Impulse response

% read measured impulse response
[imp, Fs] = audioread('PC_stereo_uus.wav',[47200,inf]);
N = length(imp);            % samples
Ts = 1/Fs;                  % sample time 
t = 0:Ts:(N-1)*Ts;          % time vector

figure();
plot(t, imp); grid on;
xlabel('Time (s)');
ylabel('Amplitude');
axis([0 0.2 -1 1]);
set(gca,'XTick',0:0.025:0.2);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
print(gcf, 'impulse', '-dpdf', '-painters'); 

%% Log squared

imp_db = 10*log10(imp.^2);

figure();
plot(t, imp_db); grid on;
xlabel('Time (s)');
ylabel('Magnitude (dB)')
axis([0 0.4 -80 0]);
set(gca,'XTick',0:0.05:0.4);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
print(gcf, 'impulse_db', '-dpdf', '-painters');

%% Waterfall plot

% spectrogram parameters
window = floor(N/8);
overlap = window/2;
nfft = N;

[s,f,t] = spectrogram(imp,window,overlap,nfft,Fs);

% real values
S = abs(s);

% matrix size
dims = size(s);  

% initiliaze matrix for loop
S_norm = zeros(dims(1), dims(2));
S_smooth = zeros(dims(1), dims(2));

% normalize
s_max = max(max(S));        % max value
for n = 1:dims(2)
    S_norm(:,n) = S(:,n)./s_max;
end

% log scale
S_db = 20*log10(S_norm);

%% 1/48 octave smoothing
for n = 1:dims(2)
    S_smooth(:,n) = smooth_spectrum(S_db(:,n),f,48);
end

%% plot

% custom colormap
map = [0 0 0.3 ; 0 0 0.4 ; 0 0 0.5 ; 0 0 0.6 ;0 0 0.7 ; 0 0 0.8 ;0 0 0.9 ; 0 0 1];

figure()
h = waterfall(t,f,S_smooth); colormap(map); 
set(h, 'LineWidth',0.05, 'LineSmoothing','on'); grid on;
view([75 40]); axis([0.1 0.75 20 20000 -82 0]); 
xlabel('Time (s)'); ylabel('Frequency (Hz)'); zlabel('Magnitude (dB)')
set(gca,'YTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'YTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
set(gca,'YMinorTick','off','MinorGridLineStyle','none'); 
set(gca,'YScale','log'); set(gco, 'Linesmoothing','on');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 24 13.5],...
        'PaperSize',[24 13.5])
print(gcf, 'waterfall', '-dpng', '-r600'); 
