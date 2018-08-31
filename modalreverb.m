%% Modal Filter Reverb
%  ELEC-E5630 - Acoustics and Audio Technology Seminar
%  Juri Lukkarila
%  2016 / 2018

close all; clearvars; clc

%% Modal Filter Reverb: number of modes

rng(0,'twister');                               % random seed
fs    = 44100;                                  % sample rate
len   = 4;                                      % length of audio (s)
modes = [128 256 512 1024 2048 4096 8192];      % number of modes

% octave band
modes_per_octave = [1 2 4 8 16 32 64 128];
f_oct = [63 125 250 500 1000 2000 4000 8000];   % middle freqs
f1    = ceil(f_oct./sqrt(2));                   % lower freq limit
f2    = floor(sqrt(2).*f_oct);                  % upper freq limit

printfigures = 1; % false = 0, true = 1;

if printfigures == 1
% plot modes per octave
figure(); loglog(f_oct, modes_per_octave, '.:','MarkerSize',14); 
axis([31.5 16000 1 128]); grid on;
xlabel('Frequency (Hz)'); ylabel('Number of modes')
set(gca,'XTick',[63 125 250 500 1000 2000 4000 8000])
set(gca,'XTickLabel',{63 125 250 500 '1k' '2k' '4k' '8k'})
set(gca,'YTick',[1 2 4 8 16 32 64 128])
set(gca,'YTickLabel',{'x' '2x' '4x' '8x' '16x' '32x' '64x' '128x'})
set(gca,'XMinorTick', 'off','YMinorTick','off','MinorGridLineStyle','none'); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],'PaperSize',[10.2 5.9])
print(gcf,'MFR2018_modes_octave', '-dpdf', '-painters');
end

if printfigures == 1
% plot octave bands
figure(); stem(f_oct, ones(1,length(f_oct)), ':r.'); hold on;
axis([31 16000 0 1.05]); grid on; set(gca,'XScale','log','YTick',[]); 
set(gca,'XTick',[63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
xlabel('Frequency (Hz)'); ylabel('Octave band')
% plot octave bands as rectangles
for n = 1:length(f_oct) 
    rectangle('Position',[f1(n) 0 f2(n)-f1(n) 1]);
end
end

%% Loop through different number of modes
for l = 1:length(modes)

% number of modes
M = modes(l);

% scale number of modes per octave
if l == 1
    modes_per_octave = [1 1 2 4 8 16 32 64];
elseif l == 2
    modes_per_octave = [1 2 4 8 16 32 64 128];
elseif l >= 3
    modes_per_octave = modes_per_octave.*2;   
end

% make sure we get the correct amount of modes
while sum(modes_per_octave) < M
    modes_per_octave(end) = modes_per_octave(end) + 1;
end

fprintf('Modes = %4d, dist. %2d %2d %3d %3d %3d %4d %4d %4d, sum = %4d\n',...
    M, modes_per_octave, sum(modes_per_octave))

% random frequencies inside each octave band
f = zeros(1,M);
index = 1;
for i = 1:length(modes_per_octave)
    for k = index:(index + modes_per_octave(i)-1)
        if k > M
            break
        end
        f(k) = (f2(i)-f1(i)).*rand + f1(i); % freq betweeen band limits
    end
    index = k + 1;
end

f = sort(round(f,1));   % frequency vector
w = 2*pi.*f./fs;        % normalized frequency

if printfigures == 1
% plot modal spacing
figure(); stem(f, ones(1,length(f)),'Marker','None', 'LineWidth',0.05);
axis([40 16000 0 1.1]);
xlabel('Frequency (Hz)'); ylabel('Mode')
set(gca,'XScale','log','YTick', [])
set(gca,'XTick',[63 125 250 500 1000 2000 4000 8000])
set(gca,'XTickLabel',{63 125 250 500 '1k' '2k' '4k' '8k'})
set(gca,'XMinorTick', 'off','YMinorTick','off','MinorGridLineStyle','none'); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR2018_modal_spacing_', int2str(M)), '-dpdf', '-painters'); 
end

%% RT

r1 = 3; % reverberation time, lowest freq (s) 

% parameters for RT curve smoothing for each mode count
r_begin       = [3 6 12 14 18 22 26];
smooth_factor = [0.2 0.2 0.1 0.05 0.04 0.02 0.01]; 

% RT curve with some ramdomness
t  = linspace(0,2,M); % time vector for exponential function.
RT = r1.*(0.1.*rand(1,M) + 0.9).*exp(-2*t);
RT(1:r_begin(l)) = r1; % force first values to chosen rt
RT = smooth(RT,smooth_factor(l),'loess');

if printfigures == 1
% plot reverberation time
figure(); semilogx(f,RT); grid on; hold on; stem(f,RT,':','Marker','None');
xlabel('Frequency (Hz)'); ylabel('RT60 (s)'); axis([31.5 16000 0 1.1*r1]);
set(gca,'XTick',[63 125 250 500 1000 2000 4000 8000])
set(gca,'XTickLabel',{63 125 250 500 '1k' '2k' '4k' '8k'})
set(gca,'XMinorTick', 'off','YMinorTick','off','MinorGridLineStyle','none'); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR2018_RT60_', int2str(M)), '-dpdf', '-painters'); 
end

%% mode parameters

% mode damping
a = log(1000)./(RT.*fs);

% mode gain
low     = round(M/25);             % smaller variation for low freq modes
high    = M - low;                 % rest of the modes
gp_low  = 0.3.*rand(1,low)+0.7;    % 0.3 - 1     
gp_high = 0.89.*rand(1,high)+0.01; % 0.01 - 0.9
% combine with decreasing curve
gp = horzcat(gp_low, gp_high).*nthroot(logspace(0,-0.5,M),2);
gp = gp ./ max(gp); % normalize

% random phases 0...2pi
phase = (2*pi).*rand(1,M);

% complex mode gain
g = gp.*exp(1i.*phase);

if printfigures == 1
% plot mode gain
figure(); stem(f,gp,'Marker','None','LineWidth',0.05); 
ylabel('mode gain'); xlabel('Frequency (Hz)'); axis([31.5 16000 0 1]);
set(gca,'XTick',[63 125 250 500 1000 2000 4000 8000], 'XScale','log')
set(gca,'XTickLabel',{63 125 250 500 '1k' '2k' '4k' '8k'}); grid on;
set(gca,'XMinorTick', 'off', 'XMinorGrid', 'off'); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR2018_gain_', int2str(M)), '-dpdf', '-painters'); 
end

%% impulse response

x       = zeros(1,len*fs); x(1) = 1;    % input
N       = length(x);                    % samples
ts      = 1/fs;                         % sample time 
xt      = 0:ts:(N-1)*ts;                % time vector
y       = zeros(1,N);                   % output
ym_prev = zeros(1,M);                   % init previous step

% filter loop
for t = 1:N                         
    ym = g.*x(t) + exp(1i.*w-a').* ym_prev; % filter
    y(t) = sum(ym);                         % sum to output
    ym_prev = ym;                           % store for previous
end

%% export audio

out = real(y);
norm = max(abs(min(out)),max(out)); % biggest value
imp = out ./ norm;                  % normalize
soundsc(imp,fs)
audiowrite(strcat('MFR2018_impulse_', int2str(M),'.wav'), imp, fs); 

%% plot impulse response

if printfigures == 1
figure(); plot(xt, imp); grid on; axis auto; xlim([0 2]);
xlabel('Time (s)'); ylabel('Amplitude');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR2018_imp_', int2str(M)), '-dpdf', '-painters'); 

% log squared
y_abs  = imp.^2;
y_norm = y_abs./max(y_abs);
ylog   = 10*log10(y_norm);

figure(); plot(xt, ylog); grid on; ylim([-100 0])
xlabel('Time (s)'); ylabel('Amplitude (dB)');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR2018_imp_log_', int2str(M)), '-dpdf', '-painters'); 
end
fprintf('... impulse done ');

%% spectrogram
if printfigures == 1
window = round(N/64);       % window length
if mod(window, 2) ~= 0      % if odd
    window = window + 1;    
end
overlap = window/2;         % window overlap 50%
nfft = 2^16;                % fft points
figure();
spectrogram(y,hamming(window),overlap,nfft,fs,'yaxis','MinThreshold',-80,'power'); 
set(gca,'YScale','log'); axis([0 len 0.03 16]); ylabel('Frequency (Hz)');
set(gca,'YTick',[0.0315 0.063 0.125 0.25 0.5 1 2 4 8 16]); xlabel('Time (s)');
set(gca,'YTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
% change colorbar size and location
c = colorbar; ax = gca; axpos = ax.Position; cpos = c.Position;
cpos(3) = 0.5*cpos(3); c.Position = cpos; ax.Position = axpos;
c.Label.String = 'Power (dB)';
set(gcf,'PaperUnits','centimeters','PaperPosition',[-0.4 0 18 9.9],'PaperSize',[16.4 9.6])
print(gcf,strcat('MFR2018_spectrogram_', int2str(M)), '-dpdf', '-r300'); 
end
fprintf('spectrogram done ');

%% Vocal sample

[vocal, fs] = audioread('vocal2.wav');
N           = length(vocal);
y           = zeros(1,N);
ym_prev     = zeros(1,M);

% filter loop
for t = 1:N                         
    ym      = g.*vocal(t) + exp(1i.*w-a').* ym_prev;
    y(t)    = sum(ym);
    ym_prev = ym;
end

% Export audio

out = real(y);
norm = max(abs(min(out)),max(out));
reverb = out ./ norm;
%soundsc(reverb,fs)
audiowrite(strcat('MFR_vocal_reverb_', int2str(M),'.wav'), reverb', fs); 

% mixed with dry signal
mix = 0.6.*vocal' + 0.4.*reverb;
%soundsc(mix,fs)
audiowrite(strcat('MFR_vocal_mix_', int2str(M),'.wav'), mix, fs); 
fprintf('vocal done ');

%% Drum beat sample

[beat, fs] = audioread('beat.wav');
N          = length(beat);
y          = zeros(1,N);
ym_prev    = zeros(1,M);

% filter loop
for t = 1:N
    ym      = g.*beat(t) + exp(1i.*w-a').* ym_prev;
    y(t)    = sum(ym); 
    ym_prev = ym;
end

% Export audio

out = real(y);
norm = max(abs(min(out)),max(out));
reverb = out ./ norm;
%soundsc(reverb,fs)
audiowrite(strcat('MFR_beat_reverb_', int2str(M),'.wav'), reverb, fs); 

% mixed with dry signal
mix = 0.6.*beat' + 0.4.*reverb;
%soundsc(mix,fs)
audiowrite(strcat('MFR_beat_mix_', int2str(M),'.wav'), mix, fs); 
fprintf('beat done ...\n');

close all;

end

fprintf('all done!\n');