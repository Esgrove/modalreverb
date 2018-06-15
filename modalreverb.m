%% Modal Filter Reverb
%  ELEC-E5630 - Acoustics and Audio Technology Seminar
%  Juri Lukkarila
%  2016

close all; clearvars; clc

%% Modal Filter Reverb: number of modes

rng(0,'twister');                               % random seed
fs = 44100;                                     % sample rate
len = 3;                                        % length of input (s)
modes = [125 250 500 1000 2000 4000];           % number of modes

% octave band
modes_per_octave = [1 2 4 8 16 32 64 128];
f_oct = [63 125 250 500 1000 2000 4000 8000];   % middle freqs
f1 = floor(f_oct./sqrt(2));                     % lower freq limit
f2 = floor(sqrt(2).*f_oct);                     % upper freq limit

% plot modes per octave
figure(); loglog(f_oct, modes_per_octave, '.:', 'MarkerSize',14); 
axis([31.5 16000 1 128]); grid on;
xlabel('Frequency (Hz)'); ylabel('Number of modes')
set(gca,'XTick',[63 125 250 500 1000 2000 4000 8000])
set(gca,'XTickLabel',{63 125 250 500 '1k' '2k' '4k' '8k'})
set(gca,'YTick',[1 2 4 8 16 32 64 128])
set(gca,'YTickLabel',{'x' '2x' '4x' '8x' '16x' '32x' '64x' '128x'})
set(gca,'XMinorTick', 'off','YMinorTick','off','MinorGridLineStyle','none'); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
print(gcf,'MFR_modes_octave', '-dpdf', '-painters');

%% Loop through different number of modes
for l = 1:6

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

fprintf('Modes = %d, dist. %d %d %d %d %d %d %d %d, sum = %d\n',...
    M, modes_per_octave, sum(modes_per_octave))

% random frequencies inside each octave band
f = zeros(1,M);
index = 1;
for i = 1:8
    for k = index:(index + modes_per_octave(i)-1)
        if k > M
            break
        end
        f(k) =  (f2(i)-f1(i)).*rand + f1(i); % freq. betweeen band limits
    end
    index = k + 1;
end

% frequency vector
f = sort(round(f,1));
w = 2*pi.*f./fs;

% plot modal spacing
figure(); stem(f, ones(1,length(f)),'Marker','None', 'LineWidth',0.05);
axis([40 16000 0 1.1]);
xlabel('Frequency (Hz)'); ylabel('Mode')
set(gca,'XScale','log','YTick', [])
set(gca,'XTick',[63 125 250 500 1000 2000 4000 8000])
set(gca,'XTickLabel',{63 125 250 500 '1k' '2k' '4k' '8k'})
set(gca,'XMinorTick', 'off','YMinorTick','off','MinorGridLineStyle','none'); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR_modal_spacing_', int2str(M)), '-dpdf', '-painters'); 

%% RT

r1 = 3;     % reverberation time, lowest freq (s) 

r_begin = [3 6 12 14 18 22];
smooth_factor = [0.2 0.2 0.1 0.05 0.04 0.02]; 

% RT curve with some ramdomness
t = linspace(0,2,M);        % time vector for exponential function.
RT = r1.*(0.1.*rand(1,M) + 0.9).*exp(-2*t);
RT(1:r_begin(l)) = r1;              % force first values to chosen number
RT = smooth(RT,smooth_factor(l),'loess');

% plot
figure(); semilogx(f,RT); grid on;
xlabel('Frequency (Hz)'); ylabel('RT60 (s)');
axis([31.5 16000 0 1.1*r1]);
set(gca,'XTick',[63 125 250 500 1000 2000 4000 8000])
set(gca,'XTickLabel',{63 125 250 500 '1k' '2k' '4k' '8k'})
set(gca,'XMinorTick', 'off','YMinorTick','off','MinorGridLineStyle','none'); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR_RT60_', int2str(M)), '-dpdf', '-painters'); 

%% mode damping
a = log(1000)./(RT.*fs);

low = round(M/25);
high = M - low;

% mode gain
gp_low = 0.3.*rand(1,low)+0.7;                              % 0.3 - 1
gp_high = 0.89.*rand(1,high)+0.01;                          % 0.01 - 0.9
gp = horzcat(gp_low, gp_high).*nthroot(logspace(0,-1,M),2);

figure(); stem(f,gp,'Marker','None', 'LineWidth',0.05); 
ylabel('mode gain'); xlabel('Frequency (Hz)'); axis([31.5 16000 0 1]);
set(gca,'XTick',[63 125 250 500 1000 2000 4000 8000], 'XScale','log')
set(gca,'XTickLabel',{63 125 250 500 '1k' '2k' '4k' '8k'})
set(gca,'XMinorTick', 'off'); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR_gain_', int2str(M)), '-dpdf', '-painters'); 

% random phases 0...2pi
phase = (2*pi).*rand(1,M);

% complex mode gain
g = gp.*exp(1i.*phase);

%% impulse response
x = zeros(1,len*fs);          
x(1) = 1;

N = length(x);              % samples
ts = 1/fs;                  % sample time 
xt = 0:ts:(N-1)*ts;         % time vector

ym = zeros(1,M);            % current values
y = zeros(1,N);             % output

% first value, ym(t-1) = 0
for m = 1:M
    ym(m) = g(m).*x(1);
end

y(1) = sum(ym);

% filter loop
for t = 2:N
    ym_prev = ym;       % store previous value
    for m = 1:M
        ym(m) = g(m).*x(t) + exp(1i*w(m)-a(m))*ym_prev(m);
    end
    y(t) = sum(ym);     % store output
end

%% export audio
miny = min(real(y));
maxy = max(real(y));
if abs(miny) >= maxy
    impulse = real(y)./abs(miny);
else
    impulse = real(y)./abs(maxy);
end
%soundsc(impulse,fs)
audiowrite(strcat('MFR_impulse_', int2str(M),'.wav'), impulse, fs); 

%% Impulse Response plot

figure(); plot(xt, impulse); grid on; axis auto; xlim([0 2]);
xlabel('Time (s)'); ylabel('Amplitude');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR_imp_', int2str(M)), '-dpdf', '-painters'); 

% log squared
y_abs = impulse.^2;
y_norm = y_abs./max(y_abs);
ylog = 10*log10(y_norm);

figure(); plot(xt, ylog); grid on; ylim([-100 0])
xlabel('Time (s)'); ylabel('Amplitude (dB)');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
print(gcf,strcat('MFR_imp_log_', int2str(M)), '-dpdf', '-painters'); 

%% spectrogram
window = round(N/100);      % divide to approx. 100 windows
if mod(window, 2) ~= 0      % if odd
    window = window + 1;    
end
overlap = window/2;
nfft = N;
figure();
spectrogram(y,window,overlap,nfft,fs,'yaxis','MinThreshold', -100, 'power'); 
set(gca,'YScale','log'); axis([0 len 0.03 16]); ylabel('Frequency (Hz)');
set(gca,'YTick',[0.035 0.063 0.125 0.25 0.5 1 2 4 8 16]); xlabel('Time (s)');
set(gca,'YTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
% change colorbar size and location
c = colorbar; ax = gca; axpos = ax.Position; cpos = c.Position;
cpos(3) = 0.5*cpos(3); c.Position = cpos; ax.Position = axpos;
c.Label.String = 'Power (dB)';
set(gcf,'PaperUnits','centimeters','PaperPosition',[-0.4 0 18 9.9],...
        'PaperSize',[16.4 9.6])
print(gcf,strcat('MFR_spectrogram_', int2str(M)), '-dpdf', '-r300'); 
%% 3D
% view([100 50]); colorbar('off');
% set(gcf,'PaperUnits','centimeters','PaperPosition',[-1.2 0 21.4 12],...
%         'PaperSize',[19.4 11.6])
% print(gcf,strcat('MFR_spectrogram_3D_', int2str(M)), '-dpdf', '-r300'); 

%% Vocal sample

[vocal, fs] = audioread('vocal2.wav');

N = length(vocal);

ym = zeros(1,M);
y = zeros(1,N);

% first value, ym(t-1) = 0
for m = 1:M
    ym(m) = g(m).*vocal(1);
end

y(1) = sum(ym);

% filter loop
for t = 2:N
    ym_prev = ym;
    for m = 1:M
        ym(m) = g(m).*vocal(t) + exp(1i*w(m)-a(m))*ym_prev(m);
    end
    y(t) = sum(ym);  
end

%% Export audio
miny = min(real(y));
maxy = max(real(y));
if abs(miny) >= maxy
    reverb = real(y)./abs(miny);
else
    reverb = real(y)./abs(maxy);
end
%soundsc(reverb,fs)
audiowrite(strcat('MFR_vocal_reverb_', int2str(M),'.wav'), reverb', fs); 

%% mixed with dry signal
mix = (0.7.*vocal)' + 0.3.*reverb;
%soundsc(real(vocal_mix),fs)
audiowrite(strcat('MFR_vocal_mix_', int2str(M),'.wav'), mix, fs); 

%% Drum beat sample

[beat, fs] = audioread('beat.wav');

N = length(beat);

ym = zeros(1,M);
y = zeros(1,N);

% first value, ym(t-1) = 0
for m = 1:M
    ym(m) = g(m).*beat(1);
end

y(1) = sum(ym);

% filter loop
for t = 2:N
    ym_prev = ym;
    for m = 1:M
        ym(m) = g(m).*beat(t) + exp(1i*w(m)-a(m))*ym_prev(m);
    end
    y(t) = sum(ym);  
end

%% Export audio
miny = min(real(y));
maxy = max(real(y));
if abs(miny) >= maxy
    reverb = real(y)./abs(miny);
else
    reverb = real(y)./abs(maxy);
end
%soundsc(reverb,fs)
audiowrite(strcat('MFR_beat_reverb_', int2str(M),'.wav'), reverb, fs); 

%% mixed with dry signal
mix = (0.7.*beat)' + 0.3.*reverb;
%soundsc(real(vocal_mix),fs)
audiowrite(strcat('MFR_beat_mix_', int2str(M),'.wav'), mix, fs); 

end

close all;