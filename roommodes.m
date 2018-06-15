% ELEC-E5630 - Acoustics and Audio Technology Seminar
% Juri Lukkarila, 2016
% Room Modes

%% Axial

L = 4;          % length
W = 6;          % width
H = 3;          % height

v = 343;        % speed of sound
    
axial_L = zeros(1,3300);
axial_W = zeros(1,3300);
axial_H = zeros(1,3300);

% calculate modes
for a = 1:length(axial_L)
    axial_L(a) = v/2*sqrt((a/L)^2);
    axial_W(a) = v/2*sqrt((a/W)^2);
    axial_H(a) = v/2*sqrt((a/H)^2);
end

% combine axial modes
axial = sort(round(horzcat(axial_L,axial_W,axial_H),1));

% count unique
uv = unique(axial);
n  = histc(axial,uv);

%% modes
figure(); stem(uv,n, 'Marker', 'none'); set(gca,'XScale','log');
set(gca,'XTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
axis([20 24000 0 3.5]);
xlabel('Frequency (Hz)')
ylabel('Modes')

%% mode distribution
figure(); stem(axial(1:200), ones(1,length(axial(1:200))),'Marker', 'none',...
    'LineWidth',0.05); set(gca,'XScale','log');
set(gca,'XTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
axis([20 3000 0 1.1]);
xlabel('Frequency (Hz)')
ylabel('Mode')
set(gca, 'YTick', []);
set(gcf,'PaperUnits','centimeters','PaperPosition',[-0.6 0 16 9],...
        'PaperSize',[14.6 8.8])
print(gcf, 'axialmodes', '-dpdf', '-painters'); 

%% Octave Band

f = [31.5 63 125 250 500 1000 2000 4000 8000 16000];
f1 = floor(f./sqrt(2));
f2 = floor(sqrt(2).*f);

f_oct = horzcat(f1',f',f2');

% count modes on octave band
modes_per_octave = zeros(1,10);

for i = 1:10
    modes_per_octave(i) = sum(axial(:)>f1(i) & axial(:)<f2(i));
end

figure(); semilogx(f, modes_per_octave);
set(gca,'XTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
set(gca, 'XMinorTick', 'off', 'MinorGridLineStyle', 'none'); 
set(gca, 'YMinorTick', 'off');
xlim([31.5 20000]); grid on;
xlabel('Frequency (Hz)')
ylabel('Number of modes')
set(gcf,'PaperUnits','centimeters','PaperPosition',[-0.6 0 16 9],...
        'PaperSize',[14.6 8.8])
print(gcf, 'modes_per_octave', '-dpdf', '-painters'); 

%% Tangential

tangential_LW = zeros(1,10000);
tangential_LH = zeros(1,10000);
tangential_WH = zeros(1,10000);

for a = 1:100
    for b = 1:100
        i = (a-1)*100 + b;
        tangential_LW(i) =  v/2*sqrt((a/L)^2 + (b/W)^2);
        tangential_LH(i) =  v/2*sqrt((a/L)^2 + (b/H)^2);
        tangential_WH(i) =  v/2*sqrt((a/W)^2 + (b/H)^2);
    end
end

tangential = sort(round(horzcat(tangential_LW,tangential_LH,tangential_WH),1));

uv = unique(tangential);
n  = histc(tangential,uv);

figure(); stem(uv,n, 'Marker', 'none'); set(gca,'XScale','log');
set(gca,'XTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
axis([20 24000 0 20]);
xlabel('Frequency (Hz)')
ylabel('Modes')

%% Oblique

oblique = v/2*sqrt((1/L)^2 + (1/W)^2 + (1/H)^2);
for a = 1:22
    for b = 1:22
        for c = 1:22
            result = v/2*sqrt((a/L)^2 + (b/W)^2 + (c/H)^2);
            oblique = [oblique, result];
        end
    end
end

oblique = sort(round(oblique(2:end),1));

uv = unique(oblique);
n  = histc(oblique,uv);

figure(); stem(uv,n, 'Marker', 'none'); set(gca,'XScale','log');
set(gca,'XTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
axis([20 24000 0 20]);
xlabel('Frequency (Hz)')
ylabel('Modes')

%% All

modes = sort(horzcat(axial,tangential,oblique));

uv = unique(modes);
n  = histc(modes,uv);

figure(); stem(uv,n, 'Marker', 'none'); set(gca,'XScale','log');
set(gca,'XTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
axis([20 24000 0 20]);
xlabel('Frequency (Hz)')
ylabel('Unique modes')

figure(); stem(modes, ones(1,length(modes)),'Marker','none','LineWidth',0.01); 
set(gca,'XScale','log', 'YTick', []);
set(gca,'XTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000]);
set(gca,'XTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
axis([20 24000 0 1.1]); xlabel('Frequency (Hz)'); ylabel('Mode')
set(gcf,'PaperUnits','centimeters','PaperPosition',[-0.6 0 16 9],...
        'PaperSize',[14.6 8.8])
print(gcf, 'all_modes', '-dpdf', '-painters'); 