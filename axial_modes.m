% ELEC-E5630 - Acoustics and Audio Technology Seminar
% Juri Lukkarila, 2016

%% Axial Room Modes

L = 4;          % length
W = 6;          % width
H = 3;          % height

v = 343;        % speed of sound
    
axial_L = zeros(1,1000);
axial_W = zeros(1,1000);
axial_H = zeros(1,1000);

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
figure(); stem(axial(1:200), ones(1,length(axial(1:200))),'Marker', 'none'); set(gca,'XScale','log');
set(gca,'XTick',[31.5 63 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{31.5 63 125 250 500 '1k' '2k' '4k' '8k' '16k'})
axis([20 3000 0 1.1]);
xlabel('Frequency (Hz)')
ylabel('Mode')
set(gca, 'YTick', []);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
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
%xlim([31.5 20000]);
axis tight; grid on;
xlabel('Frequency (Hz)')
ylabel('Number of modes')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10.7 6],...
        'PaperSize',[10.2 5.9])
print(gcf, 'modes_per_octave', '-dpdf', '-painters'); 