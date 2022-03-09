clc
clear
close all

ffts = 2^12;
fs = 5e6;
df = fs/ffts;
row = 2;
beginSS = 160000;
endSS = 400000;

%Import FLDI Data
one = importdata('C:/Users/mpohlma3/Documents/Roughness Sleeves/FLDI Data/C4Smooth800000.dat');
meansubtracted1 = one(beginSS:endSS,row) - mean(one(beginSS:endSS,row));
[PSD1,f1] = pwelch(meansubtracted1,hann(ffts),round(0.5*ffts),ffts,fs);
G1 = (PSD1*df) / var(meansubtracted1);
t = one(beginSS:endSS,1);

two = importdata('C:/Users/mpohlma3/Documents/Roughness Sleeves/FLDI Data/C4Diamond00000.dat');
meansubtracted2 = two(beginSS:endSS,row) - mean(two(beginSS:endSS,row));
[PSD2,f2] = pwelch(meansubtracted2,hann(ffts),round(0.5*ffts),ffts,fs);
G2 = (PSD2*df) / var(meansubtracted2);

three = importdata('C:/Users/mpohlma3/Documents/Roughness Sleeves/FLDI Data/C4RealL00000.dat');
meansubtracted3 = three(beginSS:endSS,row) - mean(three(beginSS:endSS,row));
[PSD3,f3] = pwelch(meansubtracted3,hann(ffts),round(0.5*ffts),ffts,fs);
G3 = (PSD3*df) / var(meansubtracted3);

%Wave Packet Analysis
up = 470e3;
down = 490e3;
B1 = bandpass(meansubtracted1,[up down],fs);
B2 = bandpass(meansubtracted2,[up down],fs);
B3 = bandpass(meansubtracted3,[up down],fs);

[env1,low1] = envelope(B1);
[env2,low2] = envelope(B2);
[env3,low3] = envelope(B3);

peaks1 = findpeaks(env1);
peaks2 = findpeaks(env2);
peaks3 = findpeaks(env3);
bins = 100;

%Plotting
% figure()
% loglog(f1,G1,f2,G2,f3,G3)
% title('Varying Roughness Geometry (Channel 4)')
% xlabel('Frequency (Hz)')
% ylabel('G_{xx}*df / \sigma^2')
% l = legend('Smooth','Diamond','Realistic Large');
% grid on
% set(gcf,'color','w');
% 
% figure()
% loglog(f1,PSD1,f2,PSD2,f3,PSD3)
% title('Non-Normalized')
% xlabel('Frequency (Hz)')
% ylabel('PSD')
% l = legend('Smooth','Diamond','Realistic Large');
% grid on
% set(gcf,'color','w');

figure()
subplot(3,1,1)
histogram(peaks1,bins,'EdgeColor','none','FaceColor','#0072BD');
xlabel('Amplitude [mV]')
subplot(3,1,2)
histogram(peaks2,bins,'EdgeColor','none','FaceColor','#D95319')
xlabel('Amplitude [mV]')
subplot(3,1,3)
histogram(peaks3,bins,'EdgeColor','none','FaceColor','#EDB120')
xlabel('Amplitude [mV]')
set(gcf,'color','w');

