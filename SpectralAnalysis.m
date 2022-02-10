%Mason Pohlman, 3-12-21, HORIZON Group
%Takes in video, positions it on a ROI, and performs a PSD spectral
%analysis on the video. Also, loads swag.

clc
clear
close all

fprintf('Loading 100kHz...');
tic;
file = '\\150.182.19.62\Horizon\Gragston\Mach 2.3\xoverd~10.1_S0001\xoverd~10.1_S0001.mat';
video = load(file);
t1=toc;
[r,c,d] = size(video.video);
fprintf('%f sec, Loading #2...',t1);

tic;
%Resizes and rotates individual frames
for i = 1:5000
  frame = video.video(:,:,i);
  tilt = imrotate(frame, 10);
  focus = tilt(72:102,85:143);
  focused(:,:,i) = focus;
end
t2=toc;

imagesc(focused(:,:,1000))

doub = im2double(focused);
[R,C,D] = size(doub);

fprintf('%f sec, Plastic surgery...',t2);
tic;
%Transforms matrix into useful form for pwelch
for j = 1:D
    m(j,:) = reshape(doub(:,:,j), [1,R*C]);
end
t3=toc;
fprintf('%f sec, Here comes the PSD...',t3);

tic;
%PSD Transformation
fs = 100000;
ffts = 2^10;
df = fs/ffts;
%m = m - mean(m);
[psd1,f1] = pwelch(m,hann(D),0,[],fs);
[psd,f] = pwelch(m,hann(D),round(0.5*ffts),ffts,fs);
G1 = (psd*df)/(var(m));
G2 = (psd.*f)/(var(m));

%Two normalizations, regular and SBLI community (*f, not df)
t4=toc;

fprintf('%f sec, Finna write video...',t4);

%All to write video
factor = 1/200; %inches/pixel ratio
x = (1:C)*factor;
y = (1:R)*factor;
writer = VideoWriter('C:\Users\mpohlma3\Documents\HORIZON\Schlieren 2.3\50.avi');
open(writer);

% figure()
% movegui([500 100]);
% %averages = mean(psd,2);
% set(gcf,'color','w');
% loglog(f,psd)
% xlabel('Frequency (Hz)')
% ylabel('PSD')
% title('Power Spectral Density of Schlieren')
% grid on

% figure()
% movegui([500 100]);
% set(gcf,'color','w');
% loglog(f1,psd1)
% xlabel('Frequency (Hz)')
% ylabel('PSD')
% title('Power Spectral Density of Schlieren')
% grid on

% figure()
% movegui([500 100]);
% G1 = mean(G1,2);
% set(gcf,'color','w');
% loglog(f,G1)
% xlabel('Frequency (Hz)')
% ylabel('PSD')
% title('Power Spectral Density of Schlieren')
% grid on

% figure()
% movegui([500 100]);
% G2 = mean(G2,2);
% set(gcf,'color','w');
% loglog(f,G2)
% xlabel('Frequency (Hz)')
% ylabel('PSD')
% title('Power Spectral Density of Schlieren')
% grid on

for i = 1:size(psd,1)
   PSDimage(:,:,i) = reshape(psd(i,:),[R,C]);
   subplot(2,1,1);
   imagesc(x,y,PSDimage(:,:,i));
   axis image
   set(gca,'XTick',[], 'YTick', [])
   grid on
   xlabel(num2str(f(i),'f = %.4f Hz'));
   frame = getframe(gcf);
   writeVideo(writer,frame);
   
   %Plotting PSD
   averages = mean(psd,2);
   subplot(2,1,2);
   set(gcf,'color','w');
   loglog(f,averages)
   xlabel('Frequency (Hz)')
   xline(f(i),'Color','r')
   ylabel('PSD')
   title('Power Spectral Density of Schlieren')
   grid on
end

mx = 8105;
figure()
imagesc(x,y,PSDimage(:,:,mx));
set(gcf,'color','w');
xlabel('Frequency (Hz)')
ylabel('PSD')
title('Peak Frequency')
grid on

hold off
close(writer)

fprintf('Done');



