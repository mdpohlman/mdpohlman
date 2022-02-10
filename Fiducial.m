%Mason Pohlman, 2-2-21, Image Stabilizer
%Code to read in TIF files, find fiducial marking, and adjust each frame

clc
clear
close all
tic

%Read in data
foldername = sprintf('C:/Users/mpohlma3/Documents/HORIZON/Other/allframes_S0001/');
a = dir('C:/Users/mpohlma3/Documents/HORIZON/Other/allframes_S0001/*.tif');
number = length(a);
pic = imread([foldername,'allframes_S0001000001.tif']);
[r,c] = size(pic);

tol = 17.85; %Intensity value that marks edge of square
horizontal = zeros(1,number);
vertical = zeros(1,number);
diffh = zeros(1,number);
diffv = zeros(1,number);
horizontal2 = zeros(1,number);
vertical2 = zeros(1,number);

%Loads each individual image
for i = 1:number
    filename = sprintf('allframes_S0001%06d.tif',i);
    file = [foldername,filename];
    picture = imread(file);
    
    ii = find(picture(10,900:c) <= tol,1,'first') + 900 - 1;
    
    jj = find(picture(1:r,990) >= tol,1,'first');
    
    horizontal(i) = ii;
    vertical(i) = jj;
end

floorh = min(horizontal);
floorv = min(vertical);

for i = 1:number
    diffh(i) = -horizontal(i) + floorh;
    diffv(i) = -vertical(i) + floorv;
    
    filename = sprintf('allframes_S0001%06d.tif',i);
    file = [foldername,filename];
    picture = imread(file);
    
    new = circshift(picture,[diffv(i),diffh(i)]);
    
    ii2 = find(new(10 + diffv(i),900:c) <= tol,1,'first') + 900 - 1;
    jj2 = find(new(1:r,990 + diffh(i)) >= tol,1,'first');
    horizontal2(i) = ii2; 
    vertical2(i) = jj2;
    
    %DANGER: Saves off stabilized images as TIF images
%     output = sprintf('C:/Users/mpohlma3/Documents/HORIZON/Other/Stabilized/pic%d.tif',i);
%     imwrite(new,output);
end

figure(1)
subplot(2,2,1)
plot(horizontal)
ylabel('Horizontal Position (Pixel #)')
subplot(2,2,3)
plot(horizontal2)
xlabel('Frame #')
ylabel('Horizontal Position (Pixel #)')

subplot(2,2,2)
plot(vertical)
ylabel('Vertical Position (Pixel #)')
subplot(2,2,4)
plot(vertical2)
xlabel('Frame #')
ylabel('Vertical Position (Pixel #)')
movegui([500,100])

toc