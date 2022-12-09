%% calculation
clear;
close all;
clc;

bbp = [0.4222, 0.4215, 0.3188, 0.2203]; % == code length
compression = [56.8404, 56.9407, 75.2879, 108.9490];
PSNR = [31.3566, 31.3564, 29.5680, 27.6620];

%% figure

figure(); %1
subplot(1,2,1);
plot(bbp,compression,'o-');
xlabel('bits per pixel');
ylabel('Compression Rate');
title('Compression Rate');

subplot(1,2,2); %2
plot(bbp,PSNR,'o-');
xlabel('bits per pixel');
ylabel('PSNR');
title('PSNR');

