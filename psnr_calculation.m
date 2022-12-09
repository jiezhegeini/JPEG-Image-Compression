clear;
close all;
clc;

im1 = imread('/Users/apple/Desktop/34240/project/good/image.bmp');
im2 = imread('/Users/apple/Desktop/34240/project/good/reconstructed_image.bmp');
PSNR = PeakSignaltoNoiseRatio(im1,im2)


%% PSNR

function PSNR = PeakSignaltoNoiseRatio(im,imq)
    [height width] = size(im);
    temp = double(im)-double(imq); 
    temp = sum(temp.^2,'all');
    MSE = temp/height/width;
    PSNR = 10*log10(255^2/MSE);
end