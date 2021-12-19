function [Red, Green, Blue]= ImageProcessing(index)
%Read images
Image = imread(['LR_Tiger_', num2str(index,'%.2d'),'.tif']);
%integer to float type 
Norm_Image=double(Image./255);

Red = Norm_Image(:,:,1);
Green = Norm_Image(:,:,2);
Blue = Norm_Image(:,:,3);

%thereshold setting by median

thresh_R=median(median(Red));
Red(Red<thresh_R) = 0;

thresh_G=median(median(Green));
Green(Green<thresh_G) = 0;

thresh_B=median(median(Blue));
Blue(Blue<thresh_B) = 0;