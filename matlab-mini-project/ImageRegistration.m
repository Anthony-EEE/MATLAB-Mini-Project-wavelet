function [Tx_RGB Ty_RGB]= ImageRegistration
% *************************************************************************
% Wavelets and Applications Course - Dr. P.L. Dragotti
% MATLAB mini-project 'Sampling Signals with Finite Rate of Innovation'
% Exercice 6
% *************************************************************************
% 
% FOR STUDENTS
%
% This function registers the set of 40 low-resolution images
% 'LR_Tiger_xx.tif' and returns the shifts for each image and each layer
% Red, Green and Blue. The shifts are calculated relatively to the first
% image 'LR_Tiger_01.tif'. Each low-resolution image is 64 x64 pixels.
%
%
% OUTPUT:   Tx_RGB: horizontal shifts, a 40x3 matrix
%           Ty_RGB: vertical shifts, a 40x3 matrix
%
% NOTE: _Tx_RGB(1,:) = Ty_RGB(1,:) = (0 0 0) by definition.
%       _Tx_RGB(20,2) is the horizontal shift of the Green layer of the
%       20th image relatively to the Green layer of the firs image.
%
%
% OUTLINE OF THE ALGORITHM:
%
% 1.The first step is to compute the continuous moments m_00, m_01 and m_10
% of each low-resolution image using the .mat file called:
% PolynomialReproduction_coef.mat. This file contains three matrices
% 'Coef_0_0', 'Coef_1_0' and 'Coef_0_1' used to calculate the continuous
% moments.
%
% 2.The second step consists in calculating the barycenters of the Red,
% Green and Blue layers of the low-resolution images.
%
% 3.By computing the difference between the barycenters of corresponding 
% layers between two images, the horizontal and vertical shifts can be 
% retrieved for each layer.
%
%
% Author:   Loic Baboulaz
% Date:     August 2006
%
% Imperial College London
% *************************************************************************


% Load the coefficients for polynomial reproduction
load('PolynomialReproduction_coef.mat','Coef_0_0','Coef_1_0','Coef_0_1');

% -------- include your code here -----------
Tx_RGB=ones(40,3);
Ty_RGB=ones(40,3);
%GENERATE Reference of the first camera or first image
[R, G ,B] = ImageProcessing(1);
dx1=0;
dy1=0;
%continous moments calculation in Red layer
m00_R1=sum(Coef_0_0.*R);
m01_R1=sum(Coef_0_1.*R);
m10_R1=sum(Coef_1_0.*R);
%(x1,y1)of Red
x_R1=m10_R1/m00_R1;
y_R1=m01_R1/m00_R1;

%continous moments calculation in Green layer
m00_G1=sum(Coef_0_0.*G);
m01_G1=sum(Coef_0_1.*G);
m10_G1=sum(Coef_1_0.*G);
%(x1,y1)of GREEN
x_G1=m10_G1/m00_G1;
y_G1=m01_G1/m00_G1;

%continous moments calculation in Blue layer
m00_B1=sum(Coef_0_0.*B);
m01_B1=sum(Coef_0_1.*B);
m10_B1=sum(Coef_1_0.*B);
%(x1,y1)of Blue channel
x_B1=m10_B1/m00_B1;
y_B1=m01_B1/m00_B1;



for i=1:40
    % sampling by imageprocessing with normalization and threshold 
    % setting to limit noise
    [R, G ,B] = ImageProcessing(i);
    
    %continous moments calculation in Red layer
    m00_R=sum(Coef_0_0.*R);
    m01_R=sum(Coef_0_1.*R);
    m10_R=sum(Coef_1_0.*R);
    % shifts (dx,dy) of Red
    x_R=m10_R/m00_R;
    y_R=m01_R/m00_R;
    dx_R=x_R-x_R1;
    dy_R=y_R-y_R1;
    %feed shifts in the Tx_RGB and Ty_RGB
    Tx_RGB(i,1)=dx_R;
    Ty_RGB(i,1)=dy_R;


    %continous moments calculation in Green layer
    m00_G=sum(Coef_0_0.*G);
    m01_G=sum(Coef_0_1.*G);
    m10_G=sum(Coef_1_0.*G);
    % shifts (dx,dy) of Green
    x_G=m10_G/m00_G;
    y_G=m01_G/m00_G;
    dx_G=x_G-x_G1;
    dy_G=y_G-y_G1;
    %feed shifts in the Tx_RGB and Ty_RGB
    Tx_RGB(i,2)=dx_G;
    Ty_RGB(i,2)=dy_G;


    %continous moments calculation in Blue layer
    m00_B=sum(Coef_0_0.*B);
    m01_B=sum(Coef_0_1.*B);
    m10_B=sum(Coef_1_0.*B);
    % shifts (dx,dy) of Blue
    x_B=m10_B/m00_B;
    y_B=m01_B/m00_B;
    dx_B=x_B-x_B1;
    dy_B=y_B-y_B1;
    %feed shifts in the Tx_RGB and Ty_RGB
    Tx_RGB(i,3)=dx_B;
    Ty_RGB(i,3)=dy_B;


end







