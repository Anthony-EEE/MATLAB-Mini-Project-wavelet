clc;clear;
% maximum degree of polynomials
Max_d=3;
% length of the signal
Len=2048;
% sampling period
T=64;
% signal range [0, 32]
Max_Sig=32;
MIn_Sig=0;
% time of sampling points
Time_res=1/64;
t=0:(Time_res):(Len-1)/T;
%numbers of shift,the first shift is [0,64],2nd shift is [65, 128]......
num_shifts=Len/T;

%define B-spline function and beta_0 is a box function
beta_0 = ones(1,T);
temp = beta_0;
for i=1:4
    phi_T = [conv(beta_0, temp)] / T;
    temp = phi_T;
end

%define Kernels that can reproduce polynomials
Kernels=zeros(num_shifts,2048);
phi = zeros(1,2048);

phi(1:length(phi_T))=phi_T;
for i=1:num_shifts   
    Kernels(i,:)=[zeros(1,(i-1)*T), phi(1: end-(i-1)*T)];
end

%dual basis exists because of b-spline rather that dB4!!!!
dualKernels=((Kernels.')*inv(Kernels*(Kernels.'))).';

%t^m
polynomials=zeros(Max_d,2048);

for degree=0:Max_d
    polynomials(degree+1,:)=t.^(degree);
end

%using formular of Cm,n=<t^m, phi>, this time phi is dual!!!
coefficients=zeros(Max_d, num_shifts);
coefficients=dualKernels*polynomials.';

%reproduce polynomials of maximum degree 3 using coefficients
repoly=Kernels.'*coefficients;
repoly=repoly';

%PLOT degree=0, original poly, reproduced poly, shafted kernel shapes
subplot(2,2,1)
plot(t,polynomials(1,:),'b')
hold on;
plot(t,repoly(1,:),"r")
hold on;
for i=1:num_shifts
    hold on;
    plot(t, Kernels(i,:).*coefficients(i,1),"k:")
end
legend("original polynomials"," poly reproduced by B-Spline","shifted kernels shapes")
title("Degree-0 polynomial")
%PLOT degree=1
subplot(2,2,2)
plot(t,polynomials(2,:),'b')
hold on;
plot(t,repoly(2,:),"r")
hold on;
for i=1:num_shifts
    hold on;
    plot(t, Kernels(i,:).*coefficients(i,2),"k:")
end
legend("original polynomials"," poly reproduced by B-Spline","shifted kernels shapes")
title("Degree-1 polynomial")
%PLOT degree=2
subplot(2,2,3)
plot(t,polynomials(3,:),'r')
hold on;
plot(t,repoly(3,:),"b")
hold on;
for i=1:num_shifts
    hold on;
    plot(t, Kernels(i,:).*coefficients(i,3),"k:")
end
legend("original polynomials"," poly reproduced by B-Spline","shifted kernels shapes")
title("Degree-2 polynomial")
%PLOT degree=3
subplot(2,2,4)
plot(t,polynomials(4,:),'r')
hold on;
plot(t,repoly(4,:),"b")
hold on;
for i=1:num_shifts
    hold on;
    plot(t, Kernels(i,:).*coefficients(i,4),"k:")
end
legend("original polynomials"," poly reproduced by B-Spline","shifted kernels shapes")
title("Degree-3 polynomial")


