clear;
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



%Daubechies:db4 filters with 4 vanishing moment can reproduce polynomial of Max_d=3, 
%and generate phi, psi, 
%and generate Kernel matrix combining with each shifted version
Kernels=zeros(num_shifts,2048);
phi = zeros(1,2048);
[phi_T, psi_T, xval]=wavefun('dB4', 6);
phi(1:length(phi_T))=phi_T;
for i=1:num_shifts   
    Kernels(i,:)=[zeros(1,(i-1)*T), phi(1: end-(i-1)*T)];
end

%t^m
polynomials=zeros(Max_d,2048);

for degree=0:Max_d
    polynomials(degree+1,:)=t.^(degree);
end

%using formular of Cm,n=<t^m, phi>
coefficients=zeros(Max_d, num_shifts);
coefficients=Kernels*polynomials.'/T;

%reproduce polynomials of maximum degree 3 using coefficients
repoly=Kernels.'*coefficients;
repoly=repoly';

%PLOT degree=0, original poly, reproduced poly, shafted kernel shapes
figure(1);
plot(t,polynomials(1,:),'b')
hold on;
plot(t,repoly(1,:),"r")
hold on;
for i=1:num_shifts
    hold on;
    plot(t, Kernels(i,:).*coefficients(i,1),"k:")
end
legend("original polynomials"," poly reproduced by dB4","shifted kernels shapes")
title("Degree-0 polynomial")
%PLOT degree=1
figure(2);
plot(t,polynomials(2,:),'b')
hold on;
plot(t,repoly(2,:),"r")
hold on;
for i=1:num_shifts
    hold on;
    plot(t, Kernels(i,:).*coefficients(i,2),"k:")
end
legend("original polynomials"," poly reproduced by dB4","shifted kernels shapes")
title("Degree-1 polynomial")
%PLOT degree=2
figure(3);
plot(t,polynomials(3,:),'b')
hold on;
plot(t,repoly(3,:),"r")
hold on;
for i=1:num_shifts
    hold on;
    plot(t, Kernels(i,:).*coefficients(i,3),"k:")
end
legend("original polynomials"," poly reproduced by dB4","shifted kernels shapes")
title("Degree-2 polynomial")
%PLOT degree=3
figure(4);
plot(t,polynomials(4,:),'b')
hold on;
plot(t,repoly(4,:),"r")
hold on;
for i=1:num_shifts
    hold on;
    plot(t, Kernels(i,:).*coefficients(i,4),"k:")
end
legend("original polynomials"," poly reproduced by dB4","shifted kernels shapes")
title("Degree-3 polynomial")


