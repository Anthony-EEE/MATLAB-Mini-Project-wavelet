clc;clear;
syms z;

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


%%first: create a stream of Diracs with K=2;
K=2;
% only dB5 with vanishing moment 5 can vanish the maximum polynomial degree 4
N=2*K+1;% N=5>2K=4
% maximum degree of polynomials
Max_d=N-1;

%dirac function with locations and amplitudes.
Dirac=zeros(1, 2048);
Dirac(1000)=10;
Dirac(2000)=20;


%second: sampling from x(t) to y[n] 
%1.Daubechies filters with dB5, where 5=N+1 vanishing moments;
Kernels=zeros(num_shifts,2048);
phi = zeros(1,2048);
[phi_T, psi_T, xval]=wavefun('dB5', 6);
phi(1:length(phi_T))=phi_T;
for i=1:num_shifts   
    Kernels(i,:)=[zeros(1,(i-1)*T), phi(1: end-(i-1)*T)];
end
% sampling using formula y[n]=<x(t), phi(t-n)>=phi'*x(t) row vector * column vector;
y_n=Kernels*Dirac.';
y_n=y_n';


% Third: Using 3 steps to reconstruct a_k and t_k 
% 1.generate moments tau with formula tau=sum(cmn * y_n); 
% in this exercise tau =sm
% so we need coefficient cmn at first. it can be drivated by cmn=<t^m,phi(t-n)>
% polynomails t^m
polynomials=zeros(Max_d,2048);
for degree=0:Max_d
    polynomials(degree+1,:)=t.^(degree);
end
coefficients=zeros(Max_d, num_shifts);
coefficients=Kernels*polynomials'/T;
tau=y_n*coefficients;
display(tau);
% apply Gaussian noise with different variance to tau
variance =[1, 1e3, 1e6, 1e8, 1e12];
SM=[];
for i=1:length(variance)
     %Gaussian noise
    Gaussion_noise = sqrt(variance(i)).*randn(1,5);
    sm = tau + Gaussion_noise;
    SM=[SM;sm];
    display(sm);
end
