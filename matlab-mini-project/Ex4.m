clc;clear;
syms z;
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


%%first: create a stream of Diracs with K=2;
K=2;
%so the maximum polynomial degree is 3, 
% only dB4 with vanishing moment length(tau)=4 can vanish the maximum polynomial
% degree N = 4-1;
N=2*K-1;

%dirac function with locations and amplitudes.
Dirac=zeros(1, 2048);
x_0=1000;
t_0=x_0/T;

x_1=2000;
t_1=x_1/T;

display([t_0; t_1]);
a_0=10;
a_1=20;
Dirac(x_0)=a_0;
Dirac(x_1)=a_1;


%second: sampling from x(t) to y[n] 
%1.Daubechies filters with dB4, where 4=N+1 vanishing moments;
Kernels=zeros(num_shifts,2048);
phi = zeros(1,2048);
[phi_T, psi_T, xval]=wavefun('dB4', 6);
phi(1:length(phi_T))=phi_T;
for i=1:num_shifts   
    Kernels(i,:)=[zeros(1,(i-1)*T), phi(1: end-(i-1)*T)];
end
% sampling using formula y[n]=<x(t), phi(t-n)>=phi'*x(t) row vector * column vector;
y_n=Kernels*Dirac.';
y_n=y_n';



% Third: Using 3 steps to reconstruct a_k and t_k 
% 1.generate tau with formula tau=sum(cmn * y_n);
% so we need coefficient cmn at first. it can be drivated by cmn=<t^m,phi(t-n)>
% polynomails t^m
polynomials=zeros(Max_d,2048);
for degree=0:Max_d
    polynomials(degree+1,:)=t.^(degree);
end
coefficients=zeros(Max_d, num_shifts);
coefficients=Kernels*polynomials.'/T;
tau=y_n*coefficients;

%2. find location t_k by Yule-Walker and annihilating filter
A = [tau(2),tau(1); tau(3), tau(2)];
B = [-tau(3); -tau(4)];%tau=[tau0,tau1,tau2,tau3]
% compute [h1;h2]
h=A\B;
%z-transform to get H(z) and to find t_k
h_0=1;
H_z=h_0+h(1)*z^(-1)+h(2)*z^(-2);
t_k=solve(H_z);

%%3. Using Vandermonde system to find weight a_k
%[1,1;t_0,t_1][a_0;a_1]=[tau0;tau1]
AA = [1,1; t_k(1), t_k(2)];
BB = [tau(1);tau(2)];
a_k= mldivide(AA,BB);
display(round(t_k, 5));
t_k=t_k.*T;
display(round(t_k));
display(round(a_k, 2));