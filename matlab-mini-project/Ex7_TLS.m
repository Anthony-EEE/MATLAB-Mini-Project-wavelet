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
N=7;% N=5,6,7,8,9>2K=4
% maximum degree of polynomials
Max_d=N-1;

%dirac function with locations 1000 and 2000 and amplitudes 10, 20.
Dirac=zeros(1, 2048);
Dirac(1000)=150;
Dirac(2000)=200;


%second: sampling from x(t) to y[n] 
%1.Daubechies filters with dB5, where 5=N+1 vanishing moments;
Kernels=zeros(num_shifts,2048);
phi = zeros(1,2048);
dB = sprintf('dB%d',N);
[phi_T, psi_T, xval]=wavefun(dB, 6);
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
variance =1e3;
% 
%Gaussian noise
Gaussion_noise = sqrt(variance).*randn(1, N);
sm = tau + Gaussion_noise;
%sm = tau;
%because we use dB5,6,7, so S is 
S=ones(N-K,K+1); 
for i=1:(N-K)
    S(i,:)=[sm(K+1+(i-1)),sm(K+(i-1)),sm(K-1+(i-1))];
end
S1=S;


%Total least-squares approach(TLS)
%SVD to find h
[U, lamda, V]=svd(S);
%h=V(:,end);
h=V(2:end,end)/V(1,end);


%z-transform to get H(z) and to find t_k

H_z=1+h(1)*z^(-1)+h(2)*z^(-2);
t_k=solve(H_z);

%%3. Using Vandermonde system to find weight a_k
%[1,1;t_0,t_1][a_0;a_1]=[tau0;tau1]
AA = [1,1; t_k(1), t_k(2)];
BB = [sm(1);sm(2)];
a_k= mldivide(AA,BB);
display(round(t_k, 8));
t_k=t_k.*T;

display(round(t_k));
display(round(a_k, 2));

figure(1);
scatter(t_k,a_k)

