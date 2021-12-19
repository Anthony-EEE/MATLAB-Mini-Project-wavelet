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
N=6;% N=5,6,7,8,9>2K=4
% maximum degree of polynomials
Max_d=N-1;

%dirac function with locations 1000 and 2000 and amplitudes 10, 20.
Dirac=zeros(1, 2048);
Dirac(1000)=10;
Dirac(2000)=20;


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
% variance =1e9;
% 
% %Gaussian noise
% Gaussion_noise = sqrt(variance).*randn(1, N);
% sm = tau + Gaussion_noise;
sm = tau;
%because we use dB5,6,7, so S is 
S=ones(N-K,K+1); 
for i=1:(N-K)
    S(i,:)=[sm(K+1+(i-1)),sm(K+(i-1)),sm(K-1+(i-1))];
end
S1=S;

%only using Total least-squares approach(TLS)
%SVD to find h
[U, lamda, V]=svd(S);


% [U, lamda, V]=svd(S);
% display(size(V))
% lamda_new=zeros(size(lamda));
% %Using Cadzow to update S
% for i=1:40
% 
%     %keep the k largest diagonal coefficients of lamda to define new lamda
%     lamda_new(1:2,1:2)=lamda(1:2,1:2);
%     %reconstruct S
%     S=U*lamda_new*V';
%     %make S be the form of Toeplitz   
%     %average diagonals and make it toeplitz
%     count =1;
%     for j = (-(N-K-1)) : 1 : (K)
%         M = mean(diag(S,j));
%         if j<0
%             col(count) = M; 
%             count = count+1;
%         elseif j==0
%             row(j+1) = M; 
%             col(count) = M;
%         else
%             row(j+1)=M; 
%         end
%     end
% 
%     col = flip(col); %flip th col vector but not flip the col vector
%     S=toeplitz(col,row);
%     S=S(1:(N-K),1:3);
%     display(S)
%     %TLS
%     [U, lamda, V]=svd(S);
% end

%SVD to find h
h=V(:,end);


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

