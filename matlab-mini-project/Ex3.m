clc;clear;
syms z;
%%tau=[tau0,tau1,tau2,tau3]
load('tau.mat');
K=2;% K diracs


% 1st step: N+1 moments which is len(tau)
% can find polynomials of maximum degree
N=length(tau)-1;


% 2nd Step: find location t_k by Yule-Walker and annihilating filter
% Yule-Walker [tau1, tau0;tau2, tau1][h1;h2]=-[tau2;tau3]
A = [tau(2),tau(1); tau(3), tau(2)];
B = [-tau(3); -tau(4)];%tau=[tau0,tau1,tau2,tau3]

% compute [h1;h2]
h=A\B;

%z-transform to get H(z) and to find t_k
h_0=1;
H_z=h_0+h(1)*z^(-1)+h(2)*z^(-2);
t_k=solve(H_z);


%%3rd Step: Using Vandermonde system to find weight a_k
%[1,1;t_0,t_1][a_0;a_1]=[tau0;tau1]
AA = [1,1; t_k(1), t_k(2)];
BB = [tau(1);tau(2)];
a_k= AA\BB;

display(t_k);
display(a_k);
