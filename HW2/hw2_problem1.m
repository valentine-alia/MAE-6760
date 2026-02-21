%% MAE 6760 Model Based Estimation
% Cornell University
% M Campbell
%
% Homework #2
% Problem #1: Least Squares Estimaiton
%   for a bridge circuit
%

%% define measurement output matrix: zk = H * x + vk
H=[1 -1;
    10 1;
   1 10];
xtrue=[2;10]; %true values
n=20; %number of measurements
rng(0); %gives repeatable random values

%% generate measurement data
sig1=10;
sig2=1;
sig3=100;
va=diag([sig1 sig2 sig3])*randn(3,n);
za=repmat(H*xtrue,[1,n])+va;

%% Part (a): LS estimate using noisy data z
s1=10;
s2=1;
s3=100;
va=diag([s1 s2 s3])*randn(3,n);
za=repmat(H*xtrue,[1,n])+va;


xhata = inv(n*(H'*H)) * H' * sum(za,2);




%% Part (b): use the same z as above, but now you know the noise
zb=za;

R = diag([s1^2 s2^2 s3^2]);
Rinv = inv(R);
xhatb = inv(n*H'*Rinv*H)*sum(H'*Rinv*zb,2);

d1=norm(xhata - xtrue)^2;
d2=norm(xhatb - xtrue)^2;


%% Part (c): Estimate LS estimator covariance

Pxhatb = inv(H' * Rinv * H)/n;

%% Part (d): Estimate covariance of x for (a) from data 

vaest = za - repmat(H*xhata,[1,n]);
Rhat = cov(vaest');
Rhatinv = inv(Rhat)
xhatd = inv(n*H'*Rhatinv*H)*sum(H'*Rhatinv*zb,2)
%Rtrueish = cov(va')



%% Part (e): LS estimator using Gaussian and uniform noise
Ue=173;
v3e=Ue*[[rand(1,n)-0.5]*2];
ve=[va(1:2,:);v3e];
ze=repmat(H*xtrue,[1,n])+ve;

