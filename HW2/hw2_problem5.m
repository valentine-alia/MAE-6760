%% MAE 6760 Model Based Estimation
% Cornell University
% M Campbell
%
% Homework #2
% Problem #5: Nonlinear Least Squares
%   for multi-sensor range localization problem
%
rng(10);

%% Problem set up with beacons and msmt generation
%Beacon locations
bA = [-10;100];
bB = [490;20];
bC = [500;40];
beacons=[bA bB bC];
%True location of the car
xtrue = [-5;2];

%Perfect measurements
RA = sqrt((bA(1) - xtrue(1))^2 + (bA(2) - xtrue(2))^2);
RB = sqrt((bB(1) - xtrue(1))^2 + (bB(2) - xtrue(2))^2);
RC = sqrt((bC(1) - xtrue(1))^2 + (bC(2) - xtrue(2))^2);

%a total of 10 measurements
n=10;


%% Part (a) - two ranging beacons (A,B,C)
%generate noisy measurements to beacons A,B,C
Rpart_a = diag([10 10 10]);
v_a = sqrtm(Rpart_a)*randn(3,n);
z_a = repmat([RA;RB;RC],[1,n]) + v_a;


%% Part (b) - two ranging beacons (A,B)
%Get noisy measurements to beacons A,B from part (a)
%simply so that the same noise is used
ii_beacons=[1 2];
Rpart_b1 = Rpart_a(ii_beacons,ii_beacons); 
z_b1 = z_a(ii_beacons,:);


%% Part (b) - two ranging beacons (B,C)
%Get noisy measurements to beacons B,C from part (a)
%simply so that the same noise is used
ii_beacons=[2 3];
Rpart_b2 = Rpart_a(ii_beacons,ii_beacons); 
z_b2 = za(ii_beacons,:);
x0=[0;0];


%% Part (c) - perfect linearization

%% Part (d) - correlated noise model
Rpart_d = [10 0 0; 0 10 9; 0 9 10];
v_d = sqrtm(Rpartd)*randn(3,n);
z_d = repmat([RA;RB;RC],[1,n]) + v_d;
