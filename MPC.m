%% MPC Controller Simulation for Quadrotor Model

clc
close all
clear all

%define parameters value
g = 9.8;
% 4 motor, 1 batre, 1 body frame, 1 hub
m = ((4*0.288) + 0.356 + 0.701 + 0.293);
ixx = 0.1328;
iyy = 0.1324;
izz = 0.2638;

%State Space Model
A= [0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 -g 0 0 0 0 0 0 0 0 0 0;
    g 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0];

B= [0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 1/ixx 0 0 0 0 0 1/ixx 0 0;
    0 0 1/iyy 0 0 0 0 0 1/iyy 0;
    0 0 0 1/izz 0 0 0 0 0 1/izz;
    0 0 0 0 1/m 0 0 0 0 0;
    0 0 0 0 0 1/m 0 0 0 0;
    1/m 0 0 0 0 0 1/m 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0];

C =[1 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 1];

D = zeros(size(C,1),size(B,2));

sys = ss(A,B,C,D);
sys = setmpcsignals(sys,'MV',[1 2 3 4],'UD',[5 6 7 8 9 10]);
sys.InputName = {'f_t','\tau_x','\tau_y','\tau_z','f_w_x','f_w_y','f_w_z','\tau_w_x','\tau_w_y','\tau_w_z'};
sys.OutputName = {'\gamma','\Theta','\Psi','x','y','z'};
sys.StateName = {'\gamma','\Theta','\Psi','\gamma''','\Theta''','\Psi''','x''','y''','z''','x','y','z'};

%MPC
%Create MPCobj
Ts = 0.01;   %sampling time
p = 10;     %prediction horizon
m = 2;      %control horizon
T = 50;
num_sim_steps = round(T/Ts);


MPCobj = mpc(sys,Ts,p,m,[]);
simOptions = mpcsimopt(MPCobj);
simOptions.PlantInitialState = [0; 0.05; 0.025; 0.01; 0.01; 0.01; 0.1; 0.1; 0.1; 2; 3; 4];

%Define weights
MPCobj.W.Output = [100, 100, 100, 10, 10, 10];
MPCobj.W.Input = [1 1 1 1];
MPCobj.W.InputRate = [0.1 0.1 0.1 0.1];

%Define constraints
MPCobj.MV = struct('Min',{-10;-10;-10;-10},'Max',{10;10;10;10},'RateMin',{-1;-1;-1;-1},'RateMax',{1;1;1;1});
MPCobj.OV = struct('Min',{-0.1;-0.1;-0.1;-1;-1;-1},'Max',{0.1;0.1;0.1;1;1;1});

%define disturbance
mod1 = tf(1,1);
indist = mod1*eye(6);
setindist(MPCobj,'model',indist);
simOptions.UnmeasuredDisturbance = wgn(num_sim_steps,6,-12);

%define measurement noise model
mod1 = tf(1,1);
noise_model = mod1*eye(6);
MPCobj.Model.Noise = noise_model;
simOptions.OutputNoise = wgn(num_sim_steps,6,-12);

%simulasi
T = 50;
num_sim_steps = round(T/Ts);
r = [0 0 0 0 0 0];  %reference
sim(MPCobj, num_sim_steps, r, simOptions);
