%% Respond checking for open loop system, controllability checking and observability checking

clc
close all
clear all

t = 0:0.001:20; 
dt = t(2) - t(1);

g = 9.8;
% 4 motor, 1 batre, 1 body frame, 1 hub
m = ((4*0.288) + 0.356 + 0.701 + 0.293);
ixx = 0.1328;
iyy = 0.1324;
izz = 0.2638;

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

B= [0 0 0 0;
    0 0 0 0;
    0 0 0 0; 
    0 1/ixx 0 0;
    0 0 1/iyy 0;
    0 0 0 1/izz;
    0 0 0 0;
    0 0 0 0;
    1/m 0 0 0;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0];

C =[1 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0; 
    0 0 1 0 0 0 0 0 0 0 0 0; 
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 1];

G = [0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 1/ixx 0 0;
     0 0 0 0 1/iyy 0;
     0 0 0 0 0 1/izz;
     1/m 0 0 0 0 0;
     0 1/m 0 0 0 0;
     0 0 1/m 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0];
 
D = zeros(6,4);
H=zeros(6,6);

% Check Controllability
ctrb_matrix = ctrb(A,B);
[baris, kolom]=size(ctrb_matrix);
if rank(ctrb_matrix) == baris
    disp('Controllable') 
end
% Check Observability
obsv_matrix = obsv(A,C);
[baris, kolom]=size(obsv_matrix);
if rank(obsv_matrix) == kolom
    disp('Observable') 
end

% Block of code for getting lqr gain and for simulation
N = 0;
Q = 100;
R = 1;
[K,S,e] = lqr(A,B,Q,R,N);
K_lqr=K;

X0(:,1)=[0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]';
X1(:,1)=[0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]';

for i = 2:length(t)
      u= .5*ones(4,1);
      u1(:,i)=-K*X1(:,i-1);
    
      X0(:,i) = X0(:,i-1)  +dt * (A*X0(:,i-1) + B*u);
      y0(:,i) = C*X0(:,i);
    
    X1(:,i) = X1(:,i-1)  +dt * ((A - B*K_lqr)*X1(:,i-1));
    y1(:,i) = C*X1(:,i);
    
end

% Block of code for getting kalman gain and for simulation
Q_k=eye(6,6);
R_k=eye(6,6);

[kest, L , P] = kalman(ss(A,[B G], C, [D H]),Q_k,R_k,0);



X_hat(:,1) = zeros(1,12);
y_hat(:,1) = C*X_hat;


X(:,1)=zeros(1,12);

%for i = 2:length(t)
 %   u= .5*ones(4,1);
    
%    X(:,i) = X(:,i-1)  +dt * (A*X(:,i-1) + B*u);
%    y(:,i) = C*X(:,i) + sqrt(R)*randn(size(C,1),1);

%    X_hat(:,i) = X_hat(:,i-1)  +dt * (A*X_hat(:,i-1) + B*u +L *(y(:,i-1)-y_hat(:,i-1)));
%    y_hat(:,i) = C*X_hat(:,i) ;
%end

% Block of code for simulate LQG
X(:,1)=[0 0 0 0 0 0 0 0 0 0 0 0]';
X_hat(:,1)=[1 1 1 1 1 1 1 1 1 1 1 1]';
y_hat(:,1) = C*X_hat;
u(:, i) = zeros(4,1);

for i = 2:length(t)
    
     u(:, i)=(-K_lqr*X_hat(:, i-1));
    
    X(:,i) = X(:,i-1)  + dt * (A*X(:,i-1) + B*(-K_lqr*X_hat(:, i-1))+ G*wgn(size(C,1),1,0));
    y(:,i) = C*X(:,i) + sqrt(R_k)*wgn(size(C,1),1,0);
     

    X_hat(:,i) = X_hat(:,i-1)  +dt*(A*X_hat(:,i-1) + B*(-K_lqr*X_hat(:, i-1)) +L *(y(:,i-1)-y_hat(:,i-1)));
    y_hat(:,i) = C*X_hat(:,i) ;
end

% Plot Control Signal

% figure(1);
% subplot(2,2,1)
% plot(t,u(1,:))
% xlabel('time')
% ylabel('ft (Newton)')
% subplot(2,2,2)
% plot(t,u(2,:))
% xlabel('time')
% ylabel('Tau_x (Newton)')
% subplot(2,2,3)
% plot(t,u(3,:))
% xlabel('time')
% ylabel('Tau_y (Newton)')
% subplot(2,2,4)
% plot(t,u(4,:))
% xlabel('time')
% ylabel('Tau_z (Newton)')

pole_lqr = eig(A-B*K)  
pole_lqe = eig(A-L*C) 
pole_lqg = [eig(A-B*K);eig(A-L*C)]  


% Plot State
figure(1);
plot(t,X(1,:),t,X_hat(1,:))
legend('actual','estimated')
xlabel('time')
ylabel('Roll (rad)')

figure(2);
plot(t,X(2,:),t,X_hat(2,:))
legend('actual','estimated')
xlabel('time')
ylabel('Pitch (rad)')

figure(3);
plot(t,X(3,:),t,X_hat(3,:))
legend('actual','estimated')
xlabel('time')
ylabel('yaw (rad)')

figure(4);
plot(t,X(1,:),t,X_hat(1,:))
legend('actual','estimated')
xlabel('time')
ylabel('linear x (m)')

figure(5);
plot(t,X(2,:),t,X_hat(2,:))
legend('actual','estimated')
xlabel('time')
ylabel('linear y (m)')

figure(6);
plot(t,X(3,:),t,X_hat(3,:))
legend('actual','estimated')
xlabel('time')
ylabel('linear z (m)')
