% Week 4: Concept check 2, Question 2
% 2022.03.9

clc
clear
close all
format long

% Spacecraft params
I1 = 100; % kg * m^2
I2 = 75;  % kg * m^2
I3 = 80;  % kg * m^2

% Initial conditions
sigma_o = [0.1, 0.2, -0.1]'; % MRP
omega_o = [3, 1, -2]';       % deg/sec
omega_o = deg2rad(omega_o);

% Control gains and limits
K = 0.11;       % Nm
P = 3 * eye(3); % Nm * s

% Simulation time params, s
dt = 0.01;
stop_time = 120;
time = 0:dt:stop_time;

% Allocate memotry
state = zeros(6,length(time));
state(:,1) = [sigma_o; omega_o];

for i = 1:length(time)-1
    state_dot = linear_mrp_cld(P, K, state(:,i));
    state(:,i+1) = state(:,i) + state_dot * dt;
end
plot(time, state(1,:)); hold on;
plot(time, state(2,:));
plot(time, state(3,:));
legend("\sigma_1","\sigma_2","\sigma_3");
title("MRP response for linear closed loop system");
grid on; xlabel("Time (s)"); ylabel("MRPs");

%%%%%%%%%%%%%%%%%%%
%    FUNCTIONS    %
%%%%%%%%%%%%%%%%%%%

% MRP closed loop response
function [state_dot] = linear_mrp_cld(P, K, state)
sigma = state(1:3);
sigma_dot = state(4:6);
sigma_ddot = - (P * sigma_dot + K * sigma);

state_dot = [sigma_dot; sigma_ddot];
end
