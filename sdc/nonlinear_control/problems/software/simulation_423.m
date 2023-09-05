% Week 4: Concept check 2, Question 3
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
omega_o = [30, 10, -20]';    % deg/sec

% Control gains and limits
K = 0.11;       % Nm
P = 3 * eye(3); % Nm * s

% Simulation time params, s
dt = 0.01;
stop_time = 150;
time = 0:dt:stop_time;
display_norm = 50;

% Precondition parameters
I = diag([I1, I2, I3]);
omega_o = deg2rad(omega_o);

% Allocate memory for states
sigmas = zeros(3, length(time));
omegas = zeros(3, length(time));
sigmas(:,1) = sigma_o;
omegas(:,1) = omega_o;

% Simulation loop
for i = 1:length(time)-1
    % Unpack states wrt. inertial
    sigma = sigmas(:,i);
    omega = omegas(:,i);
    
    % Compute state differentials
    u = control(I, omega, sigma, K, P);
    omega_dot = kinetics(I, omega, u);
    sigma_dot = kinematics(sigma, omega);
    
    % Numerical intergration
    omegas(:,i+1) = omega + omega_dot * dt;
    sigmas(:,i+1) = sigma + sigma_dot * dt;
    sigmas(:,i+1) = switch_mrp(sigmas(:,i+1));
end

fprintf("MRP tracking error norm at %d secs: %f\n",...
    display_norm, norm(sigmas(:,display_norm/dt)))

plot(time, sigmas(1,:)); hold on;
plot(time, sigmas(2,:));
plot(time, sigmas(3,:));
legend("\sigma_1","\sigma_2","\sigma_3");
title("MRP response for linear closed loop system");
grid on; xlabel("Time (s)"); ylabel("MRPs");


%%%%%%%%%%%%%%%%%%%
%    FUNCTIONS    %
%%%%%%%%%%%%%%%%%%%

% MRP differential equation
function [sigma_dot] = kinematics(sigma, omega)
M = sigma_matrix(sigma);
sigma_dot = 0.25 * M * omega;
end

% Rotational dynamics
function [omega_dot] = kinetics(I, omega, u)
omega_dot = I\(- cross(omega, I * omega) + u);
end

% Reference tracking control law
function [u] = control(I, omega, sigma, K, P)

% Compute omega_dot
temp = (4.0 * K) / (1 + sigma' * sigma) - (omega' * omega) / 2.0;
omega_dot = - P * omega - (omega * omega' + temp * eye(3)) * sigma;

% Plug omega_dot to kinetics
u = I * omega_dot + cross(omega, I * omega);
end

% Switch to MRP shadow
function [sigma] = switch_mrp(sigma)
sq = sigma' * sigma;
if(sq > 1)
    sigma = [-sigma(1)/sq, -sigma(2)/sq, -sigma(3)/sq]';
end
end

% MRP matrix
function [M] = sigma_matrix(mrp)
skew = [0, -mrp(3), mrp(2); mrp(3), 0, -mrp(1); -mrp(2), mrp(1), 0];
M = (1 - mrp' * mrp) * eye(3) + 2 * skew + 2 * (mrp * mrp');
end

