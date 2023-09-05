% Week 3: Concept check 1, Question 4
% 2022.02.25

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

% Control gains
K = 5;           % Nm
P = 10 * eye(3); % Nm * s

% Simulation time params, s
dt = 0.01;
stop_time = 120;
time = 0:dt:stop_time;
display_norm = 30;

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
    % Unpack states
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

fprintf("MRP norm at %d secs: %f\n",...
    display_norm,  norm(sigmas(:,display_norm/dt)))

plot(time, sigmas(1,:), time, sigmas(2,:), time, sigmas(3,:));
legend('\sigma_1', '\sigma_2', '\sigma_3');
grid on;

%%%%%%%%%%%%%%%%%%%
%    FUNCTIONS    %
%%%%%%%%%%%%%%%%%%%

% MRP differential equation
function [sigma_dot] = kinematics(mrp, omega)
skew = [0, -mrp(3), mrp(2); mrp(3), 0, -mrp(1); -mrp(2), mrp(1), 0];
M = (1 - mrp' * mrp) * eye(3) + 2 * skew + 2 * (mrp * mrp');
sigma_dot = 0.25 * M * omega;
end

% Rotational dynamics
function [omega_dot] = kinetics(I, omega, u)
omega_dot = I\(- cross(omega, I * omega) + u);
end

% Regulaton control law
function [u] = control(I, omega, sigma, K, P)
u = - K * sigma - P * omega + cross(omega, I * omega);
end

% Switch to MRP shadow
function [sigma] = switch_mrp(sigma)
sq = sigma' * sigma;
if(sq > 1)
    sigma = [-sigma(1)/sq, -sigma(2)/sq, -sigma(3)/sq]';
end
end
