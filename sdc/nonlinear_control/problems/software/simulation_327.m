% Week 3: Concept check 2, Question 6
% 2022.03.4

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
K = 5;                 % Nm
P = 10 * eye(3);       % Nm * s

% Simulation time params, s
dt = 0.01;
stop_time = 120;
time = 0:dt:stop_time;
display_norm = 70;

% MRP reference trajectory
f = 0.05; % rad/sec
sigmas_rn = zeros(3, length(time));
sigmas_rn(1,:) = 0.2 * sin(f * time);
sigmas_rn(2,:) = 0.3 * cos(f * time);
sigmas_rn(3,:) = -0.3 * sin(f * time);

% Derivative of reference trajectory
sigmas_dot_rn = zeros(3, length(time));
sigmas_dot_rn(1,:) = 0.2 * f * cos(f * time);
sigmas_dot_rn(2,:) = -0.3 * f * sin(f * time);
sigmas_dot_rn(3,:) = -0.3 * f * cos(f * time);

% Rate and rate of rate trajectory
omegas_rn = zeros(3, length(time));
omegas_dot_rn = zeros(3, length(time));

for i = 1:length(time)
    omegas_rn(:,i) = 4 * (sigma_matrix(sigmas_rn(:,i)) \ sigmas_dot_rn(:,i));
    
    if(i~=1)
        omegas_dot_rn(:,i) = (omegas_rn(:,i) - omegas_rn(:,i-1)) / dt;
    end
end

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
    
    % Unpack reference states
    sigma_rn = sigmas_rn(:,i);
    omega_rn = omegas_rn(:,i);
    omega_dot_rn = omegas_dot_rn(:,i);
    
    % Compute state differentials
    u = control(I, omega, omega_rn, omega_dot_rn, sigma, sigma_rn, K, P);
    omega_dot = kinetics(I, omega, u);
    sigma_dot = kinematics(sigma, omega);
    
    % Numerical intergration
    omegas(:,i+1) = omega + omega_dot * dt;
    sigmas(:,i+1) = sigma + sigma_dot * dt;
    sigmas(:,i+1) = switch_mrp(sigmas(:,i+1));
end

dcm_bn = mrp_to_dcm(sigmas(:,display_norm/dt));
dcm_nr = mrp_to_dcm(sigmas_rn(:,display_norm/dt))';
sigma_br = dcm_to_mrp(dcm_bn * dcm_nr);
fprintf("MRP tracking error norm at %d secs: %f\n",...
    display_norm, norm(sigma_br))

plot(time, sigmas(1,:), time, sigmas(2,:), time, sigmas(3,:));
legend('\sigma_1', '\sigma_2', '\sigma_3');
grid on;


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
function [u] = control(I, omega_bn, omega_rn, omega_dot_rn, ...
    sigma_bn, sigma_rn, K, P)

% Error MRP as seen from body frame: sigma_br
dcm_bn = mrp_to_dcm(sigma_bn);
dcm_nr = mrp_to_dcm(sigma_rn)';
sigma_br = dcm_to_mrp(dcm_bn * dcm_nr);

% Transform desired references to body frame
dcm_br = dcm_bn * dcm_nr;
omega_rn = dcm_br * omega_rn;  
omega_dot_rn = dcm_br * omega_dot_rn;

% Error omega as seen from body frame: omega_br
omega_br = omega_bn - omega_rn;

% Control law
u = -K * sigma_br - P * omega_br + cross(omega_bn, I * omega_bn) ...
    + I * (omega_dot_rn - cross(omega_bn, omega_rn));
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

% Rotation matrix to MRP
function [sigma] = dcm_to_mrp(dcm)
sigma = zeros(3,1);
sigma(1) = dcm(2,3) - dcm(3,2);
sigma(2) = dcm(3,1) - dcm(1,3);
sigma(3) = dcm(1,2) - dcm(2,1);

zeta = sqrt(trace(dcm) + 1);
sigma = sigma ./ (zeta^2 + 2 * zeta);
end

% MRP to rotation matrix
function [dcm] = mrp_to_dcm(mrp)
sq = mrp' * mrp;
M = [0, -mrp(3), mrp(2); mrp(3), 0, -mrp(1); -mrp(2), mrp(1), 0];
dcm = eye(3) + ((4 * (1 - sq) * eye(3) + 8 * M) * M) / (1 + sq)^2;
dcm = dcm';
end
