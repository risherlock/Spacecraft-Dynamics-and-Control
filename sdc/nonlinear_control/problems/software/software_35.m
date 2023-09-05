% Week 3: Concept check 5, Question 1 and 3
% 2022.03.6

clc
clear
close all

% Spacecraft params
I1 = 100; % kg * m^2
I2 = 75;  % kg * m^2
I3 = 80;  % kg * m^2
I = diag([I1, I2, I3]);

% Control gain
K = 5; % N * m

P = diag([22.3607,19.3649,20.000]);
xi = get_damping_ratio(P, I, K);
T = get_time_decay(P, I);
fprintf("P_1:\n"); disp(P)
fprintf("xi_1: "); disp(xi')
fprintf("T_1: "); disp(T')
fprintf("*****\n");

P = diag([0.223607,0.258199,0.25]);
xi = get_damping_ratio(P, I, K);
T = get_time_decay(P, I);
fprintf("P_2:\n"); disp(P)
fprintf("xi_2: "); disp(xi')
fprintf("T_2: "); disp(T')
fprintf("*****\n");

P = diag([50.000,43.3013,44.7214]);
xi = get_damping_ratio(P, I, K);
T = get_time_decay(P, I);
fprintf("P_3:\n"); disp(P)
fprintf("xi_3: "); disp(xi')
fprintf("T_3: "); disp(T')


%%%%%%%%%%%%%%%%%%%
%    FUNCTIONS    %
%%%%%%%%%%%%%%%%%%%

function xi = get_damping_ratio(P, I, K)
xi = zeros(3,1);
for i = 1:3
    xi(i) = P(i,i) / sqrt(K * I(i,i));
end
end


function T = get_time_decay(P, I)
T = zeros(3,1);
for i = 1:3
    T(i) = 2 * I(i,i) / P(i,i);
end
end