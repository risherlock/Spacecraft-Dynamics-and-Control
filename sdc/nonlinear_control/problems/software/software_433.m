% Week 4: Concept check 3, Question 3
% 2022.03.9

clc
clear
close all

% 4 RW configurations
gs1 = [0.267261, 0.534522, 0.801784]';
gs2 = [-0.267261, 0.534522, 0.801784]';
gs3 = [0.534522, 0.267261, 0.801784]';
gs4 = [-0.666667, 0.666667, 0.333333]';

% Required control
Lr = [0.1, 0.20, 0.4]';

% Minimum motor torques for Lr
Gs = [gs1, gs2, gs3, gs4];
us = lsqminnorm(Gs, Lr);

% Output result
fprintf("Motor torques:\n");
disp(us);
