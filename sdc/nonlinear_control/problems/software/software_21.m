% Week 2: Concept check 1 
% 2022.02.21

clc
clear
close all
format long

% Question 6
fprintf("        ~~~ Question 6 ~~~");
matrix = [1.53947, -0.0422688, -0.190629;...
    -0.0422688, 1.4759, 0.459006;...
    -0.190629, 0.459006, 1.48463];
eigen_values = eig(matrix)';
display(matrix);
display(eigen_values);

% Question 7
fprintf("        ~~~ Question 7 ~~~");
matrix = [-0.984331, -1.10006, -0.478579;...
    -1.10006, 1.03255, 0.338318;...
    -0.478579, 0.338318, 1.45178];
eigen_values = eig(matrix)';
display(matrix);
display(eigen_values);

% Question 8
fprintf("        ~~~ Question 8 ~~~");
matrix = [-2.0353, 0.296916, 0.365128;...
    0.296916, -1.10369, -0.074481;...
    -0.365128, -0.074481, -2.86101];
eigen_values = eig(matrix)';
display(matrix);
display(eigen_values);
