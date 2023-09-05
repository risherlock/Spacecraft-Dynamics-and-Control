% Fig. 8.3: Analytical Mechanics of Space Systems (4th ed.) 
% Rishav (2023.02.14)

clc
clear
close all

% X-axis vector
x = -20:0.001:20;

% a
subplot(3,2,1);
y = x.^2;
plot(x,y,'.'); hold on;
y = y+100; plot(x,y,'.');
y = y-200; plot(x,y,'.'); 
title("x^{2}+c");
xlabel("x"); ylabel("V(x)");
legend("c>0", "c=0", "c<0")

% b
subplot(3,2,2);
y = x.^2.*sin(x/2).^2;
plot(x,y,'.');
title("x^{2}+sin^2(0.5x)");
xlabel("x"); ylabel("V(x)");

% c
subplot(3,2,3);
y = (sin(x)+x).^2;
plot(x,y,'.');
title("(sin(x)+x)^2")
xlabel("x"); ylabel("V(x)");

% d
subplot(3,2,4);
y = x.^2-x.^4/256;
plot(x,y,'.');
title("x^{2}-x^4/256")
xlabel("x"); ylabel("V(x)");

% e
subplot(3,2,5);
y = x.^2.*exp(-x.^2/50);
plot(x,y,'.');
title("x^{2} exp(-x^{2}/50)")
xlabel("x"); ylabel("V(x)");

% f
subplot(3,2,6);
y = x.^3;
plot(x,y,'.');
title("x^{3}")
xlabel("x"); ylabel("V(x)");

