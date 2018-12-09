function DiffusionEquationNonLinear

clear;
close all;

% Setup constants and vectors.
D=1;
dx=0.02;
x_vec=0:dx:1;
dt=0.5*(dx^2)/(2*D);
t_vec=0:dt:1;

% Solve the system with FDM
U=FDMInitial(D,dx,dt,[0 1],[0 1],[0 0],0,[85 0.9]);
U=U';

% 3D Plot of the numerical solution.
subplot(1,2,1)
surf(x_vec,t_vec,U);
title('Numerical Solution');
xlabel('Space');
ylabel('Time');
shading interp
colorbar

% Plot the numerical solution.
subplot(1,2,2)
imagesc(t_vec,x_vec,U);
colorbar;
title('Numerical Solution');
xlabel('Space');
ylabel('Time');

% Plot the defined steady state solution.
figure
hold on
x1 = 0:0.01:0.4;
x2 = 0.4:0.01:0.6;
x3 = 0.6:0.01:1;
y1 = 0.1*x1;
y2 = -(x2.^2/2)+0.5.*x2-0.08;
y3 = -0.1.*x3+0.1;
plot(x_vec,U(end,:),x1,y1,'r',x2,y2,'r',x3,y3,'r');
title('Solution at t = 1');
xlabel('x');
ylabel('U(x)');
legend('Numerical Solution', 'Exact')
