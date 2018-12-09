% Solutions to the diffusion equation with initial condition x^2 and
% Non Homogeneous Dirichlet boundary condition
function HeatEquationZeroBCNonH
clear;
close all;

% Initial variables.
D=1; % Diffusion coefficient.
x_min=0; % Min domain
x_max=1; % Max domain
t_max=1; % Max time
xv=linspace(x_min,x_max,20); % Spatial vector 
dx=xv(2)-xv(1); % x step size
dt=0.5*(dx^2)/(2*D); % t stable step size
tv=0:dt:t_max; % Time vector

% Fourier series solution.
a0=xv;
% Series for the domain [0,1]. Requires recalculation should domain change.
series=@(x,t,n)2*((pi*n*sin(pi*n)+2*cos(pi*n)-2)/(pi^3*n^3))*sin(n*pi*x)*exp(-(n*pi)^2*t);

% FDM Solution
ic=@(x)x.^2;
fdmU=FDM(1,dx,dt,[x_min x_max],[0 t_max],[0 1],ic);
fdmU=fdmU';

% pdepe Solution
sol=pdepe(0,@pdex1pde,@pdex1ic,@pdex1bc,xv,tv);
% Extract the first solution component as u.
pdepeU=sol(:,:,1);

% Plotting 3D Graphs
figure
subplot(2,2,1)
surf(xv,tv,fdmU)
title('Numerical solution FDM.')
xlabel('Distance x')
ylabel('Time t')
shading interp

subplot(2,2,2)
surf(xv,tv,pdepeU) 
title('Numerical solution PDEPE.')
xlabel('Distance x')
ylabel('Time t')
shading interp

subplot(2,2,3)
imagesc(xv,tv,fdmU) 
title('Numerical solution FDM.')
xlabel('Distance x')
ylabel('Time t')
shading interp
colorbar

subplot(2,2,4)
imagesc(xv,tv,pdepeU) 
title('Numerical solution PDEPE.')
xlabel('Distance x')
ylabel('Time t')
shading interp
colorbar

% Plotting 2D solutions at fixed time t
figure
suptitle('Solutions at different times')
p1=subplot(2,2,1);
plot(xv,a0+seriesSum(series,xv,0,5),'-',xv,pdepeU(1,:),'-.',xv,fdmU(1,:),'*')
title('Solution at t = 0.0')
xlabel('Distance x')
ylabel('u(x,0)')
legend('Fourier series first three non zero terms','pdepe() solution','Finite Difference method')

p2=subplot(2,2,2);
plot(xv,a0+seriesSum(series,xv,0.1,5),xv,pdepeU(round(length(pdepeU)/10)*1,:),'-.',xv,fdmU(round(length(fdmU)/10)*1,:),'*')
title('Solution at t = 0.1')
xlabel('Distance x')
ylabel('u(x,0.1)')
legend('Fourier series first three non zero terms','pdepe() solution','Finite Difference method')

p3=subplot(2,2,3);
plot(xv,a0+seriesSum(series,xv,0.2,5),xv,pdepeU(round(length(pdepeU)/10)*2,:),'-.',xv,fdmU(round(length(fdmU)/10)*2,:),'*')
title('Solution at t = 0.2')
xlabel('Distance x')
ylabel('u(x,0.2)')
legend('Fourier series first three non zero terms','pdepe() solution','Finite Difference method')

p4=subplot(2,2,4);
plot(xv,a0+seriesSum(series,xv,0.3,5),xv,pdepeU(round(length(pdepeU)/10)*3,:),'-.',xv,fdmU(round(length(fdmU)/10)*3,:),'*')
title('Solution at t = 0.3')
xlabel('Distance x')
ylabel('u(x,0.3)')
legend('Fourier series first three non zero terms','pdepe() solution','Finite Difference method')
linkaxes([p1,p2,p3,p4], 'xy');

% Error Figures
figure
loglog(xv,abs((a0+seriesSum(series,xv,0,5))-fdmU(1,:)),xv,abs((a0+seriesSum(series,xv,0,5))-pdepeU(2,:)));
title('FDM vs pdepe() Error Comparison at t = 0')
xlabel('Domain x')
ylabel('Error')
legend('FDM Error', 'pdepe() Error')
xlim([x_min x_max])

% Pdepe function handlers
function [c,f,s]=pdex1pde(x,t,u,DuDx)
c=1;
f=DuDx;
s=0;

function u0=pdex1ic(x)
u0=x^2;

function [pl,ql,pr,qr]=pdex1bc(xl,ul,xr,ur,t)
pl=0;
ql=ul;
pr=0;
qr=ur-1;