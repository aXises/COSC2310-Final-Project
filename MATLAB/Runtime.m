% Comparison of Runtimes of pdepe() and FDM
function Runtime
clear;
close all;

% Initial variables.
D=1; % Diffusion coefficient.
x_min=0; % Min domain
x_max=1; % Max domain
t_max=1; % Max time
xv=linspace(x_min,x_max,50); % Spatial vector 
dx=xv(2)-xv(1); % x step size
dt=0.5*(dx^2)/(2*D); % t stable step size
ic=@(x)x.*(1-x);

% Fourier series solution.
a0=0;
% Series for the domain [0,1]. Requires recalculation should domain change.
series=@(x,t,n)2*(-(pi*n*sin(pi*n)+2*cos(pi*n)-2)/(pi^3*n^3))*sin(n*pi*x)*exp(-(n*pi)^2*t);

% Fourier Series Runtime
figure
N=1000;
timeA=zeros(N,1);
for i=1:N
    tic;
    seriesSum(series,xv,0,i);
    timeA(i)=toc;
end
loglog(1:N,timeA);
title('Runtime of Fourier Series Sum')
xlabel('Fourier Series Sum up to N terms')
ylabel('Runtime')
legend('Fourier Series Runtime')

% FDM Runtime
figure
N=100;
timeB=zeros(N,1);
for i=1:N
    x_min=0; % Min domain
    x_max=1; % Max domain
    t_max=1; % Max time
    xv=linspace(x_min,x_max,N); % Spatial vector 
    dx=xv(2)-xv(1); % x step size
    dt=0.5*(dx^2)/(2*D); % t stable step size
    tic;
    fdmU=FDM(1,dx,dt,[x_min x_max],[0 t_max],[0 0],ic);
    timeB(i)=toc;
end
loglog(1:N,timeB);
title('Runtime of Finite Difference Scheme')
xlabel('Spatial Grid Points N')
ylabel('Runtime')
legend('Finite Difference Scheme Runtime')

% pdepe() Runtime
figure
N=100;
timeC=zeros(N,1);
for i=1:N
    x_min=0; % Min domain
    x_max=1; % Max domain
    t_max=1; % Max time
    xv=linspace(x_min,x_max,N); % Spatial vector 
    dx=xv(2)-xv(1); % x step size
    dt=0.5*(dx^2)/(2*D); % t stable step size
    tv=0:dt:t_max; % Time vector
    tic;
    pdepe(0,@pdex1pde,@pdex1ic,@pdex1bc,xv,tv);
    timeC(i)=toc;
end
loglog(1:N,timeC);
title('Runtime of pdepe()')
xlabel('Spatial Grid Points N')
ylabel('Runtime')
legend('pdepe() Runtime')

% Pdepe function handlers
function [c,f,s]=pdex1pde(x,t,u,DuDx)
c=1;
f=DuDx;
s=0;

function u0=pdex1ic(x)
u0=x*(1-x);

function [pl,ql,pr,qr]=pdex1bc(xl,ul,xr,ur,t)
pl=0;
ql=ul;
pr=0;
qr=ur;