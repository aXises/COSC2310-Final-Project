% Comparison of Errors of pdepe() and FDM
function ErrorComparison
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
tv=0:dt:t_max; % Time vector

% Fourier series solution.
a0=0;
% Series for the domain [0,1]. Requires recalculation should domain change.
series=@(x,t,n)2*(-(pi*n*sin(pi*n)+2*cos(pi*n)-2)/(pi^3*n^3))*sin(n*pi*x)*exp(-(n*pi)^2*t);

% FDM Solution
ic=@(x)x.*(1-x);
fdmU=FDM(1,dx,dt,[x_min x_max],[0 t_max],[0 0],ic);
fdmU=fdmU';

% pdepe Solution
sol=pdepe(0,@pdex1pde,@pdex1ic,@pdex1bc,xv,tv);
% Extract the first solution component as u.
pdepeU=sol(:,:,1);

% Error Figures
figure
N=50;
hold on
for i=1:N
   loglog(xv,abs((seriesSum(series,xv,0,i))-fdmU(1,:)));
end
title('Comparison of FDM Error')
xlabel('Domain x')
ylabel('Error')
legend('Fourier Series Sum 1','Fourier Series Sum 2','Fourier Series Sum 3','Fourier Series Sum 4','Fourier Series Sum ...')
xlim([x_min x_max])

% Error Figures
figure
N=50;
hold on
for i=1:N
   loglog(xv,abs((seriesSum(series,xv,0,i))-pdepeU(2,:)));
end
title('Comparison of pdepe() Error')
xlabel('Domain x')
ylabel('Error')
legend('Fourier Series Sum 1','Fourier Series Sum 2','Fourier Series Sum 3','Fourier Series Sum 4','Fourier Series Sum ...')
xlim([x_min x_max])

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