% Finite Difference Method with Neumann Boundary Conditions.
% Inputs:
% D - Diffusion Coefficient
% dx - Spatial step size
% dt - time step size
% xbound - Boundaries for x
% tbound - Boundaries for t
% bc - Boundary values at the end points
% ic - Initial condition
% Outputs:
% U - Solution to the diffusion equation.
function U=FDMNeumann(D,dx,dt,xbound,tbound,bc,ic)
    x_vec=xbound(1):dx:xbound(2);
    t_vec=tbound(1):dt:tbound(2);
    U=zeros(length(x_vec),length(t_vec));
    U(1,:)=bc(1);
    U(end,:)=bc(2);
    U(:,1)=ic(xbound(1):dx:xbound(2));
    for t=1:length(t_vec)-1
        U(1,t+1)=U(1,t)+(D*dt/dx^2)*(U(2,t)-2*U(1,t)+bc(1));
        for x=2:length(x_vec)-1
            U(x,t+1)=U(x,t)+(D*dt/dx^2)*(U(x+1,t)-2*U(x,t)+U(x-1,t));
        end
        U(x+1,t+1)=U(x+1,t)+(D*dt/dx^2)*(bc(2)-2*U(x+1,t)+U(x,t));
        bc(1)=U(2,t+1);
        bc(2)=U(x,t+1);
    end
end