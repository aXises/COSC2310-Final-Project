% Finite Difference Method with Dirichlet Boundary Conditions.
% Inputs:
% D - Diffusion Coefficient
% dx - Spatial step size
% dt - time step size
% xbound - Boundaries for x
% tbound - Boundaries for t
% bc - Boundary values at the end points
% ic - Initial condition
% fd - Field bounds, where the constant input area is defined.
% Outputs:
% U - Solution to the diffusion equation.
function U=FDMInitial(D,dx,dt,xbound,tbound,bc,ic,fd)
    x_vec=xbound(1):dx:xbound(2);
    t_vec=tbound(1):dt:tbound(2);
    U=zeros(length(x_vec),length(t_vec));
    U(1,:)=bc(1);
    U(end,:)=bc(2);
    U(:,1)=ic;
    for t=1:length(t_vec)-1
        for x=2:length(x_vec)-1
            xk=x*dx;
            if (xk < 0.2 && xk > 0.1) || (xk < 0.9 && xk > 0.8)
                U(x,t+1)=U(x,t)+(D*dt/dx^2)*(U(x+1,t)-2*U(x,t)+U(x-1,t))+1*dt;
            else
                U(x,t+1)=U(x,t)+(D*dt/dx^2)*(U(x+1,t)-2*U(x,t)+U(x-1,t));
            end
       end
    end
end