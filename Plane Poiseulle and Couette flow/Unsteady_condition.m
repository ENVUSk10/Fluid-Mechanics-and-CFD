%% 2-D  UNSTEADY PLANE POISEULLE FLOW AND COUETTE FLOWow
clc
close all
clear all
%% CREATING THE MESH
domain_size=1;                          % Define the domain size
node_points=100;                        % Number of Grid Points
h=domain_size/(node_points-1);          % Size of each element of the mesh
mu = 2;                                 % Viscosity

dt = 0.00001;
dpdx = -10;
Re = 20;
alpha = 1/(Re*h*h);

u(1, 1:node_points)=1;                  % For couette Flow change this to say 1m/s            
u(node_points,1:node_points)=0;       

u_new(1, 1:node_points)=1;                                 
u_new(node_points,1:node_points)=0;    

v(1,1) = 0;
v(node_points,node_points) = 0;


if u(1, 1:node_points) == 1
    disp("This is Couette flow simulation")
else
    disp("This is Plane Poiseulle flow simulation")
end

x_dom = ((1:node_points)-1).*h;
y_dom = 1-((1:node_points)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);                             

error_mag=1;                            % Error estimation
error_req=1e-2;                         % Threshold Error
iterations=0;                           % No of iterations
error_record = 0;                       % Array to record Errors
%% SOLUTION
tStart = tic;
while error_mag>error_req;
    for i=2:node_points-1
        for j = 1: node_points
            u_new(i,j) = u(i,j) + (dt/(h*h))*(mu*(u(i-1,j) + u(i+1,j) - 2*u(i,j)) - (dpdx*h*h)); % using the finite difference method
        end
     end
    iterations=iterations+1;
    error_mag=0;
    for i=2:node_points
        for j = 2:node_points
            error_mag=error_mag+abs(u(i,j) - u_new(i,j));
            error_record(iterations) = error_mag;
        end
    end
     if(rem(iterations, 1000)) == 0                        % Plots the residual error
       figure(1);
       hold on
       subplot(1,2,1)
       semilogy(iterations, error_record(iterations), '-kx')
       xlabel('Iterations')
       ylabel('Residual Error')
       
      
       figure(1);
       hold on
      
       title(sprintf('Error at iteration = %d', iterations))
       subplot(1,2,2)
        
       plot(u, y_dom)
       xlabel('X')
       ylabel('u(y)')
    
       title(sprintf('Velocity Profile at\n time = %g',(iterations*dt)))
       figure(2)
       hold on
       contourf(X,Y,u, 21, 'LineStyle', 'none')
       colorbar
       colormap('jet')
       xlabel('x')
       ylabel('y')
       title(sprintf('Velocity Contours at\n time = %g',(iterations*dt)))
       
       
      
    end
    u=u_new;
end
tEnd = toc(tStart);
f = msgbox(sprintf('Solution Converged\nTime taken = %g', tEnd));
quiver(X,Y,u,v,2);



    
    
