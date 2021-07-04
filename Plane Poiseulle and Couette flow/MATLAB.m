%% 2-D Plane Poiseulle Flow and Couette flow
clc
close all
clear all
%% CREATING THE GRIDS
domain_size=1;                          % Define the domain size
node_points=50;                         % Number of Grid Points
h=domain_size/(node_points-1);          % Size of each element of the mesh
mu = 2;                                 % Viscosity

u(1)=0;                                 
u(node_points)=0;                       % For couette Flow change this to say 1m/s            

u_new(1)=0;
u_new(node_points)=0;

dpdx = 1;                               % Pressure Gradient, positive value is favourable pressure gradient here

error_mag=1;                            % Error estimation
error_req=1e-6;                         % Threshold Error
iterations=0;                           % No of iterations
error_record = 0;                       % Array to record Errors
%% SOLUTION
while error_mag>error_req;
    for i=2:node_points-1
        u_new(i) = (u(i-1) + u(i+1) + (h*h*dpdx/mu))/2; % using the finite difference method
        
    end
    iterations=iterations+1;
    error_mag=0;
    for i=2:node_points-1;
        error_mag=error_mag+abs(u(i) - u_new(i));
        error_record(iterations) = error_mag;
    end
     if(rem(iterations, 100)) == 0                        % Plots the residual error
       figure(1);
       hold on
       semilogy(iterations, error_record(iterations), '-ko')
       xlabel('Iterations')
       ylabel('Residual Error')
       
    end
    u=u_new;
end
%% VISUALIZATION
plot(u,linspace(0, domain_size, node_points))               % Velocity profile plot


    








    
    
