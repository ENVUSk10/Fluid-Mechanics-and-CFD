%% Stokes' First problem
clc
close all
clear all
%% CREATING THE MESH
domain_size=1;                          % Define the domain size
node_points=51;                         % Number of Grid Points
h=domain_size/(node_points-1);          % Size of each element of the mesh
mu = 5;                                 % Viscosity

dt = 0.00001;
nu = 1e-2;                              % Kinematic Viscosity

u(1)=1;                                 % The Velocity is 1m/s
u(node_points)=0;                       % The velocity at the bottomost Layer is 0          

u_new(1)=1;
u_new(node_points)=0;

u_transient = 0;


                             

error_mag=1;                            % Error estimation
error_req=1e-5;                         % Threshold Error
iterations=0;                           % No of iterations
error_record = 0;                       % Array to record Errors
%% SOLUTION
while error_mag>error_req;

    for i=2:node_points-1
        u_new(i) = u(i) + (dt/(h*h))*(nu)*(u(i-1) + u(i+1) - 2*u(i)); % using the finite difference method
        u_transient(iterations+1, 1:node_points) = u_new; % Records velocity w.r.t iterations
    end
    iterations=iterations+1;
    error_mag=0;
    for i=2:node_points-1;
        error_mag=error_mag+abs(u(i) - u_new(i));
        error_record(iterations) = error_mag;
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
        
       plot(u, linspace(0, domain_size, node_points))
       xlabel('X')
       ylabel('u(y)')
    
       title(sprintf('Velocity Profile at\n time = %g',(iterations*dt)))
       
       
       
      
    end
    u=u_new;
end
f = msgbox('Solution Converged');
%% Time Visualization 
time=(1:iterations).*dt;
semilogy(time, error_record)
xlabel('Time(seconds)')
ylabel('Residual Error')
title('Residual Error over time')


plot(u, linspace(0, domain_size, node_points))
xlabel('X')
ylabel('u(y)')
