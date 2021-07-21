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
dpdx = 10;
Re = 20;
alpha = 1/(Re*h*h);

u(1)=0;                                 
u(node_points)=1;                       % For couette Flow change this to say 1m/s            

u_new(1)=0;
u_new(node_points)=1;

u_transient = 0;

if u(node_points) == 1
    sprintf("This is Couette flow simulation")
else
    sprintf("This is Plane Poiseulle flow simulation")
end

                             

error_mag=1;                            % Error estimation
error_req=1e-5;                         % Threshold Error
iterations=0;                           % No of iterations
error_record = 0;                       % Array to record Errors
%% SOLUTION
tStart = tic;
while error_mag>error_req;
    for i=2:node_points-1
        u_new(i) = u(i) + (dt/(h*h))*(mu*(u(i-1) + u(i+1) - 2*u(i)) - (dpdx*h*h)); % using the finite difference method
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
tEnd = toc(tStart);
f = msgbox(sprintf('Solution Converged\nTime taken = %g', tEnd));
%% Time Visualization 
time=(1:iterations).*dt;
semilogy(time, error_record)
xlabel('Time(seconds)')
ylabel('Residual Error')
title('Residual Error over time')


plot(u, linspace(0, domain_size, node_points))
xlabel('X')
ylabel('u(y)')




    
    