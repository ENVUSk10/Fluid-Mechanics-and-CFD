%% D1Q2 1-D DIFFUSION
clc 
clear
close
%% DEFINIG THE LATTICE AND PARAMETERS
domain_length = 1;                      
node_points = 101;
h = domain_length/(node_points-1);

dt = 1;                                
dx = 1;
alpha = 0.25;                          
D = 1;
source = 0;
timestep = 30000;
Beta = (dt*D)/(dx*dx);
w = 1/(alpha/(dt*Beta) + 0.5);
temp1 = 100; % Left side Temperature
temp2 = 60;  % Right side Temperature
%% INITIALIZATION
T(node_points) = 0;
f1 = 0.5.*(T);
f2 = 0.5.*(T);
x = linspace(0,domain_length,node_points);
%% SOLUTION
tStart = tic;
for time = 1:timestep
    % Collision
    for i = 1:node_points
      T(i) = f1(i) + f2(i);
      feq1(i) = 0.5.*(T(i));
      feq2(i) = 0.5.*(T(i));
      % source = q/(T(i)*cp)
      f1(i) = (1-w)*f1(i) + w*feq1(i) + source;
      f2(i) = (1-w)*f2(i) + w*feq2(i) + source;
    end
    % Streaming
    f3 = f1;
    for i = 2: node_points-1 
        f1(i) = f3(i-1);
    end
    for i = 2:node_points-1
        f2(i) = f2(i+1);
    end
    f1(1) = temp1 - f2(1); 
    f2(node_points) = temp2 - f1(node_points);
    if rem(time,100) == 0
        figure(1);
        hold on
        plot(x,T)
    end
end
tEnd = toc(tStart);
f = msgbox(sprintf('Solution Converged\nTime taken = %g', tEnd));
