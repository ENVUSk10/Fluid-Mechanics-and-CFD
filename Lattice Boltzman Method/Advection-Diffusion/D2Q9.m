%% D2Q9 2-D  Advection Diffusion
clc
close 
clear 
%% CREATING THE LATTICE
domain_length = 1;                      
grid_points = 101;
h = domain_length/(grid_points-1);

dt = 1;                                
dx = 1;
alpha = 0.25;                          
D = 1;
source = 0;
timestep = 9000;
Beta = (dt*D)/(dx*dx);
omega = 1/(3*alpha/(dt*Beta) + 0.5);
t_left = 100; % Left side Temperature
t_right = 300;  % Right side Temperature
t_top = 400;  % Top side Temperature
t_bottom  = 200;  % Bottom side Temperature

f = zeros(grid_points,grid_points,9);
w = zeros(9,1);
w(1) = 4/9;
w(2:5) = 1/9;
w(6:9) = 1/36;

u = 1;
v = 0.4;
ck  = 1;
cx = [0 1 0 -1 0 1 -1 -1 1];
cy = [0 0 1 0 -1 1 1 -1 -1];

x=((1:grid_points)-1).*h;
y=1-((1:grid_points)-1).*h;
[X, Y] = meshgrid(x,y);
%% SOLVING BOLTZMAN SOLUTION
tStart = tic;
for time = 1: timestep
    
    T = f(:,:,1);
    for k = 2:9
        T = T + f(:,:,k);
    end
    for i = 1:9
        feq(:,:,i) = w(i).*T.*(1+3*(cx(i)*u + cy(i)*v)/ck);
    end
    for k = 1:9
        feq = w(k)*T;
        f(:,:,k) = omega*feq + (1-omega)*f(:,:,k);
    end
    
    % West to East Streaming
    fhelp = f(:,:,2);
    for i = 1:grid_points
        for j = 2:grid_points
            f(i,j,2) = fhelp(i,j-1);
        end
    end
    
    % North to south Streamimg
    fhelp = f(:,:,3);
    for i = 2:grid_points-1
        for j = 1:grid_points
            f(i,j,3) = fhelp(i+1,j);
        end
    end
    
    % East to west Streaming
    for i = 1:grid_points
        for j = 1:grid_points-1
            f(i,j,4) = f(i,j+1,4);
        end
    end
    
    % North to south Streaming
    fhelp = f(:,:,5);
    for i = 2:grid_points
        for j = 1:grid_points
            f(i,j,5) = fhelp(i-1,j);
        end
    end
    
     % North-east Streaming
     fhelp = f(:,:,6);
     for i = 1:grid_points-1
        for j = 2:grid_points
            f(i,j,6) = fhelp(i+1,j-1);
        end
     end
    
    % North-west Streaming
    for i = 1:grid_points-1
        for j = 1:grid_points-1
            f(i,j,7) = f(i+1,j+1,7);
        end
    end
    
    % South-west Streaming
    fhelp = f(:,:,8);
    for i = 2:grid_points
        for j = 1:grid_points-1
            f(i,j,8) = fhelp(i-1,j+1);
        end
    end
    
    % South-east Streaming
     fhelp = f(:,:,9);
     for i = 2:grid_points
        for j = 2:grid_points
            f(i,j,9) = fhelp(i-1,j-1);
        end
     end
    
    f(:,1,2) = w(2)*t_left + w(4)*t_left - f(:,1,4); % Left boundary
    f(:,1,6) = w(6)*t_left + w(8)*t_left - f(:,1,8); % Left boundary
    f(:,1,9) = w(9)*t_left + w(7)*t_left - f(:,1,7); % Left boundary
    
    f(:,grid_points,4) = w(4)*t_right + w(2)*t_right - f(:,grid_points,2); % Right boundary
    f(:,grid_points,8) = w(8)*t_right + w(6)*t_right - f(:,grid_points,6); % Right boundary
    f(:,grid_points,7) = w(7)*t_right + w(9)*t_right - f(:,grid_points,9); % Right boundary
    
    f(1,:,5) = w(5)*t_top + w(3)*t_top - f(1,:,3); % Top boundary 
    f(1,:,8) = w(8)*t_top + w(6)*t_top - f(1,:,6); % Top boundary
    f(1,:,9) = w(9)*t_top + w(7)*t_top - f(1,:,7); % Top boundary
    
    f(grid_points,:,7) = w(7)*t_bottom + w(9)*t_bottom - f(grid_points,:,9); % Bottom boundary
    f(grid_points,:,3) = w(3)*t_bottom + w(5)*t_bottom - f(grid_points,:,5); % Bottom boundary 
    f(grid_points,:,6) = w(6)*t_bottom + w(8)*t_bottom - f(grid_points,:,8); % Bottom boundary 
    
    if rem(time,500) == 0
        figure(1);
        contourf(X, Y, T,40,'Linestyle','none');
        title(sprintf('Temperature Contours after %d timestep', time))
        colorbar
        colormap('jet')
    end
end
tEnd = toc(tStart);
f = msgbox(sprintf('Solution Converged\nTime taken = %g', tEnd));










