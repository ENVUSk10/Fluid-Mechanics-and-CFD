clc
close 
clear 
%%
domain_length = 1;                      
node_points = 101;
h = domain_length/(node_points-1);

dt = 1;                                
dx = 1;
alpha = 0.25;                          
D = 1;
source = 0;
timestep = 400;
Beta = (dt*D)/(dx*dx);
omega = 1/(3*alpha/(dt*Beta) + 0.5);
twall = 100; % Left side Temperature
t_right = 300;  % Right side Temperature
t_top = 400;
t_bottom  = 200;

f = zeros(node_points,node_points,9);
w = zeros(9,1);
w(1) = 4/9;
w(2:5) = 1/9;
w(6:9) = 1/36;
x = linspace(0,domain_length,node_points);


x = linspace(0,domain_length,node_points);
y = linspace(0,domain_length,node_points);
[X, Y] = meshgrid(x,y);
%%
for time = 1: timestep
    
    rho = f(:,:,1);
    for k = 2:9
        rho = rho + f(:,:,k);
    end
    for k = 1:9
        feq = w(k)*rho;
        f(:,:,k) = omega*feq + (1-omega)*f(:,:,k);
    end
    
    
    % West to East Streaming
    fhelp = f(:,:,2);
    for i = 1:node_points
        for j = 2:node_points
            f(i,j,2) = fhelp(i,j-1);
        end
    end
    
    % North to south Streamimg
    fhelp = f(:,:,3);
    for i = 2:node_points-1
        for j = 1:node_points
            f(i,j,3) = fhelp(i+1,j);
        end
    end
    
    % East to west Streaming
    for i = 1:node_points
        for j = 1:node_points-1
            f(i,j,4) = f(i,j+1,4);
        end
    end
    
    % North to south Streaming
    fhelp = f(:,:,5);
    for i = 2:node_points
        for j = 1:node_points
            f(i,j,5) = fhelp(i-1,j);
        end
    end
    
     % North-east Streaming
     fhelp = f(:,:,6);
     for i = 1:node_points-1
        for j = 2:node_points
            f(i,j,6) = fhelp(i+1,j-1);
        end
     end
    
    % North-west Streaming
    for i = 1:node_points-1
        for j = 1:node_points-1
            f(i,j,7) = f(i+1,j+1,7);
        end
    end
    
    % South-west Streaming
    fhelp = f(:,:,8);
    for i = 2:node_points
        for j = 1:node_points-1
            f(i,j,8) = fhelp(i-1,j+1);
        end
    end
    
    % South-east Streaming
     fhelp = f(:,:,9);
     for i = 2:node_points
        for j = 2:node_points
            f(i,j,9) = fhelp(i-1,j-1);
        end
     end
     
    f(:,1,2) = w(2)*twall + w(4)*twall - f(:,1,4); % Left boundary, T = 1.0.
    f(:,1,6) = w(6)*twall + w(8)*twall - f(:,1,8); % Left boundary, T = 1.0.
    f(:,1,9) = w(9)*twall + w(7)*twall - f(:,1,7); % Left boundary, T = 1.0.
    
    for k = 2:9
        f(1,:,k) = f(2,:,k); % Bottom boundary, adiabatic.
    end
    
    if rem(time,100) == 0
        figure(1);
        contourf(X, Y, rho,'edgecolor','none');
        colorbar
       colormap('jet')
    end
end
rho = f(:,:,1);
for k = 2:9
    rho = rho + f(:,:,k);
end










