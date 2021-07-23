clear all
close all
clc
%% DEFINING THE MESH
L = 0.5;
B = 0.05;
M = 51;
N = 101;
dx = L/(M-1);
dy = B/(N-1);
g = 9.81;
L = 0.5;
Tw = 100;
Tinf = 40;
Tmean = (Tw + Tinf)/2;
beta = 1/(Tmean + 273);
pr = 0.69;
nu = 2e-5;
dt = 5*10^-2;
alpha = nu/pr;
nstep = 2000;
final_iteration = 100;

for i = 1: M
    for j = 1: N
        x(i,j) = (i-1)*dx;
        y(i,j) = (j - 1)*dy;
    end
end
%% Initializing the variables
T = Tinf*ones(M,N);
u = zeros(M,N);
v = zeros(M,N);
p = zeros(M,N);


%% Solving the governing equations
error_req = 1e-6;
tStart = tic;
for step = 1: final_iteration
        
        T_old = T;
        u_old = u;
        v_old = v;
        
        
        u(1, 2:N)= 0;
        u(2:M, 1)=0;
        u(1,1)=0;
        u(:, N)=0;
        
        v(1, 2:N)= 0;
        v(2:M, 1)=0;
        v(1,1)=0;
        
        T(1, 1) = Tw;
        T(2:M,1)= Tw;
        T(:,N)= Tinf;
   
    
    error_mag_T = 1;
    iterations_T = 0;
    while error_mag_T > error_req
        guess = T;
        for i = 2 : M
            for j = 2 : N - 1
       
                aE = 1/dt + u_old(i,j)/dx - v_old(i,j)/dy + 2*alpha/dy^2;
                aW = v_old(i,j)/dy - alpha/dy^2;
                aN = -alpha/dy^2;
                aS = -u_old(i,j)/dx;
                ae = T_old(i,j)/dt;
                dT = 1;
            
                T(i,j) = dT*(ae - aW*T(i, j+1) - aN*T(i, j-1) - aS*T(i-1,j))/aE + (1 - dT)*T_old(i,j);
            end
        end
        iterations_T = iterations_T + 1;
       error_mag_T = (sum(sum(abs(guess - T))));
                         % Plots the residual error
     
       
    end
    error_mag_u = 1;
    iterations_u = 0;
    while error_mag_u > error_req
        guess = u;
        for i = 2 : M
            for j = 2 : N - 1
               
            
                aE= 1/dt + u_old(i,j)/dx - v_old(i,j)/dy + 2*nu/dy^2;
                aW=v_old(i,j)/dy - nu/dy^2;
                aN=-nu/dy^2;
                aS= -u_old(i,j)/dx;
                an=u_old(i, j)/dt + g*beta*(T(i,j) - Tinf);
                dU = 1;
            
                u(i,j) = dU*(an - aW*u(i, j+1) - aN*u(i, j-1) - aS*u(i-1,j))/aE + (1 - dU)*u_old(i,j);
            end
        end
        error_mag_u = (sum(sum(abs(guess - u))));
        
    end   
  
    error_mag_v = 1;
    iterations_v = 0;
    
    while error_mag_v > error_req 
        guess = v;
        for i = 2: M
            for j = 2: N
                aE= 1/dy;
                aW= -1/dy;
                aN= (u(i-1,j) - u(i,j))/dx;
                dV = 1;
                v(i,j) = dV*(aN - aW*v(i, j-1))/aE + (1 - dV)*v_old(i, j);
            end
        end
        error_mag_v = (sum(sum(abs(guess - v))));
        
    
    end    
    p(:, 1) = 0;
    for j = 2: N
        p(:, j) = p(:, j-1) + u(:, j)*dy;
    end
    figure(1)
    temp = 1;
     while temp < M
        hold on
        plot(((u(temp,:) + (temp-1)*dx)), y(temp, :));
        temp = temp + 4;
     end
     figure(2)
     contourf(x,y,T, 21, 'LineStyle', 'none')
     colorbar
     colormap('jet')
     %F(step) = getframe;
     if(rem(step, 50)) == 0
         disp([num2str(step) 'iteration done = ']) 
    end
   
end
tEnd = toc(tStart);
f = msgbox(sprintf('Solution Converged\nTime taken = %g seconds', tEnd));
%% VISUALIZATION
% movie(F)