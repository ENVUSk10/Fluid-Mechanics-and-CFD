% Variation of The temperature in convection-diffusion problems for different
% Peclet number using Upwind scheme
clc
clear all
close all
%% DEFINIG THE MESH
dom_length = 2;
node_points = 101;
dx = dom_length/(node_points - 1);


U = 0.02;
alpha = 0.01;

pe = [0, 0.5, 1, 2.5, 10];


error_req = 1e-6;
T_record = 0
%% SOLVING USING FINITE VOLUME METHOD (UPWIND SCHEME)
for j = 1: 5
    T(1) = 0;
    T(node_points) = 1;

    T_new(1) = 0;
    T_new(node_points) = 1;

    Pe = pe(j);
    ap = Pe + (2/dx);
    ae = 1/dx;
    aw = Pe + (1/dx);
    error_mag = 1;
    iterations = 0;
    while error_mag > error_req    
        for i = 2: node_points-1
            T_new(i)=(aw*T(i-1) + ae*T(i+1))/(ap);
        end
        error_mag = 0;
        for i = 2: node_points
            error_mag = error_mag + abs(T(i) - T_new(i));
        end
        iterations = iterations + 1;
        T = T_new;     
 
    end
T_record(j, 1:node_points) = T;
end
%% VISUALIZATION
x = ((1:node_points)-1).*dx;
 for i = 1: 5
     figure(1)
     txt = ['Pe = ',num2str(pe(i))];
     hold on
     plot(x, T_record(i, :), 'DisplayName', txt)
     xlabel('Domain Length \rightarrow')
     ylabel('Temperature \rightarrow')
     legend show
 end
        
        
        
        