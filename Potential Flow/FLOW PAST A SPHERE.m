%% FLOW PAST A SPHERE
clc
clear all
close all
%% DEFINING THE GEOMETRY
V = 20;
G = 50;

a = 1 ;
c = -a*2;
b = a*2;
timestep = a*50;
[X, Y] = meshgrid([c:(b-c)/timestep:b],[c:(b-c)/timestep:b]');


for i = 1: timestep+1
  for j = 1 : timestep+1
    if(((X(i,j)*X(i,j)) + (Y(i,j)*Y(i,j))) < a)
      X(i,j) = 0;
      Y(i,j) = 0;
    end
  end
end

rho = sqrt(X.^2+Y.^2);
theta = atan2(Y,X);

% Creation of the streamline function
z = V.*sin(theta).*rho.*(1-(a^2./(rho.^2)));


timestep = 100;
r = ones(1,timestep+1)*a;
t = [0:2*pi/timestep:2*pi];

Xcircle = r.*cos(t);
Ycircle = r.*sin(t);
%% PLOTTING THE STREAMLINES
figure(1);
contourf(X,Y,z,25);
colorbar;
colormap('jet')
hold on;
fill(Xcircle, Ycircle, 'k');
title('Stream Lines');
xlabel('x \rightarrow');
ylabel('y \rightarrow');
axis square;
%% PLOTTING THE VELOCITY VECTORS
x = [-a*2:a/(1.75*3):a*2];
[x] = meshgrid(x);
y = x';

for i = 1: length(x)
  for j = 1 : length(y)
    if(((x(i,j)*x(i,j)) + (y(i,j)*y(i,j))) < a)
      x(i,j) = 0;
      y(i,j) = 0;
    end
  end
end

r=sqrt(x.^2+y.^2);
theta=atan2(y,x);
ur=V*cos(theta).*(1-a^2./(r.^2));
ut=-V*sin(theta).*(1+a^2./(r.^2))+G./(2*pi*r);
u=ur.*cos(theta)-ut.*sin(theta);
v=ur.*sin(theta)+ut.*cos(theta);

figure(2);
contour(X, Y, z, 25)
hold on
fill(Xcircle, Ycircle, 'k')
hold on
quiver(x,y,u,v);
hold on
title('Speed Vectors')
xlabel('x \rightarrow');
ylabel('y \rightarrow');
axis square;
grid off;

%% PLOTTING THE PRESSURE DISTRIBUTION
% Analytical solution for Cp distribution
t=0:.1:2*pi;
cp = 1 - 4*sin(t).^2 + 2* G / (pi*a*V) *sin(t) - (2* G/ (pi*a*V) )^2 ;

% Non lifting solution
cp_sim = 1 - 4*sin(t).^2 ;

% Lift from Kutta Jokowaski theorm
L = - 1.225*V*G;

L = strcat('Kutta Joukowski Lift: ',num2str(L),' [N]');

% plot cp distribution
figure(3);
plot(t,cp,t,cp_sim,'--r');
axis([0 2*pi min(cp) max(cp_sim)]);
title(strvcat('Pressure coefficient around the surface (standard air density)',L));
xlabel('Theta (angle with horizontal)')
ylabel('C_p')
legend('Lifting solution','Symmetrical solution');
grid on;

% Plot cp distribition on cylinder
figure(4);
cpx = cp.*cos(t);
cpy = cp.*sin(t);
fill(Xcircle, Ycircle, 'y');
hold on;
scale = 2; % used to control length of Cp in quiver plot
quiver(a*cos(t),a*sin(t),cpx,cpy,scale);
title('Pressure coefficient around the cylinder surface (standard air density)');
xlabel('x \rightarrow');
ylabel('y \rightarrow');
axis square;





















