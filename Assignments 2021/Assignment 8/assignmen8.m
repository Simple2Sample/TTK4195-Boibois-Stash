clear;

%% Parameters
l1 = 0.125;
lc1 = 0.063;
m1 = 0.02;
m2 = 0.3;
J1 = 4.7e-5;
J2 = 3.2e-5;
gravity = 9.81;

d1 = m1*lc1^2 + J1 + m2*l1^2 + J2;
d2 = J2;
gamma = gravity*(m1*lc1 + m2*l1);

%% Original system
% x = [q1; q1_dot; q2; q2_dot]
xdot = @(x,u) [x(2); 
               -1/(d1-d2)*(gamma*cos(x(1)) + u);
               x(4);
               1/(d1-d2)*(gamma*cos(x(1)) + d1*u/d2)];

%% Linearized system

A = [0 1 0 0;
     0 0 1 0;
     0 0 0 1;
     0 0 0 0];
B = [0; 0; 0; 1];

z = @(x) [d1*x(1) + d2*x(3);
          d1*x(2) + d2*x(4);
          -gamma*cos(x(1));
          gamma*sin(x(1))*x(2)];

%% References
q1d = @(t) -pi/2 + 0.3*sin(2*t);
q1d_dot = @(t) 0.6*cos(2*t);
q1d_dotdot = @(t) -1.2*sin(2*t);

q2d = @(t) pi/2 + 0.3*sin(2*t);
q2d_dot = @(t) 0.6*cos(2*t);
q2d_dotdot = @(t) -1.2*sin(2*t);

% Returns [zd1; zd2; zd3; zd4; zd4dot]
zd = @(x,xd) [d1*xd(1) + d2*x(3);
              d1*xd(2) + d2*x(4);
              -gamma*cos(xd(1));
              gamma*sin(xd(1))*xd(2);
              gamma*cos(xd(1))*xd(2)^2 + gamma*sin(xd(1))*xd(3)];

%% Control
Q = diag([100 10 1 1]);
R = 1;

K = lqr(A,B,Q,R);         

% ew_place = [-1,-2,-3,-4];
% K = place(A,B,ew_place);

v = @(z,zd) zd(5) - K*(z - zd(1:4));
u = @(x,v) -(d1-d2)/(gamma*sin(x(1)))*(v - gamma*cos(x(1))*x(2)^2 ...
    + gamma^2*sin(x(1))*cos(x(1))/(d1-d2));
%% Simulation

% Euler integration
h = 0.005;
T0 = 0;
T = 10;
N = (T-T0)/h;
t = linspace(T0, T, N);

x = zeros(4,N);
x0 = [-pi/2; 0; 0; 0];
x(:,1) = x0;

u_vec = zeros(1,N);

xref_1 = [q1d(t); q1d_dot(t); q1d_dotdot(t)];
xref_2 = [q2d(t); q2d_dot(t); q2d_dotdot(t)];

for n = 1:N-1
    % Get the current states (x = [q1; q1_dot; q2; q2_dot])
    x_n = x(:,n);
    
    % Compute the transformed states z from the state x
    z_n = z(x_n);
    
    % Get the reference for x
    xref_n = xref_1(:,n);
    
    % Compute the transformed reference for z
    zref_n = zd(x_n, xref_n);
    
    % Find the input signal for the lin. system (v = z4ref_dot - K(z-zref))
    v_n = v(z_n, zref_n);
    
    % Find u from v
    u_n = u(x_n, v_n);
    u_vec(:,n) = u_n;
    
    % Perform Euler integration
    x(:,n+1) = x_n + h*xdot(x_n, u_n);
end

%% Plots
figure();
plot(t, xref_1(1,:), "--", "DisplayName", "Desired", "LineWidth", 2);
hold on;
plot(t, x(1,:), "DisplayName", "Actual", "LineWidth", 2);
% plot(t, xref_2(1,:), "DisplayName", "Desired 2", "LineWidth", 2);
legend();

figure();
plot(t, u_vec, "DisplayName", "Input", "LineWidth", 2);
legend();
