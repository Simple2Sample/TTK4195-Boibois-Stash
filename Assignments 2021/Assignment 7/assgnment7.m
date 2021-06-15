clear;

%% Definitions

global M C g

M = @(x) 1 - 0.5*cos(x(1));
C = @(x) 0.5*sin(x(1))*x(2);
g = @(x) 4*sin(x(1));

Kp = 10;
Kd = 10;

K = 10;
Lambda = 10;

theta_d = @(t) pi + 0.3*sin(2*t);
theta_d_dot = @(t) 0.6*cos(2*t);
theta_d_dotdot = @(t) -1.2*sin(2*t);

%% Simulation

% Euler integration
h = 0.005;
T0 = 0;
T = 10;
N = (T-T0)/h;
x0 = [pi; 1];
t = linspace(T0, T, N);

x = zeros(2,N);
xref = [theta_d(t); theta_d_dot(t); theta_d_dotdot(t)];
x(:,1) = x0;

for n = 1:N-1    
%     u_n = u_PD(x(:,n), xref(:,n), Kp, Kd);
%     u_n = u_FBL(x(:,n), xref(:,n), Kp, Kd);
    u_n = u_PASSIVE(x(:,n), xref(:,n), Kp, Kd, K, Lambda);
    
%     x(:,n+1) = x(:,n) + h*x_dot(x(:,n), u_n);
    x(:,n+1) = x(:,n) + h*x_dot(x(:,n), u_n);
end
%%
figure();
plot(t, xref(1,:), "DisplayName", "Desired", "LineWidth", 2);
hold on;
plot(t, x(1,:), "DisplayName", "Actual", "LineWidth", 2);
legend();

function u = u_PD(x, xd, Kp, Kd)
    u = -Kp*(x(1)-xd(1)) - Kd*(x(2)-xd(2));
end

function u = u_FBL(x,xd,Kp,Kd)
    global M C g
    u = M(x)*xd(3) + C(x)*x(2) + g(x) -Kp*(x(1)-xd(1)) -Kd*(x(2)-xd(2));
end

function u = u_PASSIVE(x, xd, Kp, Kd, K, Lambda)
    global M C g
    a = xd(3) - Lambda*(x(2) - xd(2));
    b = xd(2) - Lambda*(x(1) - xd(1));
    r = x(2) - xd(2) + Lambda*(x(1) - xd(1));
    
    u = M(x)*a + C(x)*b + g(x) - K*r;
end

function dx = x_dot(x,u)
    global M C g
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = 1/M(x)*(-C(x)*x(2) - g(x) + u);
end