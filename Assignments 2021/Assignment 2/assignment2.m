
syms t real

x = 100*cos(t) - 150*cos(2*t);
y = 100*sin(t) - 150*sin(2*t);
z = 5*sin(t);

xdot = diff(x,t);
ydot = diff(y,t);
zdot = diff(z,t);

r = [x,y,z];
rdot = diff(r,t);
rddot = diff(rdot,t);

%% Task 1
d = vpa(int(sqrt(xdot^2 + ydot^2 + zdot^2), t, 0, 2*pi));

%% Task 2
tau = rdot / norm(rdot);
n = diff(tau, t) / norm(diff(tau, t)) * norm(rdot);
b = cross(tau, n);

R = simplify([tau', n', b']);

%% Task 3
R_func = matlabFunction(R);
[eigvecs, eigvals] = eig(R_func(pi/2));

%% Task 4
taudot = diff(tau, t);
ndot = diff(n, t);
bdot = diff(b, t);
omega = 0.5*[cross(tau, taudot) + cross(b, bdot) + cross(b, bdot)];

%% Task 5
Omega = 10;
rod_length = 1;

psi = pi/3 * sin(Omega*t);
ABe = [rod_length*cos(psi), rod_length*sin(psi), 0];
ABedot = diff(ABe, t);

vel = rdot + cross(omega, ABe) + ABedot*[tau; n; b];

epsilon = [0;0;0]; % Placeholder for angular acceleration which I don't 
                   % know how to find

ABeddot = diff(ABedot, t);
acc = rddot + cross(epsilon, AB) + cross(omega, cross(omega, ABe)) ...
    + ABeddot + 2*cross(omega, ABedot*[tau; n; b]);



