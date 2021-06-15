%% Parameters

syms x theta dx dtheta real
syms m1 m2 l k g real


%% Task 2

xdotdot = 1/(m1+m2-m2*(cos(theta))^2) * (m2*l*sin(theta)*dtheta ...
    + m2*l*g*cos(theta)*sin(theta) + k*x);

thetadotdot = -cos(theta)/(l*(m1+m2-m2*cos(theta)))*(m2*l*sin(theta) ...
    + m2*l*cos(theta)*sin(theta) + k*x) - g*sin(theta);

q = [x; dx; theta; dtheta];
f = [dx; xdotdot; dtheta; thetadotdot];

dfdq = [diff(f(1),x), diff(f(1),dx), diff(f(1),theta), diff(f(1),dtheta);
    diff(f(2),x), diff(f(2),dx), diff(f(2),theta), diff(f(2),dtheta);
    diff(f(3),x), diff(f(3),dx), diff(f(3),theta), diff(f(3),dtheta);
    diff(f(4),x), diff(f(4),dx), diff(f(4),theta), diff(f(4),dtheta)];
    
A = subs(dfdq, {x,dx,theta,dtheta}, {0,0,pi,0});
B = [0;1;0;0];

CO = [B, A*B, A^2*B, A^3*B];

%% Task 3
Q = diag([10, 10, 10, 10]);
R = 1;

A_num = double(subs(A, {m1, m2, k, l, g}, {1, 0.5, 25, 0.8, 9.81}));


K = lqr(test,B,Q,R);



