syms t real

omega_a = [1; 0; 0];
omega_b = [cos(t); sin(t); 0];
omega_c = [(cos(2*t))^3; (sin(2*t))^3; 0];

rho = 0.5;

T0 = 0;
T = 3;

%% Task 1

v_a = rho*omega_a;
v_b = rho*omega_b;
v_c = rho*omega_c;

d_a = vpa(int(norm(v_a), t, T0, T));
d_b = vpa(int(norm(v_b), t, T0, T));
d_c = vpa(int(norm(v_c), t, T0, T));

%% Task 2

% Skew symmetric matrix of a vector
S = @(x) [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

% Euler integration
h = 0.01;
N = (T-T0)/h;
R0 = eye(3);
t = linspace(T0, T, N);

% a)
R_a = zeros(3, 3, N);
omega_a_body = zeros(3,N);
R_a(:,:,1) = R0;
S_a = S(omega_a);

for n = 2:N
    R_a(:,:,n) = R_a(:,:,n-1) + h*S_a*R_a(:,:,n-1);
    omega_a_body(:,n) = R_a(:,:,n)'*omega_a;
end

%b)
R_b = zeros(3, 3, N);
omega_b_body = zeros(3,N);
R_b(:,:,1) = R0;
S_b = matlabFunction(S(omega_b));
omega_b = matlabFunction(omega_b);

for n = 2:N
    R_b(:,:,n) = R_b(:,:,n-1) + h*S_b(t(n))*R_b(:,:,n-1);
    omega_b_body(:,n) = R_b(:,:,n)'*omega_b(t(n));
end

%c)
R_c = zeros(3, 3, N);
omega_c_body = zeros(3,N);
R_c(:,:,1) = R0;
S_c = matlabFunction(S(omega_c));
omega_c = matlabFunction(omega_c);

for n = 2:N
    R_c(:,:,n) = R_c(:,:,n-1) + h*S_c(t(n))*R_c(:,:,n-1);
    omega_c_body(:,n) = R_c(:,:,n)'*omega_c(t(n));
end

%% Plots

%a)

figure(1)
subplot(3,1,1);
plot(t, omega_a_body(1,:));
title("\omega_a in body frame");
subplot(3,1,2);
plot(t, omega_a_body(2,:));
subplot(3,1,3);
plot(t, omega_a_body(3,:));

%b)

figure(2)
subplot(3,1,1);
plot(t, omega_b_body(1,:));
title("\omega_b in body frame");
subplot(3,1,2);
plot(t, omega_b_body(2,:));
subplot(3,1,3);
plot(t, omega_b_body(3,:));

%c)

figure(3)
subplot(3,1,1);
plot(t, omega_c_body(1,:));
title("\omega_c in body frame");
subplot(3,1,2);
plot(t, omega_c_body(2,:));
subplot(3,1,3);
plot(t, omega_c_body(3,:));
