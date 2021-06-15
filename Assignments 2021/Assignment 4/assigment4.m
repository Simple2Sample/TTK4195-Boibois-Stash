
%% Parameters for normal DH-convention
q1 = 0;
q2 = 0;
q3 = -pi/3;
q4 = 0;

rad2deg = 180/pi;
theta20 = 0.1192;
d40 = 1.813;
thetas = [q1, q2 + theta20, q3 + pi/2-theta20, 0]*rad2deg;
ds = [2.202, 0, 0, d40 + q4];
as = [0, 1.4, 0.104, 0];
alphas = [pi/2, 0, pi/2, -pi/2]*rad2deg;

%% Parameters for modified DH-convention
q1 = 0;
q2 = 0;
q3 = -pi/3;
q4 = 0;

rad2deg = 180/pi;
theta20 = 0.1192;
d40 = 1.813;
thetas = [q1, q2 + theta20, q3 + pi/2-theta20, 0]*rad2deg;
ds = [2.202, 0, 0, d40 + q4];
as = [0, 0, 1.4, 0.104];
alphas = [0, pi/2, 0, pi/2]*rad2deg;

%% Frames for standard DH convention
figure()

T01 = DH_trans_matrix(thetas(1), ds(1), as(1), alphas(1));
T12 = DH_trans_matrix(thetas(2), ds(2), as(2), alphas(2));
T23 = DH_trans_matrix(thetas(3), ds(3), as(3), alphas(3));
T34 = DH_trans_matrix(thetas(4), ds(4), as(4), alphas(4));

T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;

plot_frame(eye(4));
plot_frame(T01);
plot_frame(T02);
plot_frame(T03);
plot_frame(T04);

%% Frames for modified DH-convention
figure()

T01 = mod_DH_trans_matrix(thetas(1), ds(1), as(1), alphas(1));
T12 = mod_DH_trans_matrix(thetas(2), ds(2), as(2), alphas(2));
T23 = mod_DH_trans_matrix(thetas(3), ds(3), as(3), alphas(3));
T34 = mod_DH_trans_matrix(thetas(4), ds(4), as(4), alphas(4));

T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;

plot_frame(eye(4));
plot_frame(T01);
plot_frame(T02);
plot_frame(T03);
plot_frame(T04);

%% DH parameter to transformation matrix functions
function T = DH_trans_matrix(theta, d, a, alpha)
    [Rz, Tz, Tx, Rx] = DH_sub_matrices(theta, d, a, alpha);
    T = Rz*Tz*Tx*Rx;
end

function T = mod_DH_trans_matrix(theta, d, a, alpha)
    [Rz, Tz, Tx, Rx] = DH_sub_matrices(theta, d, a, alpha);
    T = Rx*Tx*Rz*Tz;
end

function [Rz, Tz, Tx, Rx] = DH_sub_matrices(theta, d, a, alpha)
    % Rotation around z-axis
    Rz = [rotz(theta), zeros(3,1); zeros(1,3), 1];
    % Translation along z-axis
    Tz = eye(4);
    Tz(3,4) = d;
    % Translation around x-axis
    Tx = eye(4);
    Tx(1,4) = a;
    % Rotation around x-axis
    Rx = [rotx(alpha), zeros(3,1); zeros(1,3), 1];
end

%% Plotting function
function plot_frame(T)
    origin = [0; 0; 0; 1]; % homogeneous coordinates
    x_basis = [1; 0; 0; 1];
    y_basis = [0; 1; 0; 1];
    z_basis = [0; 0; 1; 1];
    
    origin_new_frame_hom = T*origin;
    origin_new_frame = origin_new_frame_hom(1:3)/origin_new_frame_hom(4);
    
    x_basis_new_frame_hom = T*x_basis;
    x_basis_new_frame = x_basis_new_frame_hom(1:3)/x_basis_new_frame_hom(4);
    y_basis_new_frame_hom = T*y_basis;
    y_basis_new_frame = y_basis_new_frame_hom(1:3)/y_basis_new_frame_hom(4);
    z_basis_new_frame_hom = T*z_basis;
    z_basis_new_frame = z_basis_new_frame_hom(1:3)/z_basis_new_frame_hom(4);
    
    quiver3(origin_new_frame(1), origin_new_frame(2), origin_new_frame(3),...
        x_basis_new_frame(1)-origin_new_frame(1), ...
        x_basis_new_frame(2)-origin_new_frame(2), ...
        x_basis_new_frame(3)-origin_new_frame(3), ...
        "linewidth", 1.5, "color", "red");
    
    hold on;
    
    quiver3(origin_new_frame(1), origin_new_frame(2), origin_new_frame(3),...
        y_basis_new_frame(1)-origin_new_frame(1), ...
        y_basis_new_frame(2)-origin_new_frame(2), ...
        y_basis_new_frame(3)-origin_new_frame(3), ...
        "linewidth", 1.5, "color", "green");
    
    quiver3(origin_new_frame(1), origin_new_frame(2), origin_new_frame(3),...
        z_basis_new_frame(1)-origin_new_frame(1), ...
        z_basis_new_frame(2)-origin_new_frame(2), ...
        z_basis_new_frame(3)-origin_new_frame(3), ...
        "linewidth", 1.5, "color", "blue");
    
    xlim([-1 5]);
    ylim([-1 5]);
    zlim([-1 5]);
    xlabel("x");
    ylabel("y");
    zlabel("z");
end