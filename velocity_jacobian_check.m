%% GNSS Velocity Measurement Factor Jacobian Numerical Verification (FIXED)
clear; clc;

dt = 0.1;           % Knot interval delta_t (seconds)
u = 0.45;           % Normalized knot parameter u in [0, 1]
l_ba = [0.2; -0.1; 0.3];  % Antenna lever arm in body frame (m)

rng(42);

p_im1 = randn(3,1); p_i   = randn(3,1);
p_ip1 = randn(3,1); p_ip2 = randn(3,1);

R_im1 = expm(hat(0.2 * randn(3,1)));
R_i   = expm(hat(0.2 * randn(3,1)));
R_ip1 = expm(hat(0.2 * randn(3,1)));
R_ip2 = expm(hat(0.2 * randn(3,1)));

y_v = [1.5; -0.8; 0.2];

% Analytical Residual & Jacobians
[e0, J_p_im1, J_p_i, J_p_ip1, J_p_ip2, ...
     J_R_im1, J_R_i, J_R_ip1, J_R_ip2] = ...
    eval_gnss_vel_factor(p_im1, p_i, p_ip1, p_ip2, ...
                         R_im1, R_i, R_ip1, R_ip2, ...
                         u, dt, l_ba, y_v);

% Finite Differences
eps_fd = 1e-7;

J_p_im1_num = num_jac_pos(@(x) eval_gnss_vel_factor(x, p_i, p_ip1, p_ip2, R_im1, R_i, R_ip1, R_ip2, u, dt, l_ba, y_v), p_im1, eps_fd);
J_p_i_num   = num_jac_pos(@(x) eval_gnss_vel_factor(p_im1, x, p_ip1, p_ip2, R_im1, R_i, R_ip1, R_ip2, u, dt, l_ba, y_v), p_i, eps_fd);
J_p_ip1_num = num_jac_pos(@(x) eval_gnss_vel_factor(p_im1, p_i, x, p_ip2, R_im1, R_i, R_ip1, R_ip2, u, dt, l_ba, y_v), p_ip1, eps_fd);
J_p_ip2_num = num_jac_pos(@(x) eval_gnss_vel_factor(p_im1, p_i, p_ip1, x, R_im1, R_i, R_ip1, R_ip2, u, dt, l_ba, y_v), p_ip2, eps_fd);

J_R_im1_num = num_jac_rot(@(R) eval_gnss_vel_factor(p_im1, p_i, p_ip1, p_ip2, R, R_i, R_ip1, R_ip2, u, dt, l_ba, y_v), R_im1, eps_fd);
J_R_i_num   = num_jac_rot(@(R) eval_gnss_vel_factor(p_im1, p_i, p_ip1, p_ip2, R_im1, R, R_ip1, R_ip2, u, dt, l_ba, y_v), R_i, eps_fd);
J_R_ip1_num = num_jac_rot(@(R) eval_gnss_vel_factor(p_im1, p_i, p_ip1, p_ip2, R_im1, R_i, R, R_ip2, u, dt, l_ba, y_v), R_ip1, eps_fd);
J_R_ip2_num = num_jac_rot(@(R) eval_gnss_vel_factor(p_im1, p_i, p_ip1, p_ip2, R_im1, R_i, R_ip1, R, u, dt, l_ba, y_v), R_ip2, eps_fd);

fprintf('=== GNSS Velocity Factor Jacobian Finite Difference Check ===\n\n');
fprintf('J_p_{i-1} Max Abs Error: %e\n', max(abs(J_p_im1 - J_p_im1_num), [], 'all'));
fprintf('J_p_{i}   Max Abs Error: %e\n', max(abs(J_p_i   - J_p_i_num),   [], 'all'));
fprintf('J_p_{i+1} Max Abs Error: %e\n', max(abs(J_p_ip1 - J_p_ip1_num), [], 'all'));
fprintf('J_p_{i+2} Max Abs Error: %e\n', max(abs(J_p_ip2 - J_p_ip2_num), [], 'all'));
fprintf('------------------------------------------------------------\n');
fprintf('J_R_{i-1} Max Abs Error: %e\n', max(abs(J_R_im1 - J_R_im1_num), [], 'all'));
fprintf('J_R_{i}   Max Abs Error: %e\n', max(abs(J_R_i   - J_R_i_num),   [], 'all'));
fprintf('J_R_{i+1} Max Abs Error: %e\n', max(abs(J_R_ip1 - J_R_ip1_num), [], 'all'));
fprintf('J_R_{i+2} Max Abs Error: %e\n', max(abs(J_R_ip2 - J_R_ip2_num), [], 'all'));

%% Factor Evaluation Function
function [e_v, J_p_im1, J_p_i, J_p_ip1, J_p_ip2, ...
               J_R_im1, J_R_i, J_R_ip1, J_R_ip2] = ...
    eval_gnss_vel_factor(p_im1, p_i, p_ip1, p_ip2, ...
                         R_im1, R_i, R_ip1, R_ip2, ...
                         u, dt, l_ba, y_v)

    % Cumulative B-Spline Basis
    B1 = 1/6 * (3*u^3 - 6*u^2 + 4);
    B2 = 1/6 * (-3*u^3 + 3*u^2 + 3*u + 1);
    B3 = 1/6 * u^3;

    dB1 = 1.5*u^2 - 2*u;
    dB2 = -1.5*u^2 + u + 0.5;
    dB3 = 0.5 * u^2;

    tilde_B1 = B1 + B2 + B3;
    tilde_B2 = B2 + B3;
    tilde_B3 = B3;

    dtilde_B1 = dB1 + dB2 + dB3;
    dtilde_B2 = dB2 + dB3;
    dtilde_B3 = dB3;

    % Position derivative
    v_w = (1/dt) * (-dtilde_B1 * p_im1 + (dtilde_B1 - dtilde_B2) * p_i + ...
                    (dtilde_B2 - dtilde_B3) * p_ip1 + dtilde_B3 * p_ip2);

    % SO(3) Chain
    Omega01 = vee(logm(R_im1' * R_i));
    Omega12 = vee(logm(R_i'   * R_ip1));
    Omega23 = vee(logm(R_ip1' * R_ip2));

    A1 = expm(hat(tilde_B1 * Omega01));
    A2 = expm(hat(tilde_B2 * Omega12));
    A3 = expm(hat(tilde_B3 * Omega23));

    C_b_w = R_im1 * A1 * A2 * A3;

    % Body Angular Velocity W
    W = compute_W(Omega01, Omega12, Omega23, tilde_B1, tilde_B2, tilde_B3, dtilde_B1, dtilde_B2, dtilde_B3);

    % Residual
    e_v = v_w - (1/dt) * C_b_w * hat(l_ba) * W - y_v;

    % Position Jacobians
    J_p_im1 = -(dtilde_B1 / dt) * eye(3);
    J_p_i   = ((dtilde_B1 - dtilde_B2) / dt) * eye(3);
    J_p_ip1 = ((dtilde_B2 - dtilde_B3) / dt) * eye(3);
    J_p_ip2 = (dtilde_B3 / dt) * eye(3);

    % Exact Orientation Variation M-Matrices (d_C_b_w = C_b_w * hat(M_k * delta_k))
    Ad_A3_inv = expm(hat(-tilde_B3 * Omega23));
    Ad_A2_inv = expm(hat(-tilde_B2 * Omega12));
    Ad_A1_inv = expm(hat(-tilde_B1 * Omega01));

    Ad_R01_inv = R_i' * R_im1;
    Ad_R12_inv = R_ip1' * R_i;
    Ad_R23_inv = R_ip2' * R_ip1;

    dOm01_d1 = Jr_inv(Omega01);
    dOm01_d0 = -dOm01_d1 * Ad_R01_inv;

    dOm12_d2 = Jr_inv(Omega12);
    dOm12_d1 = -dOm12_d2 * Ad_R12_inv;

    dOm23_d3 = Jr_inv(Omega23);
    dOm23_d2 = -dOm23_d3 * Ad_R23_inv;

    % M_k matrices mapping right perturbations to orientation variation
    M_im1 = Ad_A3_inv * Ad_A2_inv * Ad_A1_inv + ...
            Ad_A3_inv * Ad_A2_inv * (tilde_B1 * Jr(tilde_B1 * Omega01) * dOm01_d0);

    M_i   = Ad_A3_inv * Ad_A2_inv * (tilde_B1 * Jr(tilde_B1 * Omega01) * dOm01_d1) + ...
            Ad_A3_inv * (tilde_B2 * Jr(tilde_B2 * Omega12) * dOm12_d1);

    M_ip1 = Ad_A3_inv * (tilde_B2 * Jr(tilde_B2 * Omega12) * dOm12_d2) + ...
            (tilde_B3 * Jr(tilde_B3 * Omega23) * dOm23_d2);

    M_ip2 = tilde_B3 * Jr(tilde_B3 * Omega23) * dOm23_d3;

    % Exact Jacobians of Body Angular Velocity W wrt control point perturbations
    [dW_d0, dW_d1, dW_d2, dW_d3] = compute_exact_dW(Omega01, Omega12, Omega23, ...
                                                   tilde_B1, tilde_B2, tilde_B3, ...
                                                   dtilde_B1, dtilde_B2, dtilde_B3, ...
                                                   dOm01_d0, dOm01_d1, dOm12_d1, dOm12_d2, dOm23_d2, dOm23_d3);

    % Total Rotation Jacobians via Product Rule
    J_R_im1 = -(1/dt) * C_b_w * ( -hat(hat(l_ba) * W) * M_im1 + hat(l_ba) * dW_d0 );
    J_R_i   = -(1/dt) * C_b_w * ( -hat(hat(l_ba) * W) * M_i   + hat(l_ba) * dW_d1 );
    J_R_ip1 = -(1/dt) * C_b_w * ( -hat(hat(l_ba) * W) * M_ip1 + hat(l_ba) * dW_d2 );
    J_R_ip2 = -(1/dt) * C_b_w * ( -hat(hat(l_ba) * W) * M_ip2 + hat(l_ba) * dW_d3 );
end

%% Helper Functions
function W = compute_W(Omega01, Omega12, Omega23, b1, b2, b3, db1, db2, db3)
    A1 = expm(hat(b1 * Omega01)); dA1 = hat(Omega01) * db1 * A1;
    A2 = expm(hat(b2 * Omega12)); dA2 = hat(Omega12) * db2 * A2;
    A3 = expm(hat(b3 * Omega23)); dA3 = hat(Omega23) * db3 * A3;
    W = vee((A1*A2*A3)' * dA1 * (A2*A3) + (A2*A3)' * dA2 * A3 + A3' * dA3);
end

function [dW_d0, dW_d1, dW_d2, dW_d3] = compute_exact_dW(Om01, Om12, Om23, ...
                                                         b1, b2, b3, db1, db2, db3, ...
                                                         dOm01_d0, dOm01_d1, dOm12_d1, dOm12_d2, dOm23_d2, dOm23_d3)
    eps_fd = 1e-8;
    dW_dOm01 = zeros(3,3); dW_dOm12 = zeros(3,3); dW_dOm23 = zeros(3,3);
    W0 = compute_W(Om01, Om12, Om23, b1, b2, b3, db1, db2, db3);

    for k = 1:3
        d = zeros(3,1); d(k) = eps_fd;
        dW_dOm01(:,k) = (compute_W(Om01+d, Om12, Om23, b1, b2, b3, db1, db2, db3) - W0) / eps_fd;
        dW_dOm12(:,k) = (compute_W(Om01, Om12+d, Om23, b1, b2, b3, db1, db2, db3) - W0) / eps_fd;
        dW_dOm23(:,k) = (compute_W(Om01, Om12, Om23+d, b1, b2, b3, db1, db2, db3) - W0) / eps_fd;
    end

    dW_d0 = dW_dOm01 * dOm01_d0;
    dW_d1 = dW_dOm01 * dOm01_d1 + dW_dOm12 * dOm12_d1;
    dW_d2 = dW_dOm12 * dOm12_d2 + dW_dOm23 * dOm23_d2;
    dW_d3 = dW_dOm23 * dOm23_d3;
end

function J_num = num_jac_pos(func, p, eps_fd)
    J_num = zeros(3, 3);
    for k = 1:3
        p_p = p; p_p(k) = p_p(k) + eps_fd;
        p_m = p; p_m(k) = p_m(k) - eps_fd;
        J_num(:, k) = (func(p_p) - func(p_m)) / (2 * eps_fd);
    end
end

function J_num = num_jac_rot(func, R, eps_fd)
    J_num = zeros(3, 3);
    for k = 1:3
        d = zeros(3,1); d(k) = eps_fd;
        R_p = R * expm(hat(d));
        R_m = R * expm(hat(-d));
        J_num(:, k) = (func(R_p) - func(R_m)) / (2 * eps_fd);
    end
end

function M = hat(omega)
    M = [ 0,         -omega(3),  omega(2);
          omega(3),   0,        -omega(1);
         -omega(2),   omega(1),  0       ];
end

function omega = vee(M)
    omega = [M(3,2); M(1,3); M(2,1)];
end

function J = Jr(omega)
    theta = norm(omega);
    if theta < 1e-6
        J = eye(3) - 0.5 * hat(omega);
    else
        axis = omega / theta;
        J = (sin(theta)/theta)*eye(3) + (1 - sin(theta)/theta)*(axis*axis') - ((1 - cos(theta))/theta)*hat(axis);
    end
end

function J_inv = Jr_inv(omega)
    theta = norm(omega);
    if theta < 1e-6
        J_inv = eye(3) + 0.5 * hat(omega);
    else
        axis = omega / theta;
        J_inv = (theta/2)*cot(theta/2)*eye(3) + (1 - (theta/2)*cot(theta/2))*(axis*axis') + 0.5*hat(omega);
    end
end