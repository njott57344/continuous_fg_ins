%% Numerical Verification of SO(3) x R^3 B-Spline Position Jacobians
clear; clc;

%% 1. Problem Setup & Random Initialization
rng(42); % Fixed seed for reproducibility

% Random initial orientation control points R_w_b_{k} in SO(3)
R0 = expm_so3(rand(3,1)); % R_{b, i-1}^w
R1 = expm_so3(rand(3,1)); % R_{b, i}^w
R2 = expm_so3(rand(3,1)); % R_{b, i+1}^w
R3 = expm_so3(rand(3,1)); % R_{b, i+2}^w

% Random initial position control points p_{k} in R^3
p0 = rand(3,1) * 10; % p_{i-1}
p1 = rand(3,1) * 10; % p_i
p2 = rand(3,1) * 10; % p_{i+1}
p3 = rand(3,1) * 10; % p_{i+2}

% Sensor lever arm in body frame & measurement
l_ba_b = [0.2; -0.1; 0.5];
y_p    = rand(3,1) * 10;

% Spline time parameter u in [0, 1] and cumulative B-spline basis weights
u = 0.45;
B1 = 0.5 * (1 - u)^3;           % Cumulative B-spline weights
B2 = (3*u^3 - 6*u^2 + 4) / 6;
B3 = u^3 / 6;

% -------------------------------------------------------------------------
%% 2. Evaluate Nominal Trajectory & Orientation Chain
% Relative step vectors Omega_{k, k+1} in R^3
Omega01 = logm_so3(R0' * R1);
Omega12 = logm_so3(R1' * R2);
Omega23 = logm_so3(R2' * R3);

% Scaled step rotations A_j in SO(3)
A1 = expm_so3(B1 * Omega01);
A2 = expm_so3(B2 * Omega12);
A3 = expm_so3(B3 * Omega23);

% Complete Orientation at spline time u
C_w_b = R0 * A1 * A2 * A3;

% Complete Position at spline time u
p_u = (1 - B1)*p0 + (B1 - B2)*p1 + (B2 - B3)*p2 + B3*p3;

% -------------------------------------------------------------------------
%% 3. Analytical Jacobians (Rotation & Translation)

% Lie Algebra Jacobians for relative step angles
Jl_01 = Jl(Omega01); Jl_inv_01 = inv(Jl_01); Jr_inv_01 = inv(Jr(Omega01)); Jr_01_B1 = Jr(B1 * Omega01);
Jl_12 = Jl(Omega12); Jl_inv_12 = inv(Jl_12); Jr_inv_12 = inv(Jr(Omega12)); Jr_02_B2 = Jr(B2 * Omega12);
Jl_23 = Jl(Omega23); Jl_inv_23 = inv(Jl_23); Jr_inv_23 = inv(Jr(Omega23)); Jr_03_B3 = Jr(B3 * Omega23);

% --- Rotation Jacobians (R0, R1, R2, R3) ---
M_i_minus_1 = (A1 * A2 * A3)' - (A2 * A3)' * (B1 * Jr_01_B1 * Jl_inv_01);
J_rot_0 = -C_w_b * skew(l_ba_b) * M_i_minus_1;

M_i = (A2 * A3)' * (B1 * Jr_01_B1 * Jr_inv_01) - A3' * (B2 * Jr_02_B2 * Jl_inv_12);
J_rot_1 = -C_w_b * skew(l_ba_b) * M_i;

M_i_plus_1 = A3' * (B2 * Jr_02_B2 * Jr_inv_12) - (B3 * Jr_03_B3 * Jl_inv_23);
J_rot_2 = -C_w_b * skew(l_ba_b) * M_i_plus_1;

M_i_plus_2 = B3 * Jr_03_B3 * Jr_inv_23;
J_rot_3 = -C_w_b * skew(l_ba_b) * M_i_plus_2;

% --- Translation Jacobians (p0, p1, p2, p3) ---
J_trans_0 = (1 - B1) * eye(3);
J_trans_1 = (B1 - B2) * eye(3);
J_trans_2 = (B2 - B3) * eye(3);
J_trans_3 = B3 * eye(3);

% -------------------------------------------------------------------------
%% 4. Numerical Jacobians via Central Differences
eps = 1e-7;

% Storage for numerical Jacobians
J_num_rot   = cell(1, 4);
J_num_trans = cell(1, 4);
for i = 1:4
    J_num_rot{i}   = zeros(3, 3);
    J_num_trans{i} = zeros(3, 3);
end

% Function handle for residual evaluation
res_fn = @(r0, r1, r2, r3, P0, P1, P2, P3) evaluate_position_residual( ...
    r0, r1, r2, r3, P0, P1, P2, P3, B1, B2, B3, l_ba_b, y_p);

for k = 1:3
    d_plus  = zeros(3, 1); d_plus(k)  =  eps;
    d_minus = zeros(3, 1); d_minus(k) = -eps;
    
    % --- Rotation Perturbations ---
    J_num_rot{1}(:, k) = (res_fn(R0*expm_so3(d_plus), R1, R2, R3, p0, p1, p2, p3) - res_fn(R0*expm_so3(d_minus), R1, R2, R3, p0, p1, p2, p3)) / (2*eps);
    J_num_rot{2}(:, k) = (res_fn(R0, R1*expm_so3(d_plus), R2, R3, p0, p1, p2, p3) - res_fn(R0, R1*expm_so3(d_minus), R2, R3, p0, p1, p2, p3)) / (2*eps);
    J_num_rot{3}(:, k) = (res_fn(R0, R1, R2*expm_so3(d_plus), R3, p0, p1, p2, p3) - res_fn(R0, R1, R2*expm_so3(d_minus), R3, p0, p1, p2, p3)) / (2*eps);
    J_num_rot{4}(:, k) = (res_fn(R0, R1, R2, R3*expm_so3(d_plus), p0, p1, p2, p3) - res_fn(R0, R1, R2, R3*expm_so3(d_minus), p0, p1, p2, p3)) / (2*eps);
    
    % --- Translation Perturbations ---
    J_num_trans{1}(:, k) = (res_fn(R0, R1, R2, R3, p0+d_plus, p1, p2, p3) - res_fn(R0, R1, R2, R3, p0+d_minus, p1, p2, p3)) / (2*eps);
    J_num_trans{2}(:, k) = (res_fn(R0, R1, R2, R3, p0, p1+d_plus, p2, p3) - res_fn(R0, R1, R2, R3, p0, p1+d_minus, p2, p3)) / (2*eps);
    J_num_trans{3}(:, k) = (res_fn(R0, R1, R2, R3, p0, p1, p2+d_plus, p3) - res_fn(R0, R1, R2, R3, p0, p1, p2+d_minus, p3)) / (2*eps);
    J_num_trans{4}(:, k) = (res_fn(R0, R1, R2, R3, p0, p1, p2, p3+d_plus) - res_fn(R0, R1, R2, R3, p0, p1, p2, p3+d_minus)) / (2*eps);
end

% -------------------------------------------------------------------------
%% 5. Display Verification Results
J_ana_rot   = {J_rot_0, J_rot_1, J_rot_2, J_rot_3};
J_ana_trans = {J_trans_0, J_trans_1, J_trans_2, J_trans_3};
rot_names   = {'R_{b,i-1}^w (R0)', 'R_{b,i}^w (R1)', 'R_{b,i+1}^w (R2)', 'R_{b,i+2}^w (R3)'};
trans_names = {'p_{i-1} (p0)', 'p_i (p1)', 'p_{i+1} (p2)', 'p_{i+2} (p3)'};

max_errors = zeros(1, 8);

fprintf('====================================================\n');
fprintf('             ROTATION CONTROL POINTS                \n');
fprintf('====================================================\n');
for i = 1:4
    err = max(max(abs(J_ana_rot{i} - J_num_rot{i})));
    max_errors(i) = err;
    fprintf('Max Error [%s]: %.2e\n', rot_names{i}, err);
end

fprintf('\n====================================================\n');
fprintf('           TRANSLATION CONTROL POINTS               \n');
fprintf('====================================================\n');
for i = 1:4
    err = max(max(abs(J_ana_trans{i} - J_num_trans{i})));
    max_errors(i+4) = err;
    fprintf('Max Error [%s]: %.2e\n', trans_names{i}, err);
end

fprintf('====================================================\n');
if max(max_errors) < 1e-5
    fprintf('VERIFICATION SUCCESSFUL: All 8 analytical Jacobians are correct!\n');
else
    fprintf('VERIFICATION FAILED: Mismatch detected.\n');
end

% =========================================================================
%% Helper Functions
% =========================================================================

function e_p = evaluate_position_residual(R0, R1, R2, R3, p0, p1, p2, p3, B1, B2, B3, l_ba_b, y_p)
    O01 = logm_so3(R0' * R1);
    O12 = logm_so3(R1' * R2);
    O23 = logm_so3(R2' * R3);
    
    A1_loc = expm_so3(B1 * O01);
    A2_loc = expm_so3(B2 * O12);
    A3_loc = expm_so3(B3 * O23);
    
    C = R0 * A1_loc * A2_loc * A3_loc;
    p = (1 - B1)*p0 + (B1 - B2)*p1 + (B2 - B3)*p2 + B3*p3;
    
    e_p = p + C * l_ba_b - y_p;
end

function S = skew(v)
    S = [    0, -v(3),  v(2);
          v(3),     0, -v(1);
         -v(2),  v(1),     0];
end

function v = vee(S)
    v = [S(3,2); S(1,3); S(2,1)];
end

function R = expm_so3(w)
    theta = norm(w);
    if theta < 1e-10
        R = eye(3) + skew(w);
    else
        k = w / theta;
        R = eye(3) + sin(theta)*skew(k) + (1 - cos(theta))*(skew(k)^2);
    end
end

function w = logm_so3(R)
    cos_theta = (trace(R) - 1) / 2;
    cos_theta = min(max(cos_theta, -1), 1);
    theta = acos(cos_theta);
    if theta < 1e-10
        w = vee(R - eye(3)) / 2;
    else
        w = theta / (2 * sin(theta)) * vee(R - R');
    end
end

function Jl_mat = Jl(w)
    theta = norm(w);
    if theta < 1e-7
        Jl_mat = eye(3) + 0.5 * skew(w);
    else
        k = w / theta;
        K_sk = skew(k);
        Jl_mat = sin(theta)/theta * eye(3) + ...
                 (1 - sin(theta)/theta) * (k * k') + ...
                 (1 - cos(theta))/theta * K_sk;
    end
end

function Jr_mat = Jr(w)
    Jr_mat = Jl(-w);
end