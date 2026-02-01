clc;
clear;

%% ================= USER PARAMETERS =================
Vd = 800;          % DC link voltage [V]
A  = 400;          % Desired fundamental amplitude [V]
T  = 0.02;         % Fundamental period (50 Hz)
E  = 1e-5;         % Small edge time for switching sequence

max_iter = 100;
tol      = 1e-6;
lambda   = 0.6;    % Damping factor (0 < lambda <= 1)

min_sep  = 1*pi/180;   % Minimum angular separation (1 degree)

%% =============== INITIAL GUESS (rad) ================
alpha = [15; 30; 45] * pi/180;

%% =============== NEWTON–RAPHSON LOOP =================
for iter = 1:max_iter

    a1 = alpha(1);
    a2 = alpha(2);
    a3 = alpha(3);

    %% -------- Normalized SHE Equations --------
    F = zeros(3,1);
    F(1) =  cos(a1) - cos(a2) + cos(a3) - A*pi/(2*Vd); % Fundamental
    F(2) = (1/5)*(cos(5*a1) - cos(5*a2) + cos(5*a3)); % 5th harmonic
    F(3) = (1/7)*(cos(7*a1) - cos(7*a2) + cos(7*a3)); % 7th harmonic

    %% -------- Jacobian Matrix --------
    J = zeros(3,3);

    J(1,1) = -sin(a1);
    J(1,2) =  sin(a2);
    J(1,3) = -sin(a3);

    J(2,1) = -sin(5*a1);
    J(2,2) =  sin(5*a2);
    J(2,3) = -sin(5*a3);

    J(3,1) = -sin(7*a1);
    J(3,2) =  sin(7*a2);
    J(3,3) = -sin(7*a3);

    %% -------- Newton Update (Damped) --------
    delta = J \ F;
    alpha = alpha - lambda * delta;

    %% -------- Enforce Physical Constraints --------
    alpha(1) = max(alpha(1), 0);
    alpha(2) = max(alpha(2), alpha(1) + min_sep);
    alpha(3) = max(alpha(3), alpha(2) + min_sep);
    alpha(3) = min(alpha(3), pi/2);

    %% -------- Convergence Check --------
    if norm(delta) < tol
        break;
    end

    %% -------- Divergence Protection --------
    if any(isnan(alpha)) || any(isinf(alpha))
        error('Newton-Raphson diverged');
    end
end

if iter == max_iter
    warning('Newton-Raphson did not converge');
end

%% ================== RESULTS ==================
alpha = sort(alpha);            % Final sorted angles
alpha_deg = rad2deg(alpha);

disp('Switching angles (Phase A) [deg]:');
disp(alpha_deg);

%% ================== HARMONIC CHECK ==================
a1 = alpha(1); a2 = alpha(2); a3 = alpha(3);

V1 = (2*Vd/pi)*(cos(a1) - cos(a2) + cos(a3));
V5 = (2*Vd/(5*pi))*(cos(5*a1) - cos(5*a2) + cos(5*a3));
V7 = (2*Vd/(7*pi))*(cos(7*a1) - cos(7*a2) + cos(7*a3));

fprintf('\nHarmonic Check:\n');
fprintf('Fundamental V1 = %.2f V\n', V1);
fprintf('5th Harmonic V5 = %.6f V\n', V5);
fprintf('7th Harmonic V7 = %.6f V\n', V7);

%% ================== 3-PHASE ANGLES ==================
alpha_A = alpha_deg;
alpha_B = mod(alpha_A + 120, 360);
alpha_C = mod(alpha_A + 240, 360);

disp('Phase A angles (deg):'); disp(alpha_A);
disp('Phase B angles (deg):'); disp(alpha_B);
disp('Phase C angles (deg):'); disp(alpha_C);

%% ================== TIME SEQUENCE ==================
t = alpha_A / 360 * T;
Th = T / 2;

time_seq = [ ...
    0, ...
    t(1), t(1)+E, ...
    t(2), t(2)+E, ...
    t(3), t(3)+E, ...
    Th-t(3), Th-t(3)+E, ...
    Th-t(2), Th-t(2)+E, ...
    Th-t(1), Th-t(1)+E, ...
    T ...
];

output_seq = [ ...
    0 0 1 1 0 0 1 1 0 0 1 1 0 0 ...
];

disp('Time sequence (s):');
disp(time_seq);
disp('Output sequence:');
disp(output_seq);
