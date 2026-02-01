import numpy as np
import warnings

# ================= USER PARAMETERS =================
Vd = 800          # DC link voltage [V]
A  = 400          # Desired fundamental amplitude [V]
T  = 0.02         # Fundamental period (50 Hz)
E  = 1e-5         # Small edge time for switching sequence

max_iter = 100
tol      = 1e-6
lambda_damp = 0.6    # Damping factor (0 < lambda <= 1)

min_sep  = 1 * np.pi / 180   # Minimum angular separation (1 degree)

# =============== INITIAL GUESS (rad) ================
alpha = np.array([15, 30, 45]) * np.pi / 180

# =============== NEWTON–RAPHSON LOOP =================
for iter in range(max_iter):

    a1 = alpha[0]
    a2 = alpha[1]
    a3 = alpha[2]

    # -------- Normalized SHE Equations --------
    F = np.zeros(3)
    F[0] = np.cos(a1) - np.cos(a2) + np.cos(a3) - A * np.pi / (2 * Vd)  # Fundamental
    F[1] = (1/5) * (np.cos(5*a1) - np.cos(5*a2) + np.cos(5*a3))  # 5th harmonic
    F[2] = (1/7) * (np.cos(7*a1) - np.cos(7*a2) + np.cos(7*a3))  # 7th harmonic

    # -------- Jacobian Matrix --------
    J = np.zeros((3, 3))

    J[0, 0] = -np.sin(a1)
    J[0, 1] = np.sin(a2)
    J[0, 2] = -np.sin(a3)

    J[1, 0] = -np.sin(5*a1)
    J[1, 1] = np.sin(5*a2)
    J[1, 2] = -np.sin(5*a3)

    J[2, 0] = -np.sin(7*a1)
    J[2, 1] = np.sin(7*a2)
    J[2, 2] = -np.sin(7*a3)

    # -------- Newton Update (Damped) --------
    try:
        delta = np.linalg.solve(J, F)
    except np.linalg.LinAlgError:
        print("Warning: Singular Jacobian matrix")
        break
    
    alpha = alpha - lambda_damp * delta

    # -------- Enforce Physical Constraints --------
    alpha[0] = max(alpha[0], 0)
    alpha[1] = max(alpha[1], alpha[0] + min_sep)
    alpha[2] = max(alpha[2], alpha[1] + min_sep)
    alpha[2] = min(alpha[2], np.pi/2)

    # -------- Convergence Check --------
    if np.linalg.norm(delta) < tol:
        break

    # -------- Divergence Protection --------
    if np.any(np.isnan(alpha)) or np.any(np.isinf(alpha)):
        raise ValueError('Newton-Raphson diverged')

if iter == max_iter - 1:
    warnings.warn('Newton-Raphson did not converge')

# ================== RESULTS ==================
alpha = np.sort(alpha)            # Final sorted angles
alpha_deg = np.degrees(alpha)

print('Switching angles (Phase A) [deg]:')
print(alpha_deg)

# ================== HARMONIC CHECK ==================
a1 = alpha[0]
a2 = alpha[1]
a3 = alpha[2]

V1 = (2*Vd/np.pi) * (np.cos(a1) - np.cos(a2) + np.cos(a3))
V5 = (2*Vd/(5*np.pi)) * (np.cos(5*a1) - np.cos(5*a2) + np.cos(5*a3))
V7 = (2*Vd/(7*np.pi)) * (np.cos(7*a1) - np.cos(7*a2) + np.cos(7*a3))

print('\nHarmonic Check:')
print(f'Fundamental V1 = {V1:.2f} V')
print(f'5th Harmonic V5 = {V5:.6f} V')
print(f'7th Harmonic V7 = {V7:.6f} V')

# ================== 3-PHASE ANGLES ==================
alpha_A = alpha_deg
alpha_B = np.mod(alpha_A + 120, 360)
alpha_C = np.mod(alpha_A + 240, 360)

print('Phase A angles (deg):')
print(alpha_A)
print('Phase B angles (deg):')
print(alpha_B)
print('Phase C angles (deg):')
print(alpha_C)

# ================== TIME SEQUENCE ==================
t = alpha_A / 360 * T
Th = T / 2

time_seq = np.array([
    0,
    t[0], t[0] + E,
    t[1], t[1] + E,
    t[2], t[2] + E,
    Th - t[2], Th - t[2] + E,
    Th - t[1], Th - t[1] + E,
    Th - t[0], Th - t[0] + E,
    T
])

output_seq = np.array([
    0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0
])

print('Time sequence (s):')
print(time_seq)
print('Output sequence:')
print(output_seq)
