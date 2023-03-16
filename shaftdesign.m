clear all
n = ; %design factor
S_e = ; %
S_ut = ; %ultimate strength
K_f =; %fatigue stress concentration factor for bending
K_fs = ; %fatigue stress concentration factor for torsion
M_a = ; %alternating moment
T_a = ; %alternating torque
M_m = ;%midrange moment
T_m = ;midrange torque

% diameter as a function of above parameters
d = (16n/pi (1/S_e * sqrt(4 * (K_f * M_a)^2 + 3(K_fs T_a)^2) + 1/S_ut * sqrt(4 * (K_f * M_m)^2 + 3 * (K_fs * T_m)^2)))^(1/3);