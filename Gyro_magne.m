function [gyro_sol, x_k_plus, P_k_plus] = Gyro_magne(gyro, magne, i, gyro_init, x_kmins1, P_kmins1)
% Inputs:
%   gyro                Gyro data at current epoch(turn rate)
%   magne           	Magnetometer data at current epoch
%   i                   Current epoch
%   gyro_init           Initialisation of gyro solution at current epoch
%   x_kmins1            Initialisation of state vector
%   P_kmins1            Initialisation of error covariance
% Outputs:
%   gyro_sol            Corrected heading solution
%   x_k_plus            Updated state vector
%   P_k_plus            Updated error covariance


% Constants
deg_to_rad = 0.01745329252; % Degrees to radians conversion factor
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor

% trasition matrix
tao_s = 0.5;
S_rg = (1 * deg_to_rad)^2;
S_bgd = 1E-8;
Phi_kmins1 = [1, tao_s;
                0, 1];

% noise covariance matrix
Q_kmins1 = [S_rg * tao_s + 1/3 * S_bgd * tao_s^3, 1/2 * S_bgd * tao_s^2;
            1/2 * S_bgd * tao_s^2, S_bgd * tao_s];
 
% propagate the state and error covariance
x_k_mins = Phi_kmins1 * x_kmins1;
P_k_mins = Phi_kmins1 * P_kmins1 * Phi_kmins1.' + Q_kmins1;

% measurement matrix
H_k = [-1, 0];

% measurement noise covariance
tao_m = 3E-6;
R_k = tao_m^2;

% Kalman gain matrix
K_k = P_k_mins * H_k.' * inv(H_k * P_k_mins * H_k.' + R_k);

% measurement innovation vector
gyro_sol = gyro_init * deg_to_rad + tao_s * gyro(i);
magne_sol = magne(i+1) * deg_to_rad;
deltaz_k_mins = magne_sol - gyro_sol + x_k_mins(1);

% update the state vector and error covariance
x_k_plus = x_k_mins + K_k * deltaz_k_mins;
P_k_plus = (eye(2) - K_k * H_k) * P_k_mins;

% correct gyro solution
gyro_sol = gyro_sol - x_k_plus(1);
gyro_sol = gyro_sol * rad_to_deg;

end