function [gyro_sol, x_k_plus, P_k_plus] = Gyro_magne(gyro, magne, i, gyro_init, x_kmins1, P_kmins1)

% Constants
deg_to_rad = 0.01745329252; % Degrees to radians conversion factor
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor
c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]);
R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity

%heading
tao_s = 0.5;
S_rg = (1 * deg_to_rad)^2;
S_bgd = 1E-8;
Phi_kmins1 = [1, tao_s;
                0, 1];
Q_kmins1 = [S_rg * tao_s + 1/3 * S_bgd * tao_s^3, 1/2 * S_bgd * tao_s^2;
            1/2 * S_bgd * tao_s^2, S_bgd * tao_s];
        
x_k_mins = Phi_kmins1 * x_kmins1;
P_k_mins = Phi_kmins1 * P_kmins1 * Phi_kmins1.' + Q_kmins1;

H_k = [-1, 0];

tao_m = 3E-6;
R_k = tao_m^2;

K_k = P_k_mins * H_k.' * inv(H_k * P_k_mins * H_k.' + R_k);

gyro_sol = gyro_init * deg_to_rad + tao_s * gyro(i);
magne_sol = magne(i+1) * deg_to_rad;
deltaz_k_mins = magne_sol - gyro_sol + x_k_mins(1);

x_k_plus = x_k_mins + K_k * deltaz_k_mins;
P_k_plus = (eye(2) - K_k * H_k) * P_k_mins;

gyro_sol = gyro_sol - x_k_plus(1);
gyro_sol = gyro_sol * rad_to_deg;

end