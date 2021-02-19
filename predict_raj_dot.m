function r_aj_dot = predict_raj_dot(u_aj,r_aj,v_ej,r_ej,r_ea,v_ea)
% Inputs:
%   u_aj       line-of-sight unit vector
%   r_aj       user position
%   v_ej       Cartesian ECEF velocity of satellite j
%   r_ej       Cartesian ECEF position of satellite j
%   r_ea       predicted Cartesian ECEF user position
%   v_ea       predicted Cartesian ECEF user velocity
% Outputs:
%   r_aj_dot   range rates predicted from the user position

c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s

% skew symmetric matrix
A = [0, -omega_ie, 0;
    omega_ie, 0, 0;
    0, 0, 0];
% Sagnac effect compensation matrix
C = [1, omega_ie*r_aj/c, 0;
    -omega_ie*r_aj/c, 1, 0;
    0, 0, 1];

% Predict the range rates from the user position to each satellite
r_aj_dot = u_aj.'*(C*(v_ej.' + A*r_ej.') - (v_ea + A*r_ea));

