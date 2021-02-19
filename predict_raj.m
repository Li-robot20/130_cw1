function r_aj_new = predict_raj(r_aj, r_ej, r_ea) 
% Inputs:
%   r_aj       user position
%   r_ej       Cartesian ECEF position of satellite j
%   r_ea       predicted Cartesian ECEF user position
% Outputs:
%   r_aj_new   ranges predicted from the user position

c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s

% Sagnac effect compensation matrix
C = [1, omega_ie*r_aj/c, 0;
    -omega_ie*r_aj/c, 1, 0;
    0, 0, 1];

% Predict the ranges from the approximate user position to each satellite
r_aj_new = sqrt((C*r_ej.' - r_ea).' * (C*r_ej.' - r_ea));

end