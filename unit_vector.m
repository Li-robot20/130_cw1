function u_aj = unit_vector(r_aj, r_ej, r_ea)
% Inputs:
%   r_aj       user position
%   r_ej       Cartesian ECEF position of satellite j
%   r_ea       predicted Cartesian ECEF user position
% Outputs:
%   u_aj       line-of-sight unit vector

c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s

% Sagnac effect compensation matrix
C = [1, omega_ie*r_aj/c, 0;
    -omega_ie*r_aj/c, 1, 0;
    0, 0, 1];

% Compute the line-of-sight unit vector from the user position to each satellite
u_aj = (C*r_ej.' - r_ea)/r_aj;

end