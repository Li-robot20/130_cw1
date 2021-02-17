function raj_dot = calculate_raj_dot(raj, rej, rea, vej, vea, uaj)
omega_ie = 7.292115*10^-5;
omega = Skew_symmetric([0, 0, omega_ie]);
c = 299792458;
C = [1, omega_ie*raj/c, 0;
    -omega_ie*raj/c, 1, 0;
    0, 0, 1];
raj_dot = uaj.' * (C * (vej + omega * rej) - (vea + omega * rea));
end