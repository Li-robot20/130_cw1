function r_aj = calculate_raj(r_aj, r_ej, r_ea)
omega_ie = 7.292115*10^-5;
c = 299792458;
C = [1, omega_ie*r_aj/c, 0;
    -omega_ie*r_aj/c, 1, 0;
    0, 0, 1];
r_aj = sqrt((C*r_ej.' - r_ea).' * (C*r_ej.' - r_ea));
end