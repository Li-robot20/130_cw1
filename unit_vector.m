function u_aj = unit_vector(r_aj, r_ej, r_ea)
omega_ie = 7.292115*10^-5;
c = 299792458;
C = [1, omega_ie*r_aj/c, 0;
    -omega_ie*r_aj/c, 1, 0;
    0, 0, 1];
u_aj = (C*r_ej.' - r_ea)/r_aj;
end